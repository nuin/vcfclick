"""Annotation service backed by DuckDB.

Holds STATIC PUBLIC reference data: gene coordinates (RefSeq GFF3) and
ClinVar significance. Transcript / exon / CDS hierarchy lands in
annotations/transcripts.py (Phase 2). All tiers ship under the same
OSS license as the engine.

Architectural role: the two stores (ClickHouse for sample data, DuckDB
for reference data) are intentionally separated so the MCP server
composes across them at query time. Reference data updates monthly
(ClinVar) without touching the sample store; sample data grows without
touching the reference store. This separation is *operationally* useful
for everyone, and would also be the basis for a security boundary if
the hosted tier ever takes on regulated workloads.

DuckDB is the right tool here for three reasons:
  1. Embedded — no extra process, no network hop from the MCP server.
  2. Vectorised — gene overlap queries and ClinVar joins are fast.
  3. Public data — the curated DuckDB file can be shipped as a
     downloadable artefact alongside the OSS package.

Data source: RefSeq GFF3 + ClinVar VCF, GRCh38 only. Loaders live
under annotations/loaders/ (not yet written) and are idempotent batch
jobs that run on a cron, not part of the query path.
"""

from __future__ import annotations

import os
from bisect import bisect_right
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import duckdb

DUCKDB_PATH = Path(__file__).parent / "annotations.duckdb"


def _store_path() -> Path:
    """The annotation store path. `VCFCLICK_ANNOTATIONS_DB` overrides the
    bundled default so a custom/shared store (or a test store reachable by
    subprocesses) can be used."""
    override = os.environ.get("VCFCLICK_ANNOTATIONS_DB")
    return Path(override) if override else DUCKDB_PATH


SCHEMA_DDL = """
CREATE TABLE IF NOT EXISTS refseq_genes (
    gene_symbol  VARCHAR PRIMARY KEY,    -- HGNC symbol, e.g. 'BRCA1'
    chrom        VARCHAR NOT NULL,       -- 'chr17' style
    start_pos    UINTEGER NOT NULL,      -- 1-based, inclusive
    end_pos      UINTEGER NOT NULL,      -- 1-based, inclusive
    strand       VARCHAR,                -- '+' or '-'
    refseq_id    VARCHAR,                -- NCBI Gene ID
    description  VARCHAR
);

CREATE INDEX IF NOT EXISTS idx_genes_range
    ON refseq_genes (chrom, start_pos, end_pos);

CREATE TABLE IF NOT EXISTS clinvar_variants (
    chrom         VARCHAR NOT NULL,
    pos           UINTEGER NOT NULL,
    ref           VARCHAR NOT NULL,
    alt           VARCHAR NOT NULL,
    clin_sig      VARCHAR,           -- 'Pathogenic', 'Likely_pathogenic', ...
    review_status VARCHAR,           -- ClinVar gold-stars equivalent
    clinvar_id    VARCHAR,           -- VCV accession
    condition     VARCHAR,           -- semicolon-joined trait names
    PRIMARY KEY (chrom, pos, ref, alt)
);

CREATE TABLE IF NOT EXISTS gnomad_af (
    chrom     VARCHAR NOT NULL,
    pos       UINTEGER NOT NULL,
    ref       VARCHAR NOT NULL,
    alt       VARCHAR NOT NULL,
    af        DOUBLE,                -- overall population allele frequency
    af_grpmax DOUBLE,                -- highest genetic-ancestry-group AF (popmax)
    PRIMARY KEY (chrom, pos, ref, alt)
);
"""


def get_connection() -> duckdb.DuckDBPyConnection:
    """Open (and initialise on first use) the DuckDB annotation store."""
    conn = duckdb.connect(str(_store_path()))
    conn.execute(SCHEMA_DDL)
    from annotations.transcripts import TRANSCRIPTS_SCHEMA_DDL

    conn.execute(TRANSCRIPTS_SCHEMA_DDL)
    return conn


@dataclass(frozen=True)
class GeneRange:
    gene_symbol: str
    chrom: str
    start_pos: int
    end_pos: int
    strand: str | None


@dataclass(frozen=True)
class ClinVarRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    clin_sig: str | None
    review_status: str | None
    clinvar_id: str | None
    condition: str | None


@dataclass(frozen=True)
class GnomadAF:
    chrom: str
    pos: int
    ref: str
    alt: str
    af: float | None
    af_grpmax: float | None

    @property
    def popmax(self) -> float | None:
        """The rarity-relevant frequency: the highest-ancestry-group AF
        (popmax) when present, else the overall AF."""
        return self.af_grpmax if self.af_grpmax is not None else self.af


def position_for_gene(symbol: str) -> GeneRange | None:
    """Translate a gene symbol to GRCh38 coordinates.

    Used by the LLM to convert questions like "calls in BRCA1" into a
    range filter on the ClickHouse genotypes table. The match ignores case
    (the GENCODE loader stores a symbol as HGNC spells it, so `C9orf72` is
    stored with its lower-case letters); a symbol stored in exactly the
    spelling asked for wins over one that differs only in case.
    """
    conn = get_connection()
    try:
        row = conn.execute(
            """
            SELECT gene_symbol, chrom, start_pos, end_pos, strand
            FROM refseq_genes
            WHERE upper(gene_symbol) = upper(?)
            ORDER BY (gene_symbol = ?) DESC, gene_symbol
            LIMIT 1
            """,
            [symbol, symbol],
        ).fetchone()
    finally:
        conn.close()
    return GeneRange(*row) if row else None


def gene_at(chrom: str, pos: int) -> list[GeneRange]:
    """All genes overlapping a single position. Multiple results possible
    (overlapping transcripts, antisense genes). `chrom` may be spelt as the
    VCF spells it (`17`, `chr17`, `X`, `chr23`): every alias is matched
    (storage.contigs.chrom_aliases)."""
    from storage.contigs import chrom_aliases

    names = list(chrom_aliases(chrom))
    marks = ", ".join("?" for _ in names)
    conn = get_connection()
    try:
        rows = conn.execute(
            f"""
            SELECT gene_symbol, chrom, start_pos, end_pos, strand
            FROM refseq_genes
            WHERE chrom IN ({marks}) AND ? BETWEEN start_pos AND end_pos
            ORDER BY start_pos
            """,
            [*names, pos],
        ).fetchall()
    finally:
        conn.close()
    return [GeneRange(*r) for r in rows]


class GeneIndex:
    """Every gene interval, in memory, for a caller that looks up many
    positions. `gene_at` opens the store and runs a query per position
    (about 0.1 ms each once the store is open, a few ms if it is opened
    every time); `GeneIndex.at` answers the same question from a sorted
    array in about a microsecond. The table is read once, when the index
    is built."""

    def __init__(self, conn: duckdb.DuckDBPyConnection | None = None) -> None:
        own = conn is None
        if conn is None:
            conn = get_connection()
        try:
            rows = conn.execute(
                "SELECT gene_symbol, chrom, start_pos, end_pos, strand "
                "FROM refseq_genes ORDER BY chrom, start_pos"
            ).fetchall()
        finally:
            if own:
                conn.close()
        # keyed by the canonical contig name (storage.contigs.canonical_chrom),
        # so a VCF that says `17` or `X` finds the store's `chr17` and `chrX`
        from storage.contigs import canonical_chrom

        self._canonical = canonical_chrom
        self._genes: dict[str, list[GeneRange]] = {}
        for row in rows:
            self._genes.setdefault(canonical_chrom(row[1]), []).append(GeneRange(*row))
        for genes in self._genes.values():  # two spellings in one store: re-sort
            genes.sort(key=lambda g: g.start_pos)
        # per chromosome: start of each gene, and the furthest end reached
        # by it or any gene before it (so a backwards scan can stop)
        self._starts = {c: [g.start_pos for g in gs] for c, gs in self._genes.items()}
        self._reach: dict[str, list[int]] = {}
        for chrom, genes in self._genes.items():
            far, reach = 0, []
            for g in genes:
                far = max(far, g.end_pos)
                reach.append(far)
            self._reach[chrom] = reach

    def at(self, chrom: str, pos: int) -> list[GeneRange]:
        """Same genes as `gene_at(chrom, pos)`, ordered by start. `chrom` in
        any spelling of the contig (`7`, `chr7`, `X`, `23`)."""
        chrom = self._canonical(chrom)
        genes = self._genes.get(chrom)
        if not genes:
            return []
        reach = self._reach[chrom]
        hits = []
        i = bisect_right(self._starts[chrom], pos) - 1  # last gene starting <= pos
        while i >= 0 and reach[i] >= pos:
            if genes[i].end_pos >= pos:
                hits.append(genes[i])
            i -= 1
        hits.reverse()
        return hits


# The point lookups below cast the position parameter to the column's type:
# DuckDB binds a Python int as BIGINT and then casts the UINTEGER `pos` column
# to match, which stops it skipping row groups by their min/max. Cast the
# parameter instead and a lookup reads one row group, not the whole table (4.5
# ms against 0.3 ms on 5 million rows). TRY_CAST keeps "no such position"
# (negative, or past 2^32) a miss, not an error.


def clinvar_lookup(chrom: str, pos: int, ref: str, alt: str) -> ClinVarRecord | None:
    """Look up ClinVar significance for a specific allele.

    Returns None if the variant is not in ClinVar — which the caller
    should distinguish from "benign" in user-facing output.
    """
    conn = get_connection()
    row = conn.execute(
        """
        SELECT chrom, pos, ref, alt, clin_sig, review_status,
               clinvar_id, condition
        FROM clinvar_variants
        WHERE chrom = ? AND pos = TRY_CAST(? AS UINTEGER) AND ref = ? AND alt = ?
        """,
        [chrom, pos, ref, alt],
    ).fetchone()
    return ClinVarRecord(*row) if row else None


def gnomad_af(chrom: str, pos: int, ref: str, alt: str) -> GnomadAF | None:
    """Look up the gnomAD allele frequency for a specific allele.

    Returns None if the variant is absent from the loaded gnomAD slice —
    the caller should treat absence as "rare", not as AF 0, since the
    loaded slice may not cover the locus.
    """
    conn = get_connection()
    row = conn.execute(
        "SELECT chrom, pos, ref, alt, af, af_grpmax FROM gnomad_af "
        "WHERE chrom = ? AND pos = TRY_CAST(? AS UINTEGER) AND ref = ? AND alt = ?",
        [chrom, pos, ref, alt],
    ).fetchone()
    return GnomadAF(*row) if row else None


def clinvar_many(
    alleles: Iterable[tuple[str, int, str, str]],
    conn: duckdb.DuckDBPyConnection | None = None,
) -> dict[tuple[str, int, str, str], ClinVarRecord]:
    """`clinvar_lookup` for many alleles in one query, the ClinVar
    counterpart of `gnomad_af_many`: {(chrom, pos, ref, alt): ClinVarRecord}
    for alleles present, none for those absent. When an allele has several
    rows the first one returned wins, as in `clinvar_lookup`."""
    keys = list(dict.fromkeys((c, int(p), r, a) for c, p, r, a in alleles))
    if not keys:
        return {}
    own = conn is None
    if conn is None:
        conn = get_connection()
    try:
        if conn.execute("SELECT 1 FROM clinvar_variants LIMIT 1").fetchone() is None:
            return {}
        import pyarrow as pa

        batch = pa.table(
            {
                "chrom": pa.array([k[0] for k in keys], pa.string()),
                "pos": pa.array([k[1] for k in keys], pa.int64()),
                "ref": pa.array([k[2] for k in keys], pa.string()),
                "alt": pa.array([k[3] for k in keys], pa.string()),
            }
        )
        name = "_vcfclick_clinvar_alleles"
        conn.register(name, batch)
        try:
            rows = conn.execute(
                "SELECT v.chrom, v.pos, v.ref, v.alt, v.clin_sig, v.review_status, "
                "v.clinvar_id, v.condition "
                f"FROM clinvar_variants v JOIN {name} c ON v.chrom = c.chrom "
                "AND v.pos = c.pos AND v.ref = c.ref AND v.alt = c.alt"
            ).fetchall()
        finally:
            conn.unregister(name)
    finally:
        if own:
            conn.close()
    found: dict[tuple[str, int, str, str], ClinVarRecord] = {}
    for row in rows:
        found.setdefault((row[0], int(row[1]), row[2], row[3]), ClinVarRecord(*row))
    return found


def gnomad_af_many(
    alleles: Iterable[tuple[str, int, str, str]],
    conn: duckdb.DuckDBPyConnection | None = None,
) -> dict[tuple[str, int, str, str], GnomadAF]:
    """`gnomad_af` for many alleles in one query: {(chrom, pos, ref, alt):
    GnomadAF} for those present in the loaded slice. An allele absent from
    the slice has no entry, exactly as `gnomad_af` returns None for it.

    A caller that filters every candidate of a trio (thousands of alleles)
    must use this, not `gnomad_af` in a loop: a lookup costs 0.1 to 1.5 ms
    on an open connection, growing with the table, and about 5 ms when it
    has to open the store, so 9,000 candidates took 6 s to 40 s. Here the
    alleles are joined to the table as one Arrow batch (a few ms for 9,000
    against a million-row table). `conn` is the caller's open connection; with
    none, the store is opened and closed here."""
    keys = list(dict.fromkeys((c, int(p), r, a) for c, p, r, a in alleles))
    if not keys:
        return {}
    own = conn is None
    if conn is None:
        conn = get_connection()
    try:
        # the common case: no gnomAD loaded, nothing to look up
        if conn.execute("SELECT 1 FROM gnomad_af LIMIT 1").fetchone() is None:
            return {}
        import pyarrow as pa

        batch = pa.table(
            {
                "chrom": pa.array([k[0] for k in keys], pa.string()),
                "pos": pa.array([k[1] for k in keys], pa.int64()),
                "ref": pa.array([k[2] for k in keys], pa.string()),
                "alt": pa.array([k[3] for k in keys], pa.string()),
            }
        )
        name = "_vcfclick_gnomad_alleles"
        conn.register(name, batch)
        try:
            rows = conn.execute(
                "SELECT g.chrom, g.pos, g.ref, g.alt, g.af, g.af_grpmax "
                f"FROM gnomad_af g JOIN {name} c ON g.chrom = c.chrom "
                "AND g.pos = c.pos AND g.ref = c.ref AND g.alt = c.alt"
            ).fetchall()
        finally:
            conn.unregister(name)
    finally:
        if own:
            conn.close()
    found: dict[tuple[str, int, str, str], GnomadAF] = {}
    for row in rows:
        found.setdefault((row[0], int(row[1]), row[2], row[3]), GnomadAF(*row))
    return found
