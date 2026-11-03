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
    gene_symbol  VARCHAR NOT NULL,       -- HGNC symbol, e.g. 'BRCA1'
    chrom        VARCHAR NOT NULL,       -- 'chr17' style
    start_pos    UINTEGER NOT NULL,      -- 1-based, inclusive
    end_pos      UINTEGER NOT NULL,      -- 1-based, inclusive
    strand       VARCHAR,                -- '+' or '-'
    refseq_id    VARCHAR,                -- Ensembl gene id (GENCODE gene_id)
    description  VARCHAR,
    -- one locus per symbol and contig: a pseudo-autosomal gene (SHOX) has a
    -- chrX row and a chrY row. Stores made before this key have
    -- gene_symbol alone, and a replacing `annotations load` rebuilds them.
    PRIMARY KEY (gene_symbol, chrom)
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

-- The chromosomes gnomad_af has records on, written by the gnomAD loader so
-- a reader does not scan gnomad_af for them. Empty in a store whose rows did
-- not come through the loader: gnomad_chromosomes() then scans.
CREATE TABLE IF NOT EXISTS gnomad_chroms (
    chrom     VARCHAR PRIMARY KEY
);

-- Optional (`vcfclick annotations load-hgnc`): what else a gene is called.
-- Each row maps one token to the current approved HGNC symbol, which is the
-- gene_symbol of refseq_genes. Empty unless loaded.
CREATE TABLE IF NOT EXISTS gene_aliases (
    alias   VARCHAR NOT NULL,        -- the token, as HGNC lists it
    kind    VARCHAR NOT NULL,        -- approved | previous | alias | hgnc_id | entrez
    symbol  VARCHAR NOT NULL         -- the current approved symbol
);
"""


class AnnotationStoreUnavailable(RuntimeError):
    """The annotation store exists but cannot be opened: another process is
    writing it (DuckDB lets one writer or many readers hold a file), or the
    file cannot be read. Distinct from an empty or missing store, which reads
    as "nothing loaded"."""


def _all_ddl() -> str:
    from annotations.transcripts import TRANSCRIPTS_SCHEMA_DDL

    return SCHEMA_DDL + TRANSCRIPTS_SCHEMA_DDL


def ensure_schema(conn: duckdb.DuckDBPyConnection) -> None:
    """Create every annotation table and index the store lacks. Writers only:
    the loaders, through `get_connection`."""
    conn.execute(_all_ddl())


def _table_ddl() -> dict[str, str]:
    """{table: its CREATE TABLE statement}, from the schema DDL."""
    out = {}
    for stmt in _all_ddl().split(";"):
        body = "\n".join(
            line.split("--", 1)[0] for line in stmt.strip().splitlines()
        ).strip()
        head = "CREATE TABLE IF NOT EXISTS "
        if body.startswith(head):
            out[body[len(head) :].split()[0]] = body
    return out


def _unavailable(
    path: Path, error: duckdb.Error, mode: str
) -> AnnotationStoreUnavailable:
    import re

    text = str(error)
    pid = re.search(r"PID (\d+)", text)
    if pid:
        why = (
            f"another process (PID {pid.group(1)}) has it open for writing, "
            "probably a `vcfclick annotations load*`; retry when it finishes"
        )
    else:
        why = text.strip()
    return AnnotationStoreUnavailable(
        f"the annotation store {path} cannot be opened {mode}: {why}"
    )


def get_connection() -> duckdb.DuckDBPyConnection:
    """Open the annotation store for writing and create any table it lacks.

    For loaders and other writers only. DuckDB gives a writer the file to
    itself, so every lookup uses `read_connection` instead."""
    path = _store_path()
    try:
        conn = duckdb.connect(str(path))
    except duckdb.IOException as e:
        raise _unavailable(path, e, "for writing") from e
    ensure_schema(conn)
    return conn


def read_connection() -> duckdb.DuckDBPyConnection:
    """Open the annotation store for lookups: read-only, no DDL, so any number
    of processes read at once and a read-only file still answers.

    A missing store is an empty in-memory one (nothing loaded, no file made).
    A table an older store lacks is an empty temporary table on this
    connection only. Raises AnnotationStoreUnavailable when another process
    is writing the store or the file cannot be read."""
    path = _store_path()
    if not path.exists():
        conn = duckdb.connect(":memory:")
        ensure_schema(conn)
        return conn
    try:
        conn = duckdb.connect(str(path), read_only=True)
    except duckdb.ConnectionException:
        # this process already holds the store open for writing (a loader,
        # or a test seeding it): DuckDB shares that instance, not a second one
        conn = duckdb.connect(str(path))
    except duckdb.IOException as e:
        raise _unavailable(path, e, "for reading") from e
    present = {
        r[0]
        for r in conn.execute(
            "SELECT table_name FROM duckdb_tables() "
            "WHERE database_name = current_database() AND schema_name = 'main'"
        ).fetchall()
    }
    for table, ddl in _table_ddl().items():
        if table not in present:
            conn.execute(
                ddl.replace("CREATE TABLE IF NOT EXISTS", "CREATE TEMP TABLE", 1)
            )
    return conn


def refresh_gnomad_chroms(conn, *, new_rows: str | None = None) -> None:
    """Bring gnomad_chroms up to date after a gnomAD load. With `new_rows`
    (a relation holding the rows just appended) and a table already filled,
    only their chromosomes are added; otherwise the set is rebuilt from
    gnomad_af (after a replace, or the first load into an older store)."""
    filled = conn.execute("SELECT 1 FROM gnomad_chroms LIMIT 1").fetchone()
    if new_rows is not None and filled is not None:
        conn.execute(
            f"INSERT INTO gnomad_chroms SELECT DISTINCT chrom FROM {new_rows} "
            "WHERE chrom NOT IN (SELECT chrom FROM gnomad_chroms)"
        )
        return
    conn.execute("DELETE FROM gnomad_chroms")
    conn.execute("INSERT INTO gnomad_chroms SELECT DISTINCT chrom FROM gnomad_af")


def gnomad_chromosomes() -> frozenset[str]:
    """The chromosomes the loaded gnomAD data has records on: gnomad_chroms,
    or a scan of gnomad_af when that is empty (rows written another way).
    Read-only, as every lookup (`read_connection`)."""
    conn = read_connection()
    try:
        rows = conn.execute("SELECT chrom FROM gnomad_chroms").fetchall()
        if not rows:
            rows = conn.execute("SELECT DISTINCT chrom FROM gnomad_af").fetchall()
    finally:
        conn.close()
    return frozenset(r[0] for r in rows)


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

    A token that is no stored symbol resolves as a `db trio --genes` token
    does (`resolve_gene_tokens`): an Ensembl gene id and, after `annotations
    load-hgnc`, a previous symbol, an alias, an HGNC or Entrez id. A token
    whose best match names two genes is None.
    """
    conn = read_connection()
    try:
        row = conn.execute(
            """
            SELECT gene_symbol, chrom, start_pos, end_pos, strand
            FROM refseq_genes
            WHERE upper(gene_symbol) = upper(?)
            ORDER BY (gene_symbol = ?) DESC, gene_symbol, chrom
            LIMIT 1
            """,
            [symbol, symbol],
        ).fetchone()
        if row:
            return GeneRange(*row)
        matches = resolve_gene_tokens([symbol], conn=conn).get(gene_token_key(symbol))
    finally:
        conn.close()
    if not matches:
        return None
    best = matches[0][0]
    genes = {}
    for kind, g in matches:
        if kind == best:
            genes.setdefault(g.gene_symbol, g)  # the first locus of each
    return next(iter(genes.values())) if len(genes) == 1 else None


def gene_at(chrom: str, pos: int) -> list[GeneRange]:
    """All genes overlapping a single position. Multiple results possible
    (overlapping transcripts, antisense genes). `chrom` may be spelt as the
    VCF spells it (`17`, `chr17`, `X`, `chr23`): every alias is matched
    (storage.contigs.chrom_aliases)."""
    from storage.contigs import chrom_aliases

    names = list(chrom_aliases(chrom))
    marks = ", ".join("?" for _ in names)
    conn = read_connection()
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
            conn = read_connection()
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
#
# Every allele lookup matches the contig under each of its spellings
# (storage.contigs.chrom_aliases, the canonical `chr17` first): the loaders
# store the canonical name, a cohort from an Ensembl-style VCF says `17`, and
# a store loaded before the loaders canonicalised may say either. An answer
# carries the caller's spelling, and the canonical spelling wins a tie.

_CLINVAR_COLS = "chrom, pos, ref, alt, clin_sig, review_status, clinvar_id, condition"
_GNOMAD_COLS = "chrom, pos, ref, alt, af, af_grpmax"


def _spellings(chrom: str) -> list[str]:
    from storage.contigs import canonical_chrom, chrom_aliases

    canonical = canonical_chrom(chrom)
    return list(dict.fromkeys([canonical, chrom, *chrom_aliases(chrom)]))


def _point_lookup(table: str, cols: str, chrom: str, pos: int, ref: str, alt: str):
    names = _spellings(chrom)
    marks = ", ".join("?" for _ in names)
    rank = " ".join(f"WHEN ? THEN {i}" for i in range(len(names)))
    conn = read_connection()
    try:
        row = conn.execute(
            f"SELECT {cols} FROM {table} WHERE chrom IN ({marks}) "
            "AND pos = TRY_CAST(? AS UINTEGER) AND ref = ? AND alt = ? "
            f"ORDER BY CASE chrom {rank} END LIMIT 1",
            [*names, pos, ref, alt, *names],
        ).fetchone()
    finally:
        conn.close()
    return None if row is None else (chrom, *row[1:])


def _many_lookup(
    table: str,
    cols: str,
    alleles: Iterable[tuple[str, int, str, str]],
    conn: duckdb.DuckDBPyConnection | None,
) -> dict[tuple[str, int, str, str], tuple]:
    """{caller's (chrom, pos, ref, alt): row with the caller's chrom} for each
    allele `table` holds, one join for all of them."""
    keys = list(dict.fromkeys((c, int(p), r, a) for c, p, r, a in alleles))
    if not keys:
        return {}
    own = conn is None
    if conn is None:
        conn = read_connection()
    try:
        # the common case: nothing loaded, nothing to look up
        if conn.execute(f"SELECT 1 FROM {table} LIMIT 1").fetchone() is None:
            return {}
        import pyarrow as pa

        idx, rank, chroms, poss, refs, alts = [], [], [], [], [], []
        for i, (c, p, r, a) in enumerate(keys):
            for n, name in enumerate(_spellings(c)):
                idx.append(i)
                rank.append(n)
                chroms.append(name)
                poss.append(p)
                refs.append(r)
                alts.append(a)
        batch = pa.table(
            {
                "idx": pa.array(idx, pa.int64()),
                "rank": pa.array(rank, pa.int32()),
                "chrom": pa.array(chroms, pa.string()),
                "pos": pa.array(poss, pa.int64()),
                "ref": pa.array(refs, pa.string()),
                "alt": pa.array(alts, pa.string()),
            }
        )
        name = f"_vcfclick_{table}_alleles"
        select = ", ".join(f"t.{c.strip()}" for c in cols.split(","))
        conn.register(name, batch)
        try:
            rows = conn.execute(
                f"SELECT c.idx, {select} FROM {table} t JOIN {name} c "
                "ON t.chrom = c.chrom AND t.pos = c.pos AND t.ref = c.ref "
                "AND t.alt = c.alt ORDER BY c.idx, c.rank"
            ).fetchall()
        finally:
            conn.unregister(name)
    finally:
        if own:
            conn.close()
    found: dict[tuple[str, int, str, str], tuple] = {}
    for i, *row in rows:
        key = keys[i]
        found.setdefault(key, (key[0], *row[1:]))
    return found


def clinvar_lookup(chrom: str, pos: int, ref: str, alt: str) -> ClinVarRecord | None:
    """Look up ClinVar significance for a specific allele, `chrom` in any
    spelling of the contig.

    Returns None if the variant is not in ClinVar — which the caller
    should distinguish from "benign" in user-facing output.
    """
    row = _point_lookup("clinvar_variants", _CLINVAR_COLS, chrom, pos, ref, alt)
    return ClinVarRecord(*row) if row else None


def gnomad_af(chrom: str, pos: int, ref: str, alt: str) -> GnomadAF | None:
    """Look up the gnomAD allele frequency for a specific allele, `chrom` in
    any spelling of the contig.

    Returns None if the variant is absent from the loaded gnomAD slice —
    the caller should treat absence as "rare", not as AF 0, since the
    loaded slice may not cover the locus.
    """
    row = _point_lookup("gnomad_af", _GNOMAD_COLS, chrom, pos, ref, alt)
    return GnomadAF(*row) if row else None


def clinvar_many(
    alleles: Iterable[tuple[str, int, str, str]],
    conn: duckdb.DuckDBPyConnection | None = None,
) -> dict[tuple[str, int, str, str], ClinVarRecord]:
    """`clinvar_lookup` for many alleles in one query, the ClinVar
    counterpart of `gnomad_af_many`: {(chrom, pos, ref, alt): ClinVarRecord}
    for alleles present, none for those absent, keyed as the caller spelt
    them. When an allele has several rows the first one returned wins, as in
    `clinvar_lookup`."""
    found = _many_lookup("clinvar_variants", _CLINVAR_COLS, alleles, conn)
    return {k: ClinVarRecord(*row) for k, row in found.items()}


def gnomad_af_many(
    alleles: Iterable[tuple[str, int, str, str]],
    conn: duckdb.DuckDBPyConnection | None = None,
) -> dict[tuple[str, int, str, str], GnomadAF]:
    """`gnomad_af` for many alleles in one query: {(chrom, pos, ref, alt):
    GnomadAF} for those present in the loaded slice, keyed as the caller
    spelt them. An allele absent from the slice has no entry, exactly as
    `gnomad_af` returns None for it.

    A caller that filters every candidate of a trio (thousands of alleles)
    must use this, not `gnomad_af` in a loop: a lookup costs 0.1 to 1.5 ms
    on an open connection, growing with the table, and about 5 ms when it
    has to open the store, so 9,000 candidates took 6 s to 40 s. Here the
    alleles are joined to the table as one Arrow batch (a few ms for 9,000
    against a million-row table). `conn` is the caller's open connection; with
    none, the store is opened and closed here."""
    found = _many_lookup("gnomad_af", _GNOMAD_COLS, alleles, conn)
    return {k: GnomadAF(*row) for k, row in found.items()}


# How a token found its gene, best first: a token that is a current symbol is
# that gene even when it is also another gene's previous symbol.
RESOLVED_BY: tuple[str, ...] = (
    "symbol",
    "ensembl",
    "previous",
    "alias",
    "hgnc_id",
    "entrez",
)
ALIAS_KINDS: tuple[str, ...] = ("previous", "alias", "hgnc_id", "entrez")


def gene_token_key(token: str) -> str:
    """A gene token as it is compared: upper case, an Ensembl gene id without
    its version (`ENSG00000147894.16` is `ENSG00000147894`)."""
    key = token.strip().upper()
    if key.startswith("ENSG") and "." in key:
        key = key.split(".", 1)[0]
    return key


def has_gene_aliases(conn: duckdb.DuckDBPyConnection) -> bool:
    """Is the optional alias table (`annotations load-hgnc`) filled?"""
    try:
        row = conn.execute("SELECT 1 FROM gene_aliases LIMIT 1").fetchone()
    except duckdb.Error:  # an older store without the table
        return False
    return row is not None


def gene_count(conn: duckdb.DuckDBPyConnection) -> int:
    try:
        return int(conn.execute("SELECT count(*) FROM refseq_genes").fetchone()[0])
    except duckdb.Error:
        return 0


def resolve_gene_tokens(
    tokens: Iterable[str],
    conn: duckdb.DuckDBPyConnection | None = None,
) -> dict[str, list[tuple[str, GeneRange]]]:
    """Every gene each token can mean, in one query on one connection:
    {gene_token_key(token): [(resolved_by, GeneRange), ...]}, best match first
    (RESOLVED_BY order), no entry for a token that matches nothing.

    A token matches a current symbol without regard to case, or an Ensembl
    gene id (the refseq_id column, version ignored), and, when `annotations
    load-hgnc` has filled `gene_aliases`, a previous symbol, an alias, an HGNC
    id (`HGNC:1234`) or an Entrez id. A gene with loci on two contigs (a
    pseudo-autosomal gene on chrX and chrY) is one entry per locus. The caller
    decides what an ambiguous token means (one previous symbol of two genes)."""
    keys = list(dict.fromkeys(gene_token_key(t) for t in tokens if t.strip()))
    if not keys:
        return {}
    own = conn is None
    if conn is None:
        conn = read_connection()
    try:
        import pyarrow as pa

        name = "_vcfclick_gene_tokens"
        cols = "g.gene_symbol, g.chrom, g.start_pos, g.end_pos, g.strand"
        parts = [
            f"SELECT t.tok, 'symbol' AS kind, {cols} FROM {name} t "
            "JOIN refseq_genes g ON upper(g.gene_symbol) = t.tok",
            f"SELECT t.tok, 'ensembl' AS kind, {cols} FROM {name} t "
            "JOIN refseq_genes g ON upper(split_part(g.refseq_id, '.', 1)) = t.tok",
        ]
        if has_gene_aliases(conn):
            # an approved HGNC symbol is that gene or nothing: when the gene
            # store lacks it, it does not fall through to a previous symbol or
            # an alias of another gene
            parts.append(
                f"SELECT t.tok, a.kind, {cols} FROM {name} t "
                "JOIN gene_aliases a ON upper(a.alias) = t.tok "
                "AND a.kind <> 'approved' "
                "JOIN refseq_genes g ON upper(g.gene_symbol) = upper(a.symbol) "
                "WHERE t.tok NOT IN (SELECT upper(alias) FROM gene_aliases "
                "WHERE kind = 'approved')"
            )
            # the other way round: an approved symbol the gene store lacks
            # because the gene was renamed after the GENCODE release, which
            # still stores it under a previous symbol of that gene
            parts.append(
                f"SELECT t.tok, 'previous' AS kind, {cols} FROM {name} t "
                "JOIN gene_aliases a ON upper(a.symbol) = t.tok "
                "AND a.kind = 'previous' "
                "JOIN refseq_genes g ON upper(g.gene_symbol) = upper(a.alias) "
                "WHERE t.tok IN (SELECT upper(alias) FROM gene_aliases "
                "WHERE kind = 'approved') "
                "AND t.tok NOT IN (SELECT upper(gene_symbol) FROM refseq_genes)"
            )
        conn.register(name, pa.table({"tok": pa.array(keys, pa.string())}))
        try:
            rows = conn.execute(" UNION ALL ".join(parts)).fetchall()
        finally:
            conn.unregister(name)
    finally:
        if own:
            conn.close()
    rank = {k: i for i, k in enumerate(RESOLVED_BY)}
    found: dict[str, list[tuple[str, GeneRange]]] = {}
    for tok, kind, *gene in rows:
        matches = found.setdefault(tok, [])
        if all(
            (g.gene_symbol, g.chrom, k) != (gene[0], gene[1], kind) for k, g in matches
        ):
            matches.append((kind, GeneRange(*gene)))
    for matches in found.values():
        matches.sort(
            key=lambda m: (rank.get(m[0], len(rank)), m[1].gene_symbol, m[1].chrom)
        )
    return found
