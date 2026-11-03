"""gnomAD sites VCF → DuckDB gnomad_af table.

gnomAD population allele frequencies for rarity filtering. Unlike
ClinVar, the gnomAD release is far too large to bundle or auto-download
(tens of GB per chromosome), so this loader takes a VCF the caller
supplies — a region slice, an `af-only` file, or the full per-chromosome
sites VCF. It reads `AF` (overall) and the highest genetic-ancestry-group
AF (popmax): `AF_grpmax` in v4, `AF_popmax` in v2.1.1 and v3.x.

The gnomAD sites VCFs are already decomposed (one ALT per record). A
small region can be pulled directly with tabix-over-HTTPS, e.g.:

    tabix https://storage.googleapis.com/gcp-public-data--gnomad/release/\\
4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr7.vcf.bgz chr7:117480000-117510000

Loaded via `vcfclick annotations load-gnomad <vcf>`; the MCP
`gnomad_lookup` tool and `db trio --gnomad-max-af` resolve against it.
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

import duckdb
import pyarrow as pa
from cyvcf2 import VCF

from annotations.loaders import store_transaction

log = logging.getLogger(__name__)


# Rows parsed per batch: a full per-chromosome sites VCF holds 10^7 to 10^8
# records, so the file is never held in memory whole.
BATCH_ROWS = 1_000_000

# The highest genetic-ancestry-group AF: `AF_grpmax` in gnomAD v4, `AF_popmax`
# in v2.1.1 and v3.x. The first one a record carries is used.
GROUP_MAX_KEYS = ("AF_grpmax", "AF_popmax")


def _af(info, key: str) -> float | None:
    """Read a per-allele AF INFO value as a float, or None. cyvcf2 returns
    a scalar for a decomposed record's Number=A field; tolerate a tuple."""
    v = info.get(key)
    if v is None:
        return None
    if isinstance(v, (list, tuple)):
        v = v[0] if v else None
    try:
        return None if v is None else float(v)
    except (TypeError, ValueError):
        return None


def _group_max(info) -> float | None:
    for key in GROUP_MAX_KEYS:
        v = _af(info, key)
        if v is not None:
            return v
    return None


def iter_gnomad(vcf_path: Path):
    """Stream a gnomAD sites VCF, yielding (chrom, pos, ref, alt, af,
    af_grpmax). Multi-allelic records are skipped — gnomAD sites VCFs are
    decomposed, so a multi-allelic line means a non-standard input. The
    contig is stored under its canonical name (`7` is `chr7`). Only the
    AF keys are read from INFO (a v4 record carries hundreds)."""
    from storage.contigs import canonical_chrom

    vcf = VCF(str(vcf_path))
    try:
        for rec in vcf:
            if not rec.ALT or len(rec.ALT) != 1:
                continue
            info = rec.INFO
            yield (
                canonical_chrom(rec.CHROM),
                int(rec.POS),
                rec.REF,
                rec.ALT[0],
                _af(info, "AF"),
                _group_max(info),
            )
    finally:
        vcf.close()


def _to_arrow(rows: list) -> pa.Table:
    cols = list(zip(*rows)) if rows else [[]] * 6
    return pa.table(
        {
            "chrom": pa.array(cols[0], type=pa.string()),
            "pos": pa.array(cols[1], type=pa.uint32()),
            "ref": pa.array(cols[2], type=pa.string()),
            "alt": pa.array(cols[3], type=pa.string()),
            "af": pa.array(cols[4], type=pa.float64()),
            "af_grpmax": pa.array(cols[5], type=pa.float64()),
        }
    )


def _duckdb_literal(text: str) -> str:
    """A DuckDB string literal (a backslash is literal in DuckDB)."""
    return "'" + text.replace("'", "''") + "'"


def _stage(vcf_path: Path, stage_path: Path) -> tuple[int, int]:
    """Parse the VCF into a scratch DuckDB file at `stage_path`, BATCH_ROWS
    at a time, keeping the first record of each (chrom, pos, ref, alt) in
    file order. Returns (rows kept, rows kept with no group-max AF)."""
    stage = duckdb.connect(str(stage_path))
    try:
        stage.execute(
            "CREATE TABLE raw (seq BIGINT, chrom VARCHAR, pos UINTEGER, "
            "ref VARCHAR, alt VARCHAR, af DOUBLE, af_grpmax DOUBLE)"
        )
        seq = 0
        batch: list[tuple] = []

        def flush() -> None:
            nonlocal seq, batch
            if not batch:
                return
            table = _to_arrow(batch)
            table = table.add_column(
                0, "seq", pa.array(range(seq, seq + len(batch)), pa.int64())
            )
            stage.register("gnomad_batch", table)
            try:
                stage.execute("INSERT INTO raw SELECT * FROM gnomad_batch")
            finally:
                stage.unregister("gnomad_batch")
            seq += len(batch)
            batch = []

        for row in iter_gnomad(vcf_path):
            batch.append(row)
            if len(batch) >= BATCH_ROWS:
                flush()
        flush()
        stage.execute(
            "CREATE TABLE rows AS SELECT chrom, pos, ref, alt, af, af_grpmax "
            "FROM raw QUALIFY row_number() OVER "
            "(PARTITION BY chrom, pos, ref, alt ORDER BY seq) = 1"
        )
        stage.execute("DROP TABLE raw")
        n, no_max = stage.execute(
            "SELECT count(*), count(*) FILTER (WHERE af_grpmax IS NULL) FROM rows"
        ).fetchone()
    finally:
        stage.close()
    return int(n), int(no_max)


def load(vcf_path: Path, replace: bool = False) -> int:
    """Populate the gnomad_af table from a gnomAD sites VCF.

    `replace=True` clears the table first; the default appends, so several
    per-chromosome slices can be loaded incrementally. Returns the number
    of rows written. Within one input file a duplicate (chrom, pos, ref,
    alt) key keeps the first occurrence; re-loading a key already in the
    table replaces it (latest load wins).

    The file is parsed in batches into a scratch DuckDB file first, so memory
    stays flat for a full per-chromosome VCF and the store is held only for
    the swap, one transaction: a file that fails to parse, or holds no
    records, leaves the table as it was.
    """
    with tempfile.TemporaryDirectory(prefix="vcfclick-gnomad-") as tmp:
        stage_path = Path(tmp) / "stage.duckdb"
        n, no_max = _stage(Path(vcf_path), stage_path)
        if not n:
            raise RuntimeError(f"No gnomAD variants parsed from {vcf_path}")
        with store_transaction() as conn:
            conn.execute(
                f"ATTACH {_duckdb_literal(str(stage_path))} AS stage (READ_ONLY)"
            )
            if replace:
                conn.execute("DELETE FROM gnomad_af")
            else:
                # append, replacing any colliding key already loaded
                conn.execute(
                    "DELETE FROM gnomad_af WHERE (chrom, pos, ref, alt) IN "
                    "(SELECT chrom, pos, ref, alt FROM stage.rows)"
                )
            conn.execute(
                "INSERT INTO gnomad_af (chrom, pos, ref, alt, af, af_grpmax) "
                "SELECT chrom, pos, ref, alt, af, af_grpmax FROM stage.rows"
            )

    log.info("[gnomad] loaded %s allele frequencies into gnomad_af", f"{n:,}")
    if no_max:
        log.warning(
            "[gnomad] %s of %s records carry no group-max AF (%s): their "
            "popmax is the overall AF",
            f"{no_max:,}",
            f"{n:,}",
            " or ".join(GROUP_MAX_KEYS),
        )
    return n
