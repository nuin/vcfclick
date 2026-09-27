"""Variants in a gene: shared by `vcfclick db gene` and the MCP tool
`variants_in_gene`.

The gene's coordinates come from the annotation store (GENCODE genes,
loaded with `vcfclick annotations load`); the variants from the cohort
database. Each variant carries its carrier and hom-alt counts from the
genotypes table, so "who has it" is answered in the same call.

Gene coordinates are stored as chr-prefixed names (chr17); a cohort may
use bare names (17). Both spellings are matched.
"""

from __future__ import annotations

import json
from dataclasses import asdict

DEFAULT_LIMIT = 500
COLUMNS = [
    "ingest_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "qual",
    "filter",
    "info_AF",
    "carriers",
    "hom_alt",
]
_INTS = {"pos", "carriers", "hom_alt"}
_FLOATS = {"qual", "info_AF"}


class GeneNotFound(LookupError):
    pass


def chrom_aliases(chrom: str) -> list[str]:
    """chr17 <-> 17 (and chrM <-> MT) so either naming matches."""
    bare = chrom[3:] if chrom.lower().startswith("chr") else chrom
    names = {chrom, bare, f"chr{bare}"}
    if bare in ("M", "MT"):
        names |= {"chrM", "MT", "M"}
    return sorted(names)


def gene_variants_sql(chroms: list[str], start: int, end: int, limit: int) -> str:
    """Portable SQL (chDB and DuckDB) for variants in a region plus carriers."""
    in_list = ", ".join("'" + c.replace("'", "''") + "'" for c in chroms)
    region = f"chrom IN ({in_list}) AND pos BETWEEN {int(start)} AND {int(end)}"
    return (
        "SELECT v.ingest_id AS ingest_id, v.chrom AS chrom, v.pos AS pos, v.ref AS ref, v.alt AS alt, "
        "v.qual AS qual, v.filter AS filter, v.info_AF AS info_AF, "
        "coalesce(g.carriers, 0) AS carriers, coalesce(g.hom_alt, 0) AS hom_alt "
        f"FROM (SELECT * FROM variants WHERE {region}) AS v "
        "LEFT JOIN ("
        "SELECT ingest_id, chrom, pos, ref, alt, count(*) AS carriers, "
        "sum(CASE WHEN gt = 2 THEN 1 ELSE 0 END) AS hom_alt "
        f"FROM genotypes WHERE gt > 0 AND {region} "
        "GROUP BY ingest_id, chrom, pos, ref, alt"
        ") AS g ON v.ingest_id = g.ingest_id AND v.chrom = g.chrom AND v.pos = g.pos "
        "AND v.ref = g.ref AND v.alt = g.alt "
        f"ORDER BY v.pos, v.ref, v.alt, v.ingest_id LIMIT {int(limit) + 1}"
    )


def _typed(col: str, v):
    if v is None:
        return None
    if col in _INTS:
        return int(v)
    if col in _FLOATS:
        return float(v)
    return v


def gene_variants(
    sess, symbol: str, *, flank: int = 0, limit: int = DEFAULT_LIMIT
) -> dict:
    """Variants overlapping `symbol` (+/- flank bp) with carrier counts.

    Raises GeneNotFound if the symbol isn't in the annotation store.
    """
    from annotations.db import position_for_gene

    g = position_for_gene(symbol)
    if g is None:
        raise GeneNotFound(symbol)
    start, end = max(1, int(g.start_pos) - flank), int(g.end_pos) + flank
    sql = gene_variants_sql(chrom_aliases(g.chrom), start, end, limit)
    data = json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]
    truncated = len(data) > limit
    rows = [{c: _typed(c, v) for c, v in zip(COLUMNS, r)} for r in data[:limit]]
    return {
        "gene": asdict(g),
        "region": {"chrom": g.chrom, "start": start, "end": end},
        "sql": sql,
        "columns": COLUMNS,
        "variants": rows,
        "row_count": len(rows),
        "truncated": truncated,
    }
