"""SQL-native stratification of benchmark concordance against the annotation store.

The benchmark *core* is SQL-free; this layer is deliberately the opposite —
stratification is a query. It registers the classified concordance frame in an
in-memory DuckDB, attaches the (DuckDB) annotation store read-only, and reports
recall/precision per biologically meaningful stratum — gnomAD AF bin, ClinVar
significance, gene. This is the capability hap.py's design forecloses.

Counts come out of SQL as integers; all ratios are computed in Python.

NOTE: this module intentionally uses SQL/DuckDB and is excluded from the
benchmark core's no-SQL rule (it is the annotation-query layer, not the matcher).
"""

from __future__ import annotations

from typing import Any

import duckdb

from benchmark.metrics import precision, recall

# gnomAD population-AF buckets (overall af).
_AF_CASE = (
    "CASE WHEN g.af IS NULL THEN 'novel' "
    "WHEN g.af < 0.001 THEN 'rare' "
    "WHEN g.af < 0.05 THEN 'low' "
    "ELSE 'common' END"
)


def _metrics_rows(rows: list[tuple]) -> list[dict[str, Any]]:
    """Turn (stratum, truth_tp, truth_fn, query_tp, query_fp) tuples into
    metric dicts, computing recall/precision in Python from integer counts."""
    out: list[dict[str, Any]] = []
    for stratum, ttp, tfn, qtp, qfp in rows:
        ttp, tfn, qtp, qfp = int(ttp), int(tfn), int(qtp), int(qfp)
        out.append(
            {
                "stratum": stratum,
                "truth_tp": ttp,
                "truth_fn": tfn,
                "query_tp": qtp,
                "query_fp": qfp,
                "recall": recall(ttp, tfn),
                "precision": precision(qtp, qfp),
            }
        )
    return out


def _counts_sql(bucket_expr: str, join_sql: str) -> str:
    return f"""
        WITH j AS (
            SELECT c.side, c.bd, {bucket_expr} AS stratum
            FROM conc c
            {join_sql}
            WHERE c.filter_view = 'ALL'
        )
        SELECT stratum,
            count(*) FILTER (WHERE side='truth' AND bd='TP') AS truth_tp,
            count(*) FILTER (WHERE side='truth' AND bd='FN') AS truth_fn,
            count(*) FILTER (WHERE side='query' AND bd='TP') AS query_tp,
            count(*) FILTER (WHERE side='query' AND bd='FP') AS query_fp
        FROM j
        GROUP BY stratum
        ORDER BY stratum
    """


def _run(concordance, ann_path: str, bucket_expr: str, join_sql: str) -> list[dict]:
    con = duckdb.connect()
    try:
        con.register("conc", concordance)
        con.execute(f"ATTACH '{ann_path}' AS ann (READ_ONLY)")
        rows = con.execute(_counts_sql(bucket_expr, join_sql)).fetchall()
    finally:
        con.close()
    return _metrics_rows(rows)


def stratify_by_gnomad(concordance, ann_path: str) -> list[dict]:
    """Recall/precision per gnomAD AF bin (novel / rare / low / common)."""
    join = (
        "LEFT JOIN ann.gnomad_af g ON c.chrom=g.chrom AND c.pos=g.pos "
        "AND c.ref=g.ref AND c.alt=g.alt"
    )
    return _run(concordance, ann_path, _AF_CASE, join)


def stratify_by_clinvar(concordance, ann_path: str) -> list[dict]:
    """Recall/precision per ClinVar clinical significance."""
    join = (
        "JOIN ann.clinvar_variants v ON c.chrom=v.chrom AND c.pos=v.pos "
        "AND c.ref=v.ref AND c.alt=v.alt"
    )
    return _run(concordance, ann_path, "v.clin_sig", join)


def stratify_by_gene(concordance, ann_path: str) -> list[dict]:
    """Recall/precision per gene (variant position within a gene's range)."""
    join = (
        "JOIN ann.refseq_genes gn ON c.chrom=gn.chrom "
        "AND c.pos BETWEEN gn.start_pos AND gn.end_pos"
    )
    return _run(concordance, ann_path, "gn.gene_symbol", join)


def _read_bed(path: str, stratum: str) -> list[tuple]:
    rows: list[tuple] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), stratum))
    return rows


def stratify_by_regions(concordance, region_beds: dict[str, str]) -> list[dict]:
    """Recall/precision per genome-stratification region set (e.g. low-complexity,
    segdup). Each BED is a named stratum; a variant in overlapping strata counts
    in each, and variants in no stratum fall under 'none'. BED is 0-based
    half-open, matching the variant's 0-based position (`pos - 1`)."""
    import pyarrow as pa

    regions: list[tuple] = []
    for name, bed in region_beds.items():
        regions.extend(_read_bed(bed, name))
    reg = pa.table(
        {
            "chrom": [r[0] for r in regions],
            "s": [r[1] for r in regions],
            "e": [r[2] for r in regions],
            "stratum": [r[3] for r in regions],
        }
    )
    sql = _counts_sql(
        "coalesce(r.stratum, 'none')",
        "LEFT JOIN reg r ON c.chrom=r.chrom AND (c.pos-1) >= r.s AND (c.pos-1) < r.e",
    )
    con = duckdb.connect()
    try:
        con.register("conc", concordance)
        con.register("reg", reg)
        rows = con.execute(sql).fetchall()
    finally:
        con.close()
    return _metrics_rows(rows)


AXES = {
    "gnomad": stratify_by_gnomad,
    "clinvar": stratify_by_clinvar,
    "gene": stratify_by_gene,
}

_COLS = [
    "stratum",
    "truth_tp",
    "truth_fn",
    "query_tp",
    "query_fp",
    "recall",
    "precision",
]


def write_stratified(
    concordance, ann_path: str, axes: list[str], outdir: str
) -> list[str]:
    """Write `stratified_<axis>.csv` for each requested axis; return paths written."""
    import csv
    import os

    written: list[str] = []
    for axis in axes:
        fn = AXES.get(axis)
        if fn is None:
            raise ValueError(f"unknown stratification axis: {axis!r}")
        rows = fn(concordance, ann_path)
        path = os.path.join(outdir, f"stratified_{axis}.csv")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=_COLS)
            w.writeheader()
            w.writerows(rows)
        written.append(path)
    return written
