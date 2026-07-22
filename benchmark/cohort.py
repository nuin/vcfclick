"""Multi-caller / cohort-scale benchmarking (beyond hap.py, phase B).

hap.py compares one query against one truth. This layer benchmarks many callers
against a shared truth in one pass, tags each classified row with its `caller`,
and answers cross-caller questions in SQL — per-caller recall/precision and
"which variants does caller X miss that others catch". Like `stratify_db`, this
is the query layer (DuckDB); the benchmark core stays SQL-free.
"""

from __future__ import annotations

import dataclasses

import duckdb

from benchmark.metrics import precision, recall
from benchmark.pipeline import classified_rows


def benchmark_callers(truth: str, ref: str, callers: dict[str, str], **kwargs):
    """Benchmark each `{caller: query_vcf}` against `truth`; return one combined
    Arrow frame of classified rows with a `caller` column.

    `kwargs` pass through to `classified_rows` (regions, engine, on_ref_mismatch,
    conf_containment, decompose_mnp, strict).
    """
    import pyarrow as pa

    records: list[dict] = []
    for name, query_vcf in callers.items():
        for r in classified_rows(truth, query_vcf, ref, **kwargs):
            d = dataclasses.asdict(r)
            d["caller"] = name
            records.append(d)
    return pa.Table.from_pylist(records)


def per_caller_metrics(frame) -> list[dict]:
    """Recall/precision per (caller, vtype) over the ALL filter view."""
    con = duckdb.connect()
    try:
        con.register("c", frame)
        rows = con.execute(
            """
            SELECT caller, vtype,
                count(*) FILTER (WHERE side='truth' AND bd='TP') AS truth_tp,
                count(*) FILTER (WHERE side='truth' AND bd='FN') AS truth_fn,
                count(*) FILTER (WHERE side='query' AND bd='TP') AS query_tp,
                count(*) FILTER (WHERE side='query' AND bd='FP') AS query_fp
            FROM c
            WHERE filter_view='ALL'
            GROUP BY caller, vtype
            ORDER BY caller, vtype
            """
        ).fetchall()
    finally:
        con.close()
    out: list[dict] = []
    for caller, vtype, ttp, tfn, qtp, qfp in rows:
        ttp, tfn, qtp, qfp = int(ttp), int(tfn), int(qtp), int(qfp)
        out.append(
            {
                "caller": caller,
                "vtype": vtype,
                "truth_tp": ttp,
                "truth_fn": tfn,
                "query_tp": qtp,
                "query_fp": qfp,
                "recall": recall(ttp, tfn),
                "precision": precision(qtp, qfp),
            }
        )
    return out


_HISTORY_COLS = ["label", "caller", "vtype", "recall", "precision"]


def append_run(history_path: str, label: str, per_caller: list[dict]) -> None:
    """Append a run's per-caller metrics to a history CSV, tagged with `label`
    (e.g. a pipeline version) so concordance is comparable over time."""
    import csv
    import os

    exists = os.path.exists(history_path)
    with open(history_path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_HISTORY_COLS, extrasaction="ignore")
        if not exists:
            w.writeheader()
        for row in per_caller:
            w.writerow({"label": label, **row})


def load_history(history_path: str) -> list[dict]:
    """Read the regression-history CSV back as rows."""
    import csv

    with open(history_path, newline="") as fh:
        return list(csv.DictReader(fh))


def variants_missed_by(frame, caller: str, min_others: int = 1) -> list[tuple]:
    """Truth variants `caller` calls FN but at least `min_others` other callers
    recover (TP) — the caller's relative blind spots."""
    con = duckdb.connect()
    try:
        con.register("c", frame)
        rows = con.execute(
            """
            WITH t AS (
                SELECT caller, chrom, pos, ref, alt, bd
                FROM c WHERE side='truth' AND filter_view='ALL'
            )
            SELECT chrom, pos, ref, alt
            FROM t
            GROUP BY chrom, pos, ref, alt
            HAVING bool_or(caller = ? AND bd = 'FN')
               AND count(*) FILTER (WHERE bd = 'TP' AND caller != ?) >= ?
            ORDER BY chrom, pos
            """,
            [caller, caller, min_others],
        ).fetchall()
    finally:
        con.close()
    return [(c, int(p), r, a) for c, p, r, a in rows]
