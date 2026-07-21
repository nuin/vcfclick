"""Count classified rows two ways that must agree: a pyarrow group-by (default)
and a pure-Python `Counter` differential oracle. Both key on (filter_view, vtype)
and split by side/bd into the five counts metrics.py consumes."""

from __future__ import annotations

from collections import Counter

import pyarrow as pa
import pyarrow.compute as pc

from benchmark.constants import BD_FN, BD_FP, BD_N, BD_TP
from benchmark.model import ClassifiedRow

COUNT_KEYS = ("truth_tp", "truth_fn", "query_tp", "query_fp", "query_unk")


def _count_key(row: ClassifiedRow) -> str | None:
    """Which count a row falls into, or None if it counts nowhere (e.g. truth-N)."""
    if row.side == "truth":
        if row.bd == BD_TP:
            return "truth_tp"
        if row.bd == BD_FN:
            return "truth_fn"
    elif row.side == "query":
        if row.bd == BD_TP:
            return "query_tp"
        if row.bd == BD_FP:
            return "query_fp"
        if row.bd == BD_N:
            return "query_unk"
    return None


def aggregate_counts(rows: list[ClassifiedRow]) -> dict:
    """pyarrow group-by over (filter_view, vtype); returns {key: {count: int}}."""
    if not rows:
        return {}

    side = pa.array([r.side for r in rows], pa.string())
    bd = pa.array([r.bd for r in rows], pa.string())
    is_truth = pc.equal(side, "truth")
    is_query = pc.equal(side, "query")

    columns = {
        "filter_view": pa.array([r.filter_view for r in rows], pa.string()),
        "vtype": pa.array([r.vtype for r in rows], pa.string()),
        "truth_tp": pc.and_(is_truth, pc.equal(bd, BD_TP)),
        "truth_fn": pc.and_(is_truth, pc.equal(bd, BD_FN)),
        "query_tp": pc.and_(is_query, pc.equal(bd, BD_TP)),
        "query_fp": pc.and_(is_query, pc.equal(bd, BD_FP)),
        "query_unk": pc.and_(is_query, pc.equal(bd, BD_N)),
    }
    table = pa.table(
        {
            k: (pc.cast(v, pa.int64()) if k in COUNT_KEYS else v)
            for k, v in columns.items()
        }
    )
    grouped = table.group_by(["filter_view", "vtype"]).aggregate(
        [(k, "sum") for k in COUNT_KEYS]
    )

    fv = grouped.column("filter_view").to_pylist()
    vt = grouped.column("vtype").to_pylist()
    sums = {k: grouped.column(f"{k}_sum").to_pylist() for k in COUNT_KEYS}
    return {
        (fv[i], vt[i]): {k: sums[k][i] for k in COUNT_KEYS}
        for i in range(grouped.num_rows)
    }


def aggregate_counts_oracle(rows: list[ClassifiedRow]) -> dict:
    """Pure-Python `Counter` oracle; must equal aggregate_counts row-for-row."""
    counter: Counter = Counter()
    groups: set = set()
    for r in rows:
        gk = (r.filter_view, r.vtype)
        groups.add(gk)
        name = _count_key(r)
        if name is not None:
            counter[(gk, name)] += 1
    return {gk: {k: counter[(gk, k)] for k in COUNT_KEYS} for gk in groups}
