"""Count classified rows two ways that must agree: a pyarrow group-by (default)
and a pure-Python `Counter` differential oracle. Both key on (filter_view, vtype)
and split by side/bd into the five counts metrics.py consumes."""

from __future__ import annotations

from collections import Counter

import pyarrow as pa
import pyarrow.compute as pc

from benchmark.constants import (
    BD_FN,
    BD_FP,
    BD_N,
    BD_TP,
    BK_AM,
    BLT_HETALT,
    BLT_HET,
    BLT_HOMALT,
)
from benchmark.model import ClassifiedRow

COUNT_KEYS = ("truth_tp", "truth_fn", "query_tp", "query_fp", "query_unk")
# Extra hap.py-parity tallies: FP genotype-vs-allele split, and Ti/Tv + het/hom
# over each side's total (for the ratio columns).
EXTRA_KEYS = (
    "query_fp_gt",
    "truth_ti",
    "truth_tv",
    "query_ti",
    "query_tv",
    "truth_het",
    "truth_hom",
    "query_het",
    "query_hom",
)
ALL_KEYS = COUNT_KEYS + EXTRA_KEYS
_HET = {BLT_HET, BLT_HETALT}


def _row_keys(row: ClassifiedRow) -> list[str]:
    """All tally keys a row contributes to (a row feeds several: e.g. a truth TP
    transition heterozygote → truth_tp, truth_ti, truth_het)."""
    keys: list[str] = []
    truth = row.side == "truth"
    query = row.side == "query"
    if truth and row.bd == BD_TP:
        keys.append("truth_tp")
    elif truth and row.bd == BD_FN:
        keys.append("truth_fn")
    elif query and row.bd == BD_TP:
        keys.append("query_tp")
    elif query and row.bd == BD_FP:
        keys.append("query_fp")
        if row.bk == BK_AM:
            keys.append("query_fp_gt")
    elif query and row.bd == BD_N:
        keys.append("query_unk")
    # Ti/Tv and het/hom are over each side's total records.
    pfx = "truth" if truth else "query"
    if row.subtype == "ti":
        keys.append(f"{pfx}_ti")
    elif row.subtype == "tv":
        keys.append(f"{pfx}_tv")
    if row.blt in _HET:
        keys.append(f"{pfx}_het")
    elif row.blt == BLT_HOMALT:
        keys.append(f"{pfx}_hom")
    return keys


def aggregate_counts(rows: list[ClassifiedRow]) -> dict:
    """pyarrow group-by over (filter_view, vtype); returns {key: {tally: int}}."""
    if not rows:
        return {}

    side = pa.array([r.side for r in rows], pa.string())
    bd = pa.array([r.bd for r in rows], pa.string())
    bk = pa.array([r.bk for r in rows], pa.string())
    sub = pa.array([r.subtype for r in rows], pa.string())
    blt = pa.array([r.blt for r in rows], pa.string())
    is_truth = pc.equal(side, "truth")
    is_query = pc.equal(side, "query")
    is_het = pc.is_in(blt, value_set=pa.array(sorted(_HET)))

    bools = {
        "truth_tp": pc.and_(is_truth, pc.equal(bd, BD_TP)),
        "truth_fn": pc.and_(is_truth, pc.equal(bd, BD_FN)),
        "query_tp": pc.and_(is_query, pc.equal(bd, BD_TP)),
        "query_fp": pc.and_(is_query, pc.equal(bd, BD_FP)),
        "query_unk": pc.and_(is_query, pc.equal(bd, BD_N)),
        "query_fp_gt": pc.and_(
            pc.and_(is_query, pc.equal(bd, BD_FP)), pc.equal(bk, BK_AM)
        ),
        "truth_ti": pc.and_(is_truth, pc.equal(sub, "ti")),
        "truth_tv": pc.and_(is_truth, pc.equal(sub, "tv")),
        "query_ti": pc.and_(is_query, pc.equal(sub, "ti")),
        "query_tv": pc.and_(is_query, pc.equal(sub, "tv")),
        "truth_het": pc.and_(is_truth, is_het),
        "truth_hom": pc.and_(is_truth, pc.equal(blt, BLT_HOMALT)),
        "query_het": pc.and_(is_query, is_het),
        "query_hom": pc.and_(is_query, pc.equal(blt, BLT_HOMALT)),
    }
    table = pa.table(
        {
            "filter_view": pa.array([r.filter_view for r in rows], pa.string()),
            "vtype": pa.array([r.vtype for r in rows], pa.string()),
            **{k: pc.cast(v, pa.int64()) for k, v in bools.items()},
        }
    )
    grouped = table.group_by(["filter_view", "vtype"]).aggregate(
        [(k, "sum") for k in ALL_KEYS]
    )

    fv = grouped.column("filter_view").to_pylist()
    vt = grouped.column("vtype").to_pylist()
    sums = {k: grouped.column(f"{k}_sum").to_pylist() for k in ALL_KEYS}
    return {
        (fv[i], vt[i]): {k: sums[k][i] for k in ALL_KEYS}
        for i in range(grouped.num_rows)
    }


def aggregate_counts_oracle(rows: list[ClassifiedRow]) -> dict:
    """Pure-Python `Counter` oracle; must equal aggregate_counts row-for-row."""
    counter: Counter = Counter()
    groups: set = set()
    for r in rows:
        gk = (r.filter_view, r.vtype)
        groups.add(gk)
        for name in _row_keys(r):
            counter[(gk, name)] += 1
    return {gk: {k: counter[(gk, k)] for k in ALL_KEYS} for gk in groups}
