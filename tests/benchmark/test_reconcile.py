from __future__ import annotations

import pytest

from benchmark.constants import (
    BD_FN,
    BD_FP,
    BD_N,
    BD_TP,
    BK_AM,
    BK_GM,
    BK_NONE,
    FILTER_ALL,
    FILTER_PASS,
)
from benchmark.model import NormRecord
from benchmark.reconcile import classify, classify_haplotype


def nr(
    side,
    pos,
    ref,
    alt,
    gt,
    *,
    is_pass=True,
    other=False,
    in_conf=True,
    locus=None,
    vtype="SNP",
):
    return NormRecord(
        side=side,
        chrom="chr1",
        pos=pos,
        ref=ref,
        alt=alt,
        gt=gt,
        is_pass=is_pass,
        other_alt_present=other,
        locus_id=locus or (pos, ref, (alt,)),
        in_conf=in_conf,
        vtype=vtype,
        subtype="",
        blt="",
    )


def _by_side(rows):
    return {r.side: r for r in rows}


def test_tp_both_present_gt_equal():
    t = nr("truth", 10, "A", "G", "0/1")
    q = nr("query", 10, "A", "G", "0/1")
    rows = classify([t], [q], FILTER_ALL)
    assert len(rows) == 2
    s = _by_side(rows)
    assert s["truth"].bd == BD_TP and s["truth"].bk == BK_GM
    assert s["query"].bd == BD_TP and s["query"].bk == BK_GM


def test_fn_truth_only_in_conf():
    t = nr("truth", 10, "A", "G", "0/1", in_conf=True)
    rows = classify([t], [], FILTER_ALL)
    assert len(rows) == 1
    assert rows[0].side == "truth"
    assert rows[0].bd == BD_FN and rows[0].bk == BK_NONE


def test_fp_query_only_in_conf():
    q = nr("query", 10, "A", "G", "0/1", in_conf=True)
    rows = classify([], [q], FILTER_ALL)
    assert len(rows) == 1
    assert rows[0].side == "query"
    assert rows[0].bd == BD_FP and rows[0].bk == BK_NONE


def test_gt_mismatch_double_penalty():
    t = nr("truth", 10, "A", "G", "0/1")
    q = nr("query", 10, "A", "G", "1/1")
    rows = classify([t], [q], FILTER_ALL)
    assert len(rows) == 2
    s = _by_side(rows)
    assert s["truth"].bd == BD_FN and s["truth"].bk == BK_AM
    assert s["query"].bd == BD_FP and s["query"].bk == BK_AM
    # exactly +1 FN and +1 FP
    assert sum(1 for r in rows if r.bd == BD_FN) == 1
    assert sum(1 for r in rows if r.bd == BD_FP) == 1


def test_unk_gating_query_out_of_conf_is_N():
    q = nr("query", 10, "A", "G", "0/1", in_conf=False)
    rows = classify([], [q], FILTER_ALL)
    assert len(rows) == 1
    assert rows[0].bd == BD_N and rows[0].bk == BK_NONE


def test_truth_only_out_of_conf_is_dropped():
    t = nr("truth", 10, "A", "G", "0/1", in_conf=False)
    rows = classify([t], [], FILTER_ALL)
    assert rows == []


def test_pass_vs_all_filter_view():
    # A query call that key-matches truth but fails FILTER.
    t = nr("truth", 10, "A", "G", "0/1", in_conf=True)
    q = nr("query", 10, "A", "G", "0/1", is_pass=False, in_conf=True)

    all_rows = classify([t], [q], FILTER_ALL)
    s = _by_side(all_rows)
    assert s["truth"].bd == BD_TP and s["query"].bd == BD_TP

    pass_rows = classify([t], [q], FILTER_PASS)
    # query is absent in PASS view => truth becomes FN, no TP anywhere
    assert all(r.bd != BD_TP for r in pass_rows)
    assert len(pass_rows) == 1
    assert pass_rows[0].side == "truth" and pass_rows[0].bd == BD_FN


def test_hetalt_not_matched_to_het():
    # truth plain 0/1; query het-alt 1/2 split row (same key, other alt present)
    t = nr("truth", 10, "A", "G", "0/1", other=False)
    q = nr("query", 10, "A", "G", "0/1", other=True, locus=(10, "A", ("G", "T")))
    rows = classify([t], [q], FILTER_ALL)
    s = _by_side(rows)
    assert s["truth"].bd == BD_FN and s["truth"].bk == BK_AM
    assert s["query"].bd == BD_FP and s["query"].bk == BK_AM


def test_duplicate_key_routes_to_N_bucket_no_merge():
    # two truth records share the normalized key => both routed to BD_N
    t1 = nr("truth", 10, "A", "G", "0/1")
    t2 = nr("truth", 10, "A", "G", "1/1")
    q = nr("query", 10, "A", "G", "0/1", in_conf=True)
    rows = classify([t1, t2], [q], FILTER_ALL)
    truth_rows = [r for r in rows if r.side == "truth"]
    assert len(truth_rows) == 2
    assert all(r.bd == BD_N and r.bk == BK_NONE for r in truth_rows)
    # no TP was produced (no silent merge/cartesian); query is a plain FP
    assert all(r.bd != BD_TP for r in rows)
    query_rows = [r for r in rows if r.side == "query"]
    assert len(query_rows) == 1 and query_rows[0].bd == BD_FP


def test_classify_haplotype_unsupported():
    with pytest.raises(NotImplementedError):
        classify_haplotype([], [], FILTER_ALL)


def test_non_pass_truth_is_still_scored_in_pass_view():
    # hap.py filters only the QUERY by FILTER; truth is the gold standard and is
    # always considered. A non-PASS truth call matched by a PASS query is a TP,
    # not a spurious query FP (and truth must stay in the recall denominator).
    t = nr("truth", 10, "A", "G", "0/1", is_pass=False)
    q = nr("query", 10, "A", "G", "0/1", is_pass=True)
    rows = classify([t], [q], FILTER_PASS)
    by_side = _by_side(rows)
    assert by_side["truth"].bd == BD_TP
    assert by_side["query"].bd == BD_TP
