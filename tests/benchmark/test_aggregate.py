from __future__ import annotations

from benchmark.aggregate import EXTRA_KEYS, aggregate_counts, aggregate_counts_oracle
from benchmark.constants import (
    BD_FN,
    BD_FP,
    BD_N,
    BD_TP,
    BK_NONE,
    FILTER_ALL,
    FILTER_PASS,
    VT_INDEL,
    VT_SNP,
)
from benchmark.model import ClassifiedRow


def _pad(d: dict) -> dict:
    """Fill the extra hap.py-parity tallies with 0 (the base scenarios don't set
    subtype/blt/bk, so ti/tv/het/hom/fp_gt are all zero)."""
    return {gk: {**counts, **{k: 0 for k in EXTRA_KEYS}} for gk, counts in d.items()}


def _row(side, filter_view, vtype, bd):
    """Minimal ClassifiedRow; only side/filter_view/vtype/bd drive counts."""
    return ClassifiedRow(
        side=side,
        filter_view=filter_view,
        chrom="chr1",
        pos=1,
        ref="A",
        alt="C",
        vtype=vtype,
        subtype="",
        in_conf=True,
        bd=bd,
        bk=BK_NONE,
        blt="",
    )


# Hand-built scenario with three (filter_view, vtype) groups.
ROWS = [
    _row("truth", FILTER_ALL, VT_SNP, BD_TP),
    _row("truth", FILTER_ALL, VT_SNP, BD_FN),
    _row("query", FILTER_ALL, VT_SNP, BD_TP),
    _row("query", FILTER_ALL, VT_SNP, BD_TP),
    _row("query", FILTER_ALL, VT_SNP, BD_FP),
    _row("query", FILTER_ALL, VT_SNP, BD_N),
    _row("truth", FILTER_ALL, VT_SNP, BD_N),  # truth-N: counted nowhere
    _row("truth", FILTER_ALL, VT_INDEL, BD_TP),
    _row("query", FILTER_PASS, VT_SNP, BD_FP),
]

EXPECTED = {
    (FILTER_ALL, VT_SNP): {
        "truth_tp": 1,
        "truth_fn": 1,
        "query_tp": 2,
        "query_fp": 1,
        "query_unk": 1,
    },
    (FILTER_ALL, VT_INDEL): {
        "truth_tp": 1,
        "truth_fn": 0,
        "query_tp": 0,
        "query_fp": 0,
        "query_unk": 0,
    },
    (FILTER_PASS, VT_SNP): {
        "truth_tp": 0,
        "truth_fn": 0,
        "query_tp": 0,
        "query_fp": 1,
        "query_unk": 0,
    },
}


def test_aggregate_counts_exact():
    assert aggregate_counts(ROWS) == _pad(EXPECTED)


def test_oracle_exact():
    assert aggregate_counts_oracle(ROWS) == _pad(EXPECTED)


def test_arrow_and_oracle_agree():
    assert aggregate_counts(ROWS) == aggregate_counts_oracle(ROWS)


def test_empty_arrow():
    assert aggregate_counts([]) == {}


def test_empty_oracle():
    assert aggregate_counts_oracle([]) == {}


def test_single_truth_fn():
    rows = [_row("truth", FILTER_ALL, VT_SNP, BD_FN)]
    expected = {
        (FILTER_ALL, VT_SNP): {
            "truth_tp": 0,
            "truth_fn": 1,
            "query_tp": 0,
            "query_fp": 0,
            "query_unk": 0,
        }
    }
    assert aggregate_counts(rows) == _pad(expected)
    assert aggregate_counts_oracle(rows) == _pad(expected)


def test_extra_tallies_titv_hethom_fpgt():
    from benchmark.constants import BK_AM

    def cr(side, bd, subtype, blt, bk=BK_NONE):
        return ClassifiedRow(
            side=side,
            filter_view=FILTER_ALL,
            chrom="chr1",
            pos=1,
            ref="A",
            alt="C",
            vtype=VT_SNP,
            subtype=subtype,
            in_conf=True,
            bd=bd,
            bk=bk,
            blt=blt,
        )

    rows = [
        cr("truth", BD_TP, "ti", "het"),
        cr("truth", BD_FN, "tv", "homalt"),
        cr("query", BD_TP, "ti", "het"),
        cr("query", BD_FP, "tv", "het", bk=BK_AM),  # genotype-error FP
        cr("query", BD_FP, "ti", "homalt"),  # plain allele FP
    ]
    agg = aggregate_counts(rows)[(FILTER_ALL, VT_SNP)]
    assert agg["truth_ti"] == 1 and agg["truth_tv"] == 1
    assert agg["query_ti"] == 2 and agg["query_tv"] == 1
    assert agg["truth_het"] == 1 and agg["truth_hom"] == 1
    assert agg["query_het"] == 2 and agg["query_hom"] == 1
    assert agg["query_fp"] == 2 and agg["query_fp_gt"] == 1  # FP.gt=1, FP.al=1
    assert aggregate_counts(rows) == aggregate_counts_oracle(rows)
