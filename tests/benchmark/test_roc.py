"""ROC / PR quality-threshold sweep over classified concordance."""

from __future__ import annotations

from benchmark.constants import BD_FN, BD_FP, BD_TP, FILTER_ALL, VT_SNP
from benchmark.model import ClassifiedRow
from benchmark.roc import roc_curve


def _q(bd, qual, side="query"):
    return ClassifiedRow(
        side=side,
        filter_view=FILTER_ALL,
        chrom="chr1",
        pos=1,
        ref="A",
        alt="G",
        vtype=VT_SNP,
        subtype="",
        in_conf=True,
        bd=bd,
        bk=".",
        blt="het",
        qual=qual,
    )


def test_roc_sweep_trades_recall_for_precision():
    # 3 query TP (qual 30/20/10) + 1 query FP (qual 15) + 1 truth-only FN.
    rows = [
        _q(BD_TP, 30),
        _q(BD_TP, 20),
        _q(BD_TP, 10),
        _q(BD_FP, 15),
        _q(BD_FN, 0, side="truth"),  # never called; FN at every threshold
    ]
    curve = {r["threshold"]: r for r in roc_curve(rows, VT_SNP)}

    # At threshold 0 (keep all): 3 TP, 1 FP, 1 fixed FN -> recall 3/4, precision 3/4
    lo = curve[min(curve)]
    assert lo["tp"] == 3 and lo["fp"] == 1
    assert lo["recall"] == 0.75 and lo["precision"] == 0.75

    # At threshold 20 (drop the qual-10 TP and qual-15 FP): 2 TP, 0 FP
    t20 = curve[20.0]
    assert t20["tp"] == 2 and t20["fp"] == 0
    assert t20["precision"] == 1.0  # FP filtered out
    assert t20["recall"] == 0.5  # 2 of 4 truth recovered


def test_roc_empty_when_no_type_rows():
    assert roc_curve([], VT_SNP) == []
