"""ROC / precision-recall quality-threshold sweep over classified concordance.

Vary a minimum quality threshold on the *query*: a query call below threshold is
treated as filtered, so its truth match reverts to FN (recall drops) and a filtered
FP disappears (precision rises). Truth-only FNs are fixed at every threshold. Pure
Python over the classified rows (which now carry `qual`); no SQL, no engine.
"""

from __future__ import annotations

from benchmark.constants import BD_FN, BD_FP, BD_TP, FILTER_ALL
from benchmark.metrics import precision, recall


def roc_curve(rows: list, vtype: str) -> list[dict]:
    """Points of (threshold, tp, fp, recall, precision) sweeping the min query
    quality upward, for one variant type over the ALL filter view."""
    rows = [r for r in rows if r.vtype == vtype and r.filter_view == FILTER_ALL]
    if not rows:
        return []

    tp_quals = sorted(r.qual for r in rows if r.side == "query" and r.bd == BD_TP)
    fp_quals = sorted(r.qual for r in rows if r.side == "query" and r.bd == BD_FP)
    fixed_fn = sum(1 for r in rows if r.side == "truth" and r.bd == BD_FN)
    total_truth = fixed_fn + len(tp_quals)  # every query TP recovers one truth call

    thresholds = sorted({0.0, *tp_quals, *fp_quals})
    out: list[dict] = []
    for t in thresholds:
        tp = sum(1 for q in tp_quals if q >= t)
        fp = sum(1 for q in fp_quals if q >= t)
        out.append(
            {
                "threshold": t,
                "tp": tp,
                "fp": fp,
                "recall": recall(tp, total_truth - tp),
                "precision": precision(tp, fp),
            }
        )
    return out


def write_roc_tsv(rows: list, vtypes: list[str], path: str) -> None:
    """Write a combined ROC TSV (Type, Threshold, TP, FP, Recall, Precision)."""
    import csv

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Type", "Threshold", "TP", "FP", "Recall", "Precision"])
        for vt in vtypes:
            for p in roc_curve(rows, vt):
                w.writerow(
                    [
                        vt,
                        p["threshold"],
                        p["tp"],
                        p["fp"],
                        f"{p['recall']:.6f}",
                        f"{p['precision']:.6f}",
                    ]
                )
