"""Metric math for benchmarking: pure, integer-count inputs, divide-guarded."""

from __future__ import annotations

from dataclasses import dataclass


def recall(tp: int, fn: int) -> float:
    """Recall = TP / (TP + FN); 0.0 when no truth positives."""
    denom = tp + fn
    return tp / denom if denom else 0.0


def precision(tp: int, fp: int) -> float:
    """Precision = TP / (TP + FP); 0.0 when no query positives."""
    denom = tp + fp
    return tp / denom if denom else 0.0


def f1(precision: float, recall: float) -> float:
    """F1 = 2PR / (P + R); 0.0 when P + R == 0."""
    denom = precision + recall
    return 2 * precision * recall / denom if denom else 0.0


def frac_na(unk: int, total: int) -> float:
    """Fraction not-assessed = UNK / total; 0.0 when total == 0."""
    return unk / total if total else 0.0


def titv(ti: int, tv: int) -> float:
    """Transition/transversion ratio = Ti / Tv; 0.0 when Tv == 0."""
    return ti / tv if tv else 0.0


def het_hom(het: int, hom: int) -> float:
    """Het/hom ratio = het / hom; 0.0 when hom == 0."""
    return het / hom if hom else 0.0


@dataclass(frozen=True)
class Metrics:
    recall: float
    precision: float
    f1: float
    frac_na: float


def metrics_from_counts(
    truth_tp: int,
    truth_fn: int,
    query_tp: int,
    query_fp: int,
    query_unk: int,
) -> Metrics:
    """Assemble Metrics from integer counts (truth for recall, query for the rest)."""
    r = recall(truth_tp, truth_fn)
    p = precision(query_tp, query_fp)
    return Metrics(
        recall=r,
        precision=p,
        f1=f1(p, r),
        frac_na=frac_na(query_unk, query_tp + query_fp + query_unk),
    )
