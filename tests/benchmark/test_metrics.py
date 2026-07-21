from __future__ import annotations

from benchmark.metrics import (
    Metrics,
    f1,
    frac_na,
    het_hom,
    metrics_from_counts,
    precision,
    recall,
    titv,
)


# ---- recall = tp / (tp + fn) ----
def test_recall_known():
    assert recall(8, 2) == 0.8


def test_recall_perfect():
    assert recall(10, 0) == 1.0


def test_recall_zero_denominator():
    assert recall(0, 0) == 0.0


# ---- precision = tp / (tp + fp) ----
def test_precision_known():
    assert precision(3, 1) == 0.75


def test_precision_zero_denominator():
    assert precision(0, 0) == 0.0


# ---- f1 = 2PR / (P + R) ----
def test_f1_perfect():
    assert f1(1.0, 1.0) == 1.0


def test_f1_known():
    # P=0.75, R=0.8 -> 2*0.6/1.55
    assert f1(0.75, 0.8) == 2 * 0.75 * 0.8 / (0.75 + 0.8)


def test_f1_zero_denominator():
    assert f1(0.0, 0.0) == 0.0


# ---- frac_na = unk / total ----
def test_frac_na_known():
    assert frac_na(2, 8) == 0.25


def test_frac_na_zero_denominator():
    assert frac_na(0, 0) == 0.0


# ---- titv = ti / tv ----
def test_titv_known():
    assert titv(6, 3) == 2.0


def test_titv_zero_denominator():
    assert titv(5, 0) == 0.0


# ---- het_hom = het / hom ----
def test_het_hom_known():
    assert het_hom(9, 3) == 3.0


def test_het_hom_zero_denominator():
    assert het_hom(4, 0) == 0.0


# ---- metrics_from_counts end-to-end ----
def test_metrics_from_counts():
    m = metrics_from_counts(truth_tp=8, truth_fn=2, query_tp=9, query_fp=1, query_unk=2)
    assert isinstance(m, Metrics)
    assert m.recall == 0.8
    assert m.precision == 0.9
    assert m.f1 == 2 * 0.9 * 0.8 / (0.9 + 0.8)
    assert m.frac_na == 2 / 12


def test_metrics_from_counts_all_zero():
    m = metrics_from_counts(0, 0, 0, 0, 0)
    assert m.recall == 0.0
    assert m.precision == 0.0
    assert m.f1 == 0.0
    assert m.frac_na == 0.0
