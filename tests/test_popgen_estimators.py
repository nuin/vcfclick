"""Unit tests for popgen.estimators against hand-computed values.

Every expected number below was worked out by hand with exact fractions
(the fraction is given next to the decimal), or comes from an
independent derivation inside the test (Fu 1995 covariances), or from an
independent implementation run once outside the project (scikit-allel,
dadi — see test_popgen_cli.py for those literals).
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

import numpy as np
import pytest

from popgen import estimators as est

# ─────────────────────── harmonic numbers / constants ───────────────────


def test_harmonic_numbers():
    assert est.a1(4) == pytest.approx(11 / 6)  # 1 + 1/2 + 1/3
    assert est.a2(4) == pytest.approx(49 / 36)  # 1 + 1/4 + 1/9
    assert est.a1(2) == 1.0
    assert est.a1(1) == 0.0


def test_tajima_constants_n10():
    """Tajima (1989) constants for n = 10, exact fractions by hand."""
    c = est.tajima_constants(10)
    assert c.a1 == pytest.approx(7129 / 2520)
    assert c.a2 == pytest.approx(9778141 / 6350400)
    assert c.b1 == pytest.approx(11 / 27)
    assert c.b2 == pytest.approx(113 / 405)
    assert c.c1 == pytest.approx(10379 / 192483)
    assert c.c2 == pytest.approx(972076658 / 20583169605)
    assert c.e1 == pytest.approx(2906120 / 152467923)
    assert c.e2 == pytest.approx(7621080998720 / 1539945893952631)


# Fu (1995) exact covariances of the unfolded SFS under the standard
# neutral model: E[ξ_i] = θ/i, Var(ξ_i) = θ/i + σ_ii θ², Cov = σ_ij θ².
# The variance of any linear estimator Σ w_i ξ_i follows, so Tajima's
# c1/c2 and Zeng's Var(θ_π − θ_L) can be re-derived independently.


def _fu_a(i: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, i)), Fraction(0))


def _fu_beta(n: int, i: int) -> Fraction:
    return Fraction(2 * n, (n - i + 1) * (n - i)) * (
        _fu_a(n + 1) - _fu_a(i)
    ) - Fraction(2, n - i)


def _fu_sigma(n: int, i: int, j: int) -> Fraction:
    if i == j:
        if 2 * i < n:
            return _fu_beta(n, i + 1)
        if 2 * i == n:
            return 2 * (_fu_a(n) - _fu_a(i)) / (n - i) - Fraction(1, i * i)
        return _fu_beta(n, i) - Fraction(1, i * i)
    if i < j:
        i, j = j, i
    if i + j < n:
        return (_fu_beta(n, i + 1) - _fu_beta(n, i)) / 2
    if i + j == n:
        return (
            (_fu_a(n) - _fu_a(i)) / (n - i)
            + (_fu_a(n) - _fu_a(j)) / (n - j)
            - (_fu_beta(n, i) + _fu_beta(n, j + 1)) / 2
            - Fraction(1, i * j)
        )
    return (_fu_beta(n, j) - _fu_beta(n, j + 1)) / 2 - Fraction(1, i * j)


def _fu_variance(n: int, w: dict[int, Fraction]) -> tuple[Fraction, Fraction]:
    """(coefficient of θ, coefficient of θ²) of Var(Σ w_i ξ_i)."""
    lin = sum((w[i] ** 2 * Fraction(1, i) for i in range(1, n)), Fraction(0))
    quad = sum(
        (w[i] * w[j] * _fu_sigma(n, i, j) for i in range(1, n) for j in range(1, n)),
        Fraction(0),
    )
    return lin, quad


@pytest.mark.parametrize("n", [4, 5, 8, 13, 20])
def test_tajima_variance_matches_fu_1995(n):
    """Var(π − S/a1) = c1·θ + c2·θ² (Tajima 1989) re-derived from Fu's
    covariances."""
    a1 = _fu_a(n)
    w = {i: Fraction(2 * i * (n - i), n * (n - 1)) - 1 / a1 for i in range(1, n)}
    lin, quad = _fu_variance(n, w)
    c = est.tajima_constants(n)
    assert float(lin) == pytest.approx(c.c1, rel=1e-12)
    assert float(quad) == pytest.approx(c.c2, rel=1e-12)


@pytest.mark.parametrize("n", [4, 5, 8, 13, 20])
def test_zeng_h_variance_matches_fu_1995(n):
    """Zeng et al. (2006) eq. 11 Var(θ_π − θ_L) re-derived from Fu's
    covariances."""
    w = {
        i: Fraction(2 * i * (n - i), n * (n - 1)) - Fraction(i, n - 1)
        for i in range(1, n)
    }
    lin, quad = _fu_variance(n, w)
    bn1 = sum(Fraction(1, i * i) for i in range(1, n + 1))
    assert lin == Fraction(n - 2, 6 * (n - 1))
    expected_quad = (
        18 * n * n * (3 * n + 2) * bn1 - (88 * n**3 + 9 * n * n - 13 * n + 6)
    ) / Fraction(9 * n * (n - 1) ** 2)
    assert quad == expected_quad


# ───────────────────────────── per-site θ ───────────────────────────────


def test_pi_uses_each_sites_own_sample_size():
    k = np.array([1, 2, 0, 3])
    n = np.array([4, 6, 6, 3])
    # 2·1·3/(4·3) = 1/2 ; 2·2·4/(6·5) = 8/15 ; monomorphic ; fixed
    assert est.pi_per_site(k, n) == pytest.approx([0.5, 8 / 15, 0, 0])
    assert est.pi_per_site(k, n).sum() == pytest.approx(31 / 30)


def test_pi_zero_below_two_haplotypes():
    assert est.pi_per_site(np.array([1, 0]), np.array([1, 0])).tolist() == [0, 0]


def test_watterson_uses_per_site_a_n():
    k = np.array([1, 2, 0, 3])
    n = np.array([4, 6, 6, 3])
    # segregating: (1,4) and (2,6) → 1/a_4 + 1/a_6 = 6/11 + 60/137
    assert est.watterson_per_site(k, n).sum() == pytest.approx(1482 / 1507)


def test_watterson_equals_s_over_a_n_without_missing_data():
    k = np.array([1, 3, 0, 5, 2])
    n = np.full(5, 10)
    assert est.watterson_per_site(k, n).sum() == pytest.approx(4 / est.a1(10))


# ───────────────────────────── Tajima's D ───────────────────────────────


def test_tajima_d_worked_example():
    # n = 10, S = 5, π = 2: (2 − 5/a1) / sqrt(e1·5 + e2·20), by hand.
    assert est.tajima_d(2.0, 5, 10) == pytest.approx(0.527643382290888, rel=1e-12)


def test_tajima_d_undefined_cases():
    assert est.tajima_d(0.0, 0, 10) is None  # no segregating sites
    assert est.tajima_d(1.0, 1, 3) is None  # variance vanishes for n < 4
    assert est.tajima_d(1.0, 1, 2) is None


def test_tajima_d_zero_when_pi_equals_watterson():
    n, s = 12, 7
    assert est.tajima_d(s / est.a1(n), s, n) == pytest.approx(0.0, abs=1e-15)


# ─────────────────────────── Fay & Wu's H ───────────────────────────────


def test_normalised_h_worked_example():
    # n = 10, S = 5, θ_π = 2, θ_L = 6/5: Var = 0.402890864533457 by hand.
    assert est.fay_wu_h_normalised(2.0, 1.2, 5, 10) == pytest.approx(
        1.26036483337028, rel=1e-12
    )


def test_normalised_h_undefined():
    assert est.fay_wu_h_normalised(0.0, 0.0, 0, 10) is None
    assert est.fay_wu_h_normalised(1.0, 1.0, 1, 2) is None


def test_thetas_from_unfolded_sfs():
    # n = 4: ξ = [3, 2, 1, 1, 0] (three monomorphic-ancestral sites).
    t = est.thetas_from_sfs(np.array([3, 2, 1, 1, 0]))
    assert t.n == 4
    assert t.s == 4
    # θ_π = (2·1·3·2 + 2·2·2·1 + 2·3·1·1)/12 = (12 + 8 + 6)/12
    assert t.pi == pytest.approx(26 / 12)
    assert t.theta_w == pytest.approx(4 / (11 / 6))
    # θ_L = (1·2 + 2·1 + 3·1)/3 ; θ_H = 2(1·2 + 4·1 + 9·1)/12
    assert t.theta_l == pytest.approx(7 / 3)
    assert t.theta_h == pytest.approx(30 / 12)
    # identity θ_H = 2θ_L − θ_π, so raw H = θ_π − θ_H = 2(θ_π − θ_L)
    assert t.theta_h == pytest.approx(2 * t.theta_l - t.pi)


# ───────────────────────────── SFS tools ────────────────────────────────


def test_projection_by_hand():
    # (k, n) = (1, 4) and (3, 5) projected to m = 3:
    # (1,4): [1/4, 3/4, 0, 0]; (3,5): [0, 3/10, 6/10, 1/10]
    xi, used = est.project_sfs(np.array([1, 3]), np.array([4, 5]), 3)
    assert used == 2
    assert xi == pytest.approx([1 / 4, 21 / 20, 3 / 5, 1 / 10])


def test_projection_drops_sites_below_m():
    xi, used = est.project_sfs(np.array([1, 3, 1]), np.array([4, 5, 2]), 3)
    assert used == 2
    assert xi == pytest.approx([1 / 4, 21 / 20, 3 / 5, 1 / 10])


def test_projection_to_own_size_is_identity():
    k = np.array([0, 1, 1, 4, 7, 10])
    xi, used = est.project_sfs(k, np.full(6, 10), 10)
    assert used == 6
    assert xi.tolist() == est.sfs_counts(k, 10).tolist()


def test_projection_matches_exact_hypergeometric():
    rng = np.random.default_rng(7)
    n = rng.integers(30, 61, size=40)
    k = np.array([rng.integers(0, x + 1) for x in n])
    m = 24
    xi, _ = est.project_sfs(k, n, m)
    exact = [
        sum(
            Fraction(comb(int(kk), j) * comb(int(nn - kk), m - j), comb(int(nn), m))
            for kk, nn in zip(k, n, strict=True)
        )
        for j in range(m + 1)
    ]
    assert xi == pytest.approx([float(x) for x in exact], rel=1e-10, abs=1e-12)


def test_projection_preserves_pi():
    """π is the chance two distinct haplotypes differ, so it is unchanged
    by subsampling: Σ per-site π equals π of the projected SFS."""
    k = np.array([1, 5, 9, 2, 0, 12])
    n = np.array([14, 16, 12, 13, 15, 12])
    xi, _ = est.project_sfs(k, n, 12)
    assert est.thetas_from_sfs(xi).pi == pytest.approx(est.pi_per_site(k, n).sum())


def test_fold():
    assert est.fold(np.array([3, 2, 1, 1, 0])).tolist() == [3, 3, 1]  # n=4
    assert est.fold(np.array([1, 2, 3, 4, 5, 6])).tolist() == [7, 7, 7]  # n=5


def test_folded_thetas_equal_unfolded():
    xi = np.array([3.0, 2, 1, 1, 0])
    u = est.thetas_from_sfs(xi)
    f = est.thetas_from_folded(est.fold(xi), 4)
    assert (f.s, f.pi, f.theta_w) == pytest.approx((u.s, u.pi, u.theta_w))
    assert f.theta_l is None


# ───────────────────────────── Hudson F_ST ──────────────────────────────


def test_hudson_components_by_hand():
    num, den = est.hudson_fst_components(
        np.array([2, 1, 3]),
        np.array([4, 4, 6]),
        np.array([0, 0, 3]),
        np.array([4, 6, 6]),
    )
    # N = (p1−p2)² − p1q1/(n1−1) − p2q2/(n2−1); D = p1q2 + p2q1
    assert num == pytest.approx([1 / 6, 0, -1 / 10])
    assert den == pytest.approx([1 / 2, 1 / 4, 1 / 2])


def test_hudson_is_ratio_of_averages():
    k1, n1 = np.array([2, 1, 3]), np.array([4, 4, 6])
    k2, n2 = np.array([0, 0, 3]), np.array([4, 6, 6])
    # Σ N / Σ D = (1/6 + 0 − 1/10) / (1/2 + 1/4 + 1/2) = 4/75,
    # not the mean of per-site ratios (= 2/45).
    assert est.hudson_fst(k1, n1, k2, n2) == pytest.approx(4 / 75)


def test_hudson_equals_mean_pairwise_difference_form():
    """Bhatia's N equals (between − mean within) pairwise difference —
    the scikit-allel formulation."""
    rng = np.random.default_rng(3)
    n1, n2 = rng.integers(4, 30, 50), rng.integers(4, 30, 50)
    k1 = np.array([rng.integers(0, x + 1) for x in n1])
    k2 = np.array([rng.integers(0, x + 1) for x in n2])
    num, den = est.hudson_fst_components(k1, n1, k2, n2)
    p1, p2 = k1 / n1, k2 / n2
    between = p1 * (1 - p2) + p2 * (1 - p1)
    within = (est.pi_per_site(k1, n1) + est.pi_per_site(k2, n2)) / 2
    assert num == pytest.approx(between - within)
    assert den == pytest.approx(between)


def test_hudson_undefined_without_polymorphism():
    z = np.zeros(3, dtype=int)
    assert est.hudson_fst(z, np.full(3, 4), z, np.full(3, 4)) is None


def test_hudson_identical_samples_negative_or_zero():
    k, n = np.array([2, 3]), np.array([6, 6])
    assert est.hudson_fst(k, n, k, n) < 0  # sampling-bias correction


# ───────────────────────── heterozygosity / F ───────────────────────────


def test_heterozygosity_by_hand():
    h = est.heterozygosity(
        np.array([1, 2]), np.array([3, 4]), np.array([1, 4]), np.array([6, 8])
    )
    # Ho = (1/3 + 2/4)/2 = 5/12; He = (10/30 + 32/56)/2 = 19/42
    assert h.ho == pytest.approx(5 / 12)
    assert h.he == pytest.approx(19 / 42)
    assert h.f == pytest.approx(3 / 38)


def test_heterozygosity_no_variation():
    h = est.heterozygosity(np.zeros(2), np.full(2, 5), np.zeros(2), np.full(2, 10))
    assert (h.ho, h.he, h.f) == (0.0, 0.0, None)
