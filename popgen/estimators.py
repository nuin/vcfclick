"""Population-genetics estimators over allele counts.

Conventions used throughout:

  * `n`  number of called haplotypes (alleles) at a site in a group;
  * `k`  number of ALT (or, once polarised, derived) alleles among them;
  * an unfolded SFS `xi` of sample size m is a length m+1 vector where
    xi[j] is the number of sites with j derived alleles (xi[0] and xi[m]
    are the monomorphic classes);
  * a folded SFS of size m has length m//2 + 1, indexed by minor-allele
    count.

All formulas are the textbook ones; references are given per function
and collected in docs/POPGEN.md.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# ─────────────────────────── harmonic numbers ───────────────────────────


def a1(n: int) -> float:
    """Watterson's a_n = sum_{i=1}^{n-1} 1/i (0 for n < 2)."""
    return float(sum(1.0 / i for i in range(1, n)))


def a2(n: int) -> float:
    """sum_{i=1}^{n-1} 1/i^2 (Tajima's a2, Zeng's b_n)."""
    return float(sum(1.0 / (i * i) for i in range(1, n)))


def _a1_table(n_max: int) -> np.ndarray:
    """table[n] = a_n for n in 0..n_max (table[0] = table[1] = 0)."""
    table = np.zeros(n_max + 1)
    if n_max >= 2:
        table[2:] = np.cumsum(1.0 / np.arange(1, n_max))
    return table


# ────────────────────────────── per-site θ ──────────────────────────────


def pi_per_site(k: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Nucleotide diversity at each site with that site's own sample size:
    2k(n-k) / (n(n-1)), the probability that two distinct haplotypes drawn
    without replacement differ (Nei & Li 1979; Tajima 1983). Equal to
    Nei's unbiased expected heterozygosity. 0 where n < 2."""
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    out = np.zeros_like(n)
    ok = n >= 2
    out[ok] = 2.0 * k[ok] * (n[ok] - k[ok]) / (n[ok] * (n[ok] - 1.0))
    return out


def segregating(k: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Boolean mask: the site is polymorphic in this sample (0 < k < n)."""
    k = np.asarray(k)
    n = np.asarray(n)
    return (k > 0) & (k < n) & (n >= 2)


def watterson_per_site(k: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Watterson's θ contribution of each site, 1/a_{n_i} for a
    segregating site and 0 otherwise, so the sum over sites is θ_W with
    each site's own sample size (Watterson 1975; the per-site-a_n form
    for missing data, e.g. Ferretti et al. 2012). With no missing data
    this is exactly S / a_n."""
    k = np.asarray(k)
    n = np.asarray(n).astype(int)
    out = np.zeros(n.shape, dtype=float)
    seg = segregating(k, n)
    if seg.any():
        table = _a1_table(int(n.max()))
        out[seg] = 1.0 / table[n[seg]]
    return out


# ───────────────────────────── Tajima's D ───────────────────────────────


@dataclass(frozen=True)
class TajimaConstants:
    """Tajima (1989) constants for sample size n."""

    n: int
    a1: float
    a2: float
    b1: float
    b2: float
    c1: float
    c2: float
    e1: float
    e2: float


def tajima_constants(n: int) -> TajimaConstants:
    """a1, a2, b1, b2, c1, c2, e1, e2 for sample size n (Tajima 1989,
    eqs. 3-4 and the e1/e2 substitution under eq. 38)."""
    if n < 2:
        raise ValueError(f"Tajima constants need n >= 2, got {n}")
    h1, h2 = a1(n), a2(n)
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / h1
    c2 = b2 - (n + 2) / (h1 * n) + h2 / (h1 * h1)
    return TajimaConstants(
        n=n,
        a1=h1,
        a2=h2,
        b1=b1,
        b2=b2,
        c1=c1,
        c2=c2,
        e1=c1 / h1,
        e2=c2 / (h1 * h1 + h2),
    )


def tajima_d(pi: float, s: float, n: int) -> float | None:
    """Tajima's D = (π − S/a1) / sqrt(e1·S + e2·S·(S−1)).

    `pi` is the summed nucleotide diversity, `s` the number of segregating
    sites (may be fractional when taken from a projected SFS), `n` the
    sample size in haplotypes. None when undefined (S = 0, or n < 4 where
    the variance vanishes)."""
    if n < 4 or s <= 0:
        return None
    c = tajima_constants(n)
    var = c.e1 * s + c.e2 * s * (s - 1)
    if var <= 0:
        return None
    return (pi - s / c.a1) / math.sqrt(var)


# ─────────────────────────── Fay & Wu's H ───────────────────────────────


def fay_wu_h_normalised(pi: float, theta_l: float, s: float, n: int) -> float | None:
    """Normalised Fay & Wu's H (Zeng, Fu, Shi & Wu 2006, eq. 11):

        H = (θ_π − θ_L) / sqrt(Var(θ_π − θ_L))

        Var = (n−2)/(6(n−1)) · θ
            + [18 n² (3n+2) b_{n+1} − (88n³ + 9n² − 13n + 6)] / [9n(n−1)²] · θ²

    with θ = S/a_n and θ² = S(S−1)/(a_n² + b_n), b_n = Σ_{i<n} 1/i².
    Requires a polarised (unfolded) SFS. None when undefined."""
    if n < 3 or s <= 0:
        return None
    an, bn = a1(n), a2(n)
    bn1 = a2(n + 1)
    theta = s / an
    theta_sq = s * (s - 1) / (an * an + bn)
    var = (n - 2) / (6 * (n - 1)) * theta + (
        (18 * n * n * (3 * n + 2) * bn1 - (88 * n**3 + 9 * n * n - 13 * n + 6))
        / (9 * n * (n - 1) ** 2)
    ) * theta_sq
    if var <= 0:
        return None
    return (pi - theta_l) / math.sqrt(var)


# ────────────────────────────── SFS tools ───────────────────────────────


@dataclass(frozen=True)
class SfsThetas:
    """θ estimators from an unfolded SFS of sample size n (Fay & Wu 2000;
    Zeng et al. 2006). For a folded spectrum only s, pi and theta_w are
    meaningful; theta_l / theta_h are then None."""

    n: int
    s: float
    pi: float
    theta_w: float
    theta_l: float | None
    theta_h: float | None


def thetas_from_sfs(xi: np.ndarray) -> SfsThetas:
    """S, θ_π, θ_W, θ_L and θ_H from an unfolded spectrum (length n+1):

        θ_π = Σ 2 j (n−j) ξ_j / (n(n−1))     θ_W = S / a_n
        θ_L = Σ j ξ_j / (n−1)                θ_H = Σ 2 j² ξ_j / (n(n−1))

    with S = Σ_{j=1}^{n−1} ξ_j and the θ_L/θ_H sums over j = 1..n−1.
    """
    xi = np.asarray(xi, dtype=float)
    n = len(xi) - 1
    j = np.arange(n + 1, dtype=float)
    inner = slice(1, n)
    s = float(xi[inner].sum())
    pi = float((xi * 2 * j * (n - j)).sum() / (n * (n - 1)))
    theta_l = float((xi[inner] * j[inner]).sum() / (n - 1))
    theta_h = float((xi[inner] * 2 * j[inner] ** 2).sum() / (n * (n - 1)))
    return SfsThetas(
        n=n, s=s, pi=pi, theta_w=s / a1(n), theta_l=theta_l, theta_h=theta_h
    )


def thetas_from_folded(eta: np.ndarray, n: int) -> SfsThetas:
    """S, θ_π, θ_W from a folded SFS of sample size n."""
    eta = np.asarray(eta, dtype=float)
    i = np.arange(len(eta), dtype=float)
    s = float(eta[1:].sum())
    pi = float((eta * 2 * i * (n - i)).sum() / (n * (n - 1)))
    return SfsThetas(n=n, s=s, pi=pi, theta_w=s / a1(n), theta_l=None, theta_h=None)


def fold(xi: np.ndarray) -> np.ndarray:
    """Fold an unfolded SFS (length n+1) onto minor-allele counts
    (length n//2 + 1): eta[i] = xi[i] + xi[n−i], with the middle class
    (i = n/2, n even) counted once."""
    xi = np.asarray(xi, dtype=float)
    n = len(xi) - 1
    eta = np.zeros(n // 2 + 1)
    for i in range(n // 2 + 1):
        eta[i] = xi[i] if i == n - i else xi[i] + xi[n - i]
    return eta


def _log_factorial_table(n_max: int) -> np.ndarray:
    """table[x] = log(x!) for x in 0..n_max."""
    table = np.zeros(n_max + 1)
    if n_max >= 1:
        table[1:] = np.cumsum(np.log(np.arange(1, n_max + 1, dtype=float)))
    return table


def project_sfs(
    k: np.ndarray, n: np.ndarray, m: int, chunk: int | None = None
) -> tuple[np.ndarray, int]:
    """Hypergeometric projection of per-site counts down to m haplotypes.

    Each site with n_i >= m contributes P(j | k_i, n_i, m) =
    C(k_i, j) C(n_i − k_i, m − j) / C(n_i, m) to class j — the expected
    spectrum of a random subsample of size m (Marth et al. 2004;
    Gutenkunst et al. 2009). Sites with n_i < m are dropped.

    Returns (unfolded spectrum of length m+1, number of sites used).
    """
    k = np.asarray(k, dtype=np.int64)
    n = np.asarray(n, dtype=np.int64)
    if m < 1:
        raise ValueError(f"projection size must be >= 1, got {m}")
    keep = n >= m
    out = np.zeros(m + 1)
    if not keep.any():
        return out, 0
    if (n[keep] == m).all():  # nothing to project: plain counts, exact
        return sfs_counts(k[keep], m), int(keep.sum())
    if chunk is None:
        # Each chunk materialises chunk × (m+1) float arrays; keep that
        # near 2M cells (~16 MB each) whatever the projection size.
        chunk = max(1, 2_000_000 // (m + 1))
    pairs, weight = np.unique(
        np.stack([n[keep], k[keep]], axis=1), axis=0, return_counts=True
    )
    lf = _log_factorial_table(int(pairs[:, 0].max()))
    j = np.arange(m + 1)
    for start in range(0, len(pairs), chunk):
        nn = pairs[start : start + chunk, 0][:, None]
        kk = pairs[start : start + chunk, 1][:, None]
        valid = (j[None, :] <= kk) & (m - j[None, :] <= nn - kk)
        jj = np.where(valid, j[None, :], 0)
        mj = np.where(valid, m - j[None, :], 0)
        log_p = (
            lf[kk]
            - lf[jj]
            - lf[np.maximum(kk - jj, 0)]
            + lf[nn - kk]
            - lf[mj]
            - lf[np.maximum(nn - kk - mj, 0)]
            - (lf[nn] - lf[m] - lf[nn - m])
        )
        p = np.where(valid, np.exp(log_p), 0.0)
        out += (p * weight[start : start + chunk, None]).sum(axis=0)
    return out, int(keep.sum())


def projected_thetas(
    k: np.ndarray, n: np.ndarray, m: int
) -> tuple[SfsThetas | None, int]:
    """θ estimators of the SFS projected to m haplotypes, in closed form.

    Equal to `thetas_from_sfs(project_sfs(k, n, m)[0])` but O(1) per
    site instead of O(m): with J ~ Hypergeometric(n, k, m) at a site,

        S   = Σ_sites 1 − P(J=0) − P(J=m)
        θ_π = Σ_sites 2k(n−k)/(n(n−1))           (unchanged by projection)
        θ_L = Σ_sites (E[J] − m·P(J=m)) / (m−1)
        θ_H = Σ_sites 2(E[J²] − m²·P(J=m)) / (m(m−1))

    using E[J] = mk/n and Var(J) = m(k/n)(1−k/n)(n−m)/(n−1). Sites with
    n < m are dropped. Returns (thetas or None if no site, sites used).
    """
    if m < 2:
        raise ValueError(f"projection size must be >= 2, got {m}")
    k = np.asarray(k, dtype=np.int64)
    n = np.asarray(n, dtype=np.int64)
    keep = n >= m
    if not keep.any():
        return None, 0
    k, n = k[keep], n[keep]
    lf = _log_factorial_table(int(n.max()))

    def log_comb(a: np.ndarray, b: int) -> np.ndarray:
        ok = a >= b
        safe = np.where(ok, a, b)
        return np.where(ok, lf[safe] - lf[b] - lf[safe - b], -np.inf)

    log_total = log_comb(n, m)
    p0 = np.exp(log_comb(n - k, m) - log_total)
    pm = np.exp(log_comb(k, m) - log_total)
    nf, kf = n.astype(float), k.astype(float)
    mean = m * kf / nf
    var = np.divide(
        m * (kf / nf) * (1 - kf / nf) * (nf - m),
        nf - 1,
        out=np.zeros_like(nf),
        where=nf > 1,
    )
    s = float((1 - p0 - pm).sum())
    pi = float(pi_per_site(k, n).sum())
    theta_l = float((mean - m * pm).sum() / (m - 1))
    theta_h = float((2 * (var + mean**2 - m * m * pm)).sum() / (m * (m - 1)))
    thetas = SfsThetas(
        n=m, s=s, pi=pi, theta_w=s / a1(m), theta_l=theta_l, theta_h=theta_h
    )
    return thetas, int(keep.sum())


def sfs_counts(k: np.ndarray, n: int) -> np.ndarray:
    """Unfolded SFS of sites that all have exactly n haplotypes."""
    return np.bincount(np.asarray(k, dtype=np.int64), minlength=n + 1).astype(float)


# ───────────────────────────── Hudson F_ST ──────────────────────────────


def hudson_fst_components(
    k1: np.ndarray, n1: np.ndarray, k2: np.ndarray, n2: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Per-site numerator and denominator of Hudson's F_ST as written by
    Bhatia, Patterson, Sankararaman & Price (2013), eq. 10:

        N = (p1 − p2)² − p1(1−p1)/(n1−1) − p2(1−p2)/(n2−1)
        D = p1(1−p2) + p2(1−p1)

    Sites where either group has fewer than 2 haplotypes contribute 0/0."""
    k1, n1, k2, n2 = (np.asarray(x, dtype=float) for x in (k1, n1, k2, n2))
    ok = (n1 >= 2) & (n2 >= 2)
    num = np.zeros_like(n1)
    den = np.zeros_like(n1)
    p1 = np.divide(k1, n1, out=np.zeros_like(n1), where=ok)
    p2 = np.divide(k2, n2, out=np.zeros_like(n2), where=ok)
    num[ok] = (
        (p1[ok] - p2[ok]) ** 2
        - p1[ok] * (1 - p1[ok]) / (n1[ok] - 1)
        - p2[ok] * (1 - p2[ok]) / (n2[ok] - 1)
    )
    den[ok] = p1[ok] * (1 - p2[ok]) + p2[ok] * (1 - p1[ok])
    return num, den


def hudson_fst(
    k1: np.ndarray, n1: np.ndarray, k2: np.ndarray, n2: np.ndarray
) -> float | None:
    """Hudson F_ST as a ratio of averages, Σ N / Σ D (Bhatia et al. 2013
    recommend this over the average of per-site ratios). None when Σ D = 0
    (no site polymorphic in the pair)."""
    num, den = hudson_fst_components(k1, n1, k2, n2)
    total = float(den.sum())
    if total <= 0:
        return None
    return float(num.sum()) / total


# ───────────────────────── heterozygosity / F ───────────────────────────


@dataclass(frozen=True)
class Heterozygosity:
    """Mean observed and expected heterozygosity and F = 1 − Ho/He."""

    ho: float | None
    he: float | None
    f: float | None


def heterozygosity(
    het: np.ndarray, called: np.ndarray, k: np.ndarray, n: np.ndarray
) -> Heterozygosity:
    """Ho = mean over sites of het/called individuals; He = mean over the
    same sites of Nei's (1978) unbiased 2k(n−k)/(n(n−1)); F = 1 − Ho/He
    (a ratio of averages over sites). Sites with fewer than 2 haplotypes
    or no called individual are skipped."""
    het, called, k, n = (np.asarray(x, dtype=float) for x in (het, called, k, n))
    ok = (called > 0) & (n >= 2)
    if not ok.any():
        return Heterozygosity(None, None, None)
    ho = float((het[ok] / called[ok]).mean())
    he = float(pi_per_site(k[ok], n[ok]).mean())
    f = None if he == 0 else 1.0 - ho / he
    return Heterozygosity(ho, he, f)
