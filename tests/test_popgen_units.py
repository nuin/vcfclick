"""Unit tests for popgen helpers that need no database: chromosome
classification, window slicing, vectorised polarisation, projection
chunking."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from popgen import estimators as est
from popgen.analysis import GroupData, Prepared, is_autosome, iter_windows
from popgen.ancestral import aa_code, orientation_from_codes, polarise
from popgen.counts import Region


@pytest.mark.parametrize(
    "chrom",
    ["1", "chr1", "22", "chr22", "chr1_KI270706v1_random", "chr17_GL000205v2_random",
     "2L", "chr3R", "NC_000001.11", "NC_000022.11", "scaffold_12"],
)  # fmt: skip
def test_autosomes(chrom):
    assert is_autosome(chrom)


@pytest.mark.parametrize(
    "chrom",
    ["X", "chrX", "Y", "chrY", "chrM", "MT", "chrMT", "XY", "PAR1", "chrX_PAR1",
     "chrX_KI270880v1_alt", "chrY_KI270740v1_random", "chrUn_KI270302v1",
     "chrUn_JTFH01000001v1_decoy", "chrEBV", "hs37d5", "GL000192.1",
     "HLA-A*01:01:01:01", "NC_000023.11", "NC_000024.10", "NC_012920.1", "W", "Z"],
)  # fmt: skip
def test_non_autosomes(chrom):
    assert not is_autosome(chrom)


def test_orientation_codes_match_polarise():
    alleles = ["A", "G", "a", "g", "C", ".", None, "N", "G|||"]
    for ref, alt, aa, mode in itertools.product(
        ["A"], ["G"], alleles, ["aa", "aa-high", "ref", "none"]
    ):
        expected = polarise(ref, alt, aa, mode)
        got = int(orientation_from_codes([aa_code(ref, alt, aa)], mode)[0])
        assert got == (-1 if expected is None else expected), (aa, mode)


def _prep(chroms: dict[str, list[int]]) -> Prepared:
    pos, slices, offset = [], [], 0
    for name, ps in chroms.items():
        pos += ps
        slices.append((name, offset, offset + len(ps)))
        offset += len(ps)
    n = len(pos)
    g = GroupData(
        size=2,
        called=np.full(n, 2, np.uint16),
        k=np.ones(n, np.uint16),
        het=np.ones(n, np.uint16),
    )
    return Prepared(
        ingest_id="x",
        grouping="all",
        ancestral="none",
        missing_data_tracked=True,
        regions=[],
        stored_chroms=list(chroms),
        chrom_slices=slices,
        pos=np.array(pos, dtype=np.int64),
        orientation=np.full(n, -1, np.int8),
        groups={"all": g},
        report={},
    )


def test_window_slices_match_masks():
    rng = np.random.default_rng(5)
    chroms = {c: sorted(rng.integers(1, 10_000, 300).tolist()) for c in ("1", "2")}
    prep = _prep(chroms)
    seen = []
    for chrom, start, end, sel in iter_windows(prep, 1000, 400):
        name, a, b = next(s for s in prep.chrom_slices if s[0] == chrom)
        on = np.zeros(prep.n_sites, bool)
        on[a:b] = True
        mask = on & (prep.pos >= start) & (prep.pos <= end)
        assert np.array_equal(np.flatnonzero(mask), np.arange(sel.start, sel.stop))
        seen.append((chrom, start, end))
    # Tiling: start at 1, step 400, stop at the first window reaching the end.
    first = [w for w in seen if w[0] == "1"]
    assert first[0][1:] == (1, 1000)
    assert first[-1][2] == chroms["1"][-1]
    assert all(e < chroms["1"][-1] for _, _, e in first[:-1])


def test_window_region_on_stored_name():
    prep = _prep({"1": [5, 50, 500]})
    prep.regions = [Region("chr1", 1, 100)]
    out = list(iter_windows(prep, 100, 100))
    assert [(c, s, e, (x.start, x.stop)) for c, s, e, x in out] == [
        ("1", 1, 100, (0, 2))
    ]


def test_projection_chunking_does_not_change_result():
    rng = np.random.default_rng(11)
    n = rng.integers(40, 80, 500)
    k = np.array([rng.integers(0, x + 1) for x in n])
    whole, _ = est.project_sfs(k, n, 30, chunk=10_000)
    tiny, _ = est.project_sfs(k, n, 30, chunk=3)
    auto, _ = est.project_sfs(k, n, 30)
    assert tiny == pytest.approx(whole)
    assert auto == pytest.approx(whole)


@pytest.mark.parametrize("m", [2, 3, 7, 24, 60])
def test_closed_form_projected_thetas_equal_full_projection(m):
    rng = np.random.default_rng(m)
    n = rng.integers(m, 2 * m + 20, 400)
    n[:50] = m  # include sites already at size m
    k = np.array([rng.integers(0, x + 1) for x in n])
    k[:5], k[5:10] = 0, n[5:10]  # monomorphic both ways
    full = est.thetas_from_sfs(est.project_sfs(k, n, m)[0])
    fast, used = est.projected_thetas(k, n, m)
    assert used == 400
    for attr in ("s", "pi", "theta_w", "theta_l", "theta_h"):
        assert getattr(fast, attr) == pytest.approx(getattr(full, attr), rel=1e-10), (
            attr
        )


def test_closed_form_projected_thetas_drops_small_sites():
    t, used = est.projected_thetas(np.array([1, 1]), np.array([3, 10]), 4)
    assert used == 1
    assert est.projected_thetas(np.array([1]), np.array([3]), 4) == (None, 0)
