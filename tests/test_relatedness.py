"""Tests for `vcfclick db relatedness` — KING-robust kinship between samples.

Fixture: tests/fixtures/related.vcf (see make_related_vcf.py). F and M are
parents of full sibs C1 and C2; D duplicates C1; U1 and U2 are unrelated;
U2 has ~3% no-calls. related.ped wrongly declares U1 a child of F and M.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
FIX = Path(__file__).parent / "fixtures"
REL_VCF = FIX / "related.vcf"
REL_PED = FIX / "related.ped"


def _vc(home: Path, *args: str, ok: bool = True):
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    if ok:
        assert r.returncode == 0, f"{' '.join(args)} failed:\n{r.stderr}"
    return r


def _setup(home: Path, *, keep_reference: bool = False, ped: bool = True) -> None:
    _vc(home, "db", "create", "r")
    args = [
        "db",
        "ingest",
        "r",
        str(REL_VCF),
        "--cohort",
        "fam",
        "--ingest-id",
        "i1",
        "--serial",
    ]
    if keep_reference:
        args.append("--keep-reference")
    _vc(home, *args)
    if ped:
        _vc(home, "db", "ped", "r", str(REL_PED))


def _run(home: Path, *extra: str) -> dict:
    out = _vc(
        home, "db", "relatedness", "r", "--format", "json", "--all", *extra
    ).stdout
    results = json.loads(out)
    assert len(results) == 1, results
    return results[0]


def _pair(res: dict, a: str, b: str) -> dict:
    for p in res["pairs"]:
        if {p["sample_a"], p["sample_b"]} == {a, b}:
            return p
    raise AssertionError(f"pair {a}-{b} not reported")


# --- the estimator itself (no database) ---------------------------------------


def test_king_robust_on_a_hand_built_matrix():
    """Two identical samples -> kinship 0.5; opposite homozygotes -> negative."""
    from cli.db_relatedness import king_robust

    g = np.array(
        [
            # a  b  c
            [1, 1, 0],
            [2, 2, 0],
            [1, 1, 2],
            [0, 0, 2],
            [1, 1, 1],
        ],
        dtype=np.int8,
    )
    k = king_robust(g)
    assert k["kinship"][0, 1] == pytest.approx(0.5)
    assert k["ibs0"][0, 1] == 0
    # a vs c: sites 2 (2 vs 0) and 4 (0 vs 2) are opposite homozygotes
    assert k["ibs0"][0, 2] == 2
    assert k["kinship"][0, 2] < 0


def test_missing_calls_are_excluded_pairwise():
    """-1 = no-call: a site only counts for a pair where both samples are called."""
    from cli.db_relatedness import king_robust

    g = np.array([[1, 1], [1, -1], [2, 2]], dtype=np.int8)
    k = king_robust(g)
    assert k["n_sites"][0, 1] == 2


# --- end to end on the fixture --------------------------------------------------


def test_relationships_are_classified(vcfclick_home):
    _setup(vcfclick_home)
    res = _run(vcfclick_home)
    assert res["ingest_id"] == "i1"
    assert res["mode"] == "sparse"

    dup = _pair(res, "C1", "D")
    assert dup["kinship"] > 0.45 and dup["relationship"] == "duplicate"

    for parent in ("F", "M"):
        for kid in ("C1", "C2"):
            p = _pair(res, parent, kid)
            assert 0.2 < p["kinship"] < 0.3, p
            assert p["ibs0"] < 0.005, p
            assert p["relationship"] == "parent-child", p

    sibs = _pair(res, "C1", "C2")
    assert 0.18 < sibs["kinship"] < 0.32, sibs
    assert sibs["ibs0"] > 0.01, sibs
    assert sibs["relationship"] == "full-siblings", sibs

    for a, b in [("F", "M"), ("U1", "C1"), ("U1", "F"), ("U2", "M")]:
        p = _pair(res, a, b)
        assert abs(p["kinship"]) < 0.06, p
        assert p["relationship"] == "unrelated", p


def test_default_hides_unrelated_pairs(vcfclick_home):
    _setup(vcfclick_home, ped=False)
    out = _vc(vcfclick_home, "db", "relatedness", "r", "--format", "json").stdout
    pairs = json.loads(out)[0]["pairs"]
    assert pairs, "related pairs should be reported"
    assert all(p["relationship"] != "unrelated" for p in pairs)


def test_pedigree_check_flags_a_false_parent(vcfclick_home):
    _setup(vcfclick_home)
    res = _run(vcfclick_home)
    checks = {(c["child"], c["parent"]): c for c in res["pedigree_checks"]}
    assert checks[("C1", "F")]["verdict"] == "ok"
    assert checks[("C2", "M")]["verdict"] == "ok"
    assert checks[("U1", "F")]["verdict"] == "mismatch"
    assert checks[("U1", "M")]["verdict"] == "mismatch"


def test_keep_reference_mode_excludes_no_calls(vcfclick_home):
    """With 0/0 stored, an absent genotype is a real no-call and is skipped,
    so U2's pairs rest on fewer sites than a pair of fully called samples."""
    _setup(vcfclick_home, keep_reference=True, ped=False)
    res = _run(vcfclick_home)
    assert res["mode"] == "keep-reference"
    full = _pair(res, "F", "M")["n_sites"]
    with_nocalls = _pair(res, "U2", "M")["n_sites"]
    assert with_nocalls < full
    assert full == 3003  # 3000 Mendelian sites + 3 planted de novo


def test_table_output_and_missing_db(vcfclick_home):
    _setup(vcfclick_home)
    out = _vc(vcfclick_home, "db", "relatedness", "r").stdout
    assert "parent-child" in out and "duplicate" in out
    r = _vc(vcfclick_home, "db", "relatedness", "nope", ok=False)
    assert r.returncode != 0 and "does not exist" in r.stderr


def test_single_region_data_is_not_classified(vcfclick_home, tiny_vcf):
    """A few sites in one small region can't support kinship: unrelated people
    sharing a haplotype would look like duplicates. Pairs stay unclassified
    and the output says why, unless --force."""
    _vc(vcfclick_home, "db", "create", "t")
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "t",
        str(tiny_vcf),
        "--ingest-id",
        "i1",
        "--serial",
    )
    res = json.loads(
        _vc(vcfclick_home, "db", "relatedness", "t", "--format", "json", "--all").stdout
    )[0]
    assert res["warning"] and "genome-wide" in res["warning"]
    assert all(p["relationship"] == "insufficient-data" for p in res["pairs"])
    table = _vc(vcfclick_home, "db", "relatedness", "t").stdout
    assert "WARNING" in table and "no related pairs" in table


def test_genome_scale_fixture_has_no_warning(vcfclick_home):
    _setup(vcfclick_home, ped=False)
    res = _run(vcfclick_home)
    assert res["warning"] is None
    assert res["span_bp"] > 100_000_000


# --- chunked KING ---------------------------------------------------------------


def _fixture_matrix() -> np.ndarray:
    """related.vcf as a (sites x samples) int8 matrix, -1 for ./. (read
    from the text, no database)."""
    code = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}
    rows = []
    for line in REL_VCF.read_text().splitlines():
        if line.startswith("#"):
            continue
        rows.append([code.get(c.split(":")[0], -1) for c in line.split("\t")[9:]])
    return np.array(rows, dtype=np.int8)


def _king_single_pass(g: np.ndarray) -> dict[str, np.ndarray]:
    """The estimator as it was before chunking: one pass over every site."""
    called = (g >= 0).astype(np.float32)
    het = (g == 1).astype(np.float32)
    homref = (g == 0).astype(np.float32)
    homalt = (g == 2).astype(np.float32)
    hethet = het.T @ het
    ibs0 = homref.T @ homalt
    ibs0 = ibs0 + ibs0.T
    hgc = het.T @ called
    denom = hgc + hgc.T
    with np.errstate(divide="ignore", invalid="ignore"):
        kinship = np.where(denom > 0, (hethet - 2 * ibs0) / denom, np.nan)
    return {
        "het_pair": denom,
        "kinship": kinship,
        "ibs0": ibs0,
        "hethet": hethet,
        "n_sites": called.T @ called,
    }


@pytest.mark.parametrize("chunk", [1, 7, 100, 1000, 100_000])
def test_chunked_king_matches_single_pass_on_fixture(chunk):
    """Accumulating the products over site chunks gives the kinship of one
    pass over every site, to 1e-6, on the related.vcf genotypes."""
    from cli.relatedness_analysis import king_robust

    g = _fixture_matrix()
    assert (g == -1).any()  # U2's no-calls are in the matrix
    want = _king_single_pass(g)
    got = king_robust(g, chunk_sites=chunk)
    np.testing.assert_allclose(got["kinship"], want["kinship"], atol=1e-6)
    for key in ("het_pair", "ibs0", "hethet", "n_sites"):
        np.testing.assert_array_equal(got[key], want[key].astype(np.int64))


def test_chunked_king_bounds_peak_memory():
    """Peak memory follows the chunk, not the matrix: four float32 copies of
    a 20,000 x 200 matrix alone would be 64 MB."""
    import tracemalloc

    from cli.relatedness_analysis import king_robust

    rng = np.random.default_rng(0)
    g = rng.integers(-1, 3, size=(20_000, 200), dtype=np.int8)
    tracemalloc.start()
    king_robust(g, chunk_sites=1000)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 16_000_000, peak


# --- thinning in SQL --------------------------------------------------------------


def _exact_thin(total: int, max_sites: int) -> list[int]:
    """Indices floor(i * total / max_sites), i < max_sites, in exact integers."""
    if total <= max_sites:
        return list(range(total))
    return [i * total // max_sites for i in range(max_sites)]


@pytest.mark.parametrize(
    "total,max_sites", [(10, 3), (30, 22), (3003, 1000), (3003, 7), (97, 96)]
)
def test_thin_predicate_keeps_evenly_spread_sites(total, max_sites):
    """The SQL predicate on a site's 0-based rank keeps exactly max_sites
    sites, the ones at floor(i * total / max_sites)."""
    from cli.relatedness_analysis import _thin_keeps

    kept = [k for k in range(total) if _thin_keeps(k, total, max_sites)]
    assert kept == _exact_thin(total, max_sites)


@pytest.mark.parametrize("max_sites", ["1000", "7"])
def test_thinned_kinship_matches_the_thinned_fixture(vcfclick_home, max_sites):
    """Thinning in SQL keeps the same sites as thinning the fixture rows, so
    kinship equals the single-pass estimator on those rows."""
    _setup(vcfclick_home, keep_reference=True, ped=False)
    res = _run(vcfclick_home, "--max-sites", max_sites, "--force")
    g = _fixture_matrix()
    rows = _exact_thin(g.shape[0], int(max_sites))
    assert res["sites_used"] == len(rows) and res["n_sites"] == g.shape[0]
    want = _king_single_pass(g[rows])
    names = ["F", "M", "C1", "C2", "D", "U1", "U2"]
    for p in res["pairs"]:
        i, j = names.index(p["sample_a"]), names.index(p["sample_b"])
        assert p["n_sites"] == int(want["n_sites"][i, j])
        if p["kinship"] is not None:
            assert p["kinship"] == pytest.approx(want["kinship"][i, j], abs=5e-5)


_FETCH_COUNT = """
import json, sys
import cli.main  # noqa: F401  (registers the command modules)
import cli.relatedness_analysis as ra
from storage import get_session
fetched = []
real = ra._arrow
def spy(sess, sql):
    t = real(sess, sql)
    fetched.append(t.num_rows)
    return t
ra._arrow = spy
res = ra.relatedness(get_session("r"), "i1", max_sites=7, include_all=True,
                     min_kinship=0.0, force=True)
print("RESULT " + json.dumps({"fetched": fetched, "used": res["sites_used"]}))
"""


def test_genotype_fetch_is_restricted_to_the_thinned_sites(vcfclick_home, run_python):
    """--max-sites bounds what leaves the database: with 7 sites kept of
    3,003, at most 7 x 7 genotype rows are fetched."""
    _setup(vcfclick_home, ped=False)
    out = run_python(vcfclick_home, _FETCH_COUNT, VCFCLICK_DB_NAME="r")
    line = [ln for ln in out.splitlines() if ln.startswith("RESULT ")][-1]
    got = json.loads(line.removeprefix("RESULT "))
    assert got["used"] == 7
    assert got["fetched"] and all(n <= 7 * 7 for n in got["fetched"]), got


_SAMPLE_GUARD = """
import click
import cli.main  # noqa: F401  (registers the command modules)
import cli.relatedness_analysis as ra
from storage import get_session
ra._MAX_SAMPLES = 3
def boom(*a, **k):
    raise AssertionError("genotype matrix built before the sample-count check")
ra._genotype_matrix = boom
try:
    ra.relatedness(get_session("r"), "i1", max_sites=200000, include_all=False,
                   min_kinship=0.0442)
except click.ClickException as e:
    print("REFUSED " + e.message)
"""


@pytest.mark.parametrize("bad", ["0", "-3", "16000001"])
def test_max_sites_out_of_range_is_a_usage_error(vcfclick_home, bad):
    """--max-sites 0 divided by zero; a negative value kept no sites; above
    2^24 the float32 counts in king_robust are no longer exact."""
    _setup(vcfclick_home, ped=False)
    r = _vc(vcfclick_home, "db", "relatedness", "r", "--max-sites", bad, ok=False)
    assert r.returncode == 2, r.stderr
    assert "Traceback" not in r.stderr and "--max-sites" in r.stderr


def test_sample_cap_is_checked_before_the_matrix_is_built(vcfclick_home, run_python):
    _setup(vcfclick_home, ped=False)
    out = run_python(vcfclick_home, _SAMPLE_GUARD, VCFCLICK_DB_NAME="r")
    assert "REFUSED 7 samples" in out, out
