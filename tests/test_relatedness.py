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
