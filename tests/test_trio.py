"""Tests for `vcfclick db trio` — Mendelian inheritance-model filtering.

Fixture trio (tests/fixtures/trio.vcf.gz), ingested --keep-reference:
  chr1:100  CHILD het, FATHER 0/0,  MOTHER 0/0   AF=0.01  -> de novo
  chr1:200  CHILD homalt, FATHER het, MOTHER het AF=0.20  -> recessive
  chr1:300  CHILD het, FATHER het,  MOTHER 0/0   AF=0.05  -> dominant (paternal)
  chr1:400  CHILD het, FATHER 0/0,  MOTHER het   AF=0.05  -> dominant (maternal)
  chr1:500  CHILD het, FATHER ./.,  MOTHER 0/0   AF=0.01  -> NOT de novo
                                                            (father no-call)

The chr1:500 case is the defensibility crux: a naive "neither parent
carries" rule would call it de novo, but the father is a no-call, not
a confident 0/0, so it must be excluded.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
FIXTURES = Path(__file__).parent / "fixtures"
TRIO_VCF = FIXTURES / "trio.vcf.gz"
TRIO_PED = FIXTURES / "trio.ped"


def _vc(home: Path, *args: str, expect_failure: bool = False):
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    if expect_failure:
        assert r.returncode != 0, f"expected failure but got rc=0:\n{r.stdout}"
    else:
        assert r.returncode == 0, (
            f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
            f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
        )
    return r


def _setup(home: Path, *, keep_reference: bool = True) -> None:
    _vc(home, "db", "create", "trio")
    args = [
        "db",
        "ingest",
        "trio",
        str(TRIO_VCF),
        "--cohort",
        "fam",
        "--ingest-id",
        "fam1",
        "--serial",
    ]
    if keep_reference:
        args.append("--keep-reference")
    _vc(home, *args)
    _vc(home, "db", "ped", "trio", str(TRIO_PED))


def _trio(home: Path, *args: str, expect_failure: bool = False):
    return _vc(
        home,
        "db",
        "trio",
        "trio",
        "--proband",
        "CHILD",
        *args,
        expect_failure=expect_failure,
    )


def test_denovo_excludes_no_call_parent_site(vcfclick_home):
    """De novo must be exactly chr1:100 — NOT chr1:500, where the father
    is a no-call (./.), so 0/0 is unproven. This is the defensibility
    guarantee keep-reference enables."""
    _setup(vcfclick_home)
    out = _trio(vcfclick_home, "--category", "denovo").stdout
    assert "denovo candidates: 1" in out
    assert "chr1:100" in out
    assert "chr1:500" not in out, "no-call-parent site must not be de novo"


def test_all_counts_with_default_rarity(vcfclick_home):
    """Default --max-af 0.01: only the AF=0.01 de-novo site survives the
    rarity filter; recessive (AF 0.20) and dominant (AF 0.05) are
    filtered out."""
    _setup(vcfclick_home)
    out = _trio(vcfclick_home).stdout
    assert "denovo          1" in out
    assert "recessive       0" in out
    assert "dominant        0" in out


def test_relaxed_rarity_surfaces_recessive_and_dominant(vcfclick_home):
    """With --max-af 0.3 the recessive (1) and both dominant (2) sites
    appear."""
    _setup(vcfclick_home)
    out = _trio(vcfclick_home, "--max-af", "0.3").stdout
    assert "denovo          1" in out
    assert "recessive       1" in out
    assert "dominant        2" in out


def test_recessive_detail(vcfclick_home):
    _setup(vcfclick_home)
    out = _trio(vcfclick_home, "--category", "recessive", "--max-af", "0.3").stdout
    assert "recessive candidates: 1" in out
    assert "chr1:200" in out
    assert "proband_gt=2" in out  # hom-alt child


def test_quality_gate_filters(vcfclick_home):
    """A high --min-dp above the fixture's depths drops every candidate,
    proving the gate is wired (fixture DP is 28-40)."""
    _setup(vcfclick_home)
    out = _trio(vcfclick_home, "--category", "denovo", "--min-dp", "1000").stdout
    assert "denovo candidates: 0" in out


def test_denovo_needs_keep_reference(vcfclick_home):
    """Without --keep-reference, parents have no stored 0/0 rows, so the
    de-novo inner-join to f.gt=0 finds nothing → 0 candidates, and a
    note tells the user to re-ingest with --keep-reference."""
    _setup(vcfclick_home, keep_reference=False)
    r = _trio(vcfclick_home, "--category", "denovo")
    assert "denovo candidates: 0" in r.stdout
    assert "keep-reference" in (r.stdout + r.stderr)


def test_trio_requires_pedigree(vcfclick_home):
    """Without a loaded pedigree, trio can't resolve parents."""
    _vc(vcfclick_home, "db", "create", "trio")
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "trio",
        str(TRIO_VCF),
        "--cohort",
        "fam",
        "--ingest-id",
        "fam1",
        "--serial",
        "--keep-reference",
    )
    r = _trio(vcfclick_home, expect_failure=True)
    assert "pedigree" in (r.stdout + r.stderr).lower()
