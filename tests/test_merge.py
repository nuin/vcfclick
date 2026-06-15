"""Tests for `vcfclick merge` — joint VCF from per-sample VCFs.

The merger wraps `bcftools merge` to produce the one joint VCF that trio
analysis needs (all samples under one ingest_id). These tests lock the
vcfclick-level contract: disjoint-sample validation, the joint output's
sample set, that absent samples fill ./. (not 0/0), and that the merged
VCF is ingest-ready (stays decomposed).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
FIXTURES = Path(__file__).parent / "fixtures"
PSA = FIXTURES / "per_sample_a.vcf.gz"
PSB = FIXTURES / "per_sample_b.vcf.gz"

# Every test here shells out to bcftools (vcfclick merge wraps it). Skip
# the whole module cleanly where bcftools isn't installed rather than
# erroring.
pytestmark = pytest.mark.skipif(
    shutil.which("bcftools") is None, reason="bcftools not installed"
)


def _vc(*args: str, expect_failure: bool = False) -> subprocess.CompletedProcess:
    r = subprocess.run([VCFCLICK_BIN, *args], cwd=REPO, capture_output=True, text=True)
    if expect_failure:
        assert r.returncode != 0, f"expected failure but got rc=0:\n{r.stdout}"
    else:
        assert r.returncode == 0, (
            f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
            f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
        )
    return r


def test_merge_produces_joint_multisample_vcf(tmp_path):
    """Two single-sample VCFs merge into one VCF carrying both samples."""
    from cyvcf2 import VCF

    out = tmp_path / "joint.vcf.gz"
    _vc("merge", str(PSA), str(PSB), "-o", str(out))
    assert out.exists()

    v = VCF(str(out))
    assert sorted(v.samples) == ["SAMPLE_A", "SAMPLE_B"]


def test_merge_fills_absent_sample_with_missing_not_homref(tmp_path):
    """A sample absent at a site the other sample called must be ./.
    (cyvcf2 gt_type 3 = unknown), NOT 0/0 — a variant-only VCF doesn't
    assert hom-reference where it's silent. This is the same limitation
    behind 'candidate' (not rigorous) de-novo detection."""
    from cyvcf2 import VCF

    out = tmp_path / "joint.vcf.gz"
    _vc("merge", str(PSA), str(PSB), "-o", str(out))

    v = VCF(str(out))
    samples = v.samples
    a_idx, b_idx = samples.index("SAMPLE_A"), samples.index("SAMPLE_B")

    # Find at least one site where exactly one sample is missing (gt_type
    # 3) and confirm the absent one is missing, not hom-ref (gt_type 0).
    saw_missing = False
    for rec in v:
        a, b = rec.gt_types[a_idx], rec.gt_types[b_idx]
        if 3 in (a, b):
            saw_missing = True
            # The present sample should be a real call (het/hom-alt),
            # the absent one missing — never silently 0/0.
            assert 0 not in (a, b), (
                f"site {rec.CHROM}:{rec.POS} filled an absent sample with "
                f"hom-ref (0) instead of missing (3): a={a} b={b}"
            )
    assert saw_missing, "fixtures should yield at least one one-sample-only site"


def test_merged_vcf_ingests_cleanly(tmp_path):
    """The joint VCF must stay decomposed (no multi-allelic sites) so
    vcfclick ingest accepts it. Both samples land in the cohort."""
    out = tmp_path / "joint.vcf.gz"
    _vc("merge", str(PSA), str(PSB), "-o", str(out))

    env_home = tmp_path / "home"
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(env_home)
    env["VCFCLICK_BACKEND"] = os.environ.get("VCFCLICK_BACKEND", "chdb")

    def run(*args):
        r = subprocess.run(
            [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
        )
        assert r.returncode == 0, f"{' '.join(args)} failed:\n{r.stderr}"
        return r

    run("db", "create", "fam")
    run(
        "db",
        "ingest",
        "fam",
        str(out),
        "--cohort",
        "trio",
        "--ingest-id",
        "fam1",
        "--serial",
    )
    r = run(
        "db",
        "query",
        "fam",
        "SELECT count(DISTINCT sample_id) FROM samples FORMAT TabSeparated",
    )
    assert r.stdout.strip() == "2"


def test_merge_rejects_overlapping_sample_names(tmp_path):
    """Merging a file with itself (same sample twice) must error with a
    clear vcfclick message, not a raw bcftools failure."""
    r = _vc(
        "merge",
        str(PSA),
        str(PSA),
        "-o",
        str(tmp_path / "x.vcf.gz"),
        expect_failure=True,
    )
    assert "disjoint" in (r.stdout + r.stderr).lower()


def test_merge_requires_two_inputs(tmp_path):
    r = _vc("merge", str(PSA), "-o", str(tmp_path / "x.vcf.gz"), expect_failure=True)
    assert "two" in (r.stdout + r.stderr).lower()
