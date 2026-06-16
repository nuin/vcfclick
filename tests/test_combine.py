"""Tests for `vcfclick combine` — the GATK3 CombineVariants reimplementation.

combine unions VCF call sets that may share samples, resolves overlaps by
PRIORITY (input order), annotates set= provenance, and can keep only
consensus sites (--min-callsets). It is native (cyvcf2 read + text write),
so these tests need no bcftools.

Fixtures (callset_a.vcf, callset_b.vcf), sites and shared sample S2:

  site      A (S1,S2)        B (S2,S3)        union outcome
  chr1:100  S1=0/1 S2=0/0    —                set=callset_a; S3 absent → ./.
  chr1:200  S1=1/1 S2=0/1    S2=1/1 S3=0/1    Intersection; S2 PRIORITIZEs A → 0/1
  chr1:300  S1=0/1 S2=./.    S2=0/1 S3=0/0    Intersection; S2 missing in A → B's 0/1
  chr1:400  —                S2=0/0 S3=1/1    set=callset_b; S1 absent → ./.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from cyvcf2 import VCF

from ingest.combine import CombineError, combine_vcfs

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
FIXTURES = Path(__file__).parent / "fixtures"
A = FIXTURES / "callset_a.vcf"
B = FIXTURES / "callset_b.vcf"

# cyvcf2 gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN(./.), 3=HOM_ALT
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3


def _read(path: Path):
    """Return (samples, {pos: {"set": str, "gt": {sample: gt_type}}})."""
    v = VCF(str(path))
    samples = list(v.samples)
    sites = {}
    for rec in v:
        sites[rec.POS] = {
            "set": rec.INFO.get("set"),
            "gt": dict(zip(samples, rec.gt_types)),
        }
    return samples, sites


def test_union_sites_and_sample_order(tmp_path):
    """All sites across inputs appear; samples are the union in
    first-appearance (priority) order: S1, S2, S3."""
    out = tmp_path / "combined.vcf"
    combine_vcfs([A, B], out)

    samples, sites = _read(out)
    assert samples == ["S1", "S2", "S3"]
    assert sorted(sites) == [100, 200, 300, 400]


def test_set_provenance_annotation(tmp_path):
    """set= names the inputs each site came from; Intersection when all."""
    out = tmp_path / "combined.vcf"
    combine_vcfs([A, B], out)

    _, sites = _read(out)
    assert sites[100]["set"] == "callset_a"
    assert sites[200]["set"] == "Intersection"
    assert sites[300]["set"] == "Intersection"
    assert sites[400]["set"] == "callset_b"


def test_prioritize_takes_higher_priority_genotype(tmp_path):
    """At chr1:200 sample S2 differs between inputs (A=0/1, B=1/1). A is
    listed first (higher priority), so the output keeps A's het call, not
    B's hom-alt."""
    out = tmp_path / "combined.vcf"
    combine_vcfs([A, B], out)

    _, sites = _read(out)
    assert sites[200]["gt"]["S2"] == HET  # A's 0/1, not B's 1/1 (HOM_ALT)


def test_prioritize_skips_missing_and_uses_next_input(tmp_path):
    """At chr1:300 the higher-priority input A has S2 as ./. (no info), so
    the next input B's 0/1 fills it — PRIORITIZE means highest-priority
    *non-missing* call."""
    out = tmp_path / "combined.vcf"
    combine_vcfs([A, B], out)

    _, sites = _read(out)
    assert sites[300]["gt"]["S2"] == HET  # B's 0/1 (A was missing)


def test_absent_sample_is_missing_not_homref(tmp_path):
    """A sample not present in any input at a site is ./. (UNKNOWN), never
    silently 0/0 — combine does not assert hom-reference where it has no
    record."""
    out = tmp_path / "combined.vcf"
    combine_vcfs([A, B], out)

    _, sites = _read(out)
    assert sites[100]["gt"]["S3"] == UNKNOWN  # S3 not in A, no B record at 100
    assert sites[400]["gt"]["S1"] == UNKNOWN  # S1 not in B, no A record at 400


def test_min_callsets_keeps_only_consensus(tmp_path):
    """--min-callsets 2 keeps only sites present in both inputs."""
    out = tmp_path / "consensus.vcf"
    combine_vcfs([A, B], out, min_callsets=2)

    _, sites = _read(out)
    assert sorted(sites) == [200, 300]  # 100 (A-only) and 400 (B-only) dropped


def test_custom_set_names(tmp_path):
    out = tmp_path / "combined.vcf"
    combine_vcfs([A, B], out, names=["gatk", "dv"])

    _, sites = _read(out)
    assert sites[100]["set"] == "gatk"
    assert sites[400]["set"] == "dv"


def test_rejects_min_callsets_out_of_range(tmp_path):
    out = tmp_path / "x.vcf"
    with pytest.raises(CombineError, match="min-callsets"):
        combine_vcfs([A, B], out, min_callsets=3)


def test_requires_two_inputs(tmp_path):
    out = tmp_path / "x.vcf"
    with pytest.raises(CombineError, match="two"):
        combine_vcfs([A], out)


def test_rejects_multiallelic_input(tmp_path):
    """Like every vcfclick input, call sets must be decomposed. A
    multi-allelic record errors with a decompose hint, not a silent
    wrong merge."""
    multi = tmp_path / "multi.vcf"
    multi.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tA\tT,G\t.\tPASS\t.\tGT\t1/2\n"
    )
    out = tmp_path / "x.vcf"
    with pytest.raises(CombineError, match="multi-allelic"):
        combine_vcfs([A, multi], out)


def test_gzip_output(tmp_path):
    """A .gz output path is written gzip-compressed and reads back."""
    out = tmp_path / "combined.vcf.gz"
    combine_vcfs([A, B], out)
    _, sites = _read(out)
    assert sorted(sites) == [100, 200, 300, 400]


def test_cli_combine_smoke(tmp_path):
    out = tmp_path / "combined.vcf"
    r = subprocess.run(
        [VCFCLICK_BIN, "combine", str(A), str(B), "-o", str(out)],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, f"combine failed:\n{r.stderr}"
    assert out.exists()
    _, sites = _read(out)
    assert sites[200]["set"] == "Intersection"
