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


def test_contig_order_follows_header_not_first_seen(tmp_path):
    """The higher-priority input (contig_x) has only a chr2 record, but
    its header declares chr1 before chr2. Output must sort chr1 before
    chr2 (reference-dictionary order), not chr2-first because that's the
    first variant seen."""
    cx = FIXTURES / "contig_x.vcf"
    cy = FIXTURES / "contig_y.vcf"
    out = tmp_path / "combined.vcf"
    combine_vcfs([cx, cy], out)

    v = VCF(str(out))
    order = [(rec.CHROM, rec.POS) for rec in v]
    assert order == [("chr1", 50), ("chr2", 50)]


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="htslib bgzip/tabix not installed",
)
def test_gzip_output_is_bgzf_and_indexed(tmp_path):
    """A .gz output is BGZF (not plain gzip) and tabix-indexed, so it is
    usable by region-parallel ingest like every other vcfclick VCF."""
    out = tmp_path / "combined.vcf.gz"
    combine_vcfs([A, B], out)

    # BGZF magic: gzip header with the FLG.FEXTRA bit set (0x04).
    with open(out, "rb") as fh:
        assert fh.read(4) == b"\x1f\x8b\x08\x04"
    assert out.with_suffix(out.suffix + ".tbi").exists()

    _, sites = _read(out)
    assert sorted(sites) == [100, 200, 300, 400]


QA = FIXTURES / "qual_a.vcf"
QB = FIXTURES / "qual_b.vcf"


def _read_qual(path: Path):
    """Return (samples, {pos: {sample: {'GQ':int|None,'DP':int|None,
    'AD':(int,int)|None}}}) plus the FORMAT keys present."""
    v = VCF(str(path))
    samples = list(v.samples)
    out, fmt = {}, set()
    for rec in v:
        fmt.update(rec.FORMAT)
        gq, dp, ad = rec.format("GQ"), rec.format("DP"), rec.format("AD")
        per = {}
        for i, s in enumerate(samples):
            per[s] = {
                "GQ": None if gq is None or gq[i][0] < 0 else int(gq[i][0]),
                "DP": None if dp is None or dp[i][0] < 0 else int(dp[i][0]),
                "AD": None
                if ad is None or ad[i][0] < 0
                else tuple(int(x) for x in ad[i]),
            }
        out[rec.POS] = per
    return samples, out, fmt


def test_format_passthrough_from_priority_source(tmp_path):
    """At chr1:200 sample S2 is called in both inputs; A (higher priority)
    wins the genotype, so its GQ/DP/AD must travel too — not B's. Quality
    must describe the genotype that was actually kept."""
    out = tmp_path / "combined.vcf"
    combine_vcfs([QA, QB], out)

    _, q, fmt = _read_qual(out)
    assert {"GQ", "DP", "AD"} <= fmt
    assert q[200]["S2"] == {"GQ": 55, "DP": 22, "AD": (13, 9)}  # A's, not B's


def test_format_passthrough_falls_through_with_genotype(tmp_path):
    """At chr1:100 the higher-priority input A has S2 as ./. (no call), so
    both the genotype and its quality come from B — they stay consistent."""
    out = tmp_path / "combined.vcf"
    combine_vcfs([QA, QB], out)

    _, q, _ = _read_qual(out)
    assert q[100]["S2"] == {"GQ": 45, "DP": 25, "AD": (15, 10)}  # B's


def test_absent_sample_format_is_all_missing(tmp_path):
    """chr1:300 exists only in input A (samples S1, S2). S3 has no record
    there, so its cell is ./. with '.' for every passthrough field — never
    a borrowed quality value."""
    out = tmp_path / "combined.vcf"
    combine_vcfs([QA, QB], out)

    v = VCF(str(out))
    s3 = v.samples.index("S3")
    rec = next(r for r in v if r.POS == 300)
    assert rec.gt_types[s3] == UNKNOWN  # ./.
    assert rec.format("GQ")[s3][0] < 0  # missing, not borrowed from A
    assert rec.format("AD")[s3][0] < 0


def test_ad_with_wrong_cardinality_is_dropped(tmp_path):
    """Output records are biallelic, so AD must be exactly ref,alt. A
    single-ALT input record whose AD still carries three values (an
    improperly decomposed call set) must not write AD=5,3,2 — that depth
    references an allele not in the output. The cell drops AD to '.'
    rather than emit an uninterpretable Number=R value."""
    hdr = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    a = tmp_path / "a.vcf"
    b = tmp_path / "b.vcf"
    a.write_text(hdr + "chr1\t1\t.\tA\tC\t.\t.\t.\tGT:AD\t0/1:5,3,2\n")
    b.write_text(
        hdr.replace("S1", "S2") + "chr1\t1\t.\tA\tC\t.\t.\t.\tGT:AD\t0/1:9,4\n"
    )
    out = tmp_path / "out.vcf"
    combine_vcfs([a, b], out)

    v = VCF(str(out))
    rec = next(iter(v))
    s1, s2 = v.samples.index("S1"), v.samples.index("S2")
    assert rec.format("AD")[s1][0] < 0  # 3-value AD dropped → missing
    assert tuple(int(x) for x in rec.format("AD")[s2]) == (9, 4)  # valid 2-value kept


def test_gt_only_inputs_produce_gt_only_format(tmp_path):
    """Inputs without GQ/DP/AD yield a GT-only FORMAT — passthrough fields
    are listed only when some input actually carries them."""
    out = tmp_path / "combined.vcf"
    combine_vcfs([A, B], out)  # the GT-only fixtures

    v = VCF(str(out))
    rec = next(iter(v))
    assert rec.FORMAT == ["GT"]


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
