"""Equivalence test: `vcfclick combine` vs GATK3 `CombineVariants`.

vcfclick combine is a reimplementation of GATK3's CombineVariants (the
tool GATK4 dropped). This test pins that claim to *real* GATK3 output.

`tests/fixtures/gatk3/combine_priority.gatk3.vcf` is the actual output of
GATK 3.8-1-0 `CombineVariants -genotypeMergeOptions PRIORITIZE -priority
first,second` on the two checked-in call sets (the generating command is
recorded in the file's `##vcfclick_golden_cmd` header). GATK3 needs Java
8 + a reference and is not runnable in CI, so the golden output is frozen
here; this test asserts vcfclick reproduces it.

Scope of the comparison: the CombineVariants *semantics* — the site
union, the `set=` provenance value, and the per-sample genotypes
(including PRIORITIZE conflict resolution for the shared sample S2).
Not compared: GATK3 additionally recomputes AC/AF/AN INFO and uses its
own header text; vcfclick emits `set=` + passthrough FORMAT only. Those
are representation differences, not CombineVariants behavior.

To regenerate the golden file (needs Java 8 + the GATK 3.8 jar):

    java -jar GenomeAnalysisTK.jar -T CombineVariants -R ref.fasta \\
        -V:first first.vcf.gz -V:second second.vcf.gz \\
        -o combine_priority.gatk3.vcf \\
        -genotypeMergeOptions PRIORITIZE -priority first,second
"""

from __future__ import annotations

from pathlib import Path

from cyvcf2 import VCF

from ingest.combine import combine_vcfs

FIXTURES = Path(__file__).parent / "fixtures" / "gatk3"
FIRST = FIXTURES / "first.vcf"
SECOND = FIXTURES / "second.vcf"
GOLDEN = FIXTURES / "combine_priority.gatk3.vcf"


def _semantics(path: Path):
    """Reduce a VCF to the CombineVariants-comparable core:
    {(chrom,pos): (set_value, {sample: gt_type})}."""
    v = VCF(str(path))
    samples = list(v.samples)
    out = {}
    for rec in v:
        out[(rec.CHROM, rec.POS)] = (
            rec.INFO.get("set"),
            dict(zip(samples, (int(t) for t in rec.gt_types))),
        )
    return samples, out


def test_combine_matches_gatk3_combinevariants(tmp_path):
    """vcfclick combine reproduces real GATK 3.8 CombineVariants output:
    same sample set, same sites, same set= provenance, same genotypes
    (PRIORITIZE picks the first input's call for the shared sample)."""
    out = tmp_path / "vcfclick.vcf"
    combine_vcfs([FIRST, SECOND], out, names=["first", "second"])

    vc_samples, vc = _semantics(out)
    gatk_samples, gatk = _semantics(GOLDEN)

    assert vc_samples == gatk_samples == ["S1", "S2", "S3"]
    assert vc == gatk


def test_gatk3_golden_pins_prioritize_behaviour():
    """Guard the golden itself: the shared sample S2 is 0/1 in `first` and
    1/1 in `second` at chr1:202; PRIORITIZE (priority first,second) must
    keep first's 0/1. If this ever flips, the golden was regenerated
    wrong and the equivalence test above would be meaningless."""
    _, gatk = _semantics(GOLDEN)
    # cyvcf2 gt_types: 1 = HET (0/1), 3 = HOM_ALT (1/1)
    assert gatk[("chr1", 202)][0] == "Intersection"
    assert gatk[("chr1", 202)][1]["S2"] == 1  # HET, from `first` (not 1/1)
