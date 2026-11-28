"""Generate the sibling-segregation fixtures (deterministic).

    python tests/fixtures/make_segregation_vcf.py

Writes, next to this file:

  segregation.vcf            planted sites + background, 7 samples, hand-set INFO/AF
  segregation_cohortaf.vcf   the same records with INFO/AF rewritten to the
                             cohort's own ac/an (the GATK-style case)
  segregation.ped            one family: FA MO, S1 S2 affected, S3 S4
                             unaffected, S5 unknown; parents phenotype 0
  segregation_genes.gff3     the five comphet genes G1..G5 on chr1

Sample order FA MO S1 S2 S3 S4 S5, FORMAT GT:GQ:DP:AD, every call explicit
(including 0/0, so a sparse ingest and a --keep-reference ingest see the same
truth). A planted site is written FA MO S1 S2 S3 S4 S5; a value is a GT, or a
(GT, GQ, DP, "ref,alt") tuple for an override. The background sites carry AF
0.30 (the default --max-af drops them) and Random(20261020).

The expectation for every planted site is in EXPECTED: the models in which it
is `consistent` (or `consistent_inferred`, `compatible`, ...) with the default
options on a --keep-reference ingest. tests/test_segregate.py checks the CLI
against it and against a pure-Python reading of the VCF text.
"""

from __future__ import annotations

import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SAMPLES = ["FA", "MO", "S1", "S2", "S3", "S4", "S5"]
PLANTED_AF = 0.001
BACKGROUND_AF = 0.30
GENES = [
    ("G1", 5_000_000, 5_001_000),
    ("G2", 5_100_000, 5_101_000),
    ("G3", 5_200_000, 5_201_000),
    ("G4", 5_300_000, 5_301_000),
    ("G5", 5_400_000, 5_401_000),
]

R, H, A, N = "0/0", "0/1", "1/1", "./."


def lowq(gt: str, gq: int = 5):
    return (gt, gq, 30, None)


def ab(gt: str, ref: int, alt: int):
    return (gt, 50, ref + alt, f"{ref},{alt}")


# name -> (pos, [FA, MO, S1, S2, S3, S4, S5], options)
SITES: dict[str, tuple] = {
    # dominant (inherited)
    "dom_pos": (1_000_000, [H, R, H, H, R, R, H], {}),
    "dom_mat": (1_000_800, [R, H, H, H, R, R, R], {}),
    "dom_unaff_carrier": (1_000_100, [H, R, H, H, H, R, R], {}),
    "dom_aff_missing": (1_000_200, [H, R, H, R, R, R, R], {}),
    "dom_both_parents": (1_000_300, [H, H, H, H, R, R, R], {}),
    "dom_unaff_nocall": (1_000_400, [H, R, H, H, R, N, R], {}),
    "dom_unaff_lowq_ref": (1_000_500, [H, R, H, H, R, lowq(R), R], {}),
    "dom_aff_lowab": (1_000_600, [H, R, H, ab(H, 36, 4), R, R, R], {}),
    "dom_common": (1_000_700, [H, R, H, H, R, R, R], {"af": 0.2}),
    "dom_two_unaff": (1_001_000, [H, R, H, H, H, H, R], {}),
    "dom_unaff_lowab": (1_001_100, [H, R, H, H, R, ab(H, 36, 4), R], {}),
    # a thin (DP 4) call is support evidence, never proof: a DP-4 sibling's 0/0
    # is not reference and a DP-4 affected het is not a carrier
    "dom_unaff_lowdp": (1_001_200, [H, R, H, H, R, ("0/0", 50, 4, "4,0"), R], {}),
    "dom_aff_lowdp": (1_001_300, [H, R, H, ("0/1", 50, 4, "2,2"), R, R, R], {}),
    # the candidate must not depend on which affected child is first
    "anc_nocall_d": (1_100_000, [H, R, N, H, R, R, R], {}),
    "anc_lowab_d": (1_100_100, [H, R, ab(H, 36, 4), H, R, R, R], {}),
    "anc_lowgq_d": (1_100_200, [H, R, lowq(H), H, R, R, R], {}),
    "anc_ref_d": (1_100_300, [H, R, R, H, R, R, R], {}),
    "anc_nocall_dm": (1_100_400, [H, R, H, N, R, R, R], {}),
    "dom_unaff_dotone": (1_100_500, [H, R, H, H, "./1", R, R], {}),
    "dom_homalt_parent": (1_100_700, [A, R, H, H, R, R, R], {}),
    # de novo
    "dn_pos": (2_000_000, [R, R, H, H, R, R, R], {}),
    "dn_unaff_carrier": (2_000_100, [R, R, H, H, H, R, R], {}),
    "dn_one_aff": (2_000_200, [R, R, H, R, R, R, R], {}),
    "dn_parent_nocall": (2_000_400, [N, R, H, H, R, R, R], {}),
    "dn_parent_lowq": (2_000_500, [lowq(R), R, H, H, R, R, R], {}),
    "anc_nocall_dn": (2_100_000, [R, R, N, H, R, R, R], {}),
    # a partial call in ANOTHER member: the stored parental 0/0 cannot be trusted
    "dn_parent_half": (2_100_100, [R, R, H, H, R, R, "0/."], {}),
    # a partial call in a parent has no row at all (ingest), so it drops the site
    "dn_parent_partial": (2_100_150, ["0/.", R, H, H, R, R, R], {}),
    "dn_parent_altreads": (2_100_200, [("0/0", 99, 30, "27,3"), R, H, H, R, R, R], {}),
    "dn_unaff_altreads": (2_100_300, [R, R, H, H, ("0/0", 99, 30, "27,3"), R, R], {}),
    # recessive
    "rec_pos": (3_000_000, [H, H, A, A, R, R, R], {}),
    "rec_unaff_het_pair": (3_000_500, [H, H, A, A, H, H, R], {}),
    "rec_unaff_hom": (3_000_100, [H, H, A, A, A, R, R], {}),
    "rec_aff_het": (3_000_200, [H, H, A, H, R, R, R], {}),
    "rec_parent_ref": (3_000_300, [H, R, A, A, R, R, R], {}),
    "rec_unaff_nocall": (3_000_400, [H, H, A, A, R, N, R], {}),
    "anc_nocall_r": (3_100_000, [H, H, N, A, R, R, R], {}),
    "rec_aff_homalt_lowab": (3_100_100, [H, H, A, ab(A, 20, 20), R, R, R], {}),
    "rec_parent_lowab": (3_100_200, [ab(H, 36, 4), H, A, A, R, R, R], {}),
    # compound het: genes G1..G5, a paternal hit and a maternal hit each
    "g1_pat": (5_000_100, [H, R, H, H, R, H, R], {}),
    "g1_mat": (5_000_900, [R, H, H, H, H, R, R], {}),
    "g2_pat": (5_100_100, [H, R, H, H, H, R, R], {}),
    "g2_mat": (5_100_900, [R, H, H, H, H, R, R], {}),
    "g3_pat": (5_200_100, [H, R, H, R, R, R, R], {}),
    "g3_mat": (5_200_900, [R, H, H, H, R, R, R], {}),
    "g4_pat": (5_300_100, [H, R, H, H, R, N, R], {}),
    "g4_mat": (5_300_900, [R, H, H, H, R, H, R], {}),
    "g5_pat": (5_400_100, [H, R, H, H, R, R, R], {}),
    "g5_mat": (5_400_900, [H, R, H, H, R, R, R], {}),
}

# Sex and mitochondrial contigs: skipped, counted.
EXTRA_CONTIG_SITES = [
    ("chrX", 10_000_000, [R, R, H, H, R, R, R]),
    ("chrM", 100, [R, R, H, H, R, R, R]),
    ("chrY", 5_000_000, [R, R, H, H, R, R, R]),
]

# name -> {model: status} with the default options on a --keep-reference
# ingest; a model that does not list the site does not report it (or reports
# it under the options noted in tests/test_segregate.py).
EXPECTED: dict[str, dict[str, str]] = {
    "dom_pos": {"dominant": "consistent"},
    "dom_mat": {"dominant": "consistent"},
    "dom_unaff_nocall": {"dominant": "compatible"},
    "dom_unaff_lowq_ref": {"dominant": "compatible"},
    "dom_aff_lowab": {"dominant": "compatible"},
    "dom_unaff_lowab": {"dominant": "compatible"},
    "dom_unaff_lowdp": {"dominant": "compatible"},
    "dom_aff_lowdp": {"dominant": "compatible"},
    "anc_nocall_d": {"dominant": "compatible"},
    "anc_lowab_d": {"dominant": "compatible"},
    "anc_lowgq_d": {"dominant": "compatible"},
    "anc_nocall_dm": {"dominant": "compatible"},
    "dom_unaff_dotone": {"dominant": "compatible"},
    "dom_homalt_parent": {"dominant": "consistent"},
    "dn_pos": {"denovo": "consistent"},
    "anc_nocall_dn": {"denovo": "compatible"},
    "dn_parent_half": {"denovo": "compatible"},
    "rec_pos": {"recessive": "consistent"},
    "rec_unaff_het_pair": {"recessive": "consistent"},
    "rec_unaff_nocall": {"recessive": "compatible"},
    "anc_nocall_r": {"recessive": "compatible"},
    "rec_aff_homalt_lowab": {"recessive": "compatible"},
    "rec_parent_lowab": {"recessive": "consistent"},
}
# sites that must NOT appear in any model at the default options
ABSENT_EVERYWHERE = [
    "dom_unaff_carrier",
    "dom_aff_missing",
    "dom_both_parents",
    "dom_common",
    "dom_two_unaff",
    "anc_ref_d",
    "dn_unaff_carrier",
    "dn_one_aff",
    "dn_parent_nocall",
    "dn_parent_lowq",
    "dn_parent_altreads",
    "dn_unaff_altreads",
    "dn_parent_partial",
    "rec_unaff_hom",
    "rec_aff_het",
    "rec_parent_ref",
]
# comphet gene -> status of its (paternal, maternal) pair, None: no pair
COMPHET = {
    "G1": "consistent",
    "G2": None,
    "G3": None,
    "G4": "compatible",
    "G5": None,
}


def _call(value, default_ad: bool = True) -> str:
    if isinstance(value, tuple):
        gt, gq, dp, ad = value
        ad = ad or _default_ad(gt, dp)
        return f"{gt}:{gq}:{dp}:{ad}"
    gt = value
    if gt == N:
        return "./.:.:.:."
    return f"{gt}:50:30:{_default_ad(gt, 30)}"


def _default_ad(gt: str, dp: int) -> str:
    return {
        "0/0": f"{dp},0",
        "0/.": f"{dp},0",
        "0/1": f"{dp // 2},{dp - dp // 2}",
        "./1": f"{dp // 2},{dp - dp // 2}",
        "1/1": f"0,{dp}",
    }.get(gt, ".")


def _gt(value) -> str:
    return value[0] if isinstance(value, tuple) else value


def _alleles(value) -> tuple[int, int]:
    """(alt alleles, called alleles) of a fully called GT, (0, 0) otherwise."""
    gt = _gt(value)
    if "." in gt:
        return 0, 0
    parts = gt.split("/")
    return sum(int(p) for p in parts), len(parts)


def background(rng: random.Random, n: int = 90) -> list[tuple[int, list]]:
    """Random sites at 7,000,000 + 1000 i: genotypes drawn independently, so
    they land in every pattern at their natural rate."""
    out = []
    for i in range(n):
        gts = [rng.choice([R, R, R, H, H, A, N]) for _ in SAMPLES]
        out.append((7_000_000 + 1000 * i, gts))
    return out


def records() -> list[tuple[str, int, list, float | None]]:
    """(chrom, pos, values, af or None for the cohort's own) sorted by contig and
    position: chr1 sites first, then the extra contigs."""
    rng = random.Random(20261020)
    rows: list[tuple[str, int, list, float | None]] = []
    for pos, gts, opts in SITES.values():
        rows.append(("chr1", pos, gts, opts.get("af", PLANTED_AF)))
    for pos, gts in background(rng):
        rows.append(("chr1", pos, gts, BACKGROUND_AF))
    rows.sort(key=lambda r: r[1])
    for chrom, pos, gts in EXTRA_CONTIG_SITES:
        rows.append((chrom, pos, gts, PLANTED_AF))
    return rows


def cohort_af(values: list) -> float:
    alt = called = 0
    for v in values:
        a, c = _alleles(v)
        alt += a
        called += c
    return round(alt / called, 4) if called else 0.0


def write_vcf(path: Path, *, cohort_derived: bool) -> None:
    contigs = ["chr1", "chrM", "chrX", "chrY"]
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={c}>" for c in contigs]
    lines += [
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES),
    ]
    for chrom, pos, gts, af in records():
        value = cohort_af(gts) if cohort_derived else af
        lines.append(
            "\t".join(
                [
                    chrom,
                    str(pos),
                    ".",
                    "A",
                    "G",
                    "50",
                    "PASS",
                    f"AF={value}",
                    "GT:GQ:DP:AD",
                ]
                + [_call(g) for g in gts]
            )
        )
    path.write_text("\n".join(lines) + "\n")


PED = """\
#family individual father mother sex phenotype
fam1 FA 0 0 1 0
fam1 MO 0 0 2 0
fam1 S1 FA MO 1 2
fam1 S2 FA MO 2 2
fam1 S3 FA MO 1 1
fam1 S4 FA MO 2 1
fam1 S5 FA MO 1 0
"""


def gff() -> str:
    lines = ["##gff-version 3"]
    for symbol, start, end in GENES:
        gid = f"ENSG0000090{symbol[1:]}.1"
        attrs = f"ID={gid};gene_id={gid};gene_type=protein_coding;gene_name={symbol}"
        lines.append(
            "\t".join(
                ["chr1", "HAVANA", "gene", str(start), str(end), ".", "+", ".", attrs]
            )
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    write_vcf(HERE / "segregation.vcf", cohort_derived=False)
    write_vcf(HERE / "segregation_cohortaf.vcf", cohort_derived=True)
    (HERE / "segregation.ped").write_text(PED)
    (HERE / "segregation_genes.gff3").write_text(gff())


if __name__ == "__main__":
    main()
