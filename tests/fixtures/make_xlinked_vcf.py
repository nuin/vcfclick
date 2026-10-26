"""Generate the X-linked fixtures (deterministic).

    python tests/fixtures/make_xlinked_vcf.py

Writes, next to this file:

  xlinked.vcf               162 GRCh38 chrX sites, 10 samples
  xlinked.ped               the three families, sexes and phenotypes declared
  xlinked_nosex.ped         the same with every sex 0 (sex is inferred)
  xlinked_conflict.ped      XSON declared female (a sex_conflict warning)
  xlinked_badfather.ped     XDAD declared female (parent_sex_conflict)
  xlinked_affected.ped      XDAD and XMUM affected (female_het -> dominant)

Samples, FORMAT GT:GQ:DP:AD, GQ 50, males DP 15, females DP 30:

  XF1  XDAD XMUM XSON XDAU  males haploid-coded on non-PAR X ('1', '0')
  XF2  YDAD YMUM YSON       the same alleles, males diploid-coded ('1/1',
                            '0/0'): the stored genotypes must be identical
  XF3  ZDAD ZMUM ZKID       one non-reference X call in ZKID (sex 0): its sex
                            cannot be inferred

A planted site is written `[father|mother|son|daughter]`: a male allele is
0, 1, H (a het artefact, 0/1), P (a partial call ./1) or '.' (no call); a
female value is the alt-allele count 0, 1 or 2. Planted sites carry INFO/AF
0.001 (R1, S1 0.005), so --max-af 0.01 keeps them; the 80 background sites
carry 0.30 and only feed sex inference and the excluded counters. Their
alleles come from Random(20261020) and from nothing else, so the planted
sites never move.

A site name, its expectation and the reasoning are in SITES below; the
tests (tests/test_xlinked.py) compare the CLI with a pure-Python reading of
this file and with the expectations in EXPECTED.
"""

from __future__ import annotations

import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SAMPLES = [
    "XDAD", "XMUM", "XSON", "XDAU",
    "YDAD", "YMUM", "YSON",
    "ZDAD", "ZMUM", "ZKID",
]  # fmt: skip

PLANTED_AF = "0.001"

# name -> (pos, father, mother, son, daughter, options). Males on non-PAR X
# unless options["diploid"]: then every value is an alt-allele count.
# Options: af, gq/dp/ad overrides per role ("m_gq", "s_dp", "s_ad", ...).
SITES: dict[str, tuple] = {
    # A: hemizygous sons
    "A1": (10_000_001, 0, 1, 1, 0, {}),
    "A2": (10_100_001, 0, 1, 1, 1, {}),
    "A3": (10_200_001, ".", 1, 1, 0, {}),
    "A4": (10_300_001, 1, 1, 1, 1, {}),
    "A5": (10_400_001, 0, 0, 1, 0, {}),
    "A6": (10_500_001, 0, 1, "H", 0, {}),
    "A7": (10_600_001, 0, 2, 1, 1, {}),
    "A8": (10_700_001, 0, 1, 1, 0, {"m_gq": 5}),
    "A9": (10_800_001, 0, 1, 1, 0, {"s_dp": 6, "s_ad": "0,6"}),
    "A10": (10_900_001, 0, 1, 1, 0, {"s_ad": "10,20"}),
    "A11": (11_000_001, 0, 1, 1, 0, {"s_ad": "2,28"}),
    # B: hom-alt daughters
    "B1": (20_000_001, 1, 1, 1, 2, {}),
    "B2": (20_100_001, 0, 1, 1, 2, {}),
    "B3": (20_200_001, 1, 0, 1, 2, {}),
    "B4": (20_300_001, 1, 2, 1, 2, {}),
    # C: het daughters
    "C1": (30_000_001, 1, 0, 0, 1, {}),
    "C2": (30_100_001, 0, 1, 0, 1, {}),
    "C3": (30_200_001, 1, 0, 0, 1, {"d_ad": "27,3"}),
    "C4": (30_300_001, 1, 1, 0, 1, {}),
    # D: de novo
    "D1": (40_000_001, 0, 0, 1, 0, {}),
    "D2": (40_100_001, 0, 0, 0, 1, {}),
    "D3": (40_200_001, 1, 0, 1, 0, {}),
    "D4": (40_300_001, 0, 0, 0, 2, {}),
    "D5": (40_400_001, ".", 0, 1, 0, {}),
    # E/F/G: compound het over the seeded genes XCH1..XCH3
    "E1": (50_000_001, 1, 0, 0, 1, {}),
    "E2": (50_100_001, 0, 1, 0, 1, {}),
    "F1": (51_000_001, 0, 1, 0, 1, {}),
    "F2": (51_100_001, 0, 1, 0, 1, {}),
    "G1": (52_000_001, 1, 0, "H", 0, {}),
    "G2": (52_100_001, 0, 1, "H", 0, {}),
    # P: PAR and the edges (GRCh38 PAR1 ends 2,781,479, PAR2 starts
    # 155,701,383). P6 / P7 sit in the bands where the builds disagree.
    "P0": (1_000_001, 1, 1, 1, 0, {"diploid": True}),
    "P1": (1_100_001, 1, 1, 2, 0, {"diploid": True}),
    "P2": (2_781_479, 1, 1, 2, 0, {"diploid": True}),
    "P3": (2_781_480, 0, 1, 1, 0, {}),
    "P4": (155_701_381, 0, 1, 1, 0, {}),
    "P5": (155_701_384, 1, 1, 2, 0, {"diploid": True}),
    "P6": (2_750_000, 0, 1, 2, 0, {"diploid": True}),
    "P7": (154_950_000, 0, 1, 2, 0, {"diploid": True}),
    # R1 / S1: rarer than --max-af, commoner than --max-af-x
    "R1": (12_000_001, 0, 1, 1, 0, {"af": "0.005"}),
    "S1": (13_800_001, 1, 0, 0, 1, {"af": "0.005"}),
    # L: what the classic models get wrong on X
    "L1": (13_000_001, "H", 0, 0, 1, {}),
    "L2": (13_100_001, "H", 1, 1, 2, {}),
    "L3": (13_200_001, 0, 0, "H", 0, {}),
    "L4": (13_300_001, 0, 1, 1, 0, {"m_ad": "27,3"}),
    "L5": (13_400_001, 0, 0, 1, 0, {"m_ad": "24,6"}),
    "L6": (13_500_001, 0, 1, "P", 0, {"s_ad": "8,7"}),
}
# T: every father {0,1} x mother {0,1,2} x son {0,1} x daughter {0,1,2}
for i, (f, m, s, d) in enumerate(
    (f, m, s, d) for f in (0, 1) for m in (0, 1, 2) for s in (0, 1) for d in (0, 1, 2)
):
    SITES[f"T{i:02d}"] = (100_000_001 + 1000 * i, f, m, s, d, {})

ZKID_SITE = "A1"  # ZMUM and ZKID het there; every other Z call is 0/0

# Expected counts (keep-reference, defaults, xlinked.ped), derived from the
# pattern table by hand in docs/TRIO.md terms; the tests also derive them.
EXPECTED = {
    "XSON": {"count": 29, "recessive": 18, "denovo": 11},
    "XDAU": {"count": 28, "recessive": 3, "female_het": 14, "denovo": 11},
}

MALE_DP, FEMALE_DP, GQ = 15, 30, 50


def _male(value, *, haploid: bool, ad: str | None) -> tuple[str, str]:
    """(GT, AD) of a male call."""
    if value == ".":
        return ("." if haploid else "./."), "."
    if value == "H":
        return "0/1", ad or "8,7"
    if value == "P":
        return "./1", ad or "8,7"
    if haploid:
        return ("1", ad or "0,15") if value == 1 else ("0", ad or "15,0")
    gt = {0: "0/0", 1: "0/1", 2: "1/1"}[value]
    return gt, ad or {0: "15,0", 1: "8,7", 2: "0,15"}[value]


def _female(value: int, ad: str | None) -> tuple[str, str]:
    gt = {0: "0/0", 1: "0/1", 2: "1/1"}[value]
    return gt, ad or {0: "30,0", 1: "15,15", 2: "0,30"}[value]


def _cell(gt: str, ad: str, gq: int, dp: int) -> str:
    return f"{gt}:{gq}:{dp}:{ad}"


def _row(name: str, pos: int, f, m, s, d, opt: dict, af: str) -> list[str]:
    # males diploid in both families at a PAR site (or a GRCh37 PAR2 call)
    haploid = not opt.get("diploid")

    def male(role: str, value, haploid: bool) -> str:
        gt, ad = _male(value, haploid=haploid, ad=opt.get(f"{role}_ad"))
        dp = opt.get(f"{role}_dp", MALE_DP)
        return _cell(gt, ad, opt.get(f"{role}_gq", GQ), dp)

    def female(role: str, value) -> str:
        gt, ad = _female(value, opt.get(f"{role}_ad"))
        dp = opt.get(f"{role}_dp", FEMALE_DP)
        return _cell(gt, ad, opt.get(f"{role}_gq", GQ), dp)

    f_x, s_x = male("f", f, haploid), male("s", s, haploid)

    # the diploid-coded copy: a hemizygous ALT is written 1/1
    def diploid(v):
        return 2 if haploid and v == 1 else v

    f_y, s_y = male("f", diploid(f), False), male("s", diploid(s), False)
    m_c, d_c = female("m", m), female("d", d)
    cells = [f_x, m_c, s_x, d_c, f_y, m_c, s_y]
    # the Z family: reference everywhere but one site
    z_dad = _cell("0" if haploid else "0/0", "15,0", GQ, MALE_DP)
    z_ref = _cell("0/0", "30,0", GQ, FEMALE_DP)
    z_het = _cell("0/1", "15,15", GQ, FEMALE_DP)
    if name == ZKID_SITE:
        cells += [z_dad, z_het, z_het]
    else:
        cells += [z_dad, z_ref, z_ref]
    return ["chrX", str(pos), name, *_alleles(pos), "50", "PASS", f"AF={af}",
            "GT:GQ:DP:AD", *cells]  # fmt: skip


def _alleles(pos: int) -> tuple[str, str]:
    # planted sites alternate A>G / C>T by position
    return ("A", "G") if (pos // 100_000) % 2 == 0 else ("C", "T")


def background() -> list[tuple]:
    """80 common, Mendelian sites (father, mother, son, daughter alleles)."""
    rng = random.Random(20261020)
    out = []
    for i in range(80):
        pos = 60_000_001 + 100_000 * i
        mat = [rng.random() < 0.3, rng.random() < 0.3]
        f = int(rng.random() < 0.3)
        son = int(mat[rng.randrange(2)])
        dau = f + int(mat[rng.randrange(2)])
        out.append((f"H{i:02d}", pos, f, sum(mat), son, dau))
    return out


def records() -> list[list[str]]:
    rows = []
    for name, (pos, f, m, s, d, opt) in SITES.items():
        rows.append((pos, _row(name, pos, f, m, s, d, opt, opt.get("af", PLANTED_AF))))
    for name, pos, f, m, s, d in background():
        rows.append((pos, _row(name, pos, f, m, s, d, {}, "0.30")))
    return [r for _, r in sorted(rows, key=lambda x: x[0])]


HEADER = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chrX,length=156040895,assembly=GRCh38>",
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
]


def vcf_text(rows: list[list[str]] | None = None) -> str:
    rows = records() if rows is None else rows
    head = "\t".join(
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        + SAMPLES
    )
    return "\n".join([*HEADER, head, *("\t".join(r) for r in rows)]) + "\n"


# family, sample, father, mother, sex, phenotype
PED = [
    ("XF1", "XDAD", "0", "0", "1", "1"),
    ("XF1", "XMUM", "0", "0", "2", "1"),
    ("XF1", "XSON", "XDAD", "XMUM", "1", "2"),
    ("XF1", "XDAU", "XDAD", "XMUM", "2", "2"),
    ("XF2", "YDAD", "0", "0", "1", "1"),
    ("XF2", "YMUM", "0", "0", "2", "1"),
    ("XF2", "YSON", "YDAD", "YMUM", "1", "2"),
    ("XF3", "ZDAD", "0", "0", "1", "1"),
    ("XF3", "ZMUM", "0", "0", "2", "1"),
    ("XF3", "ZKID", "ZDAD", "ZMUM", "0", "2"),
]


def ped_text(sex: dict | None = None, pheno: dict | None = None, *, nosex=False) -> str:
    lines = []
    for fam, s, f, m, sx, ph in PED:
        if nosex:
            sx = "0"
        sx = (sex or {}).get(s, sx)
        ph = (pheno or {}).get(s, ph)
        lines.append("\t".join((fam, s, f, m, sx, ph)))
    return "\n".join(lines) + "\n"


PEDS = {
    "xlinked.ped": {},
    "xlinked_nosex.ped": {"nosex": True},
    "xlinked_conflict.ped": {"sex": {"XSON": "2"}},
    "xlinked_badfather.ped": {"sex": {"XDAD": "2"}},
    "xlinked_affected.ped": {"pheno": {"XDAD": "2", "XMUM": "2"}},
}


def main(out: Path = HERE) -> None:
    (out / "xlinked.vcf").write_text(vcf_text())
    for name, kw in PEDS.items():
        (out / name).write_text(ped_text(**kw))


if __name__ == "__main__":
    main()
