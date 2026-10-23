"""A minimal VCF writer for tests that build their own cohorts.

Shared so fixture generators stop re-implementing it. It is a plain text
writer (no pysam/cyvcf2 write path), deterministic, and has no dependency on
make_related_vcf.py, which this module must never import or change: that
generator and its committed outputs are pinned by existing tests.

    from tests.fixtures._vcfgen import Site, write_vcf, trio_sites

A `Site` is one decomposed (single-ALT) record. `gts` holds one GT string
per sample, optionally with FORMAT extras: pass `fmt="GT:AD"` and
`gts=["0/1:15,15", ...]`.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class Site:
    chrom: str
    pos: int
    ref: str
    alt: str
    gts: tuple[str, ...]
    info: str = "."
    fmt: str = "GT"
    qual: str = "50"
    extra: dict = field(default_factory=dict, compare=False)
    filter: str = "PASS"


def write_vcf(
    path: Path,
    samples: list[str],
    sites: list[Site],
    *,
    contigs: tuple[str, ...] = ("chr1",),
    info_header: tuple[str, ...] = (
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    ),
    filter_header: tuple[str, ...] = (),
) -> Path:
    """Write `sites` (already sorted) as an uncompressed VCF 4.2 at `path`."""
    fmt_header = (
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    )
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={c}>" for c in contigs]
    lines += list(filter_header) + list(info_header) + list(fmt_header)
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"])
                 + "\tINFO\tFORMAT\t" + "\t".join(samples))  # fmt: skip
    for s in sites:
        if len(s.gts) != len(samples):
            raise ValueError(
                f"{s.chrom}:{s.pos} has {len(s.gts)} GTs for {len(samples)} samples"
            )
        lines.append(
            "\t".join(
                [
                    s.chrom,
                    str(s.pos),
                    ".",
                    s.ref,
                    s.alt,
                    s.qual,
                    s.filter,
                    s.info,
                    s.fmt,
                    *s.gts,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n")
    return path


def trio_sites(
    n: int,
    *,
    seed: int = 7,
    af: str = "cohort",
    start: int = 1000,
    step: int = 1000,
    chrom: str = "chr1",
) -> list[Site]:
    """`n` private-het sites in a CHILD / FATHER / MOTHER trio.

    Site i alternates between a de novo het (child 0/1, parents 0/0) and a
    paternally inherited het (child 0/1, father 0/1, mother 0/0). Every
    call carries GQ 99 / DP 30 and balanced AD, so the quality gates pass.

    `af` sets INFO/AF:
      "cohort"  the callset's own frequency ac/an (1/6 for a de novo het,
                2/6 for an inherited one), what a GATK joint-called family
                VCF carries;
      "popular" a population-style frequency (0.0005 de novo, 0.002
                inherited), independent of the family;
      "none"    no INFO/AF.
    """
    rng = random.Random(seed)
    sites = []
    for i in range(n):
        denovo = i % 2 == 0
        if denovo:
            gts = ("0/1:15,15:99:30", "0/0:30,0:99:30", "0/0:30,0:99:30")
            cohort_af = 1 / 6
            pop_af = 0.0005
        else:
            gts = ("0/1:15,15:99:30", "0/1:15,15:99:30", "0/0:30,0:99:30")
            cohort_af = 2 / 6
            pop_af = 0.002
        if af == "cohort":
            info = f"AF={cohort_af:.6f}"
        elif af == "popular":
            # jitter so no two sites share a value
            info = f"AF={pop_af * (1 + rng.random() / 10):.6f}"
        elif af == "none":
            info = "."
        else:
            raise ValueError(af)
        sites.append(
            Site(
                chrom,
                start + i * step,
                "A",
                "G",
                gts,
                info=info,
                fmt="GT:AD:GQ:DP",
            )
        )
    return sites


def gatk_af(x: float) -> str:
    """`x` printed the way GATK (htsjdk) writes a Float: three significant
    digits in scientific notation below 0.01, three decimals from 0.01 to 1
    (so 1/6 reads 0.167 and 1/300 reads 3.333e-03)."""
    if x < 1:
        if x < 0.01:
            return f"{x:.3e}" if abs(x) >= 1e-20 else "0.00"
        return f"{x:.3f}"
    return f"{x:.2f}"


def cohort_samples(n_samples: int) -> list[str]:
    """CHILD, FATHER, MOTHER, then unrelated samples U004... up to
    `n_samples` in all."""
    if n_samples < 3:
        raise ValueError("a cohort that holds the trio has at least 3 samples")
    return TRIO_SAMPLES + [f"U{i:03d}" for i in range(4, n_samples + 1)]


def cohort_denovo_sites(
    n_samples: int,
    infos: list[str],
    *,
    start: int = 1000,
    step: int = 1000,
    chrom: str = "chr1",
) -> list[Site]:
    """One private de novo het in CHILD per entry of `infos` (that site's
    INFO string), in a cohort of `n_samples` where every other sample is
    0/0. Each site has ac_called = 1 and an_called = 2 * n_samples, so the
    cohort's own frequency is 1 / (2 * n_samples) at every site."""
    others = ("0/0:30,0:99:30",) * (n_samples - 1)
    return [
        Site(
            chrom,
            start + i * step,
            "A",
            "G",
            ("0/1:15,15:99:30", *others),
            info=info,
            fmt="GT:AD:GQ:DP",
        )
        for i, info in enumerate(infos)
    ]


def wgs_trio_sites(
    n: int,
    *,
    seed: int = 7,
    af: str = "population",
    denovo_frac: float = 0.01,
    start: int = 1000,
    step: int = 10,
    chrom: str = "chr1",
) -> list[Site]:
    """`n` sites with the mix of a real WGS trio, in CHILD / FATHER / MOTHER
    order. Most sites are inherited or carry nothing for the child:

      about `denovo_frac`  de novo het (child 0/1, parents 0/0): the candidates
      of the rest, 44%     het inherited from the father
                      44%  het inherited from the mother
                      10%  child hom-ref, one parent het (not a candidate)
                       2%  child hom-alt, both parents het (recessive)

    The site kinds are the same for a given `seed` whatever `af` is.

    `af` sets INFO/AF:
      "population"  a rare frequency (0.0001 to 0.001) on de novo sites and a
                    common one (0.05 to 0.45) on the others, as an annotated
                    VCF has: the default --max-af 0.01 keeps only the de novo
                    sites, so about `denovo_frac` of the store are candidates;
      "cohort"      the callset's own ac/an, as GATK joint calling writes it
                    (1/6 on a de novo or lone-parent het, 1/3 inherited,
                    2/3 hom-alt): the rarity gate is skipped and every
                    inherited het becomes a candidate.
    """
    rng = random.Random(seed)
    freqs = random.Random(seed + 1)  # own stream: the kinds do not depend on af
    ref_gt, het_gt, alt_gt = "0/0:30,0:99:30", "0/1:15,15:99:30", "1/1:0,30:99:30"
    # kind -> (child, father, mother, ac of 6 alleles)
    layout = {
        "denovo": (het_gt, ref_gt, ref_gt, 1),
        "pat": (het_gt, het_gt, ref_gt, 2),
        "mat": (het_gt, ref_gt, het_gt, 2),
        "parent": (ref_gt, het_gt, ref_gt, 1),
        "rec": (alt_gt, het_gt, het_gt, 4),
    }
    if af not in ("population", "cohort"):
        raise ValueError(af)
    sites = []
    for i in range(n):
        if rng.random() < denovo_frac:
            kind = "denovo"
        else:
            kind = rng.choices(
                ("pat", "mat", "parent", "rec"), weights=(44, 44, 10, 2)
            )[0]
        *gts, ac = layout[kind]
        if af == "cohort":
            freq = ac / 6
        elif kind == "denovo":
            freq = freqs.uniform(0.0001, 0.001)
        else:
            freq = freqs.uniform(0.05, 0.45)
        sites.append(
            Site(
                chrom,
                start + i * step,
                "A",
                "G",
                tuple(gts),
                info=f"AF={freq:.6g}",
                fmt="GT:AD:GQ:DP",
            )
        )
    return sites


TRIO_SAMPLES = ["CHILD", "FATHER", "MOTHER"]
TRIO_PED = (
    "FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\n"
    "FAM1\tFATHER\t0\t0\t1\t1\n"
    "FAM1\tMOTHER\t0\t0\t2\t1\n"
)
