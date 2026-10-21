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


def write_vcf(
    path: Path,
    samples: list[str],
    sites: list[Site],
    *,
    contigs: tuple[str, ...] = ("chr1",),
    info_header: tuple[str, ...] = (
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    ),
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
    lines += list(info_header) + list(fmt_header)
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
                    "PASS",
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


TRIO_SAMPLES = ["CHILD", "FATHER", "MOTHER"]
TRIO_PED = (
    "FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\n"
    "FAM1\tFATHER\t0\t0\t1\t1\n"
    "FAM1\tMOTHER\t0\t0\t2\t1\n"
)
