"""Generate tests/fixtures/related.vcf + related.ped (deterministic).

Seven samples with known relationships, 3000 biallelic SNVs, all calls
written explicitly (0/0 included) so the sparse and --keep-reference
ingests see the same truth:

  F, M      founders (unrelated to each other)
  C1, C2    full siblings, children of F and M (Mendelian transmission)
  D         exact duplicate of C1
  U1        unrelated founder
  U2        unrelated founder with ~3% no-calls (./.)

related.ped declares F/M as parents of C1, C2 -- and, deliberately wrongly,
of U1 (a pedigree error the relatedness check must flag).

Run:  python tests/fixtures/make_related_vcf.py
"""

from __future__ import annotations

import random
from pathlib import Path

HERE = Path(__file__).parent
N_SITES = 3000
rnd = random.Random(20260927)


def founder(p: float) -> tuple[int, int]:
    return (int(rnd.random() < p), int(rnd.random() < p))


def child(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return (rnd.choice(a), rnd.choice(b))


def gt(h: tuple[int, int] | None) -> str:
    if h is None:
        return "./."
    return f"{min(h)}/{max(h)}"


samples = ["F", "M", "C1", "C2", "D", "U1", "U2"]
lines = [
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "##contig=<ID=chr1,length=248956422>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
]
bases = "ACGT"
for i in range(N_SITES):
    p = rnd.uniform(0.05, 0.5)
    f, m, u1, u2 = founder(p), founder(p), founder(p), founder(p)
    c1, c2 = child(f, m), child(f, m)
    calls = [f, m, c1, c2, c1, u1, None if rnd.random() < 0.03 else u2]
    ref = rnd.choice(bases)
    alt = rnd.choice([b for b in bases if b != ref])
    pos = 1000 + i * 50_000  # ~150 Mb: genome-scale spread, as kinship needs
    lines.append(
        f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t"
        + "\t".join(gt(c) for c in calls)
    )

(HERE / "related.vcf").write_text("\n".join(lines) + "\n")
(HERE / "related.ped").write_text(
    "fam1\tF\t0\t0\t1\t1\n"
    "fam1\tM\t0\t0\t2\t1\n"
    "fam1\tC1\tF\tM\t1\t2\n"
    "fam1\tC2\tF\tM\t2\t2\n"
    "fam1\tU1\tF\tM\t1\t2\n"  # wrong on purpose: U1 is unrelated to F and M
)
print(f"wrote {N_SITES} sites x {len(samples)} samples")
