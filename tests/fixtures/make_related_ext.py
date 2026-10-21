"""Generate tests/fixtures/related_ext.{vcf,ped,truth.json} (deterministic).

An extended pedigree with two families, for the pedigree features. It does
not touch related.vcf / related.ped, which other suites pin (7 samples).
The seed is its own; the file is small on purpose (tests ingest it once per
backend), and a larger or recombining variant belongs to the feature that
needs one.

    famA (13)   PGF PGM  = paternal grandparents      MGF MGM = maternal
                UF, F    = children of PGF, PGM       M       = child of MGF, MGM
                M2       = unrelated mother           C1, C2  = children of F, M
                C2T      = declared full sibling of C2 (truly its MZ twin)
                H1       = child of F and M2 (a paternal half sibling)
                H2       = declared child of F and M2 (truly F and M's child)
    famB (5)    FB MB    = parents of K1, K2
                XB       = declared child of FB, MB (truly F and M's child)

Phenotypes in the PED: C1, H1, K1 affected; C2, UF and F unaffected; the
rest unknown. By the proband rule that makes the probands C1, C2T, H1, H2
(famA) and K1, K2, XB (famB): C2 is explicitly unaffected, UF is an
unaffected leaf, and F and M are trios with children of their own.

Every genotype is explicit (0/0 included), so sparse and --keep-reference
ingests see the same truth. Background sites are Mendelian for the TRUE
parents with INFO/AF 0.2-0.5, so the default --max-af 0.01 excludes them;
the planted sites (AF 0.001, chr1 from 200,000,000) are the only default
candidates, and truth.json says which sample each one is a candidate for.

Run:  python tests/fixtures/make_related_ext.py
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))

from tests.fixtures._vcfgen import Site, write_vcf  # noqa: E402

SEED = 20261021
N_BACKGROUND = 60
PLANTED_START = 200_000_000

SAMPLES = [
    "PGF", "PGM", "MGF", "MGM", "M2", "UF", "F", "M", "C1", "C2", "C2T",
    "H1", "H2", "FB", "MB", "K1", "K2", "XB",
]  # fmt: skip

# sample -> (family, declared father, declared mother, sex, phenotype code)
# phenotype: 1 unaffected, 2 affected, 0 unknown (PED convention)
DECLARED = {
    "PGF": ("famA", "0", "0", 1, 0),
    "PGM": ("famA", "0", "0", 2, 0),
    "MGF": ("famA", "0", "0", 1, 0),
    "MGM": ("famA", "0", "0", 2, 0),
    "M2": ("famA", "0", "0", 2, 0),
    "UF": ("famA", "PGF", "PGM", 1, 1),
    "F": ("famA", "PGF", "PGM", 1, 1),
    "M": ("famA", "MGF", "MGM", 2, 0),
    "C1": ("famA", "F", "M", 1, 2),
    "C2": ("famA", "F", "M", 2, 1),
    "C2T": ("famA", "F", "M", 2, 0),
    "H1": ("famA", "F", "M2", 1, 2),
    "H2": ("famA", "F", "M2", 2, 0),
    "FB": ("famB", "0", "0", 1, 0),
    "MB": ("famB", "0", "0", 2, 0),
    "K1": ("famB", "FB", "MB", 1, 2),
    "K2": ("famB", "FB", "MB", 2, 0),
    "XB": ("famB", "FB", "MB", 1, 0),
}

# Who the parents REALLY are in the genotypes (they differ from DECLARED for
# the two mislabelled samples). C2T is a copy of C2, handled separately.
TRUE_PARENTS = {
    "UF": ("PGF", "PGM"),
    "F": ("PGF", "PGM"),
    "M": ("MGF", "MGM"),
    "C1": ("F", "M"),
    "C2": ("F", "M"),
    "H1": ("F", "M2"),
    "H2": ("F", "M"),
    "K1": ("FB", "MB"),
    "K2": ("FB", "MB"),
    "XB": ("F", "M"),
}

GT = {0: "0/0:40,0:99:40", 1: "0/1:20,20:99:40", 2: "1/1:0,40:99:40"}

# Planted sites: (id, {sample: alt count}, the model(s) it is a candidate
# for and the probands it is a candidate for, under the default gates and
# a `--keep-reference` ingest). Every sample not listed is 0/0.
PLANTED = [
    ("denovo_C1", {"C1": 1}, {"denovo": ["C1"]}),
    ("denovo_K1", {"K1": 1}, {"denovo": ["K1"]}),
    # inherited from PGF through F to H1 (F is also a carrier, so F's own
    # trio sees a dominant site): H1 is the proband, F a non-proband trio
    ("dominant_H1", {"PGF": 1, "F": 1, "H1": 1}, {"dominant": ["H1"]}),
    ("recessive_K1", {"FB": 1, "MB": 1, "K1": 2}, {"recessive": ["K1"]}),
]
# The same sites seen by the trios that are not probands (--all-trios).
PLANTED_FOR_NON_PROBANDS = {"dominant_H1": {"dominant": ["F"]}}


def declared_ped() -> str:
    lines = []
    for s in SAMPLES:
        fam, f, m, sex, pheno = DECLARED[s]
        lines.append(f"{fam}\t{s}\t{f}\t{m}\t{sex}\t{pheno}")
    return "\n".join(lines) + "\n"


def build() -> tuple[list[Site], dict]:
    rnd = random.Random(SEED)

    def founder(p: float) -> tuple[int, int]:
        return (int(rnd.random() < p), int(rnd.random() < p))

    def child(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
        return (rnd.choice(a), rnd.choice(b))

    def gt(h: tuple[int, int]) -> str:
        return GT[sum(h)]

    sites: list[Site] = []
    for i in range(N_BACKGROUND):
        p = rnd.uniform(0.05, 0.5)
        hap: dict[str, tuple[int, int]] = {}
        for s in SAMPLES:
            if s in TRUE_PARENTS:
                continue
            if s != "C2T":
                hap[s] = founder(p)
        for s in ("UF", "F", "M", "C1", "C2", "H1", "H2", "K1", "K2", "XB"):
            f, m = TRUE_PARENTS[s]
            hap[s] = child(hap[f], hap[m])
        calls = {s: gt(hap[s]) for s in hap}
        # C2T: an MZ twin of C2 -- the same genotypes, 0.2% of calls flipped
        twin = sum(hap["C2"])
        if rnd.random() < 0.002:
            twin = (twin + 1) % 3
        calls["C2T"] = GT[twin]
        sites.append(
            Site(
                "chr1",
                1000 + i * 3_000_000,
                rnd.choice("ACG"),
                "T",
                tuple(calls[s] for s in SAMPLES),
                info=f"AF={rnd.uniform(0.2, 0.5):.3f}",
                fmt="GT:AD:GQ:DP",
            )
        )
    planted = {}
    for k, (name, carriers, models) in enumerate(PLANTED):
        pos = PLANTED_START + k * 100_000
        sites.append(
            Site(
                "chr1",
                pos,
                "A",
                "G",
                tuple(GT[carriers.get(s, 0)] for s in SAMPLES),
                info="AF=0.001",
                fmt="GT:AD:GQ:DP",
            )
        )
        planted[name] = {"pos": pos, "carriers": carriers, "candidate_for": models}
    truth = {
        "samples": SAMPLES,
        "families": {
            fam: [s for s in SAMPLES if DECLARED[s][0] == fam]
            for fam in ("famA", "famB")
        },
        "declared_parents": {
            s: [DECLARED[s][1], DECLARED[s][2]]
            for s in SAMPLES
            if DECLARED[s][1] != "0"
        },
        "true_parents": {s: list(p) for s, p in TRUE_PARENTS.items()},
        "mz_twin_of": {"C2T": "C2"},
        "affected": [s for s in SAMPLES if DECLARED[s][4] == 2],
        "unaffected": [s for s in SAMPLES if DECLARED[s][4] == 1],
        # by the proband rule (an analysable child that is not explicitly
        # unaffected and either has no children or is affected)
        "probands": {
            "famA": ["C1", "C2T", "H1", "H2"],
            "famB": ["K1", "K2", "XB"],
        },
        "trios_not_probands": {
            "C2": "explicitly unaffected",
            "UF": "explicitly unaffected",
            "F": "explicitly unaffected, and has children",
            "M": "has children and is not affected",
        },
        "siblings": {
            "C1": {
                "C2": "full",
                "C2T": "full",
                "H1": "half_paternal",
                "H2": "half_paternal",
            },
            "H1": {
                "C1": "half_paternal",
                "C2": "half_paternal",
                "C2T": "half_paternal",
                "H2": "full",
            },
            "K1": {"K2": "full", "XB": "full"},
            "F": {"UF": "full"},
        },
        "generation": {
            **{s: 0 for s in ("PGF", "PGM", "MGF", "MGM", "M2", "FB", "MB")},
            **{s: 1 for s in ("UF", "F", "M", "K1", "K2", "XB")},
            **{s: 2 for s in ("C1", "C2", "C2T", "H1", "H2")},
        },
        "planted": planted,
        "candidates_for_non_probands": PLANTED_FOR_NON_PROBANDS,
    }
    return sites, truth


def write(out: Path = HERE) -> None:
    sites, truth = build()
    sites.sort(key=lambda s: (s.chrom, s.pos))
    write_vcf(out / "related_ext.vcf", SAMPLES, sites)
    (out / "related_ext.ped").write_text(declared_ped())
    (out / "related_ext.truth.json").write_text(json.dumps(truth, indent=2) + "\n")


if __name__ == "__main__":
    write()
    print(f"wrote {N_BACKGROUND + len(PLANTED)} sites x {len(SAMPLES)} samples")
