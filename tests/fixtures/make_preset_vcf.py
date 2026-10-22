"""Generate the preset fixtures (deterministic, hand-planted, no randomness).

    python tests/fixtures/make_preset_vcf.py

Writes, next to this file:

  presets.vcf.gz (+ .tbi)      one trio, CHILD / FATHER / MOTHER, with sites
                               planted exactly on the boundary of each built-in
                               preset's gates
  presets_nogq.vcf.gz (+ .tbi) the same sites with GT only: no GQ, DP or AD
  presets.ped                  the trio
  presets_truth.json           what each site is and which built-in keeps it,
                               written by hand below, not computed by vcfclick

related.vcf cannot test a gate: it carries no GQ, DP or INFO/AF, and the
gates are NULL-tolerant, so its counts are the same at --min-gq 99 --min-dp 60
--max-af 0.0001. A test written against it would pass vacuously; these sites
are what makes a gate bite. Ingest with --keep-reference (a de novo needs the
parents' stored 0/0).

The boundaries (GQ 29/30, DP 19/20, het AB 0.29/0.31 and 0.69/0.71, INFO/AF
0.0009/0.0011, a gnomAD popmax either side of 0.0001) sit around the numbers
of `strict-denovo`; the other built-ins' values are all looser. gnomAD is not a
VCF field: the test seeds the two popmax values in `GNOMAD` into the annotation
store, and `GENES` as one gene over the compound-het pair.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))

from tests.fixtures._vcfgen import Site, write_vcf  # noqa: E402

SAMPLES = ["CHILD", "FATHER", "MOTHER"]
PED = (
    "FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\n"
    "FAM1\tFATHER\t0\t0\t1\t1\n"
    "FAM1\tMOTHER\t0\t0\t2\t1\n"
)
FILTERS = ('##FILTER=<ID=LowQual,Description="Low quality">',)


def call(gt: str, ref: int, alt: int, gq: int | str = 99, dp: int | str | None = None):
    """`GT:AD:GQ:DP`. DP defaults to the read total; a '.' is a missing value."""
    depth = ref + alt if dp is None else dp
    return f"{gt}:{ref},{alt}:{gq}:{depth}"


HOM_REF = call("0/0", 40, 0)
HET = call("0/1", 20, 20, dp=40)

# (name, position, child, father, mother, info, filter). A site is "de novo" when
# both parents are 0/0 and the child is het.
DE_NOVO = [
    ("dn_ok", 1000, HET, HOM_REF, HOM_REF, "AF=0.0005", "PASS"),
    (
        "dn_gq29",
        2000,
        call("0/1", 20, 20, gq=29),
        HOM_REF,
        HOM_REF,
        "AF=0.0005",
        "PASS",
    ),
    (
        "dn_gq30",
        2100,
        call("0/1", 20, 20, gq=30),
        HOM_REF,
        HOM_REF,
        "AF=0.0005",
        "PASS",
    ),
    ("dn_dp19", 3000, call("0/1", 10, 9), HOM_REF, HOM_REF, "AF=0.0005", "PASS"),
    ("dn_dp20", 3100, call("0/1", 10, 10), HOM_REF, HOM_REF, "AF=0.0005", "PASS"),
    ("dn_ab029", 4000, call("0/1", 71, 29), HOM_REF, HOM_REF, "AF=0.0005", "PASS"),
    ("dn_ab031", 4100, call("0/1", 69, 31), HOM_REF, HOM_REF, "AF=0.0005", "PASS"),
    ("dn_ab069", 5000, call("0/1", 31, 69), HOM_REF, HOM_REF, "AF=0.0005", "PASS"),
    ("dn_ab071", 5100, call("0/1", 29, 71), HOM_REF, HOM_REF, "AF=0.0005", "PASS"),
    ("dn_af0009", 6000, HET, HOM_REF, HOM_REF, "AF=0.0009", "PASS"),
    ("dn_af0011", 6100, HET, HOM_REF, HOM_REF, "AF=0.0011", "PASS"),
    # the parent's stored 0/0 has 4 alt reads of 40: not proof of absence
    (
        "dn_parent_alt",
        7000,
        HET,
        call("0/0", 36, 4, dp=40),
        HOM_REF,
        "AF=0.0005",
        "PASS",
    ),
    # ... and one alt read is noise
    (
        "dn_parent_alt_low",
        7100,
        HET,
        call("0/0", 39, 1, dp=40),
        HOM_REF,
        "AF=0.0005",
        "PASS",
    ),
    ("dn_nonpass", 8000, HET, HOM_REF, HOM_REF, "AF=0.0005", "LowQual"),
    # the father's depth is missing: tolerated by default, a fail when strict
    (
        "dn_null_dp_parent",
        9000,
        HET,
        call("0/0", 40, 0, dp="."),
        HOM_REF,
        "AF=0.0005",
        "PASS",
    ),
    ("dn_gnomad_rare", 10000, HET, HOM_REF, HOM_REF, "AF=0.0005", "PASS"),
    ("dn_gnomad_common", 10100, HET, HOM_REF, HOM_REF, "AF=0.0005", "PASS"),
]
HOM_ALT = call("1/1", 0, 40)
INHERITED = [
    # recessive: child hom-alt, both parents het
    ("rec_ok", 11000, HOM_ALT, HET, HET, "AF=0.002", "PASS"),
    # dominant: child het, father het, mother 0/0
    ("dom_pat", 12000, HET, HET, HOM_REF, "AF=0.0005", "PASS"),
    # a compound-het pair in GENES: one het from each parent
    ("ch_pat", 13100, HET, HET, HOM_REF, "AF=0.002", "PASS"),
    ("ch_mat", 13500, HET, HOM_REF, HET, "AF=0.002", "PASS"),
]
# 30 sites nobody is a candidate at (the child is 0/0, one parent is het), so
# INFO/AF is present on enough sites for the cohort-AF probe to judge it.
FILLERS = [
    (
        f"fill_{i:02d}",
        20000 + 100 * i,
        HOM_REF,
        HET,
        HOM_REF,
        f"AF={0.03 + 0.013 * i:.3f}",
        "PASS",
    )
    for i in range(30)
]

GNOMAD = {  # position -> (af, popmax); every other candidate is absent from gnomAD
    10000: (0.00003, 0.00005),
    10100: (0.0003, 0.0005),
}
GENES = [("G1", "chr1", 13000, 14000)]

# What each built-in keeps, by hand. Every built-in is run with its own
# category list and its own numbers (cli/builtin_presets.py).
EXPECT = {
    "broad-first-pass": {
        # looser than every boundary and AF <= 0.01: all seventeen de novo sites
        "denovo": [n for n, *_ in DE_NOVO],
        "recessive": ["rec_ok"],
        "dominant": ["dom_pat", "ch_pat", "ch_mat"],
        "comphet": ["G1"],
    },
    "strict-denovo": {
        "denovo": [
            "dn_ok",
            "dn_gq30",
            "dn_dp20",
            "dn_ab031",
            "dn_ab069",
            "dn_af0009",
            "dn_parent_alt_low",
            "dn_gnomad_rare",
        ],
    },  # fmt: skip
    "rare-recessive": {"recessive": ["rec_ok"], "comphet": ["G1"]},
    "rare-dominant": {
        # GQ 25, DP 12, AB 0.25-0.75, AF 0.001, gnomAD 0.0001, no PASS or parent screen
        "denovo": [
            n for n, *_ in DE_NOVO if n not in ("dn_af0011", "dn_gnomad_common")
        ],
        "dominant": ["dom_pat"],
    },
}


def sites(*, with_quality: bool = True) -> list[Site]:
    rows = [*DE_NOVO, *INHERITED, *FILLERS]
    out = []
    for _, pos, child, father, mother, info, flt in rows:
        gts = (child, father, mother)
        fmt = "GT:AD:GQ:DP"
        if not with_quality:
            gts = tuple(g.split(":")[0] for g in gts)
            fmt = "GT"
        out.append(Site("chr1", pos, "A", "G", gts, info=info, fmt=fmt, filter=flt))
    return out


def _bgzip(path: Path) -> None:
    gz = path.with_name(path.name + ".gz")
    with open(gz, "wb") as fh:
        subprocess.run(["bgzip", "-c", str(path)], stdout=fh, check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(gz)], check=True)
    path.unlink()


def main() -> None:
    for name, quality in (("presets", True), ("presets_nogq", False)):
        plain = write_vcf(
            HERE / f"{name}.vcf", SAMPLES, sites(with_quality=quality),
            filter_header=FILTERS,
        )  # fmt: skip
        _bgzip(plain)
    (HERE / "presets.ped").write_text(PED)
    truth = {
        "positions": {n: pos for n, pos, *_ in [*DE_NOVO, *INHERITED]},
        "gnomad_popmax": {str(pos): v[1] for pos, v in GNOMAD.items()},
        "gene": {"symbol": "G1", "chrom": "chr1", "start": 13000, "end": 14000},
        "kept": EXPECT,
        "note": "hand-written; see make_preset_vcf.py for the reason of each site",
    }
    (HERE / "presets_truth.json").write_text(json.dumps(truth, indent=2) + "\n")
    print(f"wrote {len(sites())} sites x {len(SAMPLES)} samples")


if __name__ == "__main__":
    main()
