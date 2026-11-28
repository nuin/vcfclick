"""Timing gate for `db segregate` on a synthetic whole-genome-style sibship.

The candidate sites are the union of every affected child's alt rows, and the
pivot query reads them three times on chDB (which inlines the `cand` CTE into
each use). This gate fails if the default call gets slower than a multiple of
one full scan of the same store's genotypes on the same machine, the budget
rule of tests/test_trio_timing.py (a slow CI runner scales both sides).

The store: `VCFCLICK_TIMING_SITES` sites (20,000 by default) in a six-sample
family, FA MO S1 S2 (affected) S3 S4 (unaffected), keep-reference, INFO/AF the
callset's own (so the rarity gate is skipped and every inherited het is a
candidate, the heavy case). The mix of a real family: about 40 percent paternal
hets and 40 percent maternal hets carried by both affected children, 10
percent carried by one affected child, 5 percent a het sibling, 3 percent
recessive-shaped, 2 percent de novo.

Measured on Apple silicon, 200,000 sites (1.2 million genotype rows), best of
two, `db segregate --format json` / text against a full scan of genotypes:

                    scan     json     text
    DuckDB          0.4 s    5.2 s    5.7 s    (13x / 14x the scan)
    chDB            0.8 s    4.0 s    2.8 s    (5x / 4x the scan)

The budget is 35x the scan. A real WGS store (10^8 genotype rows) is not
measured: this file bounds the ratio at the sizes CI can afford, and
`VCFCLICK_TIMING_SITES=200000` is the size to run before a release.
"""

from __future__ import annotations

import json
import random

import pytest

from tests.conftest import run_cli
from tests.fixtures._vcfgen import Site, write_vcf
from tests.test_trio_timing import (
    SITE_STEP,
    SYNTHETIC_SITES,
    _baseline,
    _best_of_two,
    _synthetic_budget,
)

pytestmark = pytest.mark.timing

BACKENDS = ("chdb", "duckdb")
SAMPLES = ["FA", "MO", "S1", "S2", "S3", "S4"]
PED = """\
fam1 FA 0 0 1 0
fam1 MO 0 0 2 0
fam1 S1 FA MO 1 2
fam1 S2 FA MO 2 2
fam1 S3 FA MO 1 1
fam1 S4 FA MO 2 1
"""
REF, HET, ALT = "0/0:30,0:99:30", "0/1:15,15:99:30", "1/1:0,30:99:30"

# kind -> (weight, [FA, MO, S1, S2, S3, S4])
KINDS = {
    "pat": (40, [HET, REF, HET, HET, REF, REF]),
    "mat": (40, [REF, HET, HET, HET, REF, REF]),
    "one": (10, [HET, REF, HET, REF, REF, REF]),
    "sib": (5, [HET, REF, HET, HET, HET, REF]),
    "rec": (3, [HET, HET, ALT, ALT, REF, HET]),
    "dn": (2, [REF, REF, HET, HET, REF, REF]),
}


def sibship_sites(n: int, seed: int = 7) -> list[Site]:
    rng = random.Random(seed)
    names, weights = zip(*[(k, w) for k, (w, _) in KINDS.items()], strict=True)
    sites = []
    for i in range(n):
        kind = rng.choices(names, weights)[0]
        gts = KINDS[kind][1]
        alt = sum({REF: 0, HET: 1, ALT: 2}[g] for g in gts)
        sites.append(
            Site(
                "chr1", 1000 + i * SITE_STEP, "A", "G", tuple(gts),
                info=f"AF={alt / 12:.6f}", fmt="GT:AD:GQ:DP",
            )
        )  # fmt: skip
    return sites


def _store(backend, home, tmp_path, monkeypatch):
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    vcf = write_vcf(tmp_path / "sib.vcf", SAMPLES, sibship_sites(SYNTHETIC_SITES))
    ped = tmp_path / "sib.ped"
    ped.write_text(PED)
    run_cli(home, backend, "db", "create", "t")
    run_cli(
        home, backend, "db", "ingest", "t", str(vcf),
        "--ingest-id", "g1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(home, backend, "db", "ped", "t", str(ped))


@pytest.mark.parametrize("backend", BACKENDS)
def test_default_db_segregate_on_a_synthetic_wgs_sibship(
    backend, vcfclick_home, tmp_path, monkeypatch
):
    _store(backend, vcfclick_home, tmp_path, monkeypatch)
    budget = _synthetic_budget(SYNTHETIC_SITES, _baseline(vcfclick_home, backend))
    took, result = _best_of_two(
        vcfclick_home, backend, "db", "segregate", "t", "--format", "json"
    )
    doc = json.loads(result.stdout)["families"][0]
    assert doc["gates"]["af_gate"] == "skipped_cohort_derived"
    # the heavy case really is heavy: most sites are candidates of a model
    dominant = doc["models"]["dominant"]
    assert dominant["count"] > 0.7 * SYNTHETIC_SITES and dominant["truncated"]
    assert doc["models"]["denovo"]["count"] > 0
    assert took < budget, (
        f"`db segregate --format json` on {SYNTHETIC_SITES} sites ({backend}) took "
        f"{took:.1f}s, budget {budget:.1f}s (35x a full scan of genotypes). A "
        "second pass over the genotype table, or a query per candidate, costs "
        "this much; keep the pivot to one scan per model."
    )
    text_took, _ = _best_of_two(vcfclick_home, backend, "db", "segregate", "t")
    assert text_took < budget, (
        f"text summary took {text_took:.1f}s (budget {budget:.1f}s)"
    )
