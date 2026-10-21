"""Timing gate for a default `db trio`.

The desktop apps run `db trio --format json` (category all) on every
proband selection, and the family program adds work to that call (parentage,
lineage, X models). This gate fails if the default call gets slower than a
stated budget, so a later step has to say so when it spends time.

Budgets are wall-clock seconds for the whole CLI process (interpreter start,
store open, every model, JSON), best of two runs, on the fixtures below.
`BASELINE` is what the machine this was written on measured (Apple silicon,
warm cache); the budget is several times that so a slow CI runner does not
flake, and tight enough that a second pass over the genotype table or a
per-candidate round trip shows up. Rerun this file after each step that
adds a default-on scan.

Not covered by the committed fixtures: a WGS-sized store. Set
VCFCLICK_TIMING_SITES=N (e.g. 200000) to build a synthetic N-site trio and
time it too, with VCFCLICK_TIMING_BUDGET_S as the limit.

Where the time goes today: on related.vcf (3,003 sites, 7 samples), about
2.6 of the 3.8 s is `comphet_genes` opening an annotation-store connection
for each of ~740 candidates, whether or not genes are loaded.
"""

from __future__ import annotations

import os
import time

import pytest

from tests.conftest import FIXTURES, run_cli
from tests.fixtures._vcfgen import (
    TRIO_PED,
    TRIO_SAMPLES,
    trio_sites,
    write_vcf,
)

BACKENDS = ("chdb", "duckdb")
GIAB = FIXTURES / "giab"

# name -> (vcf, `db ped` args, `db trio` args, baseline seconds, budget seconds)
CASES = {
    "related": (
        FIXTURES / "related.vcf",
        ("--proband", "C1", "--father", "F", "--mother", "M"),
        ("--proband", "C1"),
        3.8,
        15.0,
    ),
    "giab_denovo": (
        GIAB / "denovo_trio.vcf.gz",
        (str(GIAB / "denovo_trio.ped"),),
        ("--proband", "HG002"),
        0.3,
        5.0,
    ),
    "giab_cftr": (
        GIAB / "cftr_trio.vcf.gz",
        (str(GIAB / "cftr_trio.ped"),),
        ("--proband", "HG002"),
        0.3,
        5.0,
    ),
}


def _best_of_two(home, backend, *args) -> float:
    best = float("inf")
    for _ in range(2):
        t = time.perf_counter()
        run_cli(home, backend, *args)
        best = min(best, time.perf_counter() - t)
    return best


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("name", list(CASES))
def test_default_db_trio_stays_within_budget(name, backend, vcfclick_home, monkeypatch):
    vcf, ped_args, trio_args, baseline, budget = CASES[name]
    # no annotation store: time the command, not whatever a developer loaded
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(vcfclick_home / "ann.duckdb"))
    run_cli(vcfclick_home, backend, "db", "create", "t")
    run_cli(
        vcfclick_home, backend, "db", "ingest", "t", str(vcf),
        "--ingest-id", "g1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(vcfclick_home, backend, "db", "ped", "t", *ped_args)

    # the apps' call (JSON, category all) and the terminal default (text)
    for fmt in ("json", "text"):
        took = _best_of_two(
            vcfclick_home, backend, "db", "trio", "t", *trio_args, "--format", fmt
        )
        assert took < budget, (
            f"`db trio --format {fmt}` on {name} ({backend}) took {took:.1f}s, "
            f"budget {budget:g}s (baseline {baseline}s). A step that adds a "
            "default-on scan must make it lazy or raise this budget on purpose."
        )


@pytest.mark.skipif(
    not os.environ.get("VCFCLICK_TIMING_SITES"),
    reason="set VCFCLICK_TIMING_SITES=N to time a synthetic N-site trio",
)
@pytest.mark.parametrize("backend", BACKENDS)
def test_default_db_trio_on_a_synthetic_store(backend, vcfclick_home, tmp_path):
    n = int(os.environ["VCFCLICK_TIMING_SITES"])
    budget = float(os.environ.get("VCFCLICK_TIMING_BUDGET_S", "120"))
    vcf = write_vcf(
        tmp_path / "big.vcf", TRIO_SAMPLES, trio_sites(n, af="popular", step=10)
    )
    ped = tmp_path / "trio.ped"
    ped.write_text(TRIO_PED)
    run_cli(vcfclick_home, backend, "db", "create", "t")
    run_cli(
        vcfclick_home, backend, "db", "ingest", "t", str(vcf),
        "--ingest-id", "g1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(vcfclick_home, backend, "db", "ped", "t", str(ped))
    took = _best_of_two(
        vcfclick_home,
        backend,
        "db",
        "trio",
        "t",
        "--proband",
        "CHILD",
        "--format",
        "json",
    )
    assert took < budget, f"{n} sites ({backend}): {took:.1f}s, budget {budget:.0f}s"
