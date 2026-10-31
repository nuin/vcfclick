"""Timing gate for a default `db trio`.

The desktop apps run `db trio --format json` (category all) on every
proband selection, and every model added to that call (parentage, lineage,
X-linked, ...) spends time there. This gate fails if the default call gets
slower than a stated budget, so a change has to say so when it spends time.

Budgets are wall-clock seconds for the whole CLI process (interpreter start,
store open, every model, JSON), best of two runs. The numbers below were
measured on Apple silicon with a warm cache; each budget is several times the
measurement so a slow CI runner does not flake, and still far below what a
second pass over the genotype table or a per-candidate round trip costs.
Rerun this file after each step that adds a default-on scan.

Committed fixtures (`CASES`): related.vcf (3,003 sites, 7 samples) takes
0.2 s on DuckDB and 0.3 s on chDB, the GIAB trios 0.2 s and 0.3 s. Budget 5 s.

Synthetic WGS-style trio (`test_default_db_trio_on_a_synthetic_wgs_store`):
`wgs_trio_sites` builds the mix of a real trio, in which most sites are
inherited hets or carry nothing for the child and about 1% are de novo. The
annotation store is its own file with 2,000 filler genes and one gene every
3 kb across the sites, so the compound-het grouping and the per-candidate
lookups do real work. It is run twice, because the default AF gate decides
how many candidates there are:

  population  INFO/AF is a population frequency, the rarity gate applies,
              and about 1% of the sites are candidates;
  cohort      INFO/AF is the callset's own (what GATK writes for a family),
              the gate is skipped, and every inherited het is a candidate
              (87% of the sites), so the compound-het gene grouping looks up
              hundreds of thousands of positions.

Measured, `db trio --format json`, best of two:

                      20,000 sites        200,000 sites
                      DuckDB   chDB       DuckDB   chDB
    population          0.5 s   0.6 s       1.1 s   1.4 s
    cohort              1.5 s   1.6 s       3.2 s   3.0 s

Ingesting the 200,000-site store takes about 8 s. The budget is
`SYNTHETIC_BASE_S` plus `SYNTHETIC_PER_SITE_S` per site: 5.6 s at the default
20,000 sites and 11 s at 200,000, about 3.5 times the slowest figure above.
`VCFCLICK_TIMING_SITES=N` sets the number of sites and
`VCFCLICK_TIMING_BUDGET_S` replaces the budget.

`--gnomad-max-af 0.01` (`test_gnomad_filter_on_a_synthetic_wgs_store`) is
timed on the cohort store, where nearly every inherited het is a candidate
that the filter has to look up, as a JSON document and as the text summary.
`none` has no gnomAD loaded (the case that was measured when the filter
looked candidates up one by one); `slice` has gnomAD for every other site, a
quarter of them common, so the batched join runs and drops sites. JSON /
text, best of two:

                      20,000 sites                200,000 sites
                      DuckDB        chDB          DuckDB        chDB
    none              1.5 / 0.3 s   1.6 / 0.4 s   3.3 / 1.0 s   3.2 / 0.8 s
    slice             2.4 / 0.5 s   2.5 / 0.6 s   5.6 / 1.4 s   5.5 / 1.3 s

The JSON document with a slice costs more than without one because it
renders up to 1,000 candidates per model, each with a gene and a gnomAD
lookup (about 1 ms apiece, 3,000 of them at the largest size); the filter
itself, one join for all candidates of a model, takes 0.1 to 0.2 s per
170,000 candidates.

What it has caught: the cohort row at 20,000 sites took 90 s while each
candidate opened the annotation store (4 ms apiece) and the rows of the
document did the same (7 ms per lookup, 1,600 lookups). The first is fixed by
`annotations.db.GeneIndex`, which reads the gene table once, the second by
keeping one store connection open while the document is built. The gnomAD
filter then did the same per candidate (0.65 ms on a held connection, 4.6 ms
when not): on a 5,000-site cohort store JSON took 6.5 s and the text listing
40 s against 1.2 s and 0.4 s now, which is 4 and 25 minutes at 200,000 sites.
It is fixed by `annotations.gnomad_af_many`, which looks every candidate up
in one query. A lookup that costs 0.1 ms per candidate would add 17 s at
200,000 sites and show up only there, which is why that size is worth running
before a release.
"""

from __future__ import annotations

import json
import os
import random
import time

import pytest

from tests.conftest import FIXTURES, run_cli
from tests.fixtures._vcfgen import (
    TRIO_PED,
    TRIO_SAMPLES,
    wgs_trio_sites,
    write_vcf,
)

# Wall-clock gates need one machine, not every CI cell: the test matrix
# deselects them (-m "not timing") and the `timing` job runs them. A plain
# `pytest tests/` still runs them locally.
pytestmark = pytest.mark.timing

BACKENDS = ("chdb", "duckdb")
GIAB = FIXTURES / "giab"

# name -> (vcf, `db ped` args, `db trio` args, measured seconds, budget seconds)
CASES = {
    "related": (
        FIXTURES / "related.vcf",
        ("--proband", "C1", "--father", "F", "--mother", "M"),
        ("--proband", "C1"),
        0.3,
        5.0,
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

SYNTHETIC_SITES = int(os.environ.get("VCFCLICK_TIMING_SITES", "20000"))
SYNTHETIC_BASE_S = 5.0  # interpreter start, store open, rendering a document
SYNTHETIC_PER_SITE_S = 0.00003  # measured: 0.00001 on the cohort rows
SITE_STEP = 10  # base pairs between synthetic sites


def _best_of_two(home, backend, *args):
    """(best wall-clock seconds, the last run's result) of two runs."""
    best, result = float("inf"), None
    for _ in range(2):
        t = time.perf_counter()
        result = run_cli(home, backend, *args)
        best = min(best, time.perf_counter() - t)
    return best, result


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
        took, _ = _best_of_two(
            vcfclick_home, backend, "db", "trio", "t", *trio_args, "--format", fmt
        )
        assert took < budget, (
            f"`db trio --format {fmt}` on {name} ({backend}) took {took:.1f}s, "
            f"budget {budget:g}s (measured {baseline}s). A step that adds a "
            "default-on scan must make it lazy or raise this budget on purpose."
        )


def _seed_genes(n_sites: int) -> int:
    """Fill the (env-redirected) annotation store with genes: chr1 genes of
    2 kb every 3 kb across the synthetic sites, so the compound-het grouping
    finds genes, plus filler genes on chr2 so the table has a real size.
    DuckDB in-process, never a chDB session. Returns the gene count."""
    from annotations.db import get_connection

    last = 1000 + n_sites * SITE_STEP
    genes = [
        (f"G{i}", "chr1", s, s + 1999)
        for i, s in enumerate(range(1000, last + 3000, 3000))
    ]
    rng = random.Random(11)
    for i in range(2000):
        start = rng.randint(1, 5_000_000)
        genes.append((f"F{i}", "chr2", start, start + rng.randint(500, 50_000)))
    conn = get_connection()
    conn.executemany(
        "INSERT INTO refseq_genes VALUES (?, ?, ?, ?, '+', 'id', 'protein_coding')",
        genes,
    )
    conn.close()
    return len(genes)


# Wall-clock seconds on a shared CI runner are 3-6x what they are on the
# machine the constants above were measured on, so the synthetic-store
# budgets are relative: a multiple of a full scan of the same store's
# genotypes on the same machine, run just before (interpreter start plus
# one pass over the table). Measured: db trio is 9-10x that baseline on
# both backends, 16x with a gnomAD slice loaded. The regressions this file
# exists for (one query or store open per candidate, a second pass) were
# 40-250x.
SYNTHETIC_RATIO = 35.0


def _baseline(home, backend) -> float:
    """Best of two `db query` full scans of the genotypes table."""
    took, _ = _best_of_two(
        home, backend, "db", "query", "t",
        "SELECT count(*) AS n, sum(gt) AS s FROM genotypes",
    )  # fmt: skip
    return took


def _synthetic_budget(n: int, baseline: float | None = None) -> float:
    if os.environ.get("VCFCLICK_TIMING_BUDGET_S"):
        return float(os.environ["VCFCLICK_TIMING_BUDGET_S"])
    if baseline is not None:
        return SYNTHETIC_RATIO * baseline
    return SYNTHETIC_BASE_S + SYNTHETIC_PER_SITE_S * n


def _wgs_store(af, backend, home, tmp_path, monkeypatch):
    """Ingest `SYNTHETIC_SITES` synthetic trio sites into db `t` and point
    the annotation store at a file of its own with genes seeded. Returns
    (sites, genes)."""
    n = SYNTHETIC_SITES
    # an annotation store of its own: never the repo's, which another process
    # may hold open, and the same genes on every machine
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    n_genes = _seed_genes(n)

    vcf = write_vcf(
        tmp_path / "wgs.vcf",
        TRIO_SAMPLES,
        wgs_trio_sites(n, af=af, step=SITE_STEP),
    )
    ped = tmp_path / "trio.ped"
    ped.write_text(TRIO_PED)
    run_cli(home, backend, "db", "create", "t")
    run_cli(
        home, backend, "db", "ingest", "t", str(vcf),
        "--ingest-id", "g1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(home, backend, "db", "ped", "t", str(ped))
    return n, n_genes


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("af", ("population", "cohort"))
def test_default_db_trio_on_a_synthetic_wgs_store(
    af, backend, vcfclick_home, tmp_path, monkeypatch
):
    n, n_genes = _wgs_store(af, backend, vcfclick_home, tmp_path, monkeypatch)
    budget = _synthetic_budget(n, _baseline(vcfclick_home, backend))

    took, result = _best_of_two(
        vcfclick_home,
        backend,
        "db", "trio", "t", "--proband", "CHILD", "--format", "json",
    )  # fmt: skip
    doc = json.loads(result.stdout)
    counts = {m: doc["models"][m]["count"] for m in doc["models"]}
    print(
        f"\n[timing] {n} sites, {n_genes} genes, af={af}, {backend}: "
        f"{took:.2f}s (budget {budget:.1f}s) {counts} af_gate="
        f"{doc['gates']['af_gate']}"
    )

    # the scenario is what its name says, so the time means what it claims
    if af == "population":
        assert doc["gates"]["af_gate"] == "applied"
        assert 0.005 * n <= counts["denovo"] <= 0.02 * n  # about 1% candidates
        assert counts["dominant"] == 0 and counts["comphet"] == 0
    else:
        assert doc["gates"]["af_gate"] == "skipped_cohort_derived"
        assert counts["dominant"] > 0.4 * n  # every inherited het
        assert counts["comphet"] > 0  # the gene grouping ran on all of them
    assert took < budget, (
        f"`db trio` on {n} synthetic sites ({af} AF, {backend}) took {took:.1f}s, "
        f"budget {budget:.1f}s. A per-candidate query or a second pass over the "
        "genotype table shows up here."
    )


def _seed_gnomad(n_sites: int) -> int:
    """Load gnomAD popmax for every other synthetic site: common (0.2) for
    one in four of all sites, rare (0.001) for another one in four, nothing
    for the rest. Returns the number of rows."""
    from annotations.db import get_connection

    rows = [
        ("chr1", 1000 + i * SITE_STEP, "A", "G", 0.001, 0.2 if i % 4 == 0 else 0.001)
        for i in range(0, n_sites, 2)
    ]
    conn = get_connection()
    conn.executemany("INSERT INTO gnomad_af VALUES (?, ?, ?, ?, ?, ?)", rows)
    conn.close()
    return len(rows)


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("gnomad", ("none", "slice"))
def test_gnomad_filter_on_a_synthetic_wgs_store(
    gnomad, backend, vcfclick_home, tmp_path, monkeypatch
):
    """`--gnomad-max-af` on the cohort-AF store, where the AF gate is
    skipped and nearly every inherited het is a candidate that the filter has
    to look up. The text listing and the JSON document are both timed: the
    first once spent 4.6 ms per candidate (25 minutes for a 200,000-site
    family), the second 0.65 ms (4 minutes), because each candidate went to
    the annotation store on its own. `none` is no gnomAD loaded, the case
    that was measured; `slice` has gnomAD for half the sites, so the batched
    join runs and drops the common ones."""
    n, n_genes = _wgs_store("cohort", backend, vcfclick_home, tmp_path, monkeypatch)
    n_rows = _seed_gnomad(n) if gnomad == "slice" else 0
    budget = _synthetic_budget(n, _baseline(vcfclick_home, backend))
    trio = ("db", "trio", "t", "--proband", "CHILD")

    base = json.loads(run_cli(vcfclick_home, backend, *trio, "--format", "json").stdout)
    base_counts = {m: base["models"][m]["count"] for m in base["models"]}
    assert base["gates"]["af_gate"] == "skipped_cohort_derived"
    assert base_counts["dominant"] > 0.4 * n

    for fmt in ("json", "text"):
        took, result = _best_of_two(
            vcfclick_home, backend, *trio, "--gnomad-max-af", "0.01", "--format", fmt
        )
        print(
            f"\n[timing] {n} sites, {n_genes} genes, {n_rows} gnomAD rows, "
            f"--gnomad-max-af 0.01 {fmt}, {backend}: {took:.2f}s "
            f"(budget {budget:.1f}s)"
        )
        assert took < budget, (
            f"`db trio --gnomad-max-af 0.01 --format {fmt}` on {n} synthetic "
            f"sites ({gnomad} gnomAD, {backend}) took {took:.1f}s, budget "
            f"{budget:.1f}s. A gnomAD lookup per candidate shows up here."
        )
        if fmt == "json":
            doc = json.loads(result.stdout)
            counts = {m: doc["models"][m]["count"] for m in doc["models"]}
        else:  # `all` lists a count per model, comphet in genes
            counts = {
                line.split()[0]: int(line.split()[1])
                for line in result.stdout.splitlines()[1:]
            }
            # no chrX in this store: the text has no xlinked line, the JSON a 0
            counts.setdefault("xlinked", base_counts.get("xlinked", 0))
        # the filter did what the scenario says: nothing without gnomAD, and
        # with the slice the sites it calls common (one in four) are gone
        if gnomad == "none":
            assert counts == base_counts
        else:
            assert 0.6 * base_counts["dominant"] < counts["dominant"]
            assert counts["dominant"] < 0.9 * base_counts["dominant"]
            assert counts["recessive"] < base_counts["recessive"]


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_preset_adds_a_bounded_cost_on_a_synthetic_wgs_store(
    backend, vcfclick_home, tmp_path, monkeypatch
):
    """`--preset` adds the effective-gate report: one pass over the three
    samples' rows, one count over the variants and, for the GQ / DP / AD gates,
    one grouped count per selected model over its candidate join (three
    distinct joins for all four models, compound het sharing dominant's), and no
    query per candidate. The candidate counts alone take 31 ms on DuckDB and
    145 ms on chDB at 200,000 sites (`db trio`'s own counts of the same models
    take 44 and 159 ms), so the preset call is about 0.07 s (DuckDB) and 0.2 s
    (chDB) over the plain one. Timed against the same call without a preset on
    the same store, and against the full-scan baseline the other gates here
    use."""
    n, _ = _wgs_store("population", backend, vcfclick_home, tmp_path, monkeypatch)
    budget = _synthetic_budget(n, _baseline(vcfclick_home, backend))
    trio = ("db", "trio", "t", "--proband", "CHILD", "--format", "json")

    plain, _ = _best_of_two(vcfclick_home, backend, *trio)
    took, result = _best_of_two(
        vcfclick_home, backend, *trio, "--preset", "broad-first-pass"
    )
    doc = json.loads(result.stdout)
    print(
        f"\n[timing] {n} sites, --preset broad-first-pass, {backend}: {took:.2f}s "
        f"(plain {plain:.2f}s, budget {budget:.1f}s)"
    )
    # the scenario is what it says: the report read every row of the trio
    assert doc["preset"]["name"] == "broad-first-pass"
    assert doc["effective"]["rows"] == {"proband": n, "father": n, "mother": n}
    assert took < budget, (
        f"`db trio --preset` on {n} synthetic sites ({backend}) took {took:.1f}s, "
        f"budget {budget:.1f}s. The effective-gate report is one pass over the "
        "samples' rows; a query per candidate or a second pass shows up here."
    )
    assert took < 3 * plain + 1.0, (took, plain)


def _bed(path, n_intervals: int, n_sites: int):
    """`n_intervals` BED intervals spread evenly over the synthetic sites
    (chr1, 1000 + i * SITE_STEP), each covering half of its stretch."""
    span = n_sites * SITE_STEP
    step = span // n_intervals
    lines = [
        f"chr1\t{1000 + k * step - 1}\t{1000 + k * step + step // 2}\n"
        for k in range(n_intervals)
    ]
    path.write_text("".join(lines))
    return path


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_scope_of_500_intervals_on_a_synthetic_wgs_store(
    backend, vcfclick_home, tmp_path, monkeypatch
):
    """The scope's SQL pushdown at each of its limits (cli/scope.py):
    SCOPE_SQL_EXACT_MAX (64) merged intervals as exact BETWEENs,
    SCOPE_SQL_MAX_INTERVALS (500) as the bucket superset, and one interval
    more, where the SQL carries no scope and the Python backstop alone filters
    the rows; each against the same call without a scope. On the cohort-AF
    store, where nearly every inherited het is a candidate, so the backstop has
    the most rows to filter. Each interval covers half of its stretch, so the
    scoped runs keep about half the candidates.

    Measured, `db trio --format json`, best of two, Apple silicon:

                       plain   64 (sql)  500 (sql_buckets)  501 (python)
      20,000 DuckDB    0.31 s   0.31 s       0.33 s            0.33 s
             chDB      0.41 s   0.52 s       0.44 s            0.43 s
      200,000 DuckDB   1.02 s   0.98 s       1.09 s            1.24 s
              chDB     0.96 s   1.13 s       1.09 s            1.18 s

    500 intervals as BETWEENs (the first form) cost 2.7-2.9 s on chDB at
    200,000 sites, 1.7 s over plain: chDB evaluates every BETWEEN on every row
    of each of the three joined sides. That is why 65-500 intervals go as a
    bucket set, an IN over integers that costs nothing measurable.

    With `--preset broad-first-pass` the effective-gate report counts the rows
    in scope exactly in every mode (outside exact SQL: the bucket superset,
    and per-site rows only in edge buckets, trimmed in Python). 200,000 sites:

                       plain preset  64 (sql)  500 (sql_buckets)  501 (python)
              DuckDB      1.10 s      1.09 s       1.58 s            1.73 s
              chDB        1.21 s      1.46 s       1.52 s            1.62 s
    """
    from cli.scope import SCOPE_SQL_EXACT_MAX, SCOPE_SQL_MAX_INTERVALS

    n, _ = _wgs_store("cohort", backend, vcfclick_home, tmp_path, monkeypatch)
    budget = _synthetic_budget(n, _baseline(vcfclick_home, backend))
    trio = ("db", "trio", "t", "--proband", "CHILD", "--format", "json")
    plain_s, plain = _best_of_two(vcfclick_home, backend, *trio)
    base = json.loads(plain.stdout)["models"]
    preset_s, _ = _best_of_two(
        vcfclick_home, backend, *trio,
        "--preset", "broad-first-pass", "--preset-ignore", "gnomad_max_af",
    )  # fmt: skip
    for n_intervals, applied in (
        (SCOPE_SQL_EXACT_MAX, "sql"),
        (SCOPE_SQL_MAX_INTERVALS, "sql_buckets"),
        (SCOPE_SQL_MAX_INTERVALS + 1, "python"),
    ):
        bed = _bed(tmp_path / f"r{n_intervals}.bed", n_intervals, n)
        took, result = _best_of_two(
            vcfclick_home, backend, *trio, "--region-file", str(bed)
        )
        doc = json.loads(result.stdout)
        print(
            f"\n[timing] {n} sites, {n_intervals} intervals ({applied}), {backend}: "
            f"{took:.2f}s (plain {plain_s:.2f}s, budget {budget:.1f}s)"
        )
        assert doc["scope"]["applied"] == applied
        assert doc["scope"]["intervals"] == n_intervals
        # the scenario is what it says: about half the candidates are kept
        kept = doc["models"]["dominant"]["count"] / base["dominant"]["count"]
        assert 0.3 < kept < 0.7, kept
        assert took < budget, (
            f"`db trio --region-file` with {n_intervals} intervals on {n} sites "
            f"({backend}) took {took:.1f}s, budget {budget:.1f}s"
        )
        assert took < 1.5 * plain_s + 0.5, (applied, took, plain_s)
        # with a preset, the effective-gate report counts the rows in scope
        # exactly in every mode: outside exact SQL its statements carry the
        # bucket superset and group by site, and Python trims them
        preset = ("--preset", "broad-first-pass", "--preset-ignore", "gnomad_max_af")
        p_took, p_result = _best_of_two(
            vcfclick_home, backend, *trio, *preset, "--region-file", str(bed)
        )
        print(
            f"[timing] {n} sites, {n_intervals} intervals ({applied}) with a "
            f"preset, {backend}: {p_took:.2f}s (plain preset {preset_s:.2f}s)"
        )
        assert json.loads(p_result.stdout)["scope"]["applied"] == applied
        assert p_took < budget
        assert p_took < 1.5 * preset_s + 0.5, (applied, p_took, preset_s)
