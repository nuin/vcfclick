"""One string-literal quoting rule per backend, end to end.

chDB reads two escapes inside a literal (`\\` and `''`), DuckDB reads one
(`''`). A helper that doubles backslashes on both builds, on DuckDB, a
literal that is not the stored value: `'DOM\\\\kid'` (8 chars) never equals
the 7-char sample id, so the sample vanishes from every pivot, join and
export. These tests ingest sample ids and cohort names that carry a
backslash or a quote and require every quoting path to find them on both
engines.
"""

from __future__ import annotations

import json
import re

import pytest

from tests.conftest import FIXTURES, bgzip_vcf, run_cli
from tests.fixtures._vcfgen import Site, write_vcf

BACKENDS = ("chdb", "duckdb")

AWKWARD = ("DOM\\kid", "O'Brien", "plain")

# Site: (pos, GTs for AWKWARD in order). The het calls are what an export
# must reproduce; a sample the quoting loses reads 0/0 everywhere.
SITES = [
    (100, ("0/1", "0/0", "1/1")),
    (200, ("1/1", "0/1", "0/0")),
    (300, ("0/0", "1/1", "0/1")),
]


def _awkward_vcf(tmp_path):
    sites = [Site("chr1", pos, "A", "G", gts) for pos, gts in SITES]
    plain = write_vcf(tmp_path / "awkward.vcf", list(AWKWARD), sites)
    return bgzip_vcf(plain, tmp_path / "awkward.vcf.gz")


def _ingest(home, backend, vcf, name="t", **kw):
    run_cli(home, backend, "db", "create", name)
    args = ["db", "ingest", name, str(vcf), "--ingest-id", kw.get("iid", "i1")]
    args += ["--serial", "--keep-reference"]
    if "cohort" in kw:
        args += ["--cohort", kw["cohort"]]
    run_cli(home, backend, *args)


def _query(home, backend, sql, name="t"):
    r = run_cli(home, backend, "db", "query", name, sql, "--format", "JSONCompact")
    return json.loads(r.stdout)["data"]


@pytest.mark.parametrize("backend", BACKENDS)
def test_awkward_sample_ids_are_found_by_every_quoting_path(
    backend, vcfclick_home, tmp_path, monkeypatch
):
    """genotype_pivot_sql, sample_join and str_list quote through
    sql_quote_str; the literal must equal the stored id on this backend."""
    monkeypatch.setenv("VCFCLICK_BACKEND", backend)
    from storage.sql import genotype_pivot_sql, sample_join, str_list

    _ingest(vcfclick_home, backend, _awkward_vcf(tmp_path))

    stored = _query(
        vcfclick_home,
        backend,
        "SELECT sample_id, count(*) FROM genotypes GROUP BY sample_id ORDER BY sample_id",
    )
    assert {s: int(n) for s, n in stored} == {s: 3 for s in AWKWARD}

    roles = {"a": AWKWARD[0], "b": AWKWARD[1], "c": AWKWARD[2]}
    pivot = genotype_pivot_sql("i1", roles)
    counts = _query(
        vcfclick_home,
        backend,
        f"SELECT count(a_gt), count(b_gt), count(c_gt) FROM ({pivot}) x",
    )
    assert [int(n) for n in counts[0]] == [3, 3, 3]

    joined = _query(
        vcfclick_home,
        backend,
        "SELECT count(*) FROM genotypes g "
        + sample_join("f", AWKWARD[0])
        + " "
        + sample_join("m", AWKWARD[1])
        + " WHERE g.sample_id = "
        + str_list([AWKWARD[2]]),
    )
    assert int(joined[0][0]) == 3

    listed = _query(
        vcfclick_home,
        backend,
        f"SELECT count(*) FROM genotypes WHERE sample_id IN ({str_list(AWKWARD[:2])})",
    )
    assert int(listed[0][0]) == 6


@pytest.mark.parametrize("backend", BACKENDS)
def test_export_keeps_the_calls_of_an_awkward_sample(backend, vcfclick_home, tmp_path):
    """`db export --samples` filters with the sample id in SQL; a lost
    match writes 0/0 at every site instead of the stored het/hom calls."""
    _ingest(vcfclick_home, backend, _awkward_vcf(tmp_path))
    for i, sample in enumerate(AWKWARD[:2]):
        out = tmp_path / f"out{i}.vcf"
        run_cli(
            vcfclick_home,
            backend,
            "db",
            "export",
            "t",
            "-o",
            str(out),
            "--samples",
            sample,
        )
        lines = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
        got = [ln.split("\t")[9].split(":")[0] for ln in lines]
        assert got == [gts[i] for _, gts in SITES], sample


@pytest.mark.parametrize("backend", BACKENDS)
def test_db_diff_finds_a_cohort_with_a_backslash(backend, vcfclick_home):
    tiny = FIXTURES / "tiny.vcf.gz"
    run_cli(vcfclick_home, backend, "db", "create", "d")
    for iid, cohort in (("x", "c\\ne"), ("y", "A")):
        run_cli(
            vcfclick_home,
            backend,
            "db",
            "ingest",
            "d",
            str(tiny),
            "--ingest-id",
            iid,
            "--cohort",
            cohort,
            "--serial",
        )
    r = run_cli(
        vcfclick_home,
        backend,
        "db",
        "diff",
        "d",
        "--cohort-a",
        "c\\ne",
        "--cohort-b",
        "A",
        "--format",
        "JSONCompact",
    )
    assert json.loads(r.stdout)["rows"] == 5

    # A name ending in a backslash is "unknown", never a SQL syntax error.
    r = run_cli(
        vcfclick_home,
        backend,
        "db",
        "diff",
        "d",
        "--cohort-a",
        "a\\",
        "--cohort-b",
        "A",
        ok=False,
    )
    assert re.search(r"unknown cohort\(s\): a\\\. Available: A", r.stderr), r.stderr
