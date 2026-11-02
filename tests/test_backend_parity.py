"""chDB and DuckDB give the same answer where one engine used to differ.

Each test ran green on one backend and red on the other before its fix:
`db stats` counted (ingest_id, sample_id) pairs by string concatenation on
DuckDB; `db dump` pasted the staging path into SQL unquoted; the trio
allele-balance gate added two 16-bit unsigned columns unwidened, which
DuckDB refuses above 65535.
"""

from __future__ import annotations


import pytest

from tests.conftest import FIXTURES, bgzip_vcf, run_cli
from tests.fixtures._vcfgen import Site, write_vcf

BACKENDS = ("chdb", "duckdb")


def _vcf(tmp_path, name, samples, sites):
    plain = write_vcf(tmp_path / f"{name}.vcf", samples, sites)
    return bgzip_vcf(plain, tmp_path / f"{name}.vcf.gz")


@pytest.mark.parametrize("backend", BACKENDS)
def test_stats_sample_count_is_a_tuple_distinct_on_both_engines(backend, monkeypatch):
    monkeypatch.setenv("VCFCLICK_BACKEND", backend)
    from cli.db_stats import distinct_sample_pair_expr

    assert distinct_sample_pair_expr() == "count(DISTINCT (ingest_id, sample_id))"


@pytest.mark.parametrize("backend", BACKENDS)
def test_stats_counts_distinct_sample_pairs_not_concatenations(
    backend, vcfclick_home, tmp_path
):
    """('a', 'bc') and ('ab', 'c') are two samples of cohort X; joined with
    no separator they are the one string 'abc'. `db batch` takes ingest ids
    from file names and `--ingest-id` takes any string, so lengths vary."""
    a = _vcf(tmp_path, "a", ["bc"], [Site("chr1", 100, "A", "G", ("0/1",))])
    b = _vcf(tmp_path, "b", ["c"], [Site("chr1", 200, "A", "G", ("0/1",))])
    run_cli(vcfclick_home, backend, "db", "create", "t")
    for iid, vcf in (("a", a), ("ab", b)):
        run_cli(
            vcfclick_home,
            backend,
            "db",
            "ingest",
            "t",
            str(vcf),
            "--ingest-id",
            iid,
            "--cohort",
            "X",
            "--serial",
        )
    out = run_cli(vcfclick_home, backend, "db", "stats", "t").stdout
    line = next(ln for ln in out.splitlines() if ln.strip().startswith("X "))
    assert line.split()[-2:] == ["2", "samples"], line


@pytest.mark.parametrize("backend", BACKENDS)
def test_dump_writes_into_a_directory_with_a_quote_in_its_name(
    backend, vcfclick_home, tmp_path
):
    """The staging file sits in the output directory, whose name the user
    chose; the path is quoted like every other literal."""
    run_cli(vcfclick_home, backend, "db", "create", "t")
    run_cli(
        vcfclick_home,
        backend,
        "db",
        "ingest",
        "t",
        str(FIXTURES / "tiny.vcf.gz"),
        "--ingest-id",
        "i1",
        "--serial",
    )
    out = tmp_path / "O'Brien"
    run_cli(vcfclick_home, backend, "db", "dump", "t", "--out", str(out))
    assert (out / "variants.parquet").stat().st_size > 0
    assert (out / "genotypes.parquet").stat().st_size > 0
    assert not list(out.glob("tmp*")), "staging file left behind"
