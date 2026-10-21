"""The AF rarity gate treats a site AT --max-af as rare, on both backends.

`variants.info_AF` is Float32 on chDB and REAL on DuckDB, and the gate
compares it with a double: Float32(0.001) is 0.0010000000474..., which is
greater than the double 0.001, so on chDB `info_AF <= 0.001` was false for a
site whose AF is exactly 0.001 while DuckDB kept it. The gate adds 1e-7
(storage.sql.AF_EPSILON) to the threshold, far below any real frequency
step and above float32 rounding at every AF a VCF carries, so both engines
agree. This is the one deliberate change to the legacy SQL text, made in its
own commit.
"""

from __future__ import annotations

import json

import pytest

from tests.conftest import run_cli
from tests.fixtures._vcfgen import (
    TRIO_PED,
    TRIO_SAMPLES,
    Site,
    write_vcf,
)

BACKENDS = ("chdb", "duckdb")


def _denovo(pos: int, af: str) -> Site:
    return Site(
        "chr1",
        pos,
        "A",
        "G",
        ("0/1:15,15:99:30", "0/0:30,0:99:30", "0/0:30,0:99:30"),
        info=f"AF={af}",
        fmt="GT:AD:GQ:DP",
    )


# position -> INFO/AF of a de novo site
SITES = {
    1000: "0.0005",  # clearly rare
    2000: "0.001",  # exactly at --max-af 0.001 (Float32 rounds UP)
    3000: "0.01",  # exactly at --max-af 0.01 (Float32 rounds down)
    4000: "0.002",  # clearly above 0.001
    5000: "0.02",  # clearly above 0.01
    6000: "0.3",  # exactly at --max-af 0.3
}


@pytest.fixture
def af_store(vcfclick_home, tmp_path):
    def build(backend):
        sites = [_denovo(pos, af) for pos, af in SITES.items()]
        vcf = write_vcf(tmp_path / "af.vcf", TRIO_SAMPLES, sites)
        ped = tmp_path / "trio.ped"
        ped.write_text(TRIO_PED)
        run_cli(vcfclick_home, backend, "db", "create", "t")
        run_cli(
            vcfclick_home, backend, "db", "ingest", "t", str(vcf),
            "--ingest-id", "i1", "--serial", "--keep-reference",
        )  # fmt: skip
        run_cli(vcfclick_home, backend, "db", "ped", "t", str(ped))

    return build


def _positions(home, backend, max_af: str) -> set[int]:
    r = run_cli(
        home, backend, "db", "trio", "t", "--proband", "CHILD", "--category", "denovo",
        "--max-af", max_af, "--format", "json",
    )  # fmt: skip
    doc = json.loads(r.stdout)
    return {c["pos"] for c in doc["models"]["denovo"]["candidates"]}


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_site_exactly_at_max_af_is_kept_on_both_backends(
    backend, vcfclick_home, af_store
):
    af_store(backend)
    assert _positions(vcfclick_home, backend, "0.001") == {1000, 2000}
    assert _positions(vcfclick_home, backend, "0.01") == {1000, 2000, 3000, 4000}
    assert _positions(vcfclick_home, backend, "0.3") == set(SITES)
    # still a gate: nothing above the threshold gets through
    assert _positions(vcfclick_home, backend, "0.0004") == set()
    assert _positions(vcfclick_home, backend, "0.0009") == {1000}


def test_the_gate_text_carries_the_epsilon(monkeypatch):
    from cli.trio_analysis import Gates, Trio, trio_sql
    from storage.sql import AF_EPSILON, af_gate_sql

    assert AF_EPSILON == 1e-7
    assert af_gate_sql("v.info_AF", 0.01) == (
        "(v.info_AF IS NULL OR v.info_AF <= 0.01 + 1e-07)"
    )
    sql = trio_sql(
        "denovo", Trio("i", "C", "F", "M"), Gates(20, 10, 0.001, 0.25, 0.75),
        count_only=False,
    )  # fmt: skip
    assert "(v.info_AF IS NULL OR v.info_AF <= 0.001 + 1e-07)" in sql
    # NULL AF is still not filtered
    assert "v.info_AF IS NULL" in sql
