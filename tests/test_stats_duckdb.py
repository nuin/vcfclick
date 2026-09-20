"""`vcfclick db stats` on the DuckDB backend.

The original implementation was chDB-only (system.columns, countIf, ARRAY JOIN
mapKeys). These tests pin the DuckDB port against the same routing fixture the
chDB suite uses, so both backends report the same numbers.

DuckDB is the backend bioconda ships, so `db stats` failing there is a dead
command for every conda user.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from storage.db import map_keys_from, populated_expr, typed_columns_sql

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
ROUTING_VCF = Path(__file__).parent / "fixtures" / "routing.vcf.gz"


def _vc(home: Path, *args: str) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env["VCFCLICK_BACKEND"] = "duckdb"  # the point of this suite
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    assert r.returncode == 0, (
        f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
        f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    )
    return r


def _stats(home: Path) -> str:
    _vc(home, "db", "create", "demo")
    _vc(
        home,
        "db",
        "ingest",
        "demo",
        str(ROUTING_VCF),
        "--cohort",
        "acme",
        "--ingest-id",
        "batch_a",
        "--serial",
    )
    return _vc(home, "db", "stats", "demo").stdout


def _count(out: str, label: str) -> int:
    line = next(ln for ln in out.splitlines() if ln.strip().startswith(label))
    return int(line.strip().split()[-1].replace(",", ""))


def test_stats_runs_on_duckdb_with_row_counts(vcfclick_home):
    out = _stats(vcfclick_home)
    assert _count(out, "variants") == 2
    assert _count(out, "genotypes") == 3
    assert _count(out, "samples") == 2
    assert _count(out, "ingestions") == 1


def test_stats_duckdb_cohorts_and_contigs(vcfclick_home):
    out = _stats(vcfclick_home)
    assert "cohorts:" in out and "acme" in out
    assert "contigs:" in out and "chr1" in out


def test_stats_duckdb_typed_population_and_overflow_keys(vcfclick_home):
    out = _stats(vcfclick_home)
    # Typed INFO population section, and the overflow Map keys (the fixture
    # carries lab-specific INFO tags + a MYCUSTOM FORMAT field).
    assert "typed INFO column population" in out
    assert "variants.info_extra" in out
    assert "genotypes.format_extra" in out
    # info_DP is populated on one of the two variants -> counted, not 100%.
    dp = next((ln for ln in out.splitlines() if "info_DP" in ln), "")
    assert dp, "expected info_DP in the typed population table"


def test_stats_duckdb_flag_columns_not_counted_as_all_rows(vcfclick_home):
    """Flags are `UTINYINT DEFAULT 0` in DuckDB (nullable, defaulted), so a
    naive IS NOT NULL would report every row as populated. They must be
    counted as `!= 0`, matching chDB's non-Nullable handling."""
    out = _stats(vcfclick_home)
    total = _count(out, "variants")
    for line in out.splitlines():
        if "info_SOMATIC" in line or "info_H2" in line:
            n = int(line.split()[1].replace(",", ""))
            assert n < total, f"flag column counted as all rows: {line!r}"


# --- dialect SQL generation (both backends, no live engine needed) ---------
#
# chDB's binary can't always be loaded locally, so pin the *generated SQL* for
# both dialects here. This is what protects the chDB path from the port.


@pytest.mark.parametrize("be", ["chdb", "duckdb"])
def test_typed_columns_sql_per_dialect(monkeypatch, be):
    monkeypatch.setenv("VCFCLICK_BACKEND", be)
    sql = typed_columns_sql("variants")
    if be == "duckdb":
        assert "information_schema.columns" in sql
        assert "column_default = '0'" in sql  # flag detection
    else:
        assert "system.columns" in sql and "currentDatabase()" in sql
        assert "Nullable%" in sql  # non-Nullable => flag


@pytest.mark.parametrize("be", ["chdb", "duckdb"])
def test_populated_expr_per_dialect(monkeypatch, be):
    monkeypatch.setenv("VCFCLICK_BACKEND", be)
    nullable = populated_expr("qual", False)
    flag = populated_expr("info_H2", True)
    if be == "duckdb":
        assert nullable == 'count(*) FILTER (WHERE "qual" IS NOT NULL) AS "qual"'
        assert flag == 'count(*) FILTER (WHERE "info_H2" != 0) AS "info_H2"'
    else:
        # byte-identical to the pre-port chDB expressions
        assert nullable == "countIf(`qual` IS NOT NULL) AS `qual`"
        assert flag == "countIf(`info_H2` != 0) AS `info_H2`"


@pytest.mark.parametrize("be", ["chdb", "duckdb"])
def test_map_keys_from_per_dialect(monkeypatch, be):
    monkeypatch.setenv("VCFCLICK_BACKEND", be)
    src = map_keys_from("variants", "info_extra")
    if be == "duckdb":
        assert src == "(SELECT unnest(map_keys(info_extra)) AS k FROM variants)"
    else:
        assert src == "variants ARRAY JOIN mapKeys(info_extra) AS k"


def test_dialect_helpers_reject_unsafe_identifiers(monkeypatch):
    monkeypatch.setenv("VCFCLICK_BACKEND", "duckdb")
    with pytest.raises(ValueError):
        typed_columns_sql("variants; DROP TABLE x")
    with pytest.raises(ValueError):
        map_keys_from("variants", "info_extra); --")
