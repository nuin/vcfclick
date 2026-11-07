"""Commands share a database and fail without a traceback.

DuckDB gives one process the database file at a time. The desktop apps run
`db qc` and `db relatedness` on one database at once, and the second used
to die on "Could not set lock on file", so its panel came up empty. Opening
a DuckDB database now waits for the lock (VCFCLICK_DUCKDB_LOCK_WAIT), and a
wait that runs out ends the command with one `Error:` line. `db query` also
prints an engine's SQL error as one `Error:` line on both backends.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from tests.conftest import REPO, run_cli

BACKENDS = ("chdb", "duckdb")


def _cohort(home, backend, tiny_vcf, name="lk"):
    run_cli(home, backend, "db", "create", name)
    run_cli(home, backend, "db", "ingest", name, str(tiny_vcf))
    return name


def _hold_lock(home, name, seconds):
    """A process that holds the DuckDB file of `name` open for `seconds`."""
    db_file = home / "dbs" / name / "store.duckdb"
    code = (
        "import duckdb, sys, time\n"
        f"c = duckdb.connect({str(db_file)!r})\n"
        "print('held', flush=True)\n"
        f"time.sleep({seconds})\n"
    )
    p = subprocess.Popen(
        [sys.executable, "-c", code], cwd=REPO, stdout=subprocess.PIPE, text=True
    )
    assert p.stdout.readline().strip() == "held"
    return p


@pytest.mark.parametrize("backend", BACKENDS)
def test_query_sql_error_is_one_line(vcfclick_home, backend, tiny_vcf):
    name = _cohort(vcfclick_home, backend, tiny_vcf)
    r = run_cli(
        vcfclick_home,
        backend,
        "db",
        "query",
        name,
        "SELECT nope FROM variants",
        ok=False,
    )
    assert "Traceback" not in r.stderr
    # chDB on Linux prints a jemalloc line ("Default dirty pages decay
    # period") to stderr before the error, so look for the Error line
    errors = [ln for ln in r.stderr.splitlines() if ln.startswith("Error: ")]
    assert len(errors) == 1 and "nope" in errors[0], r.stderr


def test_duckdb_command_waits_for_the_lock(vcfclick_home, tiny_vcf):
    name = _cohort(vcfclick_home, "duckdb", tiny_vcf)
    holder = _hold_lock(vcfclick_home, name, 2)
    try:
        r = run_cli(vcfclick_home, "duckdb", "db", "qc", name, "--format", "json")
    finally:
        holder.wait(timeout=30)
    assert '"sample_id"' in r.stdout


def test_duckdb_lock_wait_runs_out_cleanly(vcfclick_home, tiny_vcf, monkeypatch):
    name = _cohort(vcfclick_home, "duckdb", tiny_vcf)
    monkeypatch.setenv("VCFCLICK_DUCKDB_LOCK_WAIT", "0.5")
    holder = _hold_lock(vcfclick_home, name, 20)
    try:
        r = run_cli(
            vcfclick_home, "duckdb", "db", "qc", name, "--format", "json", ok=False
        )
    finally:
        holder.kill()
        holder.wait(timeout=30)
    assert "Traceback" not in r.stderr
    assert "Error: the database" in r.stderr
    assert "is in use by process" in r.stderr
    assert "VCFCLICK_DUCKDB_LOCK_WAIT" in r.stderr
