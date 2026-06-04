"""Atomic-ingest tests.

An ingest that fails mid-stream — multi-allelic record, malformed row,
batch flush error, Ctrl-C — must leave the database in the state it
was in before the call. No half-loaded samples / variants / genotypes
that would force the user to `vcfclick db rm` and start over.

The fixture (`tests/fixtures/multiallelic.vcf.gz`) has 3 records: two
bi-allelic flanking one multi-allelic in the middle. The serial
loader's existing in-stream check fires on the middle record. Without
rollback, the first bi-allelic record (and the samples row) would be
left in the DB; with rollback, the DB is clean.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
MULTI_VCF = Path(__file__).parent / "fixtures" / "multiallelic.vcf.gz"


def _vc(home: Path, *args: str, expect_failure: bool = False) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
    )
    if expect_failure:
        assert r.returncode != 0, f"expected failure but got rc=0:\n{r.stdout}"
    else:
        assert r.returncode == 0, (
            f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
            f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
        )
    return r


def _query(home: Path, db: str, sql: str) -> str:
    return _vc(home, "db", "query", db, sql).stdout


def test_multi_allelic_ingest_fails_with_helpful_error(vcfclick_home):
    _vc(vcfclick_home, "db", "create", "smoke")
    r = _vc(
        vcfclick_home, "db", "ingest", "smoke", str(MULTI_VCF),
        "--cohort", "demo", "--ingest-id", "batch_a", "--serial",
        expect_failure=True,
    )
    # The in-stream check's error message must surface the bcftools fix.
    assert "Multi-allelic" in r.stderr
    assert "bcftools norm" in r.stderr


def test_failed_ingest_leaves_zero_rows(vcfclick_home):
    """The atomic guarantee: after a failed ingest, all four ingestion-
    scoped tables (variants, genotypes, samples, ingestions) report 0
    rows under the failed ingest_id."""
    _vc(vcfclick_home, "db", "create", "smoke")
    _vc(
        vcfclick_home, "db", "ingest", "smoke", str(MULTI_VCF),
        "--cohort", "demo", "--ingest-id", "batch_a", "--serial",
        expect_failure=True,
    )

    for table in ("variants", "genotypes", "samples", "ingestions"):
        out = _query(
            vcfclick_home, "smoke",
            f"SELECT count() FROM {table} WHERE ingest_id = 'batch_a'",
        )
        # Output is a boxed count; "0" must be the only number on the line.
        assert "│       0 │" in out or "│ 0 │" in out, (
            f"rollback did not scrub {table}: {out!r}"
        )


def test_successful_ingest_after_failed_one(vcfclick_home, tiny_vcf):
    """After a failed ingest, a subsequent good ingest under a different
    ingest_id must work and produce the expected row counts. Proves the
    rollback didn't damage the schema or the session."""
    _vc(vcfclick_home, "db", "create", "smoke")
    # Failed ingest #1
    _vc(
        vcfclick_home, "db", "ingest", "smoke", str(MULTI_VCF),
        "--cohort", "demo", "--ingest-id", "batch_bad", "--serial",
        expect_failure=True,
    )
    # Good ingest #2
    _vc(
        vcfclick_home, "db", "ingest", "smoke", str(tiny_vcf),
        "--cohort", "demo", "--ingest-id", "batch_good", "--serial",
    )

    out = _query(vcfclick_home, "smoke", "SELECT count() FROM variants")
    assert "5" in out


def test_ingest_id_rejected_with_quotes():
    """Direct library-level test: rollback_ingest interpolates the
    ingest_id into the DELETE statement, so validate_ingest_id is the
    safety wall — quotes/semicolons/spaces are rejected."""
    from storage import validate_ingest_id

    # Allowed forms
    for ok in ("batch_a", "2026.05.31", "a-b-c", "uuid_1234", "X"):
        validate_ingest_id(ok)

    # Rejected forms — would break the DELETE if interpolated.
    for bad in ("a' OR 1=1; --", "spaces here", "semi;colon", "back`tick", ""):
        with pytest.raises(ValueError, match="invalid ingest_id"):
            validate_ingest_id(bad)
