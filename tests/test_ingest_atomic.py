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


def _vc(
    home: Path, *args: str, expect_failure: bool = False
) -> subprocess.CompletedProcess:
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
    """Return TabSeparated output so the value is just the cell content,
    not a chDB or DuckDB box-drawing rendering. Assertions then read
    the integer directly rather than matching engine-specific glyphs."""
    if "FORMAT " not in sql.upper():
        sql = sql + " FORMAT TabSeparated"
    return _vc(home, "db", "query", db, sql).stdout


def _scalar(home: Path, db: str, sql: str) -> str:
    """Run a SELECT that returns a single cell. Returns the stripped value."""
    return _query(home, db, sql).strip()


def test_multi_allelic_ingest_fails_with_helpful_error(vcfclick_home):
    _vc(vcfclick_home, "db", "create", "smoke")
    r = _vc(
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(MULTI_VCF),
        "--cohort",
        "demo",
        "--ingest-id",
        "batch_a",
        "--serial",
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
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(MULTI_VCF),
        "--cohort",
        "demo",
        "--ingest-id",
        "batch_a",
        "--serial",
        expect_failure=True,
    )

    for table in ("variants", "genotypes", "samples", "ingestions"):
        out = _query(
            vcfclick_home,
            "smoke",
            f"SELECT count() FROM {table} WHERE ingest_id = 'batch_a'",
        )
        # Output is a boxed count; "0" must be the only number on the line.
        assert out.strip() == "0", f"rollback did not scrub {table}: {out!r}"


def test_successful_ingest_after_failed_one(vcfclick_home, tiny_vcf):
    """After a failed ingest, a subsequent good ingest under a different
    ingest_id must work and produce the expected row counts. Proves the
    rollback didn't damage the schema or the session."""
    _vc(vcfclick_home, "db", "create", "smoke")
    # Failed ingest #1
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(MULTI_VCF),
        "--cohort",
        "demo",
        "--ingest-id",
        "batch_bad",
        "--serial",
        expect_failure=True,
    )
    # Good ingest #2
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(tiny_vcf),
        "--cohort",
        "demo",
        "--ingest-id",
        "batch_good",
        "--serial",
    )

    out = _query(vcfclick_home, "smoke", "SELECT count() FROM variants")
    assert out.strip() == "5"


def test_reingest_same_id_truly_replaces_prior_data(vcfclick_home):
    """Re-running ingest under the same `ingest_id` MUST delete any
    rows from the previous run that aren't in the new VCF — not just
    upsert via ReplacingMergeTree dedup on the sorting key.

    Fixtures:
      tiny.vcf.gz    — 5 variants at chr1:{100,250,500,750,900},
                       3 samples (S1, S2, S3)
      routing.vcf.gz — 2 variants at chr1:{100,200},
                       2 samples (S1, S2)
    Overlap at chr1:100 only.

    After re-ingesting routing.vcf.gz under the same ingest_id:
      - The 4 tiny-only variants (250, 500, 750, 900) must be GONE.
      - Sample S3 (only in tiny) must be GONE.
      - The 2 routing variants must be present.
    """
    from pathlib import Path as _Path

    fixtures = _Path(__file__).parent / "fixtures"
    tiny = fixtures / "tiny.vcf.gz"
    routing = fixtures / "routing.vcf.gz"

    _vc(vcfclick_home, "db", "create", "smoke")

    # First ingest — tiny (5 variants, 3 samples)
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(tiny),
        "--cohort",
        "x",
        "--ingest-id",
        "batch_a",
        "--serial",
    )
    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM variants WHERE ingest_id = 'batch_a'",
    )
    assert out.strip() == "5", f"first ingest should land 5 variants, got: {out!r}"

    # Re-ingest under the SAME id — routing (2 variants, 2 samples)
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(routing),
        "--cohort",
        "x",
        "--ingest-id",
        "batch_a",
        "--serial",
    )

    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM variants WHERE ingest_id = 'batch_a'",
    )
    assert out.strip() == "2", (
        f"re-ingest should leave exactly 2 variants (replacement, not upsert), "
        f"got: {out!r}"
    )

    # Tiny-only positions are gone
    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM variants WHERE ingest_id = 'batch_a' "
        "AND pos IN (250, 500, 750, 900)",
    )
    assert out.strip() == "0", f"tiny-only positions should be deleted, got: {out!r}"

    # Routing variants present
    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM variants WHERE ingest_id = 'batch_a' "
        "AND pos IN (100, 200)",
    )
    assert out.strip() == "2", f"routing variants missing, got: {out!r}"

    # Sample S3 (tiny-only) is gone
    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM samples WHERE ingest_id = 'batch_a' AND sample_id = 'S3'",
    )
    assert out.strip() == "0", f"S3 should be gone after replace, got: {out!r}"


def test_reingest_with_corrupt_vcf_preserves_prior_data(vcfclick_home, tmp_path):
    """If the new VCF fails to open (corrupt header / unreadable file),
    the re-ingest under an existing ingest_id MUST leave the prior data
    intact. This is the contract codex flagged in the second review:
    the previous fix called rollback BEFORE opening the new VCF, which
    silently wiped good data on every failed re-ingest. Fix moved the
    rollback inside the try block so a bad header raises before the
    delete fires.
    """
    fixtures = Path(__file__).parent / "fixtures"
    tiny = fixtures / "tiny.vcf.gz"

    _vc(vcfclick_home, "db", "create", "smoke")

    # First ingest — tiny (5 variants, 3 samples).
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(tiny),
        "--cohort",
        "x",
        "--ingest-id",
        "batch_a",
        "--serial",
    )
    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM variants WHERE ingest_id = 'batch_a'",
    )
    assert out.strip() == "5"

    # Build a "corrupt" VCF — valid bgzip header but garbage content.
    # cyvcf2 will fail when constructing the VCF reader.
    corrupt = tmp_path / "corrupt.vcf.gz"
    # Plain non-bgzip bytes; cyvcf2 opens via htslib which expects BGZF.
    corrupt.write_bytes(b"this is not a VCF\n")

    # Re-ingest under the SAME id with the corrupt VCF — must FAIL.
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(corrupt),
        "--cohort",
        "x",
        "--ingest-id",
        "batch_a",
        "--serial",
        expect_failure=True,
    )

    # Prior data MUST still be queryable. The fix moved rollback_ingest
    # inside the try block so a bad-header VCF raises before the delete
    # fires.
    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM variants WHERE ingest_id = 'batch_a'",
    )
    assert out.strip() == "5", (
        f"failed re-ingest with corrupt VCF wiped prior data — "
        f"replacement should be atomic against bad headers: {out!r}"
    )


def test_reingest_with_multi_allelic_preserves_prior_data(vcfclick_home):
    """The hard case: re-ingest a VCF that opens + classifies fine but
    fails MID-STREAM (multi-allelic record in the body) under an
    existing ingest_id MUST leave the prior data intact. Closes the
    gap codex flagged in the third review pass.

    The stage-then-commit restructure means the variant loop only
    writes Parquet files to a tempdir — no chDB writes happen until
    the full VCF parses successfully. A multi-allelic record raises
    in Phase 1, before the rollback or any chDB writes fire.
    """
    fixtures = Path(__file__).parent / "fixtures"
    tiny = fixtures / "tiny.vcf.gz"
    multi = fixtures / "multiallelic.vcf.gz"

    _vc(vcfclick_home, "db", "create", "smoke")

    # First ingest — tiny (5 variants, 3 samples).
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(tiny),
        "--cohort",
        "x",
        "--ingest-id",
        "batch_a",
        "--serial",
    )
    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM variants WHERE ingest_id = 'batch_a'",
    )
    assert out.strip() == "5"

    # Re-ingest under same id with the multi-allelic fixture. The fixture
    # has a valid header and three records: two bi-allelic flanking one
    # multi-allelic in the middle. cyvcf2 opens it fine, classify_header
    # succeeds, and the variant loop fires the multi-allelic check on
    # record 2. Under stage-then-commit, that raise happens BEFORE any
    # chDB write — prior data is untouched.
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "smoke",
        str(multi),
        "--cohort",
        "x",
        "--ingest-id",
        "batch_a",
        "--serial",
        expect_failure=True,
    )

    # Prior data MUST still be queryable.
    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM variants WHERE ingest_id = 'batch_a'",
    )
    assert out.strip() == "5", (
        f"mid-stream-failed re-ingest wiped prior data — stage-then-commit "
        f"should preserve prior rows when Phase 1 (parse) raises: {out!r}"
    )

    # And the samples — S3 from the tiny ingest is only in the tiny VCF,
    # so its presence confirms the previous samples row didn't get
    # scrubbed by a premature rollback.
    out = _query(
        vcfclick_home,
        "smoke",
        "SELECT count() FROM samples WHERE ingest_id = 'batch_a' AND sample_id = 'S3'",
    )
    assert (
        out.strip() == "1"
    ), f"sample S3 from prior ingest should survive failed re-ingest: {out!r}"


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
