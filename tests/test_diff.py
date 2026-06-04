"""End-to-end tests for `vcfclick db diff`.

Uses the two existing fixtures as two cohorts:

  cohort `case`    ← tests/fixtures/tiny.vcf.gz     (5 variants, 3 samples)
  cohort `control` ← tests/fixtures/routing.vcf.gz  (2 variants, 2 samples)

Overlap at chr1:100 A>G — both cohorts call this variant — so the diff
output has at least one row with AF populated for both cohorts and
verifiable arithmetic. The fixtures' disjoint variants exercise the
"present in one cohort only" path (AF on the absent side reports 0,
not NULL).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
TINY = Path(__file__).parent / "fixtures" / "tiny.vcf.gz"
ROUTING = Path(__file__).parent / "fixtures" / "routing.vcf.gz"


def _vc(home: Path, *args: str, expect_failure: bool = False) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    if expect_failure:
        assert r.returncode != 0, f"expected failure but got rc=0:\n{r.stdout}"
    else:
        assert r.returncode == 0, (
            f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
            f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
        )
    return r


def _setup_two_cohort_db(home: Path, db: str = "demo") -> None:
    _vc(home, "db", "create", db)
    _vc(
        home, "db", "ingest", db, str(TINY),
        "--cohort", "case", "--ingest-id", "batch_case", "--serial",
    )
    _vc(
        home, "db", "ingest", db, str(ROUTING),
        "--cohort", "control", "--ingest-id", "batch_control", "--serial",
    )


def _tsv_rows(home: Path, db: str, *flags: str) -> list[list[str]]:
    """Run db diff in TSV format and return parsed rows."""
    out = _vc(home, "db", "diff", db,
              "--cohort-a", "case", "--cohort-b", "control",
              "--format", "TSV", *flags).stdout
    return [line.split("\t") for line in out.strip().splitlines() if line.strip()]


def test_diff_reports_correct_af_for_overlap_variant(vcfclick_home):
    """chr1:100 A>G is in both fixtures:
       case (tiny): S1=0/0 S2=0/1 S3=1/1 → 3 alt out of 6  → AF 0.5
       control (routing): S1=0/1 S2=1/1 → 3 alt out of 4   → AF 0.75
       diff = 0.5 - 0.75 = -0.25
    """
    _setup_two_cohort_db(vcfclick_home)
    rows = _tsv_rows(vcfclick_home, "demo")
    overlap = [r for r in rows if r[0] == "chr1" and r[1] == "100" and r[2] == "A" and r[3] == "G"]
    assert len(overlap) == 1
    r = overlap[0]
    # Columns: chrom, pos, ref, alt, ac_a, an_a, af_a, ac_b, an_b, af_b, af_diff
    assert r[4] == "3"  # ac_a
    assert r[5] == "6"  # an_a (3 samples * 2)
    assert r[6] == "0.5"  # af_a
    assert r[7] == "3"  # ac_b
    assert r[8] == "4"  # an_b (2 samples * 2)
    assert r[9] == "0.75"  # af_b
    assert r[10] == "-0.25"  # af_diff


def test_diff_includes_variants_unique_to_each_cohort(vcfclick_home):
    """Variants only in `case` (e.g. chr1:250) must appear with af_b=0,
    not be silently dropped. Same for variants only in `control`
    (chr1:200)."""
    _setup_two_cohort_db(vcfclick_home)
    rows = _tsv_rows(vcfclick_home, "demo")
    positions = {(r[0], r[1]) for r in rows}

    # case-only positions
    for pos in ("250", "500", "750", "900"):
        assert ("chr1", pos) in positions, f"missing case-only chr1:{pos}"

    # control-only position
    assert ("chr1", "200") in positions, "missing control-only chr1:200"

    # chr1:200 must have ac_a=0
    chr1_200 = next(r for r in rows if r[0] == "chr1" and r[1] == "200")
    assert chr1_200[4] == "0"  # ac_a
    assert chr1_200[6] == "0"  # af_a

    # chr1:250 must have ac_b=0
    chr1_250 = next(r for r in rows if r[0] == "chr1" and r[1] == "250")
    assert chr1_250[7] == "0"  # ac_b
    assert chr1_250[9] == "0"  # af_b


def test_diff_sorted_by_absolute_af_difference_desc(vcfclick_home):
    """The biggest-diff variant must come first in the output."""
    _setup_two_cohort_db(vcfclick_home)
    rows = _tsv_rows(vcfclick_home, "demo")
    af_diffs = [abs(float(r[10])) for r in rows]
    assert af_diffs == sorted(af_diffs, reverse=True)


def test_diff_top_n_limits_results(vcfclick_home):
    _setup_two_cohort_db(vcfclick_home)
    rows = _tsv_rows(vcfclick_home, "demo", "--top", "3")
    assert len(rows) == 3


def test_diff_errors_on_unknown_cohort(vcfclick_home):
    _setup_two_cohort_db(vcfclick_home)
    r = _vc(
        vcfclick_home, "db", "diff", "demo",
        "--cohort-a", "case", "--cohort-b", "nope",
        expect_failure=True,
    )
    assert "unknown cohort" in r.stderr
    assert "nope" in r.stderr


def test_diff_errors_when_db_missing(vcfclick_home):
    r = _vc(
        vcfclick_home, "db", "diff", "nonexistent",
        "--cohort-a", "case", "--cohort-b", "control",
        expect_failure=True,
    )
    assert "not found" in r.stderr


def test_diff_quote_helper_escapes_single_quotes():
    """SQL-injection guard: cohort names with embedded single quotes
    must be doubled, not passed through raw."""
    from cli.db import _quote_str

    assert _quote_str("plain") == "'plain'"
    assert _quote_str("o'brien") == "'o''brien'"
    # Classic injection attempt
    assert _quote_str("' OR 1=1 --") == "''' OR 1=1 --'"
    # The quoted form, if substituted into SQL like
    #   WHERE cohort = {q}
    # always closes with a single quote at the end, so the OR is
    # inside the string literal.
