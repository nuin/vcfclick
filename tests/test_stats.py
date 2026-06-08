"""Tests for `vcfclick db stats`.

Uses the routing fixture (2 variants, 2 samples, mixed typed and
overflow fields) so the expected counts and percentages are
hand-verifiable from the VCF:

  pos=100: every reserved + 3 lab-specific INFO tags, all 5 typed
           FORMAT fields, plus MYCUSTOM
  pos=200: stripped-down (GT/GQ/DP only)

This gives a clean mix of 100% / 50% / 66.7% / 0% populations across
columns and Map keys.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

# The current db_stats implementation depends on chDB-specific SQL
# (system.columns, countIf, ARRAY JOIN mapKeys). A DuckDB port lands
# as a follow-up; skip the suite for now when the DuckDB backend is
# active so the rest of the matrix is informative.
pytestmark = pytest.mark.skipif(
    os.environ.get("VCFCLICK_BACKEND", "").lower() == "duckdb",
    reason="db stats not yet ported to DuckDB backend",
)

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
ROUTING_VCF = Path(__file__).parent / "fixtures" / "routing.vcf.gz"


def _vc(
    home: Path, *args: str, expect_failure: bool = False
) -> subprocess.CompletedProcess:
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


def _stats(home: Path, extra: list[str] | None = None, db: str = "demo") -> str:
    """Ingest the routing fixture and return `db stats` stdout."""
    _vc(home, "db", "create", db)
    _vc(
        home,
        "db",
        "ingest",
        db,
        str(ROUTING_VCF),
        "--cohort",
        "acme",
        "--ingest-id",
        "batch_a",
        "--serial",
    )
    return _vc(home, "db", "stats", db, *(extra or [])).stdout


def _table_count(out: str, label: str) -> int:
    """Pull the integer at the end of a `<label>  <n>` row in the
    counts table. Avoids brittleness on column-width formatting."""
    line = next(line for line in out.splitlines() if line.strip().startswith(label))
    return int(line.strip().split()[-1].replace(",", ""))


def test_stats_reports_row_counts(vcfclick_home):
    out = _stats(vcfclick_home)
    assert _table_count(out, "variants") == 2
    assert _table_count(out, "genotypes") == 3
    assert _table_count(out, "samples") == 2
    assert _table_count(out, "ingestions") == 1


def test_stats_reports_cohort_breakdown(vcfclick_home):
    out = _stats(vcfclick_home)
    assert "cohorts:" in out
    assert "acme" in out
    # Two samples in the routing fixture: S1, S2
    cohort_line = next(line for line in out.splitlines() if "acme" in line)
    assert "2 samples" in cohort_line


def test_stats_reports_contig_breakdown(vcfclick_home):
    out = _stats(vcfclick_home)
    assert "contigs:" in out
    chr1_line = next(
        line for line in out.splitlines() if "chr1" in line and "variants" in line
    )
    assert "2 variants" in chr1_line


def test_stats_populated_typed_info_columns(vcfclick_home):
    """info_AC/AF/AN/DP are set on every variant in the fixture → 100%."""
    out = _stats(vcfclick_home)
    for col in ("info_AC", "info_AF", "info_AN", "info_DP"):
        line = next(line for line in out.splitlines() if line.strip().startswith(col))
        assert "100.0%" in line, f"{col} should be 100% populated: {line!r}"


def test_stats_partially_populated_columns_show_50pct(vcfclick_home):
    """info_AD_ref/alt + info_SOMATIC only set on variant at pos 100 → 50%."""
    out = _stats(vcfclick_home)
    for col in ("info_AD_ref", "info_AD_alt", "info_SOMATIC"):
        line = next(line for line in out.splitlines() if line.strip().startswith(col))
        assert "50.0%" in line, f"{col} should be 50% populated: {line!r}"


def test_stats_absent_typed_columns_show_0pct(vcfclick_home):
    """DRAGEN-specific columns aren't in the routing fixture → 0%.
    Same for unused GATK columns like MQ, QD, etc."""
    out = _stats(vcfclick_home)
    for col in ("info_HAPCOMP", "info_QD", "info_FS"):
        line = next(line for line in out.splitlines() if line.strip().startswith(col))
        assert "0.0%" in line, f"{col} should be 0% populated: {line!r}"


def test_stats_lists_info_extra_overflow_keys(vcfclick_home):
    """All three lab-specific INFO tags must appear in the overflow
    section with the right per-key population (each was set on 1 of 2
    rows = 50%)."""
    out = _stats(vcfclick_home)
    assert "info_extra" in out and "overflow keys" in out
    for key in ("MYRARETAG", "COSMICID", "CSQ"):
        line = next(
            line
            for line in out.splitlines()
            if line.strip().startswith(key) and "50.0%" in line
        )
        assert "1" in line, f"{key} should have count 1: {line!r}"


def test_stats_lists_format_extra_overflow_keys(vcfclick_home):
    """MYCUSTOM appears on the genotypes at pos 100 (both samples) but
    not pos 200 (only one sample stored). 2 of 3 genotype rows → 66.7%."""
    out = _stats(vcfclick_home)
    line = next(
        line for line in out.splitlines() if line.strip().startswith("MYCUSTOM")
    )
    assert "2" in line and "66.7%" in line


def test_stats_top_n_caps_overflow_listing(vcfclick_home):
    """--top 1 should keep only the top overflow key; header reports
    '(top 1 of 3)' for info_extra."""
    out = _stats(vcfclick_home, extra=["--top", "1"])
    info_extra_header = next(
        line for line in out.splitlines() if "info_extra" in line and "overflow" in line
    )
    assert "top 1 of 3" in info_extra_header


def test_stats_errors_on_missing_db(vcfclick_home):
    r = _vc(vcfclick_home, "db", "stats", "no-such-db", expect_failure=True)
    assert "not found" in r.stderr
