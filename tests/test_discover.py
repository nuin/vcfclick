"""Tests for `vcfclick discover` — the schema-routing previewer."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
ROUTING_VCF = Path(__file__).parent / "fixtures" / "routing.vcf.gz"


def _vc(*args: str) -> str:
    env = os.environ.copy()
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
    )
    assert (
        r.returncode == 0
    ), f"`vcfclick {' '.join(args)}` failed:\n{r.stdout}\n{r.stderr}"
    return r.stdout


def test_discover_reports_typed_and_overflow_sections():
    out = _vc("discover", str(ROUTING_VCF))
    assert "INFO fields:" in out
    assert "FORMAT fields:" in out
    assert "variants.info_extra Map" in out
    assert "genotypes.format_extra Map" in out


def test_discover_lists_typed_info_fields():
    """The routing fixture's typed INFO set: AC, AD, AF, AN, DP, SOMATIC."""
    out = _vc("discover", str(ROUTING_VCF))
    for fid in ("AC", "AD", "AF", "AN", "DP", "SOMATIC"):
        assert fid in out


def test_discover_lists_overflow_info_fields():
    out = _vc("discover", str(ROUTING_VCF))
    for fid in ("COSMICID", "CSQ", "MYRARETAG"):
        assert fid in out


def test_discover_suggests_promotion_ddl_for_overflow_scalars():
    """COSMICID (String, Number=1) → info_COSMICID Nullable(String)."""
    out = _vc("discover", str(ROUTING_VCF))
    assert "info_COSMICID Nullable(String)" in out
    assert "info_MYRARETAG Nullable(String)" in out


def test_discover_skips_promotion_for_variable_length_lists():
    """CSQ has Number=. — must be flagged as un-promotable."""
    out = _vc("discover", str(ROUTING_VCF))
    # CSQ should appear with the "keep in Map" marker, not a promote hint.
    csq_line = next(line for line in out.splitlines() if line.strip().startswith("CSQ"))
    assert "keep in Map" in csq_line
    assert "promote" not in csq_line


def test_discover_lowercases_format_column_names():
    """FORMAT-derived suggestions follow the gq/dp/ad_ref convention →
    MYCUSTOM (FORMAT) should suggest the lowercase `mycustom`, not `MYCUSTOM`."""
    out = _vc("discover", str(ROUTING_VCF))
    # Find the MYCUSTOM line specifically (in the FORMAT overflow section).
    mycustom_line = next(
        line for line in out.splitlines() if "MYCUSTOM" in line and "promote:" in line
    )
    assert "mycustom Nullable(String)" in mycustom_line


def test_discover_does_not_promote_already_typed_fields():
    """Fields like AC, DP already in INFO_SCALAR must not appear with a
    promote hint — they're typed already."""
    out = _vc("discover", str(ROUTING_VCF))
    # Walk the overflow section only and assert AC/DP are not there.
    in_overflow = False
    for line in out.splitlines():
        if "info_extra Map" in line:
            in_overflow = True
            continue
        if "FORMAT fields:" in line:
            break
        if in_overflow and line.strip().startswith(
            ("AC ", "AF ", "AN ", "DP ", "SOMATIC ")
        ):
            raise AssertionError(f"typed field leaked into overflow section: {line!r}")


def test_discover_errors_on_missing_vcf():
    """A missing path should exit non-zero with a clear message — not crash."""
    env = os.environ.copy()
    r = subprocess.run(
        [VCFCLICK_BIN, "discover", "/no/such/file.vcf.gz"],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
    )
    assert r.returncode != 0
    assert (
        "no such" in (r.stdout + r.stderr).lower()
        or "does not exist" in (r.stdout + r.stderr).lower()
    )
