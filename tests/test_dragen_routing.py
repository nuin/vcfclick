"""Tests for DRAGEN-specific INFO routing.

DRAGEN (Illumina) emits a handful of record-level INFO scalars that
the routing tables now promote to typed columns:

  FractionInformativeReads (Float)
  HAPCOMP                  (Integer)
  HAPDOM                   (Float)
  DragenSnvHardQUAL        (Float, somatic SNV mode)
  DragenIndelHardQUAL      (Float, somatic indel mode)

Other DRAGEN-specific fields (lab-specific tags, per-sample ML_PROB
etc.) still land in the info_extra / format_extra Maps; this test
asserts the typed promotion happens AND that unknown DRAGEN-shaped
tags still fall through to overflow.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
DRAGEN_VCF = Path(__file__).parent / "fixtures" / "dragen.vcf.gz"


def _vc(home: Path, *args: str) -> str:
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    assert r.returncode == 0, (
        f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
        f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    )
    return r.stdout


def _tsv(home: Path, db: str, sql: str) -> list[list[str]]:
    out = _vc(home, "db", "query", db, f"{sql} FORMAT TabSeparated")
    return [line.split("\t") for line in out.strip().splitlines() if line.strip()]


def _ingest_dragen(home: Path, db: str = "dr") -> None:
    _vc(home, "db", "create", db)
    _vc(
        home,
        "db",
        "ingest",
        db,
        str(DRAGEN_VCF),
        "--cohort",
        "demo",
        "--ingest-id",
        "batch_a",
        "--serial",
    )


def test_classify_header_promotes_dragen_fields():
    """The header classifier reports the DRAGEN INFO fields as typed,
    not as overflow."""
    from cyvcf2 import VCF

    from ingest.routing import classify_header

    cls = classify_header(VCF(str(DRAGEN_VCF)))

    expected_typed = {
        "AC",
        "AF",
        "AN",
        "DP",
        "FractionInformativeReads",
        "HAPCOMP",
        "HAPDOM",
        "DragenSnvHardQUAL",
        "DragenIndelHardQUAL",
    }
    assert expected_typed.issubset(set(cls["typed_info"]))
    # Lab-specific vendor tag is NOT in the routing tables → overflow.
    assert "DRAGEN_VENDOR_TAG" in cls["extra_info"]


def test_dragen_snv_fields_populate_typed_columns(vcfclick_home):
    """pos 100 is the SNV record with DragenSnvHardQUAL set; HAPCOMP
    and HAPDOM populated; DragenIndelHardQUAL is NULL (indel-only)."""
    _ingest_dragen(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "dr",
        "SELECT info_FractionInformativeReads, info_HAPCOMP, info_HAPDOM, "
        "info_DragenSnvHardQUAL, info_DragenIndelHardQUAL "
        "FROM variants WHERE pos = 100",
    )
    assert rows == [["0.875", "3", "0.92", "87.5", "\\N"]]


def test_dragen_indel_fields_populate_typed_columns(vcfclick_home):
    """pos 200 is the indel record with DragenIndelHardQUAL set;
    DragenSnvHardQUAL is NULL."""
    _ingest_dragen(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "dr",
        "SELECT info_FractionInformativeReads, info_HAPCOMP, info_HAPDOM, "
        "info_DragenSnvHardQUAL, info_DragenIndelHardQUAL "
        "FROM variants WHERE pos = 200",
    )
    assert rows == [["0.812", "5", "0.78", "\\N", "72.3"]]


def test_lab_specific_tag_still_lands_in_info_extra(vcfclick_home):
    """A non-routed DRAGEN-shaped tag must still fall through to the
    overflow Map, not get silently dropped."""
    _ingest_dragen(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "dr",
        "SELECT info_extra['DRAGEN_VENDOR_TAG'] FROM variants WHERE pos = 100",
    )
    assert rows == [["batch_x"]]


def test_no_dragen_field_leaks_into_info_extra(vcfclick_home):
    """All five DRAGEN-routed fields should be ABSENT from info_extra
    after promotion. If they showed up there too we'd be double-storing."""
    _ingest_dragen(vcfclick_home)
    # chDB exposes `mapKeys(MAP)`; DuckDB names it `map_keys(MAP)`.
    fn = (
        "map_keys"
        if os.environ.get("VCFCLICK_BACKEND", "").lower() == "duckdb"
        else "mapKeys"
    )
    rows = _tsv(
        vcfclick_home,
        "dr",
        f"SELECT {fn}(info_extra) FROM variants WHERE pos = 100",
    )
    keys = rows[0][0].strip("[]").replace("'", "").split(",") if rows else []
    keys = [k.strip() for k in keys if k.strip()]

    for routed in (
        "FractionInformativeReads",
        "HAPCOMP",
        "HAPDOM",
        "DragenSnvHardQUAL",
        "DragenIndelHardQUAL",
    ):
        assert routed not in keys, (
            f"{routed} leaked into info_extra: should be in typed column only"
        )


def test_discover_lists_dragen_fields_as_typed():
    """`vcfclick discover` should now show DRAGEN fields in the typed
    bucket, not in the overflow bucket with a promotion hint."""
    r = subprocess.run(
        [VCFCLICK_BIN, "discover", str(DRAGEN_VCF)],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0
    # Walk the INFO section: each routed DRAGEN field must appear in the
    # `typed` line, NOT in the overflow lines.
    typed_line = next(
        line for line in r.stdout.splitlines() if "AC," in line and "DragenSnv" in line
    )
    for fid in (
        "FractionInformativeReads",
        "HAPCOMP",
        "HAPDOM",
        "DragenSnvHardQUAL",
        "DragenIndelHardQUAL",
    ):
        assert fid in typed_line, f"{fid} not shown as typed in discover output"
