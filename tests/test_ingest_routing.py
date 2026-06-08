"""Tests for the schema-routing claim: VCF 4.3 reserved + common GATK
fields land in typed columns; everything else lands in info_extra /
format_extra Maps.

These exist because the "adapts to any VCF" claim in the README is the
central pitch of this project. The smoke tests cover the lifecycle but
don't check WHERE individual fields end up.

The fixture (tests/fixtures/routing.vcf.gz) has two records:

  pos=100: every reserved INFO field (AC/AF/AN/DP/AD/SOMATIC) plus
           three lab-specific tags (MYRARETAG/COSMICID/CSQ), and every
           reserved FORMAT field (GT/GQ/DP/AD/PL) plus one custom one
           (MYCUSTOM).

  pos=200: a stripped-down record with GT/GQ/DP only — used to catch
           a real cyvcf2 quirk where `variant.gt_depths` returns -1
           when FORMAT column ordering changes mid-VCF, even though
           `variant.format('DP')` returns the right value. The fix is
           in ingest/vcf_load.py::build_genotype_rows.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
ROUTING_VCF = Path(__file__).parent / "fixtures" / "routing.vcf.gz"


def _vc(home: Path, *args: str) -> str:
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
    assert r.returncode == 0, (
        f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
        f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    )
    return r.stdout


def _tsv(home: Path, db: str, sql: str) -> list[list[str]]:
    """Run a SQL query in TabSeparated format → list of row-lists."""
    out = _vc(home, "db", "query", db, f"{sql} FORMAT TabSeparated")
    return [line.split("\t") for line in out.strip().splitlines() if line.strip()]


def _ingest_routing(home: Path, db: str = "rt") -> None:
    _vc(home, "db", "create", db)
    _vc(
        home,
        "db",
        "ingest",
        db,
        str(ROUTING_VCF),
        "--cohort",
        "demo",
        "--ingest-id",
        "batch_a",
        "--serial",
    )


# --- header classification (pure function, no DB needed) ---


def test_classify_header_splits_typed_vs_extra():
    """The ingester's header classifier reports which INFO/FORMAT fields
    will land in typed columns vs the overflow Maps."""
    from cyvcf2 import VCF

    from ingest.routing import classify_header

    cls = classify_header(VCF(str(ROUTING_VCF)))

    # Typed INFO: AC, AF, AN, AD, DP, SOMATIC (6 fields)
    assert set(cls["typed_info"]) == {"AC", "AF", "AN", "AD", "DP", "SOMATIC"}
    # Overflow INFO: the three lab-specific tags
    assert set(cls["extra_info"]) == {"MYRARETAG", "COSMICID", "CSQ"}

    # Typed FORMAT: GT, GQ, DP, AD, PL (5 fields)
    assert set(cls["typed_format"]) == {"GT", "GQ", "DP", "AD", "PL"}
    # Overflow FORMAT: just MYCUSTOM
    assert set(cls["extra_format"]) == {"MYCUSTOM"}


# --- INFO column routing ---


def test_info_scalar_fields_land_in_typed_columns(vcfclick_home):
    _ingest_routing(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "rt",
        "SELECT info_AC, info_AF, info_AN, info_DP FROM variants WHERE pos = 100",
    )
    assert rows == [["1", "0.25", "4", "50"]]


def test_info_pair_AD_lands_in_typed_columns(vcfclick_home):
    """INFO=AD (Number=R) splits to info_AD_ref / info_AD_alt."""
    _ingest_routing(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "rt",
        "SELECT info_AD_ref, info_AD_alt FROM variants WHERE pos = 100",
    )
    assert rows == [["30", "20"]]


def test_info_flag_SOMATIC_lands_as_one(vcfclick_home):
    """A Flag-typed INFO field present → 1, absent → 0 (never NULL)."""
    _ingest_routing(vcfclick_home)
    rows = _tsv(
        vcfclick_home, "rt", "SELECT pos, info_SOMATIC FROM variants ORDER BY pos"
    )
    assert rows == [["100", "1"], ["200", "0"]]


def test_unknown_info_fields_land_in_info_extra(vcfclick_home):
    """Lab-specific tags must end up in info_extra as strings."""
    _ingest_routing(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "rt",
        "SELECT info_extra['MYRARETAG'], info_extra['COSMICID'], info_extra['CSQ'] "
        "FROM variants WHERE pos = 100",
    )
    assert rows == [["hello", "COSM12345", "missense|MODERATE|BRCA1"]]


def test_info_extra_is_empty_when_no_unknown_fields(vcfclick_home):
    """A record with only reserved INFO fields → info_extra is {}."""
    _ingest_routing(vcfclick_home)
    # Backend-portable empty-map check: chDB has `length(MAP)`, DuckDB
    # has `cardinality(MAP)`. Pick at SQL build time.
    fn = (
        "cardinality"
        if os.environ.get("VCFCLICK_BACKEND", "").lower() == "duckdb"
        else "length"
    )
    rows = _tsv(
        vcfclick_home, "rt", f"SELECT {fn}(info_extra) FROM variants WHERE pos = 200"
    )
    assert rows == [["0"]]


def test_unknown_info_does_NOT_leak_into_typed_columns(vcfclick_home):
    """COSMICID and MYRARETAG must not appear as columns on `variants`."""
    _ingest_routing(vcfclick_home)
    # System-catalog tables differ between backends. chDB exposes
    # `system.columns`; DuckDB uses the SQL-standard
    # `information_schema.columns` (column `column_name`, not `name`).
    if os.environ.get("VCFCLICK_BACKEND", "").lower() == "duckdb":
        sql = (
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name = 'variants'"
        )
    else:
        sql = "SELECT name FROM system.columns WHERE table = 'variants'"
    cols = _tsv(vcfclick_home, "rt", sql)
    col_names = {row[0] for row in cols}
    for forbidden in ("info_COSMICID", "info_MYRARETAG", "info_CSQ"):
        assert forbidden not in col_names


# --- FORMAT column routing ---


def test_format_scalar_fields_land_in_typed_columns(vcfclick_home):
    _ingest_routing(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "rt",
        "SELECT sample_id, gq, dp FROM genotypes WHERE pos = 100 ORDER BY sample_id",
    )
    assert rows == [["S1", "30", "20"], ["S2", "28", "18"]]


def test_format_pair_AD_lands_in_typed_columns(vcfclick_home):
    _ingest_routing(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "rt",
        "SELECT sample_id, ad_ref, ad_alt FROM genotypes WHERE pos = 100 ORDER BY sample_id",
    )
    assert rows == [["S1", "10", "10"], ["S2", "0", "18"]]


def test_format_triple_PL_lands_in_typed_columns(vcfclick_home):
    """Regression: PL (Number=G) used to silently drop because the
    genotype builder never read it. Now routed to pl_ref_ref/ref_alt/alt_alt."""
    _ingest_routing(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "rt",
        "SELECT sample_id, pl_ref_ref, pl_ref_alt, pl_alt_alt "
        "FROM genotypes WHERE pos = 100 ORDER BY sample_id",
    )
    assert rows == [
        ["S1", "50", "0", "80"],
        ["S2", "120", "40", "0"],
    ]


def test_unknown_format_lands_in_format_extra(vcfclick_home):
    """MYCUSTOM must end up in format_extra per-sample."""
    _ingest_routing(vcfclick_home)
    rows = _tsv(
        vcfclick_home,
        "rt",
        "SELECT sample_id, format_extra['MYCUSTOM'] "
        "FROM genotypes WHERE pos = 100 ORDER BY sample_id",
    )
    assert rows == [["S1", "tag-a"], ["S2", "tag-b"]]


def test_dp_populates_when_FORMAT_order_changes_mid_vcf(vcfclick_home):
    """Regression: pos 200's FORMAT is `GT:GQ:DP` (no AD/PL), coming
    after pos 100's `GT:GQ:DP:AD:PL:MYCUSTOM`. The cyvcf2 `gt_depths`
    shortcut would return -1 here even though DP is clearly present.
    The fix uses `variant.format('DP')` instead."""
    _ingest_routing(vcfclick_home)
    rows = _tsv(
        vcfclick_home, "rt", "SELECT sample_id, gq, dp FROM genotypes WHERE pos = 200"
    )
    # S1 is 0/0 (not stored). S2 should have gq=25, dp=15.
    assert rows == [["S2", "25", "15"]]


def test_format_extra_empty_when_only_reserved_fields(vcfclick_home):
    _ingest_routing(vcfclick_home)
    fn = (
        "cardinality"
        if os.environ.get("VCFCLICK_BACKEND", "").lower() == "duckdb"
        else "length"
    )
    rows = _tsv(
        vcfclick_home,
        "rt",
        f"SELECT {fn}(format_extra) FROM genotypes WHERE pos = 200",
    )
    assert rows == [["0"]]
