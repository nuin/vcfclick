"""Tests for the pedigree foundation: PED parsing + `db ped` load.

The pedigree table maps (ingest_id, sample_id) -> family/father/mother/
sex/affected, loaded from a PED file separately from VCF ingest. Trio
analysis (a later increment) resolves a proband's parents from it.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
FIXTURES = Path(__file__).parent / "fixtures"
TRIO_VCF = FIXTURES / "trio.vcf.gz"
TRIO_PED = FIXTURES / "trio.ped"


def _vc(home: Path, *args: str, expect_failure: bool = False):
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


def _rows(home: Path, db: str, sql: str) -> list[list]:
    r = _vc(home, "db", "query", db, sql, "--format", "JSONCompact")
    return json.loads(r.stdout)["data"]


def _ingest_trio(home: Path, db: str = "fam") -> None:
    _vc(home, "db", "create", db)
    _vc(
        home,
        "db",
        "ingest",
        db,
        str(TRIO_VCF),
        "--cohort",
        "trio",
        "--ingest-id",
        "fam1",
        "--serial",
    )


# ─────────────────────── parse_ped unit ───────────────────────


def test_parse_ped_normalizes_sex_affected_and_founders():
    from ingest.pedigree import parse_ped

    rows = {r["sample_id"]: r for r in parse_ped(TRIO_PED)}
    assert set(rows) == {"CHILD", "FATHER", "MOTHER"}

    child = rows["CHILD"]
    assert child["family_id"] == "FAM1"
    assert child["father_id"] == "FATHER"
    assert child["mother_id"] == "MOTHER"
    assert child["sex"] == "male"
    assert child["affected"] == "affected"

    # Founders carry '0' parents and unaffected status.
    assert rows["FATHER"]["father_id"] == "0"
    assert rows["FATHER"]["sex"] == "male"
    assert rows["FATHER"]["affected"] == "unaffected"
    assert rows["MOTHER"]["sex"] == "female"


# ─────────────────────── db ped load ───────────────────────


def test_db_ped_loads_relationships(vcfclick_home):
    _ingest_trio(vcfclick_home)
    # ingest_id is inferred (single ingestion).
    r = _vc(vcfclick_home, "db", "ped", "fam", str(TRIO_PED))
    assert "3 individuals" in r.stdout

    rows = _rows(
        vcfclick_home,
        "fam",
        "SELECT sample_id, father_id, mother_id, sex, affected "
        "FROM pedigree ORDER BY sample_id",
    )
    by_id = {row[0]: row for row in rows}
    assert by_id["CHILD"][1:] == ["FATHER", "MOTHER", "male", "affected"]
    assert by_id["FATHER"][1:] == ["0", "0", "male", "unaffected"]


def test_db_info_shows_pedigree_count(vcfclick_home):
    _ingest_trio(vcfclick_home)
    _vc(vcfclick_home, "db", "ped", "fam", str(TRIO_PED))
    out = _vc(vcfclick_home, "db", "info", "fam").stdout
    assert "pedigree:" in out
    # The count line should show 3.
    ped_line = next(ln for ln in out.splitlines() if ln.startswith("pedigree:"))
    assert "3" in ped_line


def test_db_ped_reload_replaces_not_appends(vcfclick_home):
    _ingest_trio(vcfclick_home)
    _vc(vcfclick_home, "db", "ped", "fam", str(TRIO_PED))
    _vc(vcfclick_home, "db", "ped", "fam", str(TRIO_PED))  # again
    rows = _rows(vcfclick_home, "fam", "SELECT count(*) FROM pedigree")
    assert rows == [[3]], f"reload should replace, not double: {rows}"


def test_db_ped_rejects_unknown_sample(vcfclick_home, tmp_path):
    _ingest_trio(vcfclick_home)
    bad = tmp_path / "bad.ped"
    bad.write_text("FAM1\tGHOST\t0\t0\t1\t2\n")
    r = _vc(vcfclick_home, "db", "ped", "fam", str(bad), expect_failure=True)
    assert "GHOST" in (r.stdout + r.stderr)


def test_db_ped_round_trips_through_dump(vcfclick_home, tmp_path):
    """`db dump` must include pedigree.parquet so a bundled trio DB
    keeps its relationships."""
    _ingest_trio(vcfclick_home)
    _vc(vcfclick_home, "db", "ped", "fam", str(TRIO_PED))
    dump = tmp_path / "dump"
    _vc(vcfclick_home, "db", "dump", "fam", "--out", str(dump))
    assert (dump / "pedigree.parquet").exists()
