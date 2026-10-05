"""Population panel: parsing (unit) and `db panel` (CLI, both backends)."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from ingest.panel import PanelColumns, PanelError, normalise_sex, parse_panel

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
FIXTURES = Path(__file__).parent / "fixtures"
PANEL = FIXTURES / "popgen.panel"


# ─────────────────────────────── parsing ────────────────────────────────


def test_parse_1000g_panel_with_trailing_tabs():
    rows = {r["sample_id"]: r for r in parse_panel(PANEL)}
    assert len(rows) == 11
    assert rows["A1"] == {
        "sample_id": "A1",
        "population": "YRI",
        "super_population": "AFR",
        "sex": "female",
    }
    assert rows["C1"]["sex"] == "male"


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("male", "male"),
        ("M", "male"),
        ("1", "male"),
        ("Female", "female"),
        ("f", "female"),
        ("2", "female"),
        ("0", None),
        ("unknown", None),
        ("", None),
        (None, None),
    ],
)
def test_normalise_sex(raw, expected):
    assert normalise_sex(raw) == expected


def test_parse_generic_csv_with_named_columns(tmp_path):
    p = tmp_path / "labels.csv"
    p.write_text("IID,group,region,notes\nS1,north,EU,x\nS2,south,,y\n")
    rows = parse_panel(p, PanelColumns(population="group", super_population="region"))
    assert rows == [
        {
            "sample_id": "S1",
            "population": "north",
            "super_population": "EU",
            "sex": None,
        },
        {
            "sample_id": "S2",
            "population": "south",
            "super_population": None,
            "sex": None,
        },
    ]


def test_parse_whitespace_separated(tmp_path):
    p = tmp_path / "labels.txt"
    p.write_text("sample population sex\nS1 popA 1\nS2 popB 2\n")
    rows = parse_panel(p)
    assert [(r["sample_id"], r["population"], r["sex"]) for r in rows] == [
        ("S1", "popA", "male"),
        ("S2", "popB", "female"),
    ]


def test_parse_rejects_missing_population(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("sample\tpop\nS1\t\n")
    with pytest.raises(PanelError, match="no population"):
        parse_panel(p)


def test_parse_rejects_unknown_columns(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("who\twhere\nS1\tA\n")
    with pytest.raises(PanelError, match="--sample-col"):
        parse_panel(p)
    with pytest.raises(PanelError, match="not in the panel header"):
        parse_panel(p, PanelColumns(sample="who", population="nope"))


def test_parse_rejects_conflicting_duplicates(tmp_path):
    p = tmp_path / "dup.tsv"
    p.write_text("sample\tpop\nS1\tA\nS1\tB\n")
    with pytest.raises(PanelError, match="twice"):
        parse_panel(p)
    p.write_text("sample\tpop\nS1\tA\nS1\tA\n")
    assert len(parse_panel(p)) == 1  # identical duplicate is harmless


# ─────────────────────────────── db panel ───────────────────────────────


def _vc(home: Path, backend: str, *args: str, expect_failure: bool = False):
    env = {**os.environ, "VCFCLICK_HOME": str(home), "VCFCLICK_BACKEND": backend}
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if expect_failure:
        assert r.returncode != 0, r.stdout
    else:
        assert r.returncode == 0, f"{args}:\n{r.stdout}\n{r.stderr}"
    return r


def _rows(home: Path, backend: str, sql: str) -> list[list]:
    r = _vc(home, backend, "db", "query", "pg", sql, "--format", "JSONCompact")
    return json.loads(r.stdout)["data"]


@pytest.mark.parametrize("backend", ["duckdb", "chdb"])
def test_db_panel_loads_reports_and_reloads(vcfclick_home, popgen_vcf, backend):
    home = vcfclick_home
    _vc(home, backend, "db", "create", "pg")
    for iid in ("one", "two"):
        _vc(
            home,
            backend,
            "db",
            "ingest",
            "pg",
            str(popgen_vcf),
            "--ingest-id",
            iid,
            "--serial",
        )

    out = _vc(home, backend, "db", "panel", "pg", str(PANEL)).stdout
    # Applies to both ingestions; reports both directions of mismatch.
    assert "10 samples under ingest_id=one" in out
    assert "10 samples under ingest_id=two" in out
    assert "panel samples not in the database (1): ZZ9" in out
    assert "database samples missing from the panel (2)" in out and "U1" in out

    rows = _rows(
        home,
        backend,
        "SELECT ingest_id, sample_id, population, super_population, sex "
        "FROM populations WHERE sample_id IN ('A1', 'C1') ORDER BY ingest_id, sample_id",
    )
    assert rows == [
        ["one", "A1", "YRI", "AFR", "female"],
        ["one", "C1", "CHB", "EAS", "male"],
        ["two", "A1", "YRI", "AFR", "female"],
        ["two", "C1", "CHB", "EAS", "male"],
    ]

    # Idempotent: a re-load replaces, a relabel takes effect.
    relabel = popgen_vcf.parent / "relabel.tsv"
    relabel.write_text("sample\tpop\nA1\tESN\n")
    _vc(home, backend, "db", "panel", "pg", str(PANEL))
    _vc(home, backend, "db", "panel", "pg", str(relabel), "--ingest-id", "two")
    assert _rows(home, backend, "SELECT count(*) FROM populations") == [[20]]
    assert _rows(
        home,
        backend,
        "SELECT ingest_id, population FROM populations WHERE sample_id = 'A1' "
        "ORDER BY ingest_id",
    ) == [["one", "YRI"], ["two", "ESN"]]

    # Re-ingesting a VCF does not wipe the panel (like `pedigree`).
    _vc(
        home,
        backend,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "one",
        "--serial",
    )
    assert _rows(home, backend, "SELECT count(*) FROM populations") == [[20]]


def test_db_panel_errors(vcfclick_home, popgen_vcf, tmp_path):
    home, b = vcfclick_home, "duckdb"
    _vc(home, b, "db", "create", "pg")
    _vc(
        home, b, "db", "ingest", "pg", str(popgen_vcf), "--ingest-id", "one", "--serial"
    )
    stranger = tmp_path / "s.tsv"
    stranger.write_text("sample\tpop\nNOBODY\tX\n")
    r = _vc(home, b, "db", "panel", "pg", str(stranger), expect_failure=True)
    assert "none of the 1 panel samples" in r.stderr
    r = _vc(
        home,
        b,
        "db",
        "panel",
        "pg",
        str(PANEL),
        "--ingest-id",
        "nope",
        expect_failure=True,
    )
    assert "no samples found" in r.stderr
