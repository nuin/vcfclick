"""Population panel: parsing (unit) and `db panel` (CLI, both backends)."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from ingest.panel import PanelColumns, PanelError, normalise_sex, parse_panel
from tests.conftest import run_cli

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
    return run_cli(home, backend, *args, ok=not expect_failure)


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

    # Idempotent: a re-load of the same panel changes nothing.
    _vc(home, backend, "db", "panel", "pg", str(PANEL))
    assert _rows(home, backend, "SELECT count(*) FROM populations") == [[20]]
    # A panel REPLACES the labelling of each ingestion it applies to: a
    # one-sample panel for `two` leaves exactly that label there, so no
    # sample keeps a stale label; `one` is untouched.
    relabel = popgen_vcf.parent / "relabel.tsv"
    relabel.write_text("sample\tpop\nA1\tESN\n")
    _vc(home, backend, "db", "panel", "pg", str(relabel), "--ingest-id", "two")
    assert _rows(home, backend, "SELECT count(*) FROM populations") == [[11]]
    assert _rows(
        home,
        backend,
        "SELECT ingest_id, population FROM populations WHERE sample_id = 'A1' "
        "ORDER BY ingest_id",
    ) == [["one", "YRI"], ["two", "ESN"]]
    _vc(home, backend, "db", "panel", "pg", str(PANEL))
    assert _rows(home, backend, "SELECT count(*) FROM populations") == [[20]]

    # Re-ingesting the same VCF keeps the panel (like `pedigree`).
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


@pytest.mark.parametrize("backend", ["duckdb", "chdb"])
def test_reingest_without_samples_prunes_their_labels(
    vcfclick_home, popgen_vcf, tmp_path, backend
):
    """Re-ingesting an ingest_id from a VCF that lacks some samples removes
    those samples' panel labels; the remaining samples keep theirs."""
    if not shutil.which("bcftools"):
        pytest.skip("bcftools not on PATH")
    fewer = tmp_path / "fewer.vcf.gz"
    subprocess.run(
        ["bcftools", "view", "-s", "^A1,A2", "-Oz", "-o", str(fewer), str(popgen_vcf)],
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", str(fewer)], check=True)
    home = vcfclick_home
    _vc(home, backend, "db", "create", "pg")
    _vc(
        home,
        backend,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "i1",
        "--serial",
    )
    _vc(home, backend, "db", "panel", "pg", str(PANEL))
    r = _vc(
        home,
        backend,
        "db",
        "ingest",
        "pg",
        str(fewer),
        "--ingest-id",
        "i1",
        "--workers",
        "2",
    )
    assert "removed 2 population label(s)" in r.stderr
    assert _rows(
        home,
        backend,
        "SELECT sample_id FROM populations WHERE population = 'YRI' ORDER BY sample_id",
    ) == [["A3"], ["A4"]]
