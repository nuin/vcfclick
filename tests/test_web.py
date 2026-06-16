"""Tests for the optional `vcfclick web` UI (the `[web]` extra).

The whole module skips cleanly when fastapi/httpx aren't installed. The
combine endpoint and the write-guard need no database and run fully
in-process. The meta/query endpoints use a DuckDB cohort built via the
CLI in a subprocess — DuckDB sidesteps the chDB one-server-per-process
constraint that would otherwise clash with the rest of the suite in the
same process.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient  # noqa: E402

from vcfclick_web.app import app  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
FIX = Path(__file__).parent / "fixtures"

client = TestClient(app)

_HDR = (
    "##fileformat=VCFv4.2\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
)


def _vcf(sample_cols: str, rows: list[str]) -> str:
    head = _HDR + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_cols
    return head + "\n" + "\n".join(rows) + "\n"


# ───────────────────────── no database needed ─────────────────────────


def test_index_serves_spa():
    r = client.get("/")
    assert r.status_code == 200
    assert "vcfclick" in r.text
    # the SPA must wire the API endpoints it depends on
    assert "/api/query" in r.text and "/api/combine" in r.text


def test_query_rejects_writes():
    r = client.post("/api/query", json={"sql": "DROP TABLE variants"})
    assert r.status_code == 200
    assert "read-only" in r.json()["error"].lower()


def test_combine_endpoint_prioritizes_and_annotates_set():
    first = _vcf(
        "S1\tS2",
        ["chr1\t202\t.\tC\tA\t.\t.\t.\tGT\t1/1\t0/1"],  # S2 = 0/1 here
    )
    second = _vcf(
        "S2\tS3",
        [
            "chr1\t202\t.\tC\tA\t.\t.\t.\tGT\t1/1\t0/1",  # S2 = 1/1 here (lower priority)
            "chr1\t404\t.\tT\tG\t.\t.\t.\tGT\t0/0\t1/1",
        ],
    )
    d = client.post("/api/combine", json={"first": first, "second": second}).json()

    assert d["samples"] == ["S1", "S2", "S3"]
    by_pos = {r["pos"]: r for r in d["records"]}
    # chr1:202 is in both → Intersection; S2 keeps first's 0/1, not 1/1
    assert by_pos[202]["set"] == "Intersection"
    assert by_pos[202]["cells"]["S2"].startswith("0/1")
    # chr1:404 only in second; S1 absent → ./.
    assert by_pos[404]["set"] == "second"
    assert by_pos[404]["cells"]["S1"].startswith("./.")


# ───────────────────── DuckDB-backed meta + query ─────────────────────


@pytest.fixture
def duckdb_cohort(tmp_path, monkeypatch):
    home = tmp_path / "home"
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env["VCFCLICK_BACKEND"] = "duckdb"

    def run(*args):
        r = subprocess.run(
            [VCFCLICK_BIN, *args], env=env, capture_output=True, text=True
        )
        assert r.returncode == 0, f"{' '.join(args)} failed:\n{r.stderr}"

    run("db", "create", "webdb")
    run(
        "db",
        "ingest",
        "webdb",
        str(FIX / "trio.vcf.gz"),
        "--cohort",
        "trio",
        "--ingest-id",
        "i1",
        "--serial",
    )
    # Point the in-process app at the same DuckDB database.
    monkeypatch.setenv("VCFCLICK_HOME", str(home))
    monkeypatch.setenv("VCFCLICK_BACKEND", "duckdb")
    monkeypatch.setenv("VCFCLICK_DB_NAME", "webdb")
    return home


def test_meta_lists_tables(duckdb_cohort):
    m = client.get("/api/meta").json()
    assert m["db"] == "webdb"
    names = {t["name"] for t in m["tables"]}
    assert {"variants", "genotypes", "samples"} <= names
    # columns come from the locked Arrow schemas
    variants = next(t for t in m["tables"] if t["name"] == "variants")
    assert "chrom" in variants["columns"] and "pos" in variants["columns"]


def test_query_runs_select(duckdb_cohort):
    d = client.post(
        "/api/query", json={"sql": "SELECT count(*) AS n FROM variants"}
    ).json()
    assert "error" not in d, d
    assert d["columns"] == ["n"]
    assert int(d["rows"][0][0]) >= 1
