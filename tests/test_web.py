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

from storage.sql_guard import is_read_only  # noqa: E402
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
    head = (
        _HDR + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_cols
    )
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


# The write guard must not be fooled by comments, multiple statements,
# CTE-prefixed DML, or file-write constructs — the bypasses a localhost
# server has to resist (CSRF / DNS-rebinding can POST to it).
@pytest.mark.parametrize(
    "sql",
    [
        "DROP TABLE variants",
        "/*x*/ DROP TABLE variants",
        "-- c\nDROP TABLE variants",
        "SELECT 1; DROP TABLE variants",
        "WITH c AS (SELECT 1) DELETE FROM variants",
        "SELECT 42 INTO OUTFILE '/tmp/x'",
        "COPY (SELECT 42) TO '/tmp/x'",
        "INSERT INTO variants VALUES (1)",
        "   update samples SET x = 1",
        "ATTACH 'evil.db'",
    ],
)
def test_guard_blocks_writes_and_bypasses(sql):
    assert is_read_only(sql) is False


@pytest.mark.parametrize(
    "sql",
    [
        "SELECT * FROM variants LIMIT 5",
        "select chrom, pos from variants",
        "WITH c AS (SELECT 1 AS n) SELECT n FROM c",
        "SHOW TABLES",
        "DESCRIBE variants",
        "EXPLAIN SELECT 1",
        "SELECT count(*) FROM genotypes  -- trailing comment\n",
        # truncate()/replace() are functions, not statements — must pass
        "SELECT replace(chrom, 'chr', '') FROM variants",
        "SELECT truncate(qual) FROM variants",
    ],
)
def test_guard_allows_reads(sql):
    assert is_read_only(sql) is True


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


def test_meta_lists_live_columns_of_an_older_database(duckdb_cohort):
    """A database from before the called-count columns is not upgraded by
    the read-only web UI, so /api/meta must not advertise them."""
    from storage import get_session

    get_session("webdb").query("ALTER TABLE variants DROP COLUMN n_called", "CSV")
    m = client.get("/api/meta").json()
    variants = next(t for t in m["tables"] if t["name"] == "variants")
    assert "n_called" not in variants["columns"]
    assert "an_called" in variants["columns"] and "chrom" in variants["columns"]


def test_meta_without_a_pedigree_has_no_probands(duckdb_cohort):
    m = client.get("/api/meta").json()
    assert m["probands"] == [] and m["families"] == []


def test_meta_of_a_database_with_no_pedigree_table(trio_web_store):
    """An older database has no `pedigree` table at all: the picker is empty
    and the schema list leaves the table out, instead of /api/meta failing."""
    trio_web_store(FIX / "trio.vcf.gz")
    _query("DROP TABLE pedigree")
    m = client.get("/api/meta").json()
    assert m["probands"] == [] and m["families"] == []
    assert "pedigree" not in {t["name"] for t in m["tables"]}
    # /api/trio says why, in the words `db trio` uses
    d = client.get("/api/trio", params={"proband": "CHILD"}).json()
    assert "no pedigree loaded" in d["error"]


def _query(sql: str) -> None:
    from storage import get_session

    get_session().query(sql, "CSV")


def test_meta_probands_follow_the_proband_rule_and_families_are_served(
    trio_web_store,
):
    """Several families and generations: `probands` is the flat list by the
    one rule `db trio` and `db families` use; `families` is the structure."""
    trio_web_store(FIX / "related_ext.vcf", ped=FIX / "related_ext.ped")
    m = client.get("/api/meta").json()
    # C2, UF are unaffected, F and M have children: none is offered
    assert m["probands"] == sorted(["C1", "C2T", "H1", "H2", "K1", "K2", "XB"])
    assert [f["family_id"] for f in m["families"]] == ["famA", "famB"]
    fam_a = m["families"][0]
    assert list(fam_a) == ["family_id", "probands", "trios"]
    assert fam_a["probands"] == ["C1", "C2T", "H1", "H2"]
    trios = {t["child"]: t for t in fam_a["trios"]}
    assert list(trios["C1"]) == [
        "child", "father", "mother", "kind", "affected", "siblings",
    ]  # fmt: skip
    assert trios["C1"]["affected"] == "affected" and trios["C1"]["kind"] == "trio"
    assert {s["sample"]: s["kind"] for s in trios["C1"]["siblings"]} == {
        "C2": "full",
        "C2T": "full",
        "H1": "half_paternal",
        "H2": "half_paternal",
    }
    assert trios["C2"]["affected"] == "unaffected"  # a trio, just not a proband


def test_meta_reads_a_pedigree_the_loader_would_refuse(trio_web_store):
    trio_web_store(FIX / "related_ext.vcf", ped=FIX / "related_ext.ped")
    _query("DELETE FROM pedigree WHERE sample_id = 'PGF'")
    _query(
        "INSERT INTO pedigree (ingest_id, sample_id, family_id, father_id, "
        "mother_id, sex, affected) VALUES ('i1', 'PGF', 'famA', 'F', '0', "
        "'male', NULL)"  # PGF is F's father and F is PGF's: a cycle, and a duo
    )
    m = client.get("/api/meta").json()
    assert m["probands"] == ["C1", "C2T", "H1", "H2", "K1", "K2", "XB"]
    kinds = {t["child"]: t["kind"] for t in m["families"][0]["trios"]}
    assert kinds["PGF"] == "duo"  # listed, but not offered as a proband
    assert "PGF" not in m["probands"]


def test_meta_with_several_ingestions_falls_back_to_the_flat_list(
    trio_web_store, tmp_path
):
    trio_web_store(FIX / "related_ext.vcf", ped=FIX / "related_ext.ped")
    name = os.environ["VCFCLICK_DB_NAME"]
    subprocess.run(
        [VCFCLICK_BIN, "db", "ingest", name, str(FIX / "trio.vcf.gz"),
         "--ingest-id", "i2", "--serial"],
        env=os.environ.copy(), check=True, capture_output=True,
    )  # fmt: skip
    m = client.get("/api/meta").json()
    assert m["families"] == []
    # the older rule: every individual with two declared parents
    assert {"C1", "C2", "F", "M", "K1"} <= set(m["probands"])


def test_query_runs_select(duckdb_cohort):
    d = client.post(
        "/api/query", json={"sql": "SELECT count(*) AS n FROM variants"}
    ).json()
    assert "error" not in d, d
    assert d["columns"] == ["n"]
    assert int(d["rows"][0][0]) >= 1


def test_unexpected_errors_are_logged_not_returned(monkeypatch, caplog):
    """Raw exception text stays in the server log; the browser gets a
    generic message (and the generated SQL, for /api/nl)."""
    import vcfclick_web.app as webapp
    import vcfclick_web.llm as llm

    def boom(sql):
        raise RuntimeError("/secret/path/internal detail")

    monkeypatch.setattr(webapp, "_run_sql", boom)
    client = TestClient(app)

    out = client.post("/api/query", json={"sql": "SELECT 1"}).json()
    assert "internal detail" not in out["error"]
    assert "query execution failed" in out["error"]

    monkeypatch.setattr(llm, "generate_sql", lambda *a, **k: "SELECT 1")
    out = client.post(
        "/api/nl",
        json={"question": "how many?", "provider": "gemini", "key": "k", "model": ""},
    ).json()
    assert "internal detail" not in out["error"]
    assert out["sql"] == "SELECT 1"
    assert "internal detail" in caplog.text


# ───────────────────────── /api/trio on a family ─────────────────────────


@pytest.fixture
def trio_web_store(tmp_path, monkeypatch):
    """`make(vcf, keep_reference=True)` builds a one-ingestion DuckDB trio
    (pedigree loaded) the in-process app then reads.

    The storage layer caches one session per database NAME, not per path, so
    every store gets a name of its own and its session is dropped afterwards;
    otherwise a later test would read the first test's data."""
    import uuid

    import storage.db as sdb

    names: list[str] = []

    def make(vcf, *, keep_reference=True, ped=None) -> str:
        name = f"triodb_{uuid.uuid4().hex[:8]}"
        names.append(name)
        env = os.environ.copy()
        env["VCFCLICK_HOME"] = str(tmp_path / "triohome")
        env["VCFCLICK_BACKEND"] = "duckdb"
        env["VCFCLICK_ANNOTATIONS_DB"] = str(tmp_path / "ann.duckdb")

        def run(*args):
            r = subprocess.run(
                [VCFCLICK_BIN, *args], env=env, capture_output=True, text=True
            )
            assert r.returncode == 0, f"{' '.join(args)} failed:\n{r.stderr}"

        run("db", "create", name)
        ingest = ["db", "ingest", name, str(vcf), "--ingest-id", "i1", "--serial"]
        if keep_reference:
            ingest.append("--keep-reference")
        run(*ingest)
        run("db", "ped", name, str(ped or FIX / "trio.ped"))
        for key, value in env.items():
            if key.startswith("VCFCLICK_"):
                monkeypatch.setenv(key, value)
        monkeypatch.setenv("VCFCLICK_DB_NAME", name)
        return name

    yield make
    for name in names:
        session = sdb._sessions.pop(f"duckdb::{name}", None)
        if hasattr(session, "close"):
            session.close()


def _seed_gene(tmp_path, monkeypatch):
    """GENEX over chr1:1-100000 in the (per-test) annotation store."""
    from annotations.db import get_connection

    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    conn = get_connection()
    conn.execute(
        "INSERT INTO refseq_genes VALUES "
        "('GENEX', 'chr1', 1, 100000, '+', 'id', 'protein_coding')"
    )
    conn.close()


def test_trio_endpoint_serves_every_category(tmp_path, monkeypatch, trio_web_store):
    from cli.trio_analysis import CATEGORIES

    trio_web_store(FIX / "trio.vcf.gz")
    _seed_gene(tmp_path, monkeypatch)
    out = {}
    for category in CATEGORIES:
        d = client.get(
            "/api/trio",
            params={"proband": "CHILD", "category": category, "max_af": 0.3},
        ).json()
        assert "error" not in d, (category, d)
        assert d["father"] == "FATHER" and d["mother"] == "MOTHER"
        # five sites are too few for the cohort-AF probe to judge
        assert d["af_gate"] == "unavailable"
        assert d["row_count"] == len(d["rows"])
        out[category] = d
    # the original three: one row per variant, same columns as before
    assert out["denovo"]["columns"] == [
        "chrom", "pos", "ref", "alt", "proband_gt", "father_gt", "mother_gt", "af",
    ]  # fmt: skip
    assert [r[1] for r in out["denovo"]["rows"]] == [100]
    assert [r[1] for r in out["recessive"]["rows"]] == [200]
    assert [r[1] for r in out["dominant"]["rows"]] == [300, 400]
    # comphet: GENEX has a paternal (300) and a maternal (400) het
    comphet = out["comphet"]
    assert comphet["columns"] == ["gene", "origin", "chrom", "pos", "ref", "alt", "af"]
    assert [(r[0], r[1], r[3]) for r in comphet["rows"]] == [
        ("GENEX", "paternal", 300),
        ("GENEX", "maternal", 400),
    ]
    assert comphet["genes"] == 1
    assert "note" not in comphet and "note" not in out["denovo"]


def test_trio_endpoint_comphet_says_when_no_genes_are_loaded(
    trio_web_store,
):
    trio_web_store(FIX / "trio.vcf.gz")
    d = client.get(
        "/api/trio",
        params={"proband": "CHILD", "category": "comphet", "max_af": 0.3},
    ).json()
    assert d["rows"] == [] and "No gene annotations" in d["note"]


def test_trio_endpoint_keep_reference_note_follows_the_model(trio_web_store):
    """de novo / dominant / comphet need stored hom-ref rows; recessive
    does not and gets no note."""
    trio_web_store(FIX / "trio.vcf.gz", keep_reference=False)
    for category, noted in (
        ("denovo", True),
        ("dominant", True),
        ("comphet", True),
        ("recessive", False),
    ):
        d = client.get(
            "/api/trio",
            params={"proband": "CHILD", "category": category, "max_af": 0.3},
        ).json()
        assert "error" not in d, (category, d)
        assert ("--keep-reference" in d.get("note", "")) is noted, category


def test_trio_endpoint_skips_a_cohort_derived_af_gate_unless_asked(
    tmp_path, trio_web_store
):
    """A callset-AF trio (every private het has AF 1/6 or 1/3): with no
    max_af the gate is skipped with a note, an explicit max_af is honoured
    exactly as `db trio --max-af` does."""
    from tests.fixtures._vcfgen import TRIO_SAMPLES, trio_sites, write_vcf

    vcf = write_vcf(tmp_path / "cohort.vcf", TRIO_SAMPLES, trio_sites(30, af="cohort"))
    trio_web_store(vcf)

    auto = client.get("/api/trio", params={"proband": "CHILD"}).json()
    assert auto["af_gate"] == "skipped_cohort_derived"
    assert auto["row_count"] == 15  # every de novo site survives
    assert "own allele frequency" in auto["note"] and "skipped" in auto["note"]

    asked = client.get("/api/trio", params={"proband": "CHILD", "max_af": 0.01}).json()
    assert asked["af_gate"] == "cohort_derived_suspected"
    assert asked["row_count"] == 0  # the gate ran, as asked
    assert "applied as requested" in asked["note"]

    # the page leaves max AF blank by default and sends it only when filled
    from vcfclick_web.page import INDEX_HTML

    assert 'id="t-af" type="number" step="0.01" placeholder' in INDEX_HTML
    assert "if($('t-af').value !== '') q.max_af" in INDEX_HTML
