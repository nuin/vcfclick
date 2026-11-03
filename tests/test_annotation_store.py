"""How the annotation store is opened.

Lookups open it read-only and run no DDL, so several processes (two `db
trio` runs, the web UI, the MCP server) read it at once, and a store on a
read-only file system or with mode 444 still answers. Only a loader (and
`ensure_schema`) writes. A store another process is writing is reported by
name, never read as "nothing loaded".
"""

from __future__ import annotations

import contextlib
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


def _seed(path: Path) -> None:
    import annotations.db as adb

    conn = adb.get_connection()
    conn.execute(
        "INSERT INTO refseq_genes VALUES "
        "('GENEA', 'chr1', 100, 200, '+', 'ENSG1', 'protein_coding')"
    )
    conn.close()
    assert path.exists()


@contextlib.contextmanager
def _held(path: Path, opener: str):
    """Another process with the store open through `annotations.db.<opener>`
    until the block ends."""
    code = (
        "import sys, annotations.db as adb\n"
        f"conn = adb.{opener}()\n"
        "print('ready', flush=True)\n"
        "sys.stdin.read()\n"
        "conn.close()\n"
    )
    proc = subprocess.Popen(
        [sys.executable, "-c", code],
        cwd=REPO,
        env={**os.environ, "VCFCLICK_ANNOTATIONS_DB": str(path)},
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        line = proc.stdout.readline()
        assert line.strip() == "ready", proc.stderr.read()
        yield proc.pid
    finally:
        proc.stdin.close()
        proc.wait(timeout=60)


def test_lookups_read_while_another_process_reads(isolated_annotation_db):
    import annotations.db as adb

    _seed(isolated_annotation_db)
    with _held(isolated_annotation_db, "read_connection"):
        assert adb.position_for_gene("GENEA").chrom == "chr1"
        assert [g.gene_symbol for g in adb.gene_at("chr1", 150)] == ["GENEA"]
        assert adb.clinvar_lookup("chr1", 150, "A", "G") is None
        assert adb.gnomad_af_many([("chr1", 150, "A", "G")]) == {}


def test_lookups_answer_from_a_read_only_file(isolated_annotation_db):
    import annotations.db as adb

    _seed(isolated_annotation_db)
    mode = isolated_annotation_db.stat().st_mode
    os.chmod(isolated_annotation_db, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
    try:
        assert adb.position_for_gene("GENEA").start_pos == 100
        assert adb.resolve_gene_tokens(["genea"])["GENEA"][0][0] == "symbol"
    finally:
        os.chmod(isolated_annotation_db, mode)


def test_lookups_run_no_ddl(isolated_annotation_db):
    """An older store without the optional tables is read as it is: the
    lookups see them empty and the file is not changed."""
    import duckdb

    import annotations.db as adb

    conn = duckdb.connect(str(isolated_annotation_db))
    conn.execute(
        "CREATE TABLE refseq_genes (gene_symbol VARCHAR PRIMARY KEY, chrom "
        "VARCHAR, start_pos UINTEGER, end_pos UINTEGER, strand VARCHAR, "
        "refseq_id VARCHAR, description VARCHAR)"
    )
    conn.execute(
        "INSERT INTO refseq_genes VALUES ('GENEA','chr1',100,200,'+','E1','pc')"
    )
    conn.close()

    assert adb.position_for_gene("GENEA") is not None
    assert adb.clinvar_lookup("chr1", 150, "A", "G") is None
    assert adb.gnomad_af("chr1", 150, "A", "G") is None
    assert adb.resolve_gene_tokens(["GENEA"])["GENEA"][0][0] == "symbol"
    from annotations.transcripts import cds_regions_for_gene, exon_at

    assert cds_regions_for_gene("GENEA") == []
    assert exon_at("chr1", 150) == []

    conn = duckdb.connect(str(isolated_annotation_db), read_only=True)
    tables = {r[0] for r in conn.execute("SHOW TABLES").fetchall()}
    conn.close()
    assert tables == {"refseq_genes"}


def test_lookups_on_a_missing_store_create_nothing(isolated_annotation_db):
    import annotations.db as adb

    assert not isolated_annotation_db.exists()
    assert adb.position_for_gene("GENEA") is None
    assert adb.gene_at("chr1", 1) == []
    assert adb.clinvar_many([("chr1", 1, "A", "G")]) == {}
    assert not isolated_annotation_db.exists()


def test_ensure_schema_creates_the_tables(isolated_annotation_db):
    import duckdb

    import annotations.db as adb

    conn = duckdb.connect(str(isolated_annotation_db))
    adb.ensure_schema(conn)
    tables = {r[0] for r in conn.execute("SHOW TABLES").fetchall()}
    conn.close()
    assert {"refseq_genes", "clinvar_variants", "gnomad_af", "gene_aliases"} <= tables
    assert {"transcripts", "exons", "cds"} <= tables


def test_a_store_being_written_is_reported_with_the_holder(isolated_annotation_db):
    import annotations.db as adb

    _seed(isolated_annotation_db)
    with _held(isolated_annotation_db, "get_connection") as pid:
        with pytest.raises(adb.AnnotationStoreUnavailable) as err:
            adb.position_for_gene("GENEA")
    message = str(err.value)
    assert str(isolated_annotation_db) in message
    assert f"PID {pid}" in message


def test_trio_refuses_when_the_store_is_being_written(tmp_path, isolated_annotation_db):
    """Two processes on one store: `db trio` must not print an unannotated
    table as if nothing were loaded (it used to exit 0 with empty genes)."""
    from tests.test_trio import _setup, _trio

    home = tmp_path / "home"
    _setup(home)
    _seed(isolated_annotation_db)
    with _held(isolated_annotation_db, "get_connection") as pid:
        r = _trio(home, "--format", "json", expect_failure=True)
        r2 = _trio(home, "--genes", "GENEA", expect_failure=True)
    for res in (r, r2):
        assert "Traceback" not in res.stderr
        assert "annotation store" in res.stderr
        assert f"PID {pid}" in res.stderr
