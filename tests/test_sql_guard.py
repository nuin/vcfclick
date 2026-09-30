"""Tests for storage.sql_guard — the AST read-only check shared by the
web UI and the MCP `run_sql` tool. Pure parsing: no database needed."""

from __future__ import annotations

import pytest

from storage.sql_guard import is_read_only

DIALECTS = ["clickhouse", "duckdb"]

WRITES = [
    "",
    "   ",
    "DROP TABLE variants",
    "/*x*/ DROP TABLE variants",
    "-- c\nDROP TABLE variants",
    "SELECT 1; DROP TABLE variants",
    "SELECT 1; SELECT 2",
    "WITH c AS (SELECT 1) DELETE FROM variants",
    "WITH c AS (SELECT 1) INSERT INTO variants SELECT * FROM c",
    "SELECT 42 INTO OUTFILE '/tmp/x'",
    "SELECT * INTO t2 FROM variants",
    "COPY (SELECT 42) TO '/tmp/x'",
    "INSERT INTO variants VALUES (1)",
    "   update samples SET x = 1",
    "DELETE FROM variants WHERE 1",
    "CREATE TABLE t AS SELECT 1",
    "ALTER TABLE variants DROP COLUMN qual",
    "TRUNCATE TABLE variants",
    "ATTACH 'evil.db'",
    "SET max_threads = 1",
    "USE other",
    "PRAGMA enable_profiling",
    "OPTIMIZE TABLE variants",
    "SYSTEM SHUTDOWN",
    "EXPLAIN DELETE FROM variants",
    "EXPLAIN EXPLAIN SELECT 1",
    "SELECT FROM WHERE",  # unparseable → fail closed
    # ClickHouse SHOW accepts INTO OUTFILE, which writes a file
    "SHOW TABLES INTO OUTFILE '/tmp/x'",
    "show create table variants into outfile '/tmp/x' FORMAT TSV",
    "LOAD httpfs",
    "CREATE SECRET s (TYPE s3)",
]

READS = [
    "SELECT * FROM variants LIMIT 5",
    "select chrom, pos from variants",
    "WITH c AS (SELECT 1 AS n) SELECT n FROM c",
    "SELECT 1 UNION ALL SELECT 2",
    "SHOW TABLES",
    "SHOW CREATE TABLE variants",
    "DESCRIBE variants",
    "EXPLAIN SELECT 1",
    "SELECT count(*) FROM genotypes  -- trailing comment\n",
    "/* leading */ SELECT 1;",
    "SELECT 'DROP TABLE variants; DELETE' AS s",
    "SELECT * FROM variants WHERE info_x = 'insert into'",
    # truncate()/replace() are functions, not statements — must pass
    "SELECT replace(chrom, 'chr', '') FROM variants",
    "SELECT truncate(qual) FROM variants",
    "SELECT chrom FROM variants WHERE pos IN (SELECT pos FROM genotypes)",
]


@pytest.mark.parametrize("dialect", DIALECTS)
@pytest.mark.parametrize("sql", WRITES)
def test_blocks_writes_and_bypasses(sql, dialect):
    assert is_read_only(sql, dialect) is False


@pytest.mark.parametrize("dialect", DIALECTS)
@pytest.mark.parametrize("sql", READS)
def test_allows_reads(sql, dialect):
    assert is_read_only(sql, dialect) is True


@pytest.mark.parametrize(
    "sql",
    [
        "SELECT count(*) FROM genotypes FORMAT JSONCompact",
        "SELECT 1 SETTINGS max_threads = 2",
        "SELECT arrayJoin([1, 2]) AS x",
        "EXPLAIN SYNTAX SELECT 1",
        "SHOW TABLES LIKE '%into outfile%'",
    ],
)
def test_allows_clickhouse_reads(sql):
    assert is_read_only(sql, "clickhouse") is True


def test_blocks_clickhouse_insert_into_function():
    assert (
        is_read_only("INSERT INTO FUNCTION file('x') SELECT 1", "clickhouse") is False
    )


def test_dialect_follows_backend(monkeypatch):
    monkeypatch.setenv("VCFCLICK_BACKEND", "duckdb")
    assert is_read_only("SELECT 1") is True
    assert is_read_only("DROP TABLE variants") is False
    monkeypatch.setenv("VCFCLICK_BACKEND", "chdb")
    assert is_read_only("SELECT 1 FORMAT JSONCompact") is True


def test_mcp_run_sql_refuses_writes_without_touching_db(monkeypatch):
    from vcfclick_mcp import server

    def boom():
        raise AssertionError("get_session must not be called for a write")

    monkeypatch.setattr(server, "get_session", boom)
    out = server.run_sql("DROP TABLE variants")
    assert "read-only" in out["error"]
    assert out["sql"] == "DROP TABLE variants"
