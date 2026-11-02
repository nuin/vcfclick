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
    "INSTALL httpfs",
    "EXPORT DATABASE '/tmp/x'",
    "CHECKPOINT",
    "VACUUM",
    "BEGIN",
    "PREPARE p AS DELETE FROM variants",
    "EXECUTE p",
    "KILL QUERY WHERE 1",
    "RENAME TABLE a TO b",
    "SELECT 1 /* */ ; DROP TABLE variants",
    "SELECT 1 UNION ALL SELECT 2 INTO OUTFILE '/tmp/x'",
    "EXPLAIN SELECT 1 INTO OUTFILE '/tmp/x'",
    "EXPLAIN PLAN SELECT 1 INTO OUTFILE '/tmp/x'",
    "DESCRIBE TABLE variants INTO OUTFILE '/tmp/x'",
    # DuckDB's EXPLAIN ANALYZE executes the statement it explains
    "EXPLAIN ANALYZE INSERT INTO variants SELECT * FROM variants",
    "EXPLAIN ANALYZE DELETE FROM variants",
    "SUMMARIZE variants; DROP TABLE variants",
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


# Statements that read something other than the cohort tables: a file, a
# URL, the environment, the server's own state, or an unbounded generator.
# The guard allows a FROM source only when it is a plain table identifier
# (optionally schema-qualified with main / default).
EXTERNAL_READS = [
    "SELECT * FROM file('/etc/hosts', 'LineAsString')",
    "SELECT * FROM url('http://169.254.169.254/latest', CSV)",
    "SELECT * FROM s3('http://x/y', CSV)",
    "SELECT * FROM remote('h', db.t)",
    "SELECT * FROM mysql('h:3306', 'db', 't', 'u', 'p')",
    "SELECT * FROM executable('cat /etc/hosts', TabSeparated, 'x String')",
    "SELECT content FROM read_text('/etc/hosts')",
    "SELECT * FROM read_blob('/etc/hosts')",
    "SELECT * FROM read_csv('/tmp/x.csv')",
    "SELECT * FROM read_csv_auto('/tmp/x.csv')",
    "SELECT * FROM read_parquet('/tmp/x.parquet')",
    "SELECT * FROM read_json_auto('/tmp/x.json')",
    "SELECT * FROM sniff_csv('/tmp/x.csv')",
    "SELECT * FROM glob('/Users/*/.ssh/*')",
    "SELECT * FROM '/etc/hosts'",
    "SELECT * FROM 'x.parquet'",
    'SELECT * FROM "x.parquet"',
    "SELECT * FROM numbers(1000000000000)",
    "SELECT * FROM range(1000000000000)",
    "SELECT * FROM generate_series(1, 1000000000000)",
    "SELECT * FROM duckdb_settings()",
    "SELECT * FROM duckdb_tables()",
    "SELECT * FROM pragma_database_size()",
    # scalar file / environment functions
    "SELECT file('/etc/passwd')",
    "SELECT getenv('HOME')",
    # the same sources below a plain cohort query
    "SELECT * FROM variants WHERE chrom IN (SELECT * FROM file('/x', 'CSV'))",
    "SELECT * FROM variants v JOIN read_csv('/tmp/x.csv') r ON v.pos = r.pos",
    "WITH f AS (SELECT * FROM read_csv('/tmp/x.csv')) SELECT * FROM f",
    "SELECT 1 UNION ALL SELECT * FROM numbers(3)",
    "EXPLAIN SELECT * FROM numbers(3)",
    "EXPLAIN SELECT * FROM file('/etc/hosts', 'LineAsString')",
    "DESCRIBE SELECT * FROM read_csv('/tmp/x.csv')",
    "DESCRIBE TABLE file('/etc/hosts', 'LineAsString')",
    # another database, catalog or schema
    "SELECT * FROM store.main.variants",
    "SELECT * FROM information_schema.tables",
    "SELECT * FROM pg_catalog.pg_tables",
]

# Server introspection and per-query setting overrides.
SYSTEM_READS = [
    "SELECT * FROM system.processes",
    "SELECT query FROM system.processes",
    "SELECT * FROM system.tables",
    "DESCRIBE system.processes",
    "SHOW TABLES FROM system",
    "SHOW CREATE TABLE system.processes",
    "SHOW PROCESSLIST",
    "SELECT 1 SETTINGS readonly = 0",
    "SELECT 1 SETTINGS max_execution_time = 0",
    "SELECT 1 SETTINGS max_threads = 2",
    "SELECT count() FROM variants SETTINGS max_result_rows = 0 FORMAT JSONCompact",
    "SELECT * FROM (SELECT 1 SETTINGS readonly = 0)",
]

# Cohort reads the stricter source rule must keep accepting.
COHORT_READS = [
    "SELECT chrom, pos FROM variants WHERE pos > 100 ORDER BY pos LIMIT 10",
    "WITH hi AS (SELECT chrom, pos FROM variants WHERE qual > 30) "
    "SELECT count(*) FROM hi",
    "SELECT v.chrom, g.sample_id FROM variants v "
    "JOIN genotypes g ON v.chrom = g.chrom AND v.pos = g.pos",
    "SELECT * FROM variants v, samples s WHERE v.ingest_id = s.ingest_id",
    "SELECT * FROM variants AS v LEFT JOIN genotypes AS g USING (ingest_id, chrom, pos)",
    "SELECT * FROM variants WHERE pos IN (SELECT pos FROM genotypes)",
    'SELECT * FROM "variants" LIMIT 1',
    "SELECT * FROM (SELECT 1 AS n) t",
    "SELECT * FROM (VALUES (1), (2)) AS t(x)",
    "SELECT count(*) FROM genotypes GROUP BY sample_id",
    "DESCRIBE variants",
    "SHOW TABLES",
    "SHOW CREATE TABLE variants",
    "EXPLAIN SELECT chrom FROM variants JOIN samples USING (ingest_id)",
]


@pytest.mark.parametrize("dialect", DIALECTS)
@pytest.mark.parametrize("sql", EXTERNAL_READS)
def test_blocks_file_url_and_generator_sources(sql, dialect):
    assert is_read_only(sql, dialect) is False


@pytest.mark.parametrize("dialect", DIALECTS)
@pytest.mark.parametrize("sql", SYSTEM_READS)
def test_blocks_system_tables_and_settings(sql, dialect):
    assert is_read_only(sql, dialect) is False


@pytest.mark.parametrize("dialect", DIALECTS)
@pytest.mark.parametrize("sql", COHORT_READS)
def test_allows_cohort_reads(sql, dialect):
    assert is_read_only(sql, dialect) is True


def test_allows_default_schema_qualifier():
    assert is_read_only("SELECT * FROM main.variants", "duckdb") is True
    assert is_read_only("SELECT * FROM default.variants", "clickhouse") is True


@pytest.mark.parametrize(
    "sql",
    [
        "SELECT count(*) FROM genotypes FORMAT JSONCompact",
        "SELECT arrayJoin([1, 2]) AS x",
        "EXPLAIN SYNTAX SELECT 1",
        "SHOW TABLES LIKE '%into outfile%'",
    ],
)
def test_allows_clickhouse_reads(sql):
    assert is_read_only(sql, "clickhouse") is True


@pytest.mark.parametrize(
    "sql",
    [
        "SUMMARIZE variants",
        "SUMMARIZE SELECT chrom, pos FROM variants",
        "EXPLAIN ANALYZE SELECT 1",
        "FROM variants SELECT chrom LIMIT 3",
        "SELECT * EXCLUDE (qual) FROM variants",
        "SELECT count(*) FILTER (WHERE gt = 2) FROM genotypes GROUP BY ALL",
    ],
)
def test_allows_duckdb_reads(sql):
    assert is_read_only(sql, "duckdb") is True


@pytest.mark.parametrize(
    "sql",
    [
        "SELECT count() FROM variants FINAL",
        "SELECT quantile(0.5)(gq), countIf(gt = 1) FROM genotypes",
        "SELECT * FROM variants PREWHERE chrom = 'chr1' LIMIT 10 BY chrom",
        "SELECT x FROM (SELECT [1, 2] AS a) ARRAY JOIN a AS x",
        "SELECT arrayMap(x -> x * 2, [1, 2])",
        "EXPLAIN PIPELINE SELECT 1",
        "DESC variants",
    ],
)
def test_allows_more_clickhouse_reads(sql):
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

    monkeypatch.setattr(server, "get_readonly_session", boom)
    out = server.run_sql("DROP TABLE variants")
    assert "read-only" in out["error"]
    assert out["sql"] == "DROP TABLE variants"
    out = server.run_sql("SELECT * FROM file('/etc/hosts', 'LineAsString')")
    assert "read-only" in out["error"]
