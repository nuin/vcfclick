"""`db query --format X` prints the same bytes on DuckDB as on chDB.

chDB's own renderer is the reference. The DuckDB session renders the rows
itself, and used to differ: JSON came out as arrays with no types, the
*WithNames formats dropped the header, TSVWithNames and JSONEachRow fell
through to a box table, every float went through `.6g` (1234567.25 printed
as 1.23457e+06), CSV left strings unquoted, TSV escaped nothing, and an
unknown format name printed a box table instead of failing.

Three things cannot match and are normalised before the comparison: the
`statistics` block (elapsed time, and rows/bytes read, which only chDB
knows), `rows_before_limit_at_least` (chDB's LIMIT bookkeeping), and the
`Nullable(...)` / `LowCardinality(...)` wrappers of a type name (DuckDB
does not carry nullability or dictionary encoding in a result type). The
bare type underneath must match. Pretty, PrettyCompact and Vertical are
not compared: chDB's box tables depend on terminal width and block
boundaries.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess

import pytest

from tests.conftest import CLI_TIMEOUT_S, FIXTURES, REPO, run_cli

BACKENDS = ("chdb", "duckdb")

FORMATS = (
    "JSON",
    "JSONCompact",
    "JSONEachRow",
    "TSV",
    "TabSeparated",
    "TSVWithNames",
    "TabSeparatedWithNames",
    "CSV",
    "CSVWithNames",
)

# Real columns of a real cohort: LowCardinality and Nullable strings,
# UInt32, Float32 (REAL), a Map, and literals the old renderer mangled.
TABLE_SQL = {
    b: (
        "SELECT chrom, pos, ref, alt, qual, filter, info_AF, info_extra, "
        "CAST(1234567.25 AS DOUBLE) AS d, 'x,y' AS s "
        "FROM variants ORDER BY pos"
    )
    for b in BACKENDS
}

GENOTYPE_SQL = {
    b: (
        "SELECT sample_id, pos, gt, gq, dp, ad_ref, ad_alt, "
        "CAST(ad_alt AS DOUBLE) / (ad_ref + ad_alt) AS ab "
        "FROM genotypes ORDER BY pos, sample_id LIMIT 4"
    )
    for b in BACKENDS
}

EMPTY_SQL = {b: "SELECT chrom, pos FROM variants WHERE pos = 0" for b in BACKENDS}

# One row of values that exercise escaping and number formatting. The SQL
# differs per engine only where the dialects do (escapes, nan, typed
# empties); the values are the same.
EDGE_SQL = {
    "chdb": (
        "SELECT CAST(5 AS UInt64) AS c, CAST(-3 AS Int64) AS i, "
        "concat('a', char(9), 'b', char(10), 'c\\\\d''e\"f/g', char(1), "
        "char(13), char(0), char(8), char(12), char(127), char(226, 128, 168)) AS s, "
        "'é' AS u, CAST(1e20 AS Float64) AS big, CAST(1e21 AS Float64) AS e21, "
        "CAST(1e-6 AS Float64) AS e6, CAST(1e-7 AS Float64) AS e7, "
        "CAST(0.1 AS Float64) + CAST(0.2 AS Float64) AS p, "
        "CAST(100 AS Float64) AS hundred, CAST(-0. AS Float64) AS nz, "
        "CAST(1.5e-8 AS Float32) AS f32, CAST(3.4e38 AS Float32) AS f32b, "
        "CAST(0.33 AS Float32) AS f33, "
        "nan AS nn, inf AS ii, -inf AS ni, true AS b, false AS bf, "
        "toDate('2024-01-02') AS dt, toDateTime('2024-01-02 03:04:05') AS ts, "
        "map('k', 'v', 'x', 'y\"''\\\\') AS m, map('k', ['a', 'b''c']) AS mm, "
        "CAST(map() AS Map(String, String)) AS em, "
        "CAST([] AS Array(String)) AS ea, "
        "CAST([NULL, 1] AS Array(Nullable(Int32))) AS an, "
        "[CAST(0.1 AS Float32)] AS af, "
        "CAST(NULL AS Nullable(String)) AS ns, 'q' AS `col\tname`"
    ),
    "duckdb": (
        "SELECT CAST(5 AS UBIGINT) AS c, CAST(-3 AS BIGINT) AS i, "
        "concat('a', chr(9), 'b', chr(10), 'c\\d''e\"f/g', chr(1), "
        "chr(13), chr(0), chr(8), chr(12), chr(127), chr(8232)) AS s, "
        "'é' AS u, CAST(1e20 AS DOUBLE) AS big, CAST(1e21 AS DOUBLE) AS e21, "
        "CAST(1e-6 AS DOUBLE) AS e6, CAST(1e-7 AS DOUBLE) AS e7, "
        "CAST(0.1 AS DOUBLE) + CAST(0.2 AS DOUBLE) AS p, "
        "CAST(100 AS DOUBLE) AS hundred, CAST('-0.0' AS DOUBLE) AS nz, "
        "CAST(1.5e-8 AS FLOAT) AS f32, CAST(3.4e38 AS FLOAT) AS f32b, "
        "CAST(0.33 AS FLOAT) AS f33, "
        "CAST('nan' AS DOUBLE) AS nn, CAST('inf' AS DOUBLE) AS ii, "
        "CAST('-inf' AS DOUBLE) AS ni, true AS b, false AS bf, "
        "DATE '2024-01-02' AS dt, TIMESTAMP '2024-01-02 03:04:05' AS ts, "
        "MAP {'k': 'v', 'x': 'y\"''\\'} AS m, MAP {'k': ['a', 'b''c']} AS mm, "
        "CAST(MAP {} AS MAP(VARCHAR, VARCHAR)) AS em, "
        "CAST([] AS VARCHAR[]) AS ea, CAST([NULL, 1] AS INTEGER[]) AS an, "
        "[CAST(0.1 AS FLOAT)] AS af, "
        "CAST(NULL AS VARCHAR) AS ns, 'q' AS \"col\tname\""
    ),
}

QUERIES = {
    "table": TABLE_SQL,
    "genotypes": GENOTYPE_SQL,
    "empty": EMPTY_SQL,
    "edge": EDGE_SQL,
}

_STATISTICS = re.compile(rb',\n\n\t"statistics":\n\t\{\n(?:\t\t[^\n]*\n)*\t\}\n')
_BEFORE_LIMIT = re.compile(rb',\n\n\t"rows_before_limit_at_least": \d+')
_WRAPPER = re.compile(rb"(?:Nullable|LowCardinality)\(([^()]*)\)")


def _normalise(out: bytes) -> bytes:
    out = _STATISTICS.sub(b"\n", out)
    out = _BEFORE_LIMIT.sub(b"", out)
    prev = None
    while prev != out:
        prev, out = out, _WRAPPER.sub(rb"\1", out)
    return out


@pytest.fixture(scope="module")
def homes(tmp_path_factory):
    """One tiny.vcf.gz cohort per backend, shared by the module's tests
    (they only read it)."""
    out = {}
    for backend in BACKENDS:
        home = tmp_path_factory.mktemp(f"fmt-{backend}")
        run_cli(home, backend, "db", "create", "t")
        run_cli(
            home,
            backend,
            "db",
            "ingest",
            "t",
            str(FIXTURES / "tiny.vcf.gz"),
            "--ingest-id",
            "i1",
            "--serial",
        )
        out[backend] = home
    return out


def _query(home, backend, sql, fmt, ok=True) -> subprocess.CompletedProcess:
    """`db query` with stdout kept as bytes: run_cli decodes with universal
    newlines, which would turn a raw carriage return in CSV into a newline."""
    exe = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
    env = {**os.environ, "VCFCLICK_HOME": str(home), "VCFCLICK_BACKEND": backend}
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [exe, "db", "query", "t", sql, "--format", fmt],
        cwd=REPO,
        env=env,
        capture_output=True,
        timeout=CLI_TIMEOUT_S,
    )
    if ok:
        assert r.returncode == 0, r.stderr.decode(errors="replace")
    return r


@pytest.mark.parametrize("query", sorted(QUERIES))
@pytest.mark.parametrize("fmt", FORMATS)
def test_db_query_prints_chdb_bytes_on_duckdb(homes, fmt, query):
    sql = QUERIES[query]
    chdb = _query(homes["chdb"], "chdb", sql["chdb"], fmt).stdout
    duck = _query(homes["duckdb"], "duckdb", sql["duckdb"], fmt).stdout
    assert _normalise(duck) == _normalise(chdb)


def test_json_keeps_chdb_types_rows_and_objects(homes):
    """The desktop apps read `meta[].name/type` and `data` as objects."""
    import json

    out = _query(
        homes["duckdb"],
        "duckdb",
        "SELECT chrom, pos FROM variants ORDER BY pos",
        "JSON",
    ).stdout
    doc = json.loads(out)
    assert doc["meta"] == [
        {"name": "chrom", "type": "String"},
        {"name": "pos", "type": "UInt32"},
    ]
    assert doc["data"][0] == {"chrom": "chr1", "pos": 100}
    assert doc["rows"] == 5
    assert set(doc["statistics"]) == {"elapsed"}


@pytest.mark.parametrize("backend", BACKENDS)
def test_unknown_format_is_an_error(homes, backend):
    r = _query(homes[backend], backend, "SELECT 1 AS x", "NoSuchFormat", ok=False)
    assert r.returncode != 0
    assert b"NoSuchFormat" in r.stderr


@pytest.mark.parametrize(
    ("value", "single", "text"),
    [
        (1234567.25, False, "1234567.25"),
        (0.1 + 0.2, False, "0.30000000000000004"),
        (100.0, False, "100"),
        (-0.0, False, "-0"),
        (1e20, False, "100000000000000000000"),
        (1e21, False, "1e21"),
        (1e-6, False, "0.000001"),
        (1e-7, False, "1e-7"),
        (0.33000001311302185, True, "0.33"),
        (1.5e-8, True, "1.5e-8"),
        (3.4e38, True, "3.4e38"),
        (float("nan"), False, "nan"),
        (float("-inf"), False, "-inf"),
    ],
)
def test_float_text_is_clickhouse_notation(value, single, text):
    """Shortest round-trip digits; fixed point for decimal exponents -6..20.
    Checked against chDB on 3000 random Float64 and Float32 values."""
    from storage._duckdb_format import float_text

    assert float_text(value, single=single) == text


def test_fast_float_text_matches_the_general_path():
    """The per-column Float64 / Float32 formatters print what `float_text`
    prints, on 20,000 random values across the fixed-point and exponent
    ranges and their edges."""
    import random
    import struct

    from storage._duckdb_format import _double_text, _float32_text, float_text

    rng = random.Random(7)
    edges = [1e-4, 9.999e-5, 1e16, 9.999999e15, 1e-6, 1e21, 0.5, 100.0, -0.0]
    values = edges + [
        rng.uniform(-1, 1) * 10 ** rng.uniform(-9, 23) for _ in range(20_000)
    ]
    for v in values:
        assert _double_text(v) == float_text(v), v
        f32 = struct.unpack("f", struct.pack("f", v))[0]
        assert _float32_text(f32) == float_text(f32, single=True), f32


def test_per_column_formatters_render_every_column_kind():
    """Integer, bool, float, string and NULL columns take the per-column
    path; list, map, struct and decimal take the general one. Non-finite
    floats, `/`, U+2028, quotes and newlines are each escaped as chDB does."""
    duckdb = pytest.importorskip("duckdb")
    from storage._duckdb_format import render

    cur = duckdb.connect().execute(
        "SELECT 7::INTEGER AS i, true AS b, 'nan'::DOUBLE AS d, "
        "CAST(0.33 AS REAL) AS r, 'a/b\"c' || chr(10) || chr(8232) AS s, "
        "NULL::VARCHAR AS n, [1, NULL] AS l, MAP {'k': 1.5} AS m, "
        "{'x': 'q,'} AS st, 1.50::DECIMAL(4, 2) AS dec "
        "UNION ALL SELECT NULL, NULL, 1e21, NULL, NULL, 'x', NULL, NULL, NULL, NULL"
    )
    cols = [d[0] for d in cur.description]
    types = [d[1] for d in cur.description]
    rows = cur.fetchall()

    def out(fmt):
        return render(fmt, cols, types, rows, 0.5).decode()

    assert (
        '\t\t[7, true, null, 0.33, "a\\/b\\"c\\n\\u2028", null, [1,null], '
        '{"k":1.5}, {"x":"q,"}, 1.5],\n'
        '\t\t[null, null, 1e21, null, null, "x", null, null, null, null]\n'
    ) in out("JSONCompact")
    assert out("TSV") == (
        "7\ttrue\tnan\t0.33\ta/b\"c\\n\u2028\t\\N\t[1,NULL]\t{'k':1.5}\t('q,')\t1.5\n"
        "\\N\t\\N\t1e21\t\\N\t\\N\tx\t\\N\t\\N\t\\N\t\\N\n"
    )
    assert out("CSV") == (
        '7,true,nan,0.33,"a/b""c\n\u2028",\\N,"[1,NULL]","{\'k\':1.5}","q,",1.5\n'
        '\\N,\\N,1e21,\\N,\\N,"x",\\N,\\N,\\N,\\N\n'
    )
