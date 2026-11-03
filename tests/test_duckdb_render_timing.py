"""Timing gate for the DuckDB text renderer (`storage._duckdb_format`).

Internal callers parse whole tables out of JSONCompact on DuckDB, for
example `cli.relatedness_analysis._load`, which reads every SNV site and
every genotype row. A renderer that dispatched on each value's type made
that 5 times slower: 200,000 rows x 5 columns took 0.8 s (strings and
integers) and 1.6 s (with Float64 and Float32 columns) on Apple silicon.
With one formatter per column the same renders take 0.15 s and 0.33 s.
Each budget sits between the two measurements, about 3 times the current
one, so a slow runner passes and the old cost does not.
"""

from __future__ import annotations

import time

import pytest

pytestmark = pytest.mark.timing

ROWS = 200_000

CASES = {
    # The column types of the relatedness genotype query.
    "strings_and_ints": (
        "SELECT 'chr' || (i % 22)::VARCHAR AS chrom, (i * 7)::UINTEGER AS pos, "
        "'A' AS ref, 'C' AS alt, (i % 3)::TINYINT AS gt FROM range({n}) t(i)",
        0.5,
    ),
    "with_floats": (
        "SELECT 'chr' || (i % 22)::VARCHAR AS chrom, (i * 7)::UINTEGER AS pos, "
        "(i / 7.0)::DOUBLE AS d, (i / 3.0)::REAL AS r, (i % 3)::TINYINT AS gt "
        "FROM range({n}) t(i)",
        1.0,
    ),
}


@pytest.mark.parametrize("fmt", ["JSONCompact", "TSV"])
@pytest.mark.parametrize("case", sorted(CASES))
def test_duckdb_render_of_200k_rows_is_within_budget(case, fmt):
    duckdb = pytest.importorskip("duckdb")
    from storage._duckdb_format import render

    sql, budget = CASES[case]
    cur = duckdb.connect().execute(sql.format(n=ROWS))
    cols = [d[0] for d in cur.description]
    types = [d[1] for d in cur.description]
    rows = cur.fetchall()

    best = float("inf")
    for _ in range(3):
        started = time.perf_counter()
        render(fmt, cols, types, rows, 0.0)
        best = min(best, time.perf_counter() - started)
    assert best < budget, f"{case} {fmt}: {best:.2f} s, budget {budget} s"
