"""storage.query_json / query_rows: one decoder for JSONCompact results, the
same on both engines."""

from __future__ import annotations

import json

import pytest

_SCRIPT = """
import json
from storage import get_session, query_json, query_rows

sess = get_session()
doc = query_json(sess, "SELECT 1 AS a, 'x' AS b UNION ALL SELECT 2 AS a, 'y' AS b")
rows = query_rows(sess, "SELECT 1 AS a, 'x' AS b")
print("OUT " + json.dumps({
    "names": [m["name"] for m in doc["meta"]],
    "data": sorted(doc["data"]),
    "rows": rows,
}))
"""


@pytest.mark.parametrize("backend", ["chdb", "duckdb"])
def test_query_json_and_rows_decode_both_engines(vcfclick_home, run_python, backend):
    if backend == "chdb":
        pytest.importorskip("chdb")
    (vcfclick_home / "dbs" / "q").mkdir(parents=True)
    out = run_python(
        vcfclick_home, _SCRIPT, VCFCLICK_DB_NAME="q", VCFCLICK_BACKEND=backend
    )
    line = [ln for ln in out.splitlines() if ln.startswith("OUT ")][-1]
    got = json.loads(line.removeprefix("OUT "))
    assert got["names"] == ["a", "b"]
    assert got["data"] == [[1, "x"], [2, "y"]]
    assert got["rows"] == [[1, "x"]]


_HAS_REF = """
import json
from storage import get_session, has_reference_rows

sess = get_session("t")
seen = []
real = sess.query


class Spy:
    def query(self, sql, *a, **k):
        seen.append(sql)
        return real(sql, *a, **k)


spy = Spy()
print("OUT " + json.dumps({
    "any": has_reference_rows(spy),
    "i1": has_reference_rows(spy, "i1"),
    "other": has_reference_rows(spy, "zzz"),
    "sql": seen,
}))
"""


@pytest.mark.parametrize("backend", ["chdb", "duckdb"])
@pytest.mark.parametrize("keep_reference", [True, False])
def test_has_reference_rows_stops_at_the_first_row(
    vcfclick_home, tiny_vcf, run_python, backend, keep_reference
):
    """A yes/no question: the statement is bounded by LIMIT 1 (a count would
    read every gt = 0 row of a --keep-reference store before answering)."""
    from tests.conftest import run_cli

    if backend == "chdb":
        pytest.importorskip("chdb")
    run_cli(vcfclick_home, backend, "db", "create", "t")
    args = ["db", "ingest", "t", str(tiny_vcf), "--ingest-id", "i1", "--serial"]
    if keep_reference:
        args.append("--keep-reference")
    run_cli(vcfclick_home, backend, *args)
    out = run_python(
        vcfclick_home, _HAS_REF, VCFCLICK_DB_NAME="t", VCFCLICK_BACKEND=backend
    )
    got = json.loads([ln for ln in out.splitlines() if ln.startswith("OUT ")][-1][4:])
    assert (got["any"], got["i1"], got["other"]) == (
        keep_reference,
        keep_reference,
        False,
    )
    assert all("LIMIT 1" in s and "count" not in s.lower() for s in got["sql"])
