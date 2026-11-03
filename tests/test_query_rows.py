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
