"""The bulk genotype reads (relatedness, the export's per-window calls) come
back as Arrow, not as JSONCompact text the DuckDB session renders and the
caller parses again. The answers must not change."""

from __future__ import annotations

import json

import pytest

from tests.conftest import FIXTURES, run_cli

_RUN = """
import json

import cli.main  # registers the commands; cli.db_diff needs it loaded
import numpy as np
from storage import get_session

real = get_session("t")
via = {}


class Spy:
    def query(self, sql, fmt="CSV"):
        via.setdefault(fmt, []).append(sql)
        return real.query(sql, fmt)

    def query_arrow(self, sql):
        via.setdefault("arrow", []).append(sql)
        return real.query_arrow(sql)


from cli.relatedness_analysis import _load
from export.vcf_context import ExportRequest
from export.vcf_query import _window_calls, _windows, iter_records, prepare_export

samples, g, mode, total, span = _load(Spy(), "i1", 10_000)
json_genotypes = [s for s in via.get("JSONCompact", []) if "FROM genotypes" in s
                  and "LIMIT 1" not in s]
arrow_genotypes = [s for s in via.get("arrow", []) + via.get("ArrowStream", [])
                   if "FROM genotypes" in s]
ctx = prepare_export(real, ExportRequest(db_name="t"))
via.clear()
calls = {}
for region in _windows(real, ctx):
    calls.update({f"{region[0]}:{k[0]}:{k[1]}:{k[2]}": v
                  for k, v in _window_calls(Spy(), ctx, region).items()})
export_json = [s for s in via.get("JSONCompact", []) if "FROM genotypes" in s]
records = list(iter_records(real, ctx))
matched = sum(1 for _rec, c in records if c)
with_calls = real.query(
    "SELECT count(*) FROM (SELECT DISTINCT chrom, pos, ref, alt FROM ("
    "SELECT chrom, pos, ref, alt FROM genotypes UNION ALL "
    "SELECT chrom, pos, ref, alt FROM missing_genotypes) u) d", "CSV"
).bytes().decode().strip()
print("OUT " + json.dumps({
    "samples": samples, "sum": int(np.abs(g).sum()), "shape": list(g.shape),
    "mode": mode, "total": total, "span": span,
    "json_genotypes": len(json_genotypes), "arrow_genotypes": len(arrow_genotypes),
    "export_json": len(export_json), "calls": calls,
    "records": len(records), "matched": matched, "with_calls": int(with_calls),
}, sort_keys=True))
"""


@pytest.mark.parametrize("backend", ["chdb", "duckdb"])
def test_bulk_genotype_reads_use_arrow(vcfclick_home, run_python, backend):
    if backend == "chdb":
        pytest.importorskip("chdb")
    run_cli(vcfclick_home, backend, "db", "create", "t")
    run_cli(
        vcfclick_home, backend, "db", "ingest", "t", str(FIXTURES / "related.vcf"),
        "--ingest-id", "i1", "--serial",
    )  # fmt: skip
    out = run_python(
        vcfclick_home, _RUN, VCFCLICK_DB_NAME="t", VCFCLICK_BACKEND=backend
    )
    got = json.loads([ln for ln in out.splitlines() if ln.startswith("OUT ")][-1][4:])
    assert got["json_genotypes"] == 0 and got["arrow_genotypes"] == 1
    assert got["export_json"] == 0
    assert got["total"] > 0 and got["sum"] > 0 and got["calls"]
    # one value per call column, ints as ints (not strings)
    # (a recorded no-call is None: the export writes it ./.)
    first = next(v for c in got["calls"].values() for v in c.values() if v)
    gt, gq, dp, ad_ref, ad_alt, ft = first
    assert isinstance(gt, int)
    assert all("b'" not in key for key in got["calls"])  # str, never bytes
    # every site with stored calls or recorded no-calls found them (a site
    # of 0/0 calls has none)
    assert got["matched"] == got["with_calls"] > 0
