"""A trio analysis runs each model's statement once per run, however many
times the renderers ask for its rows, count or genes."""

from __future__ import annotations

import json

import pytest

from tests.conftest import FIXTURES, run_cli

_RUN = """
import collections
import dataclasses
import json
import cli.trio_output as out
from cli.trio_analysis import Gates, Trio, TrioAnalysis, _has_reference_rows
from storage import get_session

real = get_session("t")


class Spy:
    def __init__(self):
        self.seen = []

    def query(self, sql, *a, **k):
        self.seen.append(sql)
        return real.query(sql, *a, **k)


spy = Spy()
analysis = TrioAnalysis(
    session=spy,
    trio=Trio("i1", "XDAU", "XDAD", "XMUM"),
    gates=Gates(20, 10, 0.01, 0.25, 0.75),
    has_reference=_has_reference_rows(real),
    gnomad_max_af=None,
)
doc = out._trio_json(analysis, "all", 1000)
repeated = [s for s, n in collections.Counter(spy.seen).items() if n > 1]

rows = analysis.detail_rows("denovo")
before = len(spy.seen)
again = (analysis.detail_rows("denovo"), analysis.count("denovo"),
         analysis.comphet_genes())
after = len(spy.seen)

# a copy under other gates must not read this run's answers
looser = dataclasses.replace(analysis, gates=analysis.gates._replace(min_gq=0))
looser.detail_rows("denovo")
print("OUT " + json.dumps({
    "repeated": repeated,
    "same": again[0] == rows and again[1] == len(rows),
    "new_statements": after - before,
    "looser_ran": len(spy.seen) - after,
    "models": sorted(doc["models"]),
}))
"""


@pytest.mark.parametrize("backend", ["chdb", "duckdb"])
def test_json_run_issues_each_statement_once(vcfclick_home, run_python, backend):
    if backend == "chdb":
        pytest.importorskip("chdb")
    run_cli(vcfclick_home, backend, "db", "create", "t")
    run_cli(
        vcfclick_home, backend, "db", "ingest", "t", str(FIXTURES / "xlinked.vcf"),
        "--ingest-id", "i1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(vcfclick_home, backend, "db", "ped", "t", str(FIXTURES / "xlinked.ped"))
    out = run_python(
        vcfclick_home, _RUN, VCFCLICK_DB_NAME="t", VCFCLICK_BACKEND=backend
    )
    got = json.loads([ln for ln in out.splitlines() if ln.startswith("OUT ")][-1][4:])
    assert "comphet" in got["models"] and "denovo" in got["models"]
    assert got["repeated"] == []
    assert got["same"] is True
    assert got["new_statements"] == 0
    assert got["looser_ran"] >= 1
