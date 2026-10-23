"""The trio JSON renderer fetches genes, gnomAD and ClinVar for a model's
candidates in one batch per source. Its output must equal the per-candidate
lookups it replaced, with every source loaded and with none."""

from __future__ import annotations

import json

import pytest

from tests.conftest import FIXTURES, run_cli

# Annotations for alleles of the related.vcf fixture, written by a fresh
# interpreter: genes covering the planted de novo sites and a stretch of
# chr1, gnomAD for some alleles (one common, one rare), ClinVar for one.
_SEED = """
import sys
import duckdb
from annotations.db import SCHEMA_DDL

conn = duckdb.connect(sys.argv[1])
conn.execute(SCHEMA_DDL)
conn.executemany(
    "INSERT INTO refseq_genes VALUES (?, 'chr1', ?, ?, '+', 'id', 'protein_coding')",
    [("GA", 150_990_000, 151_010_000), ("GB", 151_490_000, 151_510_000),
     ("GC", 1, 50_000_000), ("GD", 40_000_000, 90_000_000)],
)
alleles = []
for line in open(sys.argv[2]):
    if line.startswith("#"):
        continue
    f = line.split("\t", 5)
    alleles.append((f[0], int(f[1]), f[3], f[4]))
gn = [(c, p, r, a, 0.2 if i % 3 == 0 else 0.0001, None if i % 2 else 0.3)
      for i, (c, p, r, a) in enumerate(alleles[::2])]
conn.executemany("INSERT INTO gnomad_af VALUES (?, ?, ?, ?, ?, ?)", gn)
cv = [(c, p, r, a, "Pathogenic" if i % 2 else "Benign", "criteria_provided", f"VCV{i}", "x")
      for i, (c, p, r, a) in enumerate(alleles[::5])]
conn.executemany("INSERT INTO clinvar_variants VALUES (?, ?, ?, ?, ?, ?, ?, ?)", cv)
conn.close()
"""

# Render the same trio document twice in one interpreter: as shipped, and
# with every decorator's prefetch removed so each candidate goes to the
# store. Counting the per-row lookups proves the two runs differ.
_RENDER = """
import dataclasses
import json
import cli.trio_output as out
from cli.trio_analysis import Gates, Trio, TrioAnalysis, _has_reference_rows
from storage import get_session

sess = get_session("t")
analysis = TrioAnalysis(
    session=sess,
    trio=Trio("i1", "C1", "F", "M"),
    gates=Gates(20, 10, 0.5, 0.25, 0.75),
    has_reference=_has_reference_rows(sess),
    gnomad_max_af=None,
)
lookups = []
real_enrich = out._enrich
out._enrich = lambda *a: lookups.append(a) or real_enrich(*a)
batched = out._trio_json(analysis, "all", 1000)
n_batched = len(lookups)
out._DECORATORS[:] = [
    dataclasses.replace(d, prefetch=None) for d in out._DECORATORS
]
per_row = out._trio_json(analysis, "all", 1000)
print("DOCS " + json.dumps([batched, per_row, n_batched, len(lookups) - n_batched]))
"""


@pytest.mark.parametrize("seeded", (True, False), ids=("annotated", "bare"))
def test_batched_annotations_equal_per_candidate_lookups(
    seeded, vcfclick_home, run_python, tmp_path
):
    vcf = FIXTURES / "related.vcf"
    store = tmp_path / "ann.duckdb"
    env = {"VCFCLICK_ANNOTATIONS_DB": str(store), "VCFCLICK_BACKEND": "duckdb"}
    if seeded:
        run_python(vcfclick_home, _SEED, str(store), str(vcf), **env)
    run_cli(vcfclick_home, "duckdb", "db", "create", "t")
    run_cli(
        vcfclick_home, "duckdb", "db", "ingest", "t", str(vcf),
        "--ingest-id", "i1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(
        vcfclick_home, "duckdb", "db", "ped", "t",
        "--proband", "C1", "--father", "F", "--mother", "M",
    )  # fmt: skip
    out = run_python(vcfclick_home, _RENDER, VCFCLICK_DB_NAME="t", **env)
    line = [ln for ln in out.splitlines() if ln.startswith("DOCS ")][-1]
    batched, per_row, n_batched, n_per_row = json.loads(line.removeprefix("DOCS "))
    assert batched == per_row
    assert batched["annotations"] == {
        "genes": seeded,
        "clinvar": seeded,
        "gnomad": seeded,
    }
    cands = [c for m in batched["models"].values() for c in m.get("candidates", [])]
    assert cands, "the fixture should yield candidates at --max-af 0.5"
    assert n_batched == 0  # every candidate was served by the batch
    assert n_per_row == len(cands)  # and without it, each asked the store
    if seeded:  # the comparison is only meaningful if annotations land
        assert any(c["genes"] for c in cands)
        assert any(c["gnomad_popmax"] is not None for c in cands)
        assert any(c["clinvar"] is not None for c in cands)
