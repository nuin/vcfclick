"""No model can ignore a scope.

`cli.scope.SCOPE_AWARE_BUILDERS` holds, per model, the SQL its builder runs
under a scope. This file fails when a model in `cli.trio_output.MODELS` has no
builder there, when a builder's SQL does not carry the scope predicate, or when
a model's rows outside a planted scope come back from a real run (both
backends, on tests/fixtures/scope.vcf.gz, whose sites of every model lie on
both sides of the planted regions). A step that registers a model without its
scope builder fails here, and so does one whose builder forgets the predicate:
the Python backstop would hide it from the row check, not from this one.
"""

from __future__ import annotations

import json
import os
from unittest import mock

import pytest

import cli.main  # noqa: F401  (registers the commands and models)
from cli import scope as sc
from cli.trio_analysis import Gates, Trio
from cli.trio_output import MODELS, ModelSpec
from tests.conftest import FIXTURES, run_cli

BACKENDS = ("chdb", "duckdb")
TRUTH = json.loads((FIXTURES / "scope_truth.json").read_text())
# the planted scope: chr7 up to 60 kb (the panel genes, the region edges) and
# XGENE1 on chrX; everything at 70 kb and beyond on chr7, and chrX:60 Mb, is out
REGIONS = ("chr7:1-60000", "chrX:50000000-50100000")


def planted_scope(*, exclude: bool = True) -> sc.Scope:
    return sc.build_scope(
        sc.ScopeRequest(
            regions=REGIONS, exclude_regions=("chr7:29999-29999",) if exclude else ()
        )
    )


def models_missing_scope(models: dict[str, ModelSpec]) -> list[str]:
    return [name for name in models if name not in sc.SCOPE_AWARE_BUILDERS]


def test_every_model_has_a_scope_builder():
    assert models_missing_scope(MODELS) == []
    assert set(MODELS) <= set(sc.SCOPE_AWARE_BUILDERS)


def test_the_check_is_not_vacuous():
    """A model registered without a scope builder is caught."""
    spec = ModelSpec(
        "planted", lambda t: True, lambda r: {}, lambda a: None, lambda a: None
    )
    assert models_missing_scope({**MODELS, "planted": spec}) == ["planted"]


def test_the_reports_statements_are_registered_too():
    """The strict_missing check and the effective-gate report count rows under
    the scope; they are builders like the models, so the check below covers
    them."""
    assert {"strict_dropped", "effective_gates"} <= set(sc.SCOPE_AWARE_BUILDERS)


@pytest.mark.parametrize("mode", ["sql", "sql_buckets"])
@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("name", sorted(sc.SCOPE_AWARE_BUILDERS))
def test_every_builder_carries_the_scope_predicate(name, backend, mode, monkeypatch):
    monkeypatch.setenv("VCFCLICK_BACKEND", backend)
    if mode == "sql_buckets":
        monkeypatch.setattr(sc, "SCOPE_SQL_EXACT_MAX", 0)
    trio, gates = Trio("i1", "C1", "F", "M"), Gates(20, 10, 0.01, 0.25, 0.75)
    scope = planted_scope()
    assert scope.sql_mode == mode
    statements = sc.SCOPE_AWARE_BUILDERS[name](trio, gates, scope)
    assert statements
    for sql, part, alias in statements:
        if part == "stats":  # a count: exact in every mode
            pred, site = sc.stats_predicate(alias, scope)
            assert pred in sql and (site is None or site in sql), (name, sql[:200])
            assert (site is None) is (mode == "sql"), (name, mode)
            continue
        # the exclusions are in every statement; the inclusions in all but the
        # compound-het candidate queries, which apply them after grouping
        assert sc.scope_sql(alias, scope, part="exclude") in sql, (name, sql[:200])
        included = sc.scope_sql(alias, scope).split(" AND NOT ")[0]
        assert (included in sql) is (part == "all"), (name, sql[:200])
    # no scope: none of it
    for sql, _, _ in sc.SCOPE_AWARE_BUILDERS[name](trio, gates, None):
        assert "BETWEEN 50000000 AND 50100000" not in sql
        assert "50000000" not in sql


# ───────────────────────────── planted rows, real runs ─────────────────────────────


@pytest.fixture(scope="module", params=BACKENDS)
def runs(request, tmp_path_factory):
    """`db trio --format json` without and with the planted scope (exact SQL),
    and with it on the Python path alone and as the bucket superset."""
    from annotations.loaders import gencode_genes

    backend = request.param
    base = tmp_path_factory.mktemp(f"drift-{backend}")
    home, ann = base / "home", base / "ann.duckdb"
    home.mkdir()
    env = {"VCFCLICK_ANNOTATIONS_DB": str(ann)}
    with mock.patch.dict(os.environ, env):
        gencode_genes.load(FIXTURES / "scope_genes.gff3")
        run_cli(home, backend, "db", "create", "t")
        run_cli(
            home, backend, "db", "ingest", "t", str(FIXTURES / "scope.vcf.gz"),
            "--ingest-id", "g1", "--serial", "--keep-reference",
        )  # fmt: skip
        run_cli(home, backend, "db", "ped", "t", str(FIXTURES / "scope.ped"))
        args = ("db", "trio", "t", "--proband", "CHILD", "--format", "json")
        scoped = [a for r in REGIONS for a in ("--region", r)]
        scoped += ["--exclude-region", "chr7:29999-29999"]
        plain = json.loads(run_cli(home, backend, *args).stdout)
        with_scope = json.loads(run_cli(home, backend, *args, *scoped).stdout)
    from tests.test_scope_db import _FORCED
    from tests.conftest import _run_python

    forced = {}
    for mode in ("python", "sql_buckets"):
        out = _run_python(
            home, _FORCED, mode, *args, *scoped, VCFCLICK_BACKEND=backend, **env
        )
        forced[mode] = json.loads(out[out.index("{") :])
    return plain, with_scope, forced


def _variants(model_doc) -> list[dict]:
    """Every variant dict (chrom, pos) a model object holds, at any depth."""
    out, stack = [], [model_doc]
    while stack:
        v = stack.pop()
        if isinstance(v, dict):
            if "chrom" in v and "pos" in v:
                out.append(v)
            else:
                stack.extend(v.values())
        elif isinstance(v, list):
            stack.extend(v)
    return out


@pytest.mark.parametrize("model", list(MODELS))
def test_planted_out_of_scope_rows_never_come_back(runs, model):
    plain, scoped, forced = runs
    scope = planted_scope()
    before = _variants(plain["models"][model])
    # not vacuous: the model has rows on both sides of the planted scope
    assert any(not scope.contains(v["chrom"], v["pos"]) for v in before), model
    assert any(scope.contains(v["chrom"], v["pos"]) for v in before), model
    for doc in (scoped, *forced.values()):
        after = _variants(doc["models"][model])
        assert after, model
        for v in after:
            inside = scope.contains(v["chrom"], v["pos"])
            assert inside or v.get("in_scope") is False, (model, v)
    assert scoped["scope"]["applied"] == "sql"
    for mode, doc in forced.items():
        assert doc["models"][model] == scoped["models"][model], mode
        assert doc["scope"]["applied"] == mode
