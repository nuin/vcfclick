"""The decorator chains of cli/trio_output.py.

A trio candidate passes through one ordered chain. Row decorators need only
the variant key and also serve other commands' rows (`decorate_row`);
candidate decorators need the 8-column trio row and never reach those rows.
Document decorators edit the finished trio document. A decorator may also
register a prefetch, run once per batch of rows. Pure unit tests: the only
store is the in-process DuckDB annotation store, never a chDB session.
"""

from __future__ import annotations

from types import SimpleNamespace

import duckdb
import pytest

import cli.trio_output as out

NO_ANNOTATIONS = {"genes": False, "clinvar": False, "gnomad": False}
TRIO_ROW = ["chr1", 100, "A", "G", 1, 0, 0, 0.01]  # the 8-column detail row
VARIANT_ROW = ("chr1", 100, "A", "G")  # all another command has to offer
ANNOTATION_KEYS = ["genes", "gnomad_popmax", "clinvar"]


@pytest.fixture
def chains(monkeypatch):
    """Private copies of the registries, so a test can register freely."""
    monkeypatch.setattr(out, "_DECORATORS", list(out._DECORATORS))
    monkeypatch.setattr(out, "_DOCUMENT_DECORATORS", list(out._DOCUMENT_DECORATORS))


def _render(**kwargs) -> out.Render:
    return out.Render(object(), 10, NO_ANNOTATIONS, **kwargs)


# ───────────────────────────── row decorators ─────────────────────────────


def test_shipped_chain_gives_a_candidate_and_a_foreign_row_the_same_columns():
    render = _render()
    cand = out._candidate(TRIO_ROW, "denovo", render)
    assert list(cand) == [
        "chrom", "pos", "ref", "alt",
        "proband_gt", "father_gt", "mother_gt", "af",
        *ANNOTATION_KEYS,
    ]  # fmt: skip
    row = out.decorate_row({"chrom": "chr1", "pos": 100}, VARIANT_ROW, render)
    assert list(row) == ["chrom", "pos", *ANNOTATION_KEYS]
    assert {k: cand[k] for k in ANNOTATION_KEYS} == {k: row[k] for k in ANNOTATION_KEYS}


def test_candidate_decorators_never_reach_a_foreign_row(chains):
    """A decorator registered for trio candidates with models=None applies to
    every trio model and to no row from another command."""

    @out.register_row_decorator
    def row_tag(d, key, render):
        d["row_key"] = list(key)
        return d

    @out.register_candidate_decorator
    def any_model(d, row, model, render):
        d["any_model"] = [model, len(row)]
        return d

    def only_denovo(d, row, model, render):
        d["only_denovo"] = True
        return d

    out.register_candidate_decorator(only_denovo, models=frozenset({"denovo"}))

    render = _render()
    denovo = out._candidate(TRIO_ROW, "denovo", render)
    assert list(denovo)[-3:] == ["row_key", "any_model", "only_denovo"]
    assert denovo["any_model"] == ["denovo", 8]
    assert "only_denovo" not in out._candidate(TRIO_ROW, "recessive", render)

    row = out.decorate_row({}, VARIANT_ROW, render)
    assert list(row) == [*ANNOTATION_KEYS, "row_key"]
    assert row["row_key"] == ["chr1", 100, "A", "G"]


def test_registration_order_is_key_order_across_both_kinds(chains):
    def add(name):
        def fn(d, *args):
            d[name] = True
            return d

        return fn

    out.register_row_decorator(add("first"))
    out.register_candidate_decorator(add("second"))
    out.register_row_decorator(add("third"))
    cand = out._candidate(TRIO_ROW, "denovo", _render())
    assert list(cand)[-3:] == ["first", "second", "third"]
    row = out.decorate_row({}, VARIANT_ROW, _render())
    assert list(row)[-2:] == ["first", "third"]


# ───────────────────── the per-batch prefetch hook ─────────────────────


def _lookup_decorators(calls):
    """A row decorator and a candidate decorator that do their lookup in a
    prefetch, as a later feature's batched query would. The per-row call only
    reads `render.cache`, so a row that was not prefetched is a KeyError."""

    def row_prefetch(render, keys):
        calls["row"].append(list(keys))
        for key in keys:
            render.cache[("depth", key)] = key[1] * 2

    def depth(d, key, render):
        d["depth"] = render.cache[("depth", key)]
        return d

    def cand_prefetch(render, rows):
        calls["cand"].append([r[1] for r in rows])
        for r in rows:
            render.cache[("pgt", tuple(r))] = r[4]

    def pgt(d, row, model, render):
        d["pgt"] = render.cache[("pgt", tuple(row))]
        return d

    out.register_row_decorator(depth, prefetch=row_prefetch)
    out.register_candidate_decorator(
        pgt, models=frozenset({"denovo"}), prefetch=cand_prefetch
    )


def test_prefetch_runs_once_per_rendered_model_and_rows_look_nothing_up(chains):
    calls = {"row": [], "cand": []}
    _lookup_decorators(calls)
    rows = [
        ["chr1", 100, "A", "G", 1, 0, 0, 0.01],
        ["chr1", 200, "C", "T", 1, 0, 0, None],
        ["chr1", 100, "A", "G", 2, 0, 0, 0.01],  # the same variant again
    ]
    render = _render()

    denovo = out._candidates(rows, "denovo", render)
    assert calls == {"row": [[("chr1", 100, "A", "G"), ("chr1", 200, "C", "T")]],
                     "cand": [[100, 200, 100]]}  # fmt: skip
    assert [(c["depth"], c["pgt"]) for c in denovo] == [(200, 1), (400, 1), (200, 2)]

    # another model: the row decorator's prefetch runs again for the new
    # batch, the candidate decorator's does not (it is for denovo only)
    recessive = out._candidates(rows, "recessive", render)
    assert len(calls["row"]) == 2 and len(calls["cand"]) == 1
    assert all("depth" in c and "pgt" not in c for c in recessive)

    # nothing to render, nothing to look up
    assert out._candidates([], "denovo", render) == []
    assert len(calls["row"]) == 2 and len(calls["cand"]) == 1


def test_decorate_rows_prefetches_once_for_a_foreign_command(chains):
    calls = {"row": [], "cand": []}
    _lookup_decorators(calls)
    rows = [("chr1", 100, "A", "G"), ("chr2", 5, "T", "C"), ("chr1", 100, "A", "G")]
    got = out.decorate_rows([{"i": i} for i in range(3)], rows, _render())
    assert calls == {"row": [[("chr1", 100, "A", "G"), ("chr2", 5, "T", "C")]],
                     "cand": []}  # fmt: skip
    assert [(d["i"], d["depth"]) for d in got] == [(0, 200), (1, 10), (2, 200)]
    assert "pgt" not in got[0]  # a candidate decorator never reaches these

    assert out.decorate_rows([], [], _render()) == []
    assert len(calls["row"]) == 1  # an empty batch makes no call
    with pytest.raises(ValueError, match="2 dicts for 3 rows"):
        out.decorate_rows([{}, {}], rows, _render())


# ─────────────────── annotation prefetch with row decorators ───────────────────


def test_decorate_rows_batches_the_annotation_lookups(
    isolated_annotation_db, monkeypatch
):
    from annotations.db import get_connection

    conn = get_connection()
    conn.execute(
        "INSERT INTO refseq_genes VALUES ('GA', 'chr1', 50, 500, '+', 'id', 'x')"
    )
    conn.execute("INSERT INTO gnomad_af VALUES ('chr1', 100, 'A', 'G', 0.2, 0.3)")
    conn.execute(
        "INSERT INTO clinvar_variants VALUES "
        "('chr1', 100, 'A', 'G', 'Pathogenic', 'x', 'VCV1', 'c')"
    )
    conn.close()
    rows = [VARIANT_ROW, ("chr1", 700, "C", "T")]  # the second matches nothing

    assert out._prefetch is out.prefetch_annotations  # the old name still works
    have, store = out._annotation_store()
    try:
        assert have == {"genes": True, "clinvar": True, "gnomad": True}
        per_row = [
            out.decorate_row({}, r, out.Render(object(), 1, have, store)) for r in rows
        ]

        # the annotation columns are a row decorator with its own prefetch:
        # decorate_rows alone batches them, and no row asks the store
        batched = out.Render(object(), 1, have, store)
        with monkeypatch.context() as m:
            m.setattr(out, "_enrich", pytest.fail)
            got = out.decorate_rows([{}, {}], rows, batched)
        assert {k for k in batched.cache if isinstance(k, tuple)} == set(rows)
    finally:
        store.close()
    assert got == per_row
    assert got[0] == {"genes": ["GA"], "gnomad_popmax": 0.3, "clinvar": "Pathogenic"}
    assert got[1] == {"genes": [], "gnomad_popmax": None, "clinvar": None}


# ─────────────────────────── document decorators ───────────────────────────


def _analysis() -> SimpleNamespace:
    """The few attributes `_trio_json` reads when there is no model to build."""
    return SimpleNamespace(
        trio=SimpleNamespace(
            ingest_id="i1",
            proband="C",
            father="F",
            mother="M",
            family_id="fam",
            kind="trio",
        ),  # fmt: skip
        gates=SimpleNamespace(
            min_gq=20,
            min_dp=10,
            max_af=0.01,
            min_ab=0.25,
            max_ab=0.75,
            af_gate="applied",
        ),  # fmt: skip
        has_reference=True,
        siblings=[],
        gnomad_max_af=None,
    )


def test_document_decorators_edit_the_whole_document_in_order(
    chains, isolated_annotation_db, monkeypatch
):
    monkeypatch.setattr(out, "MODELS", {})
    seen = []

    @out.register_document_decorator
    def first(doc, render):
        seen.append(("first", list(doc)))
        doc["trio"]["parentage"] = {"verdict": "ok"}  # a nested key
        doc["preset"] = None  # a top-level key
        return doc

    @out.register_document_decorator
    def second(doc, render):
        seen.append(("second", "preset" in doc))
        # the annotation store is still open for a decorator that looks up
        assert render.store.execute("SELECT count(*) FROM gnomad_af").fetchone() == (0,)
        return {**doc, "replaced": True}  # returning a new dict is allowed

    doc = out._trio_json(_analysis(), "all", 10)
    assert list(doc) == [
        "trio", "keep_reference", "gates", "annotations", "models",
        "preset", "replaced",
    ]  # fmt: skip
    assert doc["trio"]["parentage"] == {"verdict": "ok"}
    assert seen == [
        ("first", ["trio", "keep_reference", "gates", "annotations", "models"]),
        ("second", True),
    ]


def test_a_document_decorator_must_return_the_document(
    chains, isolated_annotation_db, monkeypatch
):
    monkeypatch.setattr(out, "MODELS", {})
    opened = {}

    @out.register_document_decorator
    def forgot_to_return(doc, render):
        opened["store"] = render.store
        doc["x"] = 1

    with pytest.raises(TypeError, match="forgot_to_return"):
        out._trio_json(_analysis(), "all", 10)
    with pytest.raises(duckdb.Error):  # closed on the way out
        opened["store"].execute("SELECT 1")
