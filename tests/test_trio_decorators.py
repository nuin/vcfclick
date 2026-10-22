"""The decorator chains of cli/trio_output.py.

A trio candidate passes through one ordered chain. Row decorators need only
the variant key and also serve other commands' rows (`decorate_row`);
candidate decorators need the 8-column trio row and never reach those rows.
Document decorators edit the finished trio document. Pure unit tests: the
only store is the in-process DuckDB annotation store, never a chDB session.
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


# ─────────────────── annotation prefetch with row decorators ───────────────────


def test_prefetch_serves_the_row_decorator_for_a_foreign_row(
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

        batched = out.Render(object(), 1, have, store)
        out.prefetch_annotations(batched, rows)
        assert {k for k in batched.cache if isinstance(k, tuple)} == set(rows)
        with monkeypatch.context() as m:  # a prefetched row never asks the store
            m.setattr(out, "_enrich", pytest.fail)
            got = [out.decorate_row({}, r, batched) for r in rows]
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
