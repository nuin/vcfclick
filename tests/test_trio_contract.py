"""The `db trio --format json` contract the desktop apps decode.

Swift's Codable and System.Text.Json fail differently on a missing,
retyped or null field, and a failed decode blanks the whole Trio panel. So
the contract is checked structurally, on real output from both backends:

  * every model object has an int `count` and a bool `blocked`, plus the
    additive `blocked_reason` (null, or a member of the closed enum);
  * every candidate carries a `genes` list of strings and int genotypes;
  * `annotations` stays {string: bool};
  * comphet gene entries always have `gene`; `paternal` / `maternal` are
    OPTIONAL (a model that cannot tell parent of origin omits them), and
    `tier` / `k_carried` / `k_absent` may appear on them.

`assert_trio_document` is importable so each later model's tests run the
same checks on their own output.
"""

from __future__ import annotations

import json

import pytest

from cli.trio_analysis import BLOCKED_REASONS, CATEGORIES
from tests.conftest import FIXTURES, run_cli

BACKENDS = ("chdb", "duckdb")
TOP_LEVEL = ["trio", "keep_reference", "gates", "annotations", "models"]
GATE_KEYS = ["min_gq", "min_dp", "max_af", "min_ab", "max_ab", "gnomad_max_af"]


def _is_int(v) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _is_num(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _check_variant(v: dict) -> None:
    assert isinstance(v["chrom"], str) and _is_int(v["pos"])
    assert isinstance(v["ref"], str) and isinstance(v["alt"], str)
    assert v["af"] is None or _is_num(v["af"])


def assert_trio_document(doc: dict) -> None:
    """Structural checks every trio document must pass, old keys and new."""
    assert list(doc)[: len(TOP_LEVEL)] == TOP_LEVEL
    trio = doc["trio"]
    assert {"ingest_id", "proband", "father", "mother"} <= set(trio)
    assert isinstance(doc["keep_reference"], bool)
    gates = doc["gates"]
    assert list(gates)[: len(GATE_KEYS)] == GATE_KEYS
    assert _is_int(gates["min_gq"]) and _is_int(gates["min_dp"])
    assert all(_is_num(gates[k]) for k in ("max_af", "min_ab", "max_ab"))
    assert gates["gnomad_max_af"] is None or _is_num(gates["gnomad_max_af"])
    # annotations stays {string: bool}: an app decodes it as a dictionary
    # of booleans and throws on any other value.
    assert doc["annotations"], "annotations is never empty"
    assert all(
        isinstance(k, str) and isinstance(v, bool)
        for k, v in doc["annotations"].items()
    )

    assert doc["models"], "a trio document always has at least one model"
    for name, model in doc["models"].items():
        assert _is_int(model["count"]), name
        assert isinstance(model["blocked"], bool), name
        assert isinstance(model["truncated"], bool), name
        # additive, and null exactly when the model is not blocked/unmet
        assert "blocked_reason" in model, name
        assert model["blocked_reason"] is None or (
            model["blocked_reason"] in BLOCKED_REASONS
        ), (name, model["blocked_reason"])
        if model["blocked"]:
            assert model["blocked_reason"] is not None, name
        # blocked_reason sits right after blocked: a reader scanning the
        # object sees them together.
        keys = list(model)
        assert keys[keys.index("blocked") + 1] == "blocked_reason", keys

        if "candidates" in model:
            assert _is_int(len(model["candidates"]))
            for c in model["candidates"]:
                _check_variant(c)
                assert all(
                    _is_int(c[k]) for k in ("proband_gt", "father_gt", "mother_gt")
                )
                assert isinstance(c["genes"], list)
                assert all(isinstance(g, str) for g in c["genes"])
                assert c["gnomad_popmax"] is None or _is_num(c["gnomad_popmax"])
                assert c["clinvar"] is None or isinstance(c["clinvar"], str)
        if "genes" in model:
            for g in model["genes"]:
                assert isinstance(g["gene"], str)
                for side in ("paternal", "maternal"):
                    if side in g:  # optional for every consumer
                        assert isinstance(g[side], list)
                        for v in g[side]:
                            _check_variant(v)
                for opt, ok in (
                    ("tier", lambda v: isinstance(v, (str, int))),
                    ("k_carried", lambda v: isinstance(v, (int, list))),
                    ("k_absent", lambda v: isinstance(v, (int, list))),
                ):
                    if opt in g:
                        assert ok(g[opt]), (opt, g[opt])


def test_contract_checker_rejects_what_would_blank_an_app():
    """The checker is not vacuous: each of these is a decode failure in at
    least one app."""
    good = {
        "trio": {"ingest_id": "i", "proband": "C", "father": "F", "mother": "M"},
        "keep_reference": True,
        "gates": dict.fromkeys(GATE_KEYS, 1) | {"gnomad_max_af": None},
        "annotations": {"genes": False},
        "models": {
            "recessive": {
                "count": 0,
                "blocked": False,
                "blocked_reason": None,
                "truncated": False,
                "candidates": [],
            }
        },
    }
    assert_trio_document(good)

    def broken(mutate):
        doc = json.loads(json.dumps(good))
        mutate(doc)
        with pytest.raises((AssertionError, KeyError)):
            assert_trio_document(doc)

    broken(lambda d: d["models"]["recessive"].update(count=True))
    broken(lambda d: d["models"]["recessive"].update(blocked="no"))
    broken(lambda d: d["models"]["recessive"].pop("blocked_reason"))
    broken(lambda d: d["models"]["recessive"].update(blocked_reason="because"))
    broken(lambda d: d["annotations"].update(genes="yes"))
    broken(lambda d: d["models"]["recessive"].update(blocked=True))


@pytest.fixture
def annotation_dbs(tmp_path, monkeypatch):
    """An annotation store with GENEX over chr1:1-100000, and an empty one."""
    from annotations.db import get_connection

    def make(name: str, genes: bool) -> str:
        path = tmp_path / name
        with monkeypatch.context() as m:
            m.setenv("VCFCLICK_ANNOTATIONS_DB", str(path))
            conn = get_connection()
            if genes:
                conn.execute(
                    "INSERT INTO refseq_genes VALUES "
                    "('GENEX', 'chr1', 1, 100000, '+', 'id', 'protein_coding')"
                )
            conn.close()
        return str(path)

    return make("with_genes.duckdb", True), make("no_genes.duckdb", False)


def _trio(home, backend, monkeypatch, ann, *args) -> dict:
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", ann)
    r = run_cli(
        home,
        backend,
        "db",
        "trio",
        "t",
        "--proband",
        "CHILD",
        "--format",
        "json",
        *args,
    )
    return json.loads(r.stdout)


def _store(home, backend, *, keep_reference: bool) -> None:
    run_cli(home, backend, "db", "create", "t")
    args = [
        "db", "ingest", "t", str(FIXTURES / "trio.vcf.gz"),
        "--cohort", "fam", "--ingest-id", "fam1", "--serial",
    ]  # fmt: skip
    if keep_reference:
        args.append("--keep-reference")
    run_cli(home, backend, *args)
    run_cli(home, backend, "db", "ped", "t", str(FIXTURES / "trio.ped"))


@pytest.mark.parametrize("backend", BACKENDS)
def test_trio_json_contract_and_blocked_reason(
    backend, vcfclick_home, monkeypatch, annotation_dbs
):
    with_genes, no_genes = annotation_dbs
    _store(vcfclick_home, backend, keep_reference=True)

    # Genes loaded and hom-ref rows stored: nothing is blocked.
    doc = _trio(vcfclick_home, backend, monkeypatch, with_genes, "--max-af", "0.3")
    assert_trio_document(doc)
    assert list(doc["models"]) == list(CATEGORIES)
    assert {m: v["blocked_reason"] for m, v in doc["models"].items()} == dict.fromkeys(
        CATEGORIES
    )
    assert doc["models"]["comphet"]["genes"], "GENEX has a paternal and a maternal het"
    assert doc["models"]["denovo"]["candidates"]

    # No gene annotations: comphet keeps blocked=False (its original
    # meaning is "no stored hom-ref rows") and says why it is empty.
    doc = _trio(vcfclick_home, backend, monkeypatch, no_genes, "--max-af", "0.3")
    assert_trio_document(doc)
    comphet = doc["models"]["comphet"]
    assert comphet["blocked"] is False and comphet["needs_genes"] is True
    assert comphet["blocked_reason"] == "needs_genes"
    assert doc["models"]["denovo"]["blocked_reason"] is None

    # A single category keeps the same per-model shape.
    doc = _trio(
        vcfclick_home, backend, monkeypatch, no_genes, "--category", "recessive"
    )
    assert_trio_document(doc)
    assert list(doc["models"]) == ["recessive"]


@pytest.mark.parametrize("backend", BACKENDS)
def test_sparse_ingest_blocks_with_keep_reference_reason(
    backend, vcfclick_home, monkeypatch, annotation_dbs
):
    _, no_genes = annotation_dbs
    _store(vcfclick_home, backend, keep_reference=False)
    doc = _trio(vcfclick_home, backend, monkeypatch, no_genes)
    assert_trio_document(doc)
    reasons = {m: v["blocked_reason"] for m, v in doc["models"].items()}
    assert reasons == {
        "denovo": "needs_keep_reference",
        "recessive": None,
        "dominant": "needs_keep_reference",
        # blocked wins over needs_genes: fix the first requirement first
        "comphet": "needs_keep_reference",
    }
    assert [m for m, v in doc["models"].items() if v["blocked"]] == [
        "denovo",
        "dominant",
        "comphet",
    ]


def test_blocked_reason_enum_is_closed_and_complete():
    """The vocabulary the plan fixes for every model of the program."""
    assert BLOCKED_REASONS == (
        "needs_keep_reference",
        "needs_genes",
        "needs_phase",
        "needs_parents",
        "needs_sex",
        "needs_ad",
        "needs_parent_ad",
        "no_autosomes",
        "no_x_data",
        "parent_sex_conflict",
        "build_conflict",
    )
    from cli.trio_analysis import check_blocked_reason

    assert check_blocked_reason(None) is None
    assert check_blocked_reason("needs_genes") == "needs_genes"
    with pytest.raises(ValueError):
        check_blocked_reason("needs_reference")  # the pre-rename spelling


def test_comphet_entry_omits_a_side_that_is_absent():
    """paternal/maternal are optional in the contract: the builder emits only
    the sides it is given."""
    from cli.trio_output import _comphet_entry

    one = _comphet_entry("G", {"paternal": [("chr1", 5, "A", "T", 0.01)]})
    assert list(one) == ["gene", "paternal"]
    both = _comphet_entry(
        "G",
        {
            "paternal": [("chr1", 5, "A", "T", 0.01)],
            "maternal": [("chr1", 9, "C", "G", None)],
        },
    )
    assert list(both) == ["gene", "paternal", "maternal"]
    assert both["maternal"][0]["af"] is None
