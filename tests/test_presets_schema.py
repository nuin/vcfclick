"""What a preset file may contain: parsing, validation, versioning, fingerprint.

Pure and in-process: no database, no file system. Every error code of
`cli/preset_schema.py` is reached here, and a file is data: nothing in it is
executed, and no message quotes more of it than a key name.
"""

from __future__ import annotations

import json

import pytest

from cli import preset_schema as ps


def envelope(filters=None, **extra) -> dict:
    return {
        "schema": ps.SCHEMA,
        "schema_version": 1,
        "filters": {"min_gq": 20} if filters is None else filters,
        **extra,
    }


def codes(doc, **kw) -> list[str]:
    return [p.code for p in ps.validate_document(doc, **kw).problems]


# ───────────────────────────────── names ─────────────────────────────────


@pytest.mark.parametrize("name", ["id-panel", "a", "rare.v2", "x_1", "0abc", "a" * 63])
def test_good_names(name):
    assert ps.name_problem(name) is None


@pytest.mark.parametrize(
    "name",
    [
        "",
        "Id-Panel",  # lowercase only
        "-lead",
        ".hidden",
        "a" * 64,
        "a/b",
        "../x",
        "/etc/passwd",
        "x.json",  # would be read as a path on the command line
        "con",  # Windows device names
        "nul",
        "com1",
        "lpt9",
        "con.txt",
        None,
        7,
    ],
)
def test_bad_names(name):
    assert ps.name_problem(name)


# ───────────────────────────────── parsing ─────────────────────────────────


@pytest.mark.parametrize(
    "text",
    ['{"a": NaN}', '{"a": Infinity}', '{"a": -Infinity}', '{"a": 1, "a": 2}'],
)
def test_parse_refuses_nan_infinity_and_duplicate_keys(text):
    with pytest.raises(ps.PresetError):
        ps.parse_json(text)


def test_parse_errors_name_the_position_not_the_content():
    secret = "root:x:0:0:SECRET-CONTENT"
    with pytest.raises(ps.PresetError) as e:
        ps.parse_json(secret)
    assert "SECRET" not in str(e.value)
    assert "line 1" in str(e.value)
    with pytest.raises(ps.PresetError, match="UTF-8"):
        ps.parse_json(b"\xff\xfe\x00")


def test_parse_has_a_size_cap():
    with pytest.raises(ps.PresetError, match="larger"):
        ps.parse_json(b" " * (ps.MAX_BYTES + 1))


def test_1e999_parses_to_infinity_and_validation_refuses_it():
    doc = ps.parse_json(
        '{"schema": "vcfclick.preset", "schema_version": 1, '
        '"filters": {"max_af": 1e999}}'
    )
    assert "bad_value" in codes(doc)


# ───────────────────────────────── validation ─────────────────────────────────


def test_a_minimal_preset_is_valid_and_normalised():
    checked = ps.validate_document(envelope({"max_af": 1, "category": "recessive"}))
    assert checked.ok and checked.problems == ()
    # an int in a float key is a float, a category string a list: one spelling
    # per meaning, so the fingerprint does not depend on how it was typed
    assert checked.doc["filters"] == {"category": ["recessive"], "max_af": 1.0}
    assert isinstance(checked.doc["filters"]["max_af"], float)
    assert list(checked.doc) == ["schema", "schema_version", "filters", "meta"]


@pytest.mark.parametrize(
    "doc,code",
    [
        ([], "not_a_preset"),
        ({"schema": "something-else"}, "not_a_preset"),
        ({"schema": ps.SCHEMA}, "unsupported_schema"),  # a v0-shaped file
        ({"schema": ps.SCHEMA, "schema_version": 0}, "unsupported_schema"),
        ({"schema": ps.SCHEMA, "schema_version": True}, "unsupported_schema"),
        (envelope() | {"schema_version": 2}, "newer_schema"),
        (envelope(name="Bad Name"), "bad_name"),
        (envelope(name="x.json"), "bad_name"),
        (envelope(title=3), "bad_type"),
        (envelope(applies_to="trio"), "bad_type"),
        (envelope(assay="proteome"), "bad_value"),
        (envelope(filters=[]), "bad_type"),
        (envelope({"min_gq": None}), "bad_value"),
        (envelope({"min_gq": True}), "bad_value"),
        (envelope({"min_gq": 100}), "bad_value"),
        (envelope({"min_gq": "20"}), "bad_value"),
        (envelope({"max_af": float("nan")}), "bad_value"),
        (envelope({"max_af": -0.1}), "bad_value"),
        (envelope({"strict_missing": 1}), "bad_value"),
        (envelope({"category": ["nope"]}), "bad_value"),
        (envelope({"proband": "C1"}), "never_presettable"),
        (envelope({"preset": "other"}), "never_presettable"),  # no nesting
        (envelope({"min_gp": 1}), "unknown_key"),
        (envelope(meta=[]), "bad_type"),
        (envelope(meta={"locked": "yes"}), "bad_type"),
        (envelope(meta={"approved_by": 3}), "bad_type"),
        (envelope({"min_ab": 0.8, "max_ab": 0.2}), "inconsistent"),
    ],
)
def test_every_error_code(doc, code):
    checked = ps.validate_document(doc)
    assert not checked.ok
    assert code in [p.code for p in checked.errors]


def test_a_newer_schema_says_to_upgrade():
    (problem,) = ps.validate_document(envelope() | {"schema_version": 9}).problems
    assert "newer vcfclick (schema 9)" in problem.message
    assert "upgrade" in problem.message


def test_never_presettable_says_why():
    (problem,) = ps.validate_document(envelope({"attest_affected": ["C1"]})).problems
    assert problem.code == "never_presettable"
    assert "attestation" in problem.message


def test_an_unknown_key_suggests_and_names_the_vcfclick_that_wrote_it():
    doc = envelope({"max_afx": 0.1}, meta={"created_with": "0.99.0"})
    (problem,) = ps.validate_document(doc).problems
    assert problem.code == "unknown_key" and problem.level == ps.ERROR
    assert "did you mean max_af" in problem.message
    assert "written by vcfclick 0.99.0; upgrade" in problem.message


def test_allow_unknown_downgrades_to_a_warning_and_lists_the_key():
    checked = ps.validate_document(envelope({"future_key": 1}), allow_unknown=True)
    assert checked.ok
    assert checked.unknown_keys == ("future_key",)
    assert [p.level for p in checked.problems] == [ps.WARN]
    assert checked.doc["filters"]["future_key"] == 1  # kept, not applied


def test_extension_keys_are_kept_and_ignored():
    checked = ps.validate_document(envelope({"x_lab": {"note": 1}, "min_gq": 20}))
    assert checked.ok and checked.problems == ()
    assert checked.doc["filters"]["x_lab"] == {"note": 1}
    assert ps.fingerprint(checked.doc) == ps.fingerprint(
        ps.validate_document(envelope({"min_gq": 20})).doc
    )


def test_unknown_envelope_keys_survive_a_load_and_a_save():
    doc = envelope(lab_extension={"a": [1, 2]}, title="t")
    checked = ps.validate_document(doc)
    assert checked.doc["lab_extension"] == {"a": [1, 2]}
    again = ps.validate_document(json.loads(ps.to_text(checked.doc)))
    assert again.doc == checked.doc


def test_an_unknown_assay_command_is_a_warning_not_an_error():
    checked = ps.validate_document(envelope(applies_to=["trio", "nope"]))
    assert checked.ok and [p.code for p in checked.problems] == ["unknown_command"]


# ───────────────────────────── versioning ─────────────────────────────


def test_migrate_is_the_identity_at_the_current_version():
    doc = envelope()
    assert ps.migrate(doc, ps.SCHEMA_VERSION) is doc


def test_a_migration_chain_runs_in_order(monkeypatch):
    monkeypatch.setattr(ps, "SCHEMA_VERSION", 3)
    monkeypatch.setattr(
        ps,
        "_MIGRATIONS",
        {1: lambda d: {**d, "via": [*d.get("via", []), 1]},
         2: lambda d: {**d, "via": [*d.get("via", []), 2]}},
    )  # fmt: skip
    assert ps.migrate({"a": 1}, 1)["via"] == [1, 2]
    assert ps.migrate({"a": 1}, 2)["via"] == [2]


# ───────────────────────────── fingerprint ─────────────────────────────


def fp(filters, **extra) -> str:
    checked = ps.validate_document(envelope(filters, **extra))
    assert checked.ok, checked.problems
    return ps.fingerprint(checked.doc)


def test_fingerprint_is_stable_under_order_and_description_and_meta():
    base = fp({"min_gq": 20, "category": ["recessive", "comphet"]})
    assert len(base) == 12 and set(base) <= set("0123456789abcdef")
    assert fp({"category": ["comphet", "recessive"], "min_gq": 20}) == base
    assert fp({"min_gq": 20, "category": "recessive,comphet"}) == base
    assert fp({"min_gq": 20, "category": ["recessive", "comphet"]}, title="x") == base
    assert (
        fp(
            {"min_gq": 20, "category": ["recessive", "comphet"]},
            meta={"updated": "2030-01-01T00:00:00Z", "locked": True},
        )
        == base
    )


def test_fingerprint_changes_with_a_value_and_float_spelling_does_not():
    assert fp({"min_gq": 20}) != fp({"min_gq": 21})
    assert fp({"max_af": 0.01}) == fp({"max_af": 1e-2})
    assert fp({"max_af": 1}) == fp({"max_af": 1.0})
    assert fp({"strict_missing": False}) != fp({})  # explicit is not absent


def test_text_round_trips():
    checked = ps.validate_document(envelope({"max_af": 0.0001, "pass_only": True}))
    again = ps.validate_document(ps.parse_json(ps.to_text(checked.doc)))
    assert again.doc == checked.doc
    assert ps.to_text(checked.doc).endswith("}\n")
