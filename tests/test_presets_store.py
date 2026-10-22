"""The preset store and the `vcfclick preset` commands, against a temp home.

In-process (Click's CliRunner): no command here opens a database. The home is
the test's own `VCFCLICK_HOME`, never the real one.
"""

from __future__ import annotations

import json
import os
import shutil

import pytest
from click.testing import CliRunner

from cli import preset_exchange as exchange
from cli import preset_schema as ps
from cli import presets as store
from cli.main import cli


@pytest.fixture
def run(vcfclick_home, monkeypatch):
    monkeypatch.delenv("VCFCLICK_PRESET_PATH", raising=False)

    def invoke(*args, input=None, ok=True):
        r = CliRunner().invoke(cli, ["preset", *args], input=input)
        if ok:
            assert r.exit_code == 0, (args, r.output)
        else:
            assert r.exit_code != 0, (args, r.output)
        return r

    return invoke


def files(home):
    d = home / "presets"
    return sorted(p.name for p in d.iterdir()) if d.exists() else []


def saved(home, name) -> dict:
    return json.loads((home / "presets" / f"{name}.json").read_text())


# ─────────────────────────────── save and use ───────────────────────────────


def test_save_list_show_round_trip(run, vcfclick_home):
    r = run(
        "save", "id-panel", "--from-preset", "rare-recessive", "--max-af", "0.002",
        "--description", "ID recessive",
    )  # fmt: skip
    assert "saved preset 'id-panel'" in r.output
    doc = saved(vcfclick_home, "id-panel")
    assert doc["name"] == "id-panel" and doc["description"] == "ID recessive"
    assert doc["filters"]["max_af"] == 0.002
    assert doc["filters"]["category"] == ["recessive", "comphet"]  # copied
    assert doc["meta"]["created_with"] and doc["meta"]["updated"]

    rows = json.loads(run("list", "--format", "json").stdout)["presets"]
    mine = next(p for p in rows if p["name"] == "id-panel")
    assert (mine["source"], mine["available"], mine["problems"]) == ("user", True, [])
    assert mine["categories"] == ["recessive", "comphet"]
    assert mine["revision"] is None and mine["status"] is None
    assert [p["name"] for p in rows[:4]] == [
        "broad-first-pass", "strict-denovo", "rare-recessive", "rare-dominant"
    ]  # fmt: skip
    assert {"fingerprint", "applies_to", "assay", "keys", "locked"} <= set(mine)

    shown = json.loads(run("show", "id-panel", "--resolve", "--format", "json").stdout)
    assert shown["filters"]["max_af"] == 0.002
    assert shown["flags"]["max_af"] == "--max-af"
    assert shown["fingerprint"] == mine["fingerprint"]
    text = run("show", "strict-denovo").stdout
    assert "starting point, not validated" in text
    assert "covers:" in text and "does not cover:" in text and "history:" in text


def test_save_takes_only_the_flags_typed(run, vcfclick_home):
    run("save", "tiny", "--min-gq", "25", "--strict-missing")
    assert saved(vcfclick_home, "tiny")["filters"] == {
        "min_gq": 25, "strict_missing": True
    }  # nothing inherited from the defaults  # fmt: skip
    run("save", "cats", "--category", "recessive,comphet")
    assert saved(vcfclick_home, "cats")["filters"]["category"] == [
        "recessive", "comphet"
    ]  # fmt: skip


@pytest.mark.parametrize(
    "args,message",
    [
        (["con"], "device"),
        (["x.json"], ".json"),
        (["Upper"], "lowercase"),
        (["broad-first-pass"], "built-in"),
        (["ok", "--min-gq", "500"], "min_gq"),
        (["ok", "--min-ab", "0.9", "--max-ab", "0.1"], "min_ab"),
        (["ok", "--from-preset", "nope"], "nope"),
        (["ok", "--from-preset", "rare-dominant", "--from-json", "-"], "choose one"),
    ],
)
def test_save_refuses(run, args, message, vcfclick_home):
    r = run("save", *args, ok=False, input="{}")
    assert message in r.output
    assert files(vcfclick_home) == []  # nothing half-written


def test_save_refuses_an_existing_name_without_force(run):
    run("save", "p", "--min-gq", "21")
    assert "already exists" in run("save", "p", "--min-gq", "22", ok=False).output
    run("save", "p", "--min-gq", "22", "--force")


def test_save_from_json_accepts_filters_only_and_keeps_future_envelope_keys(
    run, vcfclick_home
):
    body = {"filters": {"min_dp": 12}, "description": "d", "lab_note": {"a": 1}}
    run("save", "fromjson", "--from-json", "-", input=json.dumps(body))
    doc = saved(vcfclick_home, "fromjson")
    assert doc["filters"] == {"min_dp": 12} and doc["description"] == "d"
    assert doc["lab_note"] == {"a": 1}
    bad = run("save", "x", "--from-json", "-", input="[1]", ok=False)
    assert "filters" in bad.output


# ───────────────────────── set, unset, edit, rename, rm ─────────────────────────


def test_set_unset_rename_rm(run, vcfclick_home):
    run("save", "p", "--min-gq", "21")
    run("set", "p", "max_af=0.003", "strict_missing=true", "category=recessive,comphet")
    f = saved(vcfclick_home, "p")["filters"]
    assert f == {
        "category": ["recessive", "comphet"], "min_gq": 21, "max_af": 0.003,
        "strict_missing": True,
    }  # fmt: skip
    run("unset", "p", "min_gq")
    assert "min_gq" not in saved(vcfclick_home, "p")["filters"]
    assert "Did you mean max_af" in run("set", "p", "max_afx=1", ok=False).output
    assert "attestation" in run("set", "p", "attest_affected=1", ok=False).output
    assert "must be between" in run("set", "p", "min_gq=500", ok=False).output
    assert "expected KEY=VALUE" in run("set", "p", "min_gq", ok=False).output

    run("rename", "p", "q")
    assert (
        files(vcfclick_home) == ["q.json"] and saved(vcfclick_home, "q")["name"] == "q"
    )
    assert "built-in" in run("rename", "q", "rare-recessive", ok=False).output
    assert "already exists" in run("rename", "q", "q", ok=False).output
    run("rm", "q", "--yes")
    assert files(vcfclick_home) == []
    assert "no preset named" in run("rm", "q", "--yes", ok=False).output


def test_built_in_and_shared_presets_are_read_only_and_say_how_to_copy(
    run, vcfclick_home, tmp_path, monkeypatch
):
    for verb in (
        ["set", "strict-denovo", "min_gq=31"],
        ["rm", "strict-denovo", "--yes"],
    ):
        r = run(*verb, ok=False)
        assert "built-in" in r.output and "--from-preset strict-denovo" in r.output
    shared = tmp_path / "lab"
    shared.mkdir()
    (shared / "lab-panel.json").write_text(
        ps.to_text(ps.validate_document({
            "schema": ps.SCHEMA, "schema_version": 1, "name": "lab-panel",
            "filters": {"min_gq": 22},
        }).doc)
    )  # fmt: skip
    monkeypatch.setenv("VCFCLICK_PRESET_PATH", str(shared))
    r = run("set", "lab-panel", "min_gq=30", ok=False)
    assert "shared preset" in r.output and "--from-preset lab-panel" in r.output
    assert "lab-panel" in run("list", "--source", "shared").output
    assert run("path").output.splitlines() == [
        str(vcfclick_home / "presets"),
        str(shared),
    ]


def test_locked_blocks_changes_without_force(run, vcfclick_home):
    run("save", "p", "--min-gq", "21")
    doc = saved(vcfclick_home, "p")
    doc["meta"]["locked"] = True
    (vcfclick_home / "presets" / "p.json").write_text(json.dumps(doc))
    for verb in (
        ["set", "p", "min_gq=22"],
        ["unset", "p", "min_gq"],
        ["rename", "p", "q"],
        ["rm", "p", "--yes"],
    ):
        assert "locked" in run(*verb, ok=False).output
    run("set", "p", "min_gq=22", "--force")
    assert saved(vcfclick_home, "p")["filters"]["min_gq"] == 22
    assert saved(vcfclick_home, "p")["meta"]["locked"] is True


# ───────────────────── a name is a name, never a path ─────────────────────


def test_a_name_that_is_a_path_never_reaches_the_file_system(
    run, vcfclick_home, tmp_path
):
    """`preset set ../../victim/ok`, `rm /abs/path/ok` and the rest must leave
    every file outside the presets directory exactly as it was."""
    victim = tmp_path / "victim"
    victim.mkdir()
    run("save", "mine", "--min-gq", "25")
    shutil.copy(vcfclick_home / "presets" / "mine.json", victim / "ok.json")
    (victim / "doc.json").write_text("not a preset\n")
    before = {p.name: p.read_bytes() for p in victim.iterdir()}
    out = tmp_path / "out"
    out.mkdir()
    for target in (f"../../{victim.name}/ok", str(victim / "ok"), "../../victim/doc"):
        for verb in (
            ["set", target, "max_af=0.1"],
            ["unset", target, "max_af"],
            ["rename", target, "moved"],
            ["rm", target, "--yes"],
            ["export", target, "-o", str(out)],
            ["save", "copy", "--from-preset", target],
        ):
            assert run(*verb, ok=False).exit_code != 0, verb
    assert {p.name: p.read_bytes() for p in victim.iterdir()} == before
    assert files(vcfclick_home) == ["mine.json"] and list(out.iterdir()) == []
    assert "not a preset name" in run("rm", "../x", "--yes", ok=False).output
    # the store functions refuse too, whoever calls them
    for name in ("../x", "/abs/x", "a/b", ""):
        assert store.find_preset(name) == (None, [])
        with pytest.raises(ps.PresetError, match="not a preset name"):
            store.user_path(name)
        with pytest.raises(ps.PresetError, match="not a preset name"):
            store.get_user_preset(name)


# ─────────────────────── a file the checks cannot digest ───────────────────────


def test_a_file_with_a_400_digit_number_is_listed_not_a_crash(run, vcfclick_home):
    d = vcfclick_home / "presets"
    d.mkdir(parents=True)
    body = (
        '{"schema": "vcfclick.preset", "schema_version": 1, '
        '"filters": {"max_af": 1%s}}' % ("0" * 400)
    )
    (d / "big.json").write_text(body)
    rows = json.loads(run("list", "--format", "json").stdout)["presets"]
    big = next(p for p in rows if p["name"] == "big")
    assert big["available"] is False and big["problems"][0]["code"] == "bad_value"
    assert store.read_preset_file(d / "big.json", store.USER).has_errors
    r = run("validate", str(d / "big.json"), ok=False)
    assert isinstance(r.exception, SystemExit)  # a clean exit 1, no traceback
    assert "finite number" in r.output
    run("show", "big", ok=True)  # listed, shown, never applied
    assert "cannot be used" in str(
        pytest.raises(ps.PresetError, store.load_named, "big").value
    )


def test_a_text_field_that_is_null_or_the_wrong_type_never_breaks_the_listing(
    run, vcfclick_home
):
    """`preset list --format json` is the call the desktop apps make: a bad file in
    the directory must be listed with its problem, whatever its text fields hold."""
    d = vcfclick_home / "presets"
    d.mkdir(parents=True)
    head = '{"schema": "vcfclick.preset", "schema_version": 1, "filters": {}, '
    for name, field in (
        ("nullapp", '"applies_to": null'),  # accepted as absent
        ("intapp", '"applies_to": 5'),  # an error
        ("strapp", '"applies_to": "trio"'),
        ("nulltitle", '"title": null, "description": null, "assay": null'),
        ("inttitle", '"title": 5, "description": [1], "assay": 7'),
    ):
        (d / f"{name}.json").write_text(head + field + "}")
    for args in (["list", "--format", "json"], ["list", "--command", "trio"],
                 ["list", "--command", "trio", "--format", "json"]):  # fmt: skip
        run(*args)
    rows = {
        p["name"]: p
        for p in json.loads(run("list", "--format", "json").stdout)["presets"]
    }
    assert rows["nullapp"]["applies_to"] == [] and rows["nullapp"]["available"]
    assert rows["nulltitle"]["available"] and rows["nulltitle"]["title"] == "nulltitle"
    for bad in ("intapp", "strapp", "inttitle"):
        assert rows[bad]["available"] is False
        assert rows[bad]["problems"][0]["code"] in ("bad_type", "bad_value"), bad
        assert rows[bad]["applies_to"] == [] and isinstance(rows[bad]["title"], str)
        assert rows[bad]["assay"] is None
    for name in rows:
        if name in ("nullapp", "intapp", "strapp", "nulltitle", "inttitle"):
            run("show", name, "--format", "json")
    # a null is dropped from the normalised document, not carried into a store or export
    checked = ps.validate_document(json.loads((d / "nullapp.json").read_text()))
    assert checked.ok and "applies_to" not in checked.doc
    wrong = ps.validate_document(json.loads((d / "intapp.json").read_text()))
    assert not wrong.ok and "applies_to" not in wrong.doc


def test_an_identity_key_is_refused_however_it_is_spelt_and_wherever_it_sits(
    run, vcfclick_home, tmp_path
):
    """Case, hyphen and separator variants of an identity key are the key, and a
    key under `filters` is held to the same rule (an error in every mode, not the
    warning an unknown filter key is)."""
    base = {"schema": ps.SCHEMA, "schema_version": 1, "name": "idv",
            "filters": {"min_gq": 20}}  # fmt: skip
    for top in (
        "Patient",
        "sample-id",
        "sampleId",
        "SAMPLE_ID",
        "Proband",
        "Family ID",
    ):
        found = ps.validate_document({**base, top: "x"}, allow_unknown=True)
        assert [p.code for p in found.errors] == ["never_presettable"], top
    for key in (
        "proband_id",
        "Sample",
        "family_id",
        "mrn",
        "Family-ID",
        "DB",
        "Proband",
    ):
        for allow in (False, True):
            found = ps.validate_document(
                {**base, "filters": {"min_gq": 20, key: "x"}}, allow_unknown=allow
            )
            assert [p.code for p in found.errors] == ["never_presettable"], (key, allow)
            assert found.unknown_keys == ()
    meta = ps.validate_document({**base, "meta": {"Patient-ID": "x"}})
    assert [p.code for p in meta.errors] == ["never_presettable"]
    # what is not an identity is still an unknown filter key, and the envelope's own
    # keys and the extension namespace are untouched
    plain = ps.validate_document(
        {**base, "filters": {"min_gq": 20, "colour": 1, "x_sample": 1}},
        allow_unknown=True,
    )
    assert plain.ok and plain.unknown_keys == ("colour",)
    assert (
        ps.identity_problem("name") is None and ps.identity_problem("filters") is None
    )

    # import refuses them, even with --allow-unknown-keys, and writes nothing
    f = tmp_path / "idv.json"
    f.write_text(
        json.dumps({**base, "filters": {"min_gq": 20, "Sample": "x", "mrn": 1}})
    )
    r = run("import", str(f), "--allow-unknown-keys", ok=False)
    assert "Sample" in r.output and files(vcfclick_home) == []

    # a hand-edited file in the presets directory is listed with the error, never used
    # or exported
    d = vcfclick_home / "presets"
    d.mkdir(parents=True, exist_ok=True)
    (d / "idv.json").write_text(f.read_text())
    (listed,) = [p for p in store.list_presets("user") if p.name == "idv"]
    assert listed.has_errors and listed.problems[0].code == "never_presettable"
    assert "cannot be exported" in run("export", "idv", ok=False).output
    # and whatever a Preset object holds, an export leaves them out of the filters too
    tainted = store.Preset(
        "idv", store.USER,
        {**base, "filters": {"min_gq": 20, "Sample": "x", "family_id": "F1", "x_k": 1}},
    )  # fmt: skip
    assert exchange.export_document(tainted)["filters"] == {"min_gq": 20, "x_k": 1}


def test_an_exported_file_counts_as_the_preset_when_a_write_asks_whether_it_exists(
    run, vcfclick_home
):
    """Only `a.vcfclick-preset.json` is in the presets directory: it is the preset
    `a`, so save, import and rename must refuse to shadow it, and --force replaces
    that file and not a twin beside it."""
    d = vcfclick_home / "presets"
    run("save", "a", "--min-gq", "30")
    (d / "a.json").rename(d / "a.vcfclick-preset.json")
    assert store.user_file("a") == d / "a.vcfclick-preset.json"
    assert store.user_file("nothing") is None
    with pytest.raises(ps.PresetError, match="not a preset name"):
        store.user_file("../a")

    doc = store.new_document("a", {"min_gq": 40})
    with pytest.raises(ps.PresetError, match="already exists"):
        store.save_user_preset(doc)
    assert "already exists" in run("save", "a", "--min-gq", "25", ok=False).output
    incoming = [exchange.Incoming("a", doc, "a.json")]
    with pytest.raises(ps.PresetError, match="already exists"):
        exchange.import_presets(incoming, force=False)
    run("save", "b", "--min-gq", "20")
    assert "already exists" in run("rename", "b", "a", ok=False).output
    assert files(vcfclick_home) == ["a.vcfclick-preset.json", "b.json"]

    # --force writes to the file that is there; a twin NAME.json is never created
    assert store.save_user_preset(doc, force=True) == d / "a.vcfclick-preset.json"
    assert files(vcfclick_home) == ["a.vcfclick-preset.json", "b.json"]
    assert json.loads((d / "a.vcfclick-preset.json").read_text())["filters"] == {
        "min_gq": 40
    }
    again = store.new_document("a", {"min_gq": 50})
    (written,) = exchange.import_presets(
        [exchange.Incoming("a", again, "x")], force=True
    )
    assert written == d / "a.vcfclick-preset.json"
    assert files(vcfclick_home) == ["a.vcfclick-preset.json", "b.json"]
    assert [p.name for p in store.list_presets("user")] == [
        "a",
        "b",
    ]  # no duplicate_name
    # a failing import puts the original back into that same file
    before = (d / "a.vcfclick-preset.json").read_text()
    bad = exchange.Incoming("c", {"not": "serialisable", "x": {1, 2}}, "c.json")
    with pytest.raises(TypeError):
        exchange.import_presets([exchange.Incoming("a", doc, "a"), bad], force=True)
    assert (d / "a.vcfclick-preset.json").read_text() == before
    assert files(vcfclick_home) == ["a.vcfclick-preset.json", "b.json"]


# ─────────────────── a file that leaves and comes back unchanged ───────────────────


def test_import_refuses_a_number_that_is_not_finite_before_writing(
    run, vcfclick_home, tmp_path
):
    f = tmp_path / "inf.json"
    f.write_text(
        '{"schema": "vcfclick.preset", "schema_version": 1, "name": "inf", '
        '"filters": {"x_k": 1e999}}'
    )
    r = run("import", str(f), ok=False)
    assert "out of range" in r.output and files(vcfclick_home) == []
    assert "out of range" in run("validate", str(f), ok=False).output


def test_a_preset_carries_no_identity_through_save_import_or_export(
    run, vcfclick_home, tmp_path
):
    base = {
        "schema": ps.SCHEMA, "schema_version": 1, "name": "idp",
        "filters": {"min_gq": 20},
    }  # fmt: skip
    f = tmp_path / "idp.json"
    for extra in (
        {"proband": "NA12878"},
        {"db": "/home/me/cohort"},
        {"family": "F1"},
        {"meta": {"family": "F1"}},
    ):
        f.write_text(json.dumps({**base, **extra}))
        key = next(iter(extra if "meta" not in extra else extra["meta"]))
        if "meta" in extra:  # `save` writes its own meta and never copies one
            run("save", "idp", "--from-json", str(f))
            assert "family" not in saved(vcfclick_home, "idp")["meta"]
            run("rm", "idp", "--yes")
        else:
            assert key in run("save", "idp", "--from-json", str(f), ok=False).output
        assert key in run("import", str(f), ok=False).output
        assert key in run("validate", str(f), ok=False).output
    assert files(vcfclick_home) == []
    # a key that names no identity is still carried (forward compatibility)
    f.write_text(json.dumps({**base, "lab_note": 1}))
    run("save", "idp", "--from-json", str(f))
    assert saved(vcfclick_home, "idp")["lab_note"] == 1
    # whatever a Preset object holds, an export does not carry an identity key
    tainted = store.Preset(
        "idp", store.USER, {**saved(vcfclick_home, "idp"), "proband": "NA12878"}
    )
    exported = exchange.export_document(tainted)
    assert "proband" not in exported and exported["lab_note"] == 1


def test_an_exported_file_in_a_shared_directory_is_found_by_its_name_and_pin(
    run, vcfclick_home, tmp_path, monkeypatch
):
    """The lab git-folder workflow: `export -o DIR`, DIR on VCFCLICK_PRESET_PATH,
    `--preset NAME@FINGERPRINT`."""
    run("save", "sd2", "--from-preset", "strict-denovo")
    fingerprint = ps.fingerprint(saved(vcfclick_home, "sd2"))
    lab = tmp_path / "lab"
    lab.mkdir()
    run("export", "sd2", "-o", str(lab))
    assert [p.name for p in lab.iterdir()] == ["sd2.vcfclick-preset.json"]
    (vcfclick_home / "presets" / "sd2.json").unlink()
    monkeypatch.setenv("VCFCLICK_PRESET_PATH", str(lab))

    found = store.load_named(f"sd2@{fingerprint}")
    assert (found.preset.name, found.preset.source) == ("sd2", "shared")
    assert found.pinned and found.preset.problems == ()
    (row,) = store.list_presets("shared")
    assert row.name == "sd2" and row.problems == ()
    assert "sd2" in run("list", "--source", "shared").output

    # NAME.json next to it is the same name; the plain file is the one used
    shutil.copy(lab / "sd2.vcfclick-preset.json", lab / "sd2.json")
    assert store.find_preset("sd2")[0].path.name == "sd2.json"
    rows = store.list_presets("shared")
    assert [p.name for p in rows] == ["sd2", "sd2"]
    assert "duplicate_name" in [x.code for x in rows[1].problems]
    assert rows[0].problems == ()
    # a preset cannot be named so that its own export file reads back as another
    assert ps.name_problem("sd2.vcfclick-preset")
    assert "vcfclick-preset" in run("save", "x.vcfclick-preset", ok=False).output


def test_edit_revalidates_and_does_not_save_an_invalid_file(
    run, vcfclick_home, monkeypatch
):
    run("save", "p", "--min-gq", "21")
    original = (vcfclick_home / "presets" / "p.json").read_text()
    edited = json.loads(original)
    edited["filters"]["min_gq"] = 200
    monkeypatch.setattr("click.edit", lambda *a, **k: json.dumps(edited))
    r = run("edit", "p", ok=False)
    assert "not saved" in r.output and "min_gq" in r.output
    assert (vcfclick_home / "presets" / "p.json").read_text() == original

    edited["filters"]["min_gq"] = 23
    monkeypatch.setattr("click.edit", lambda *a, **k: json.dumps(edited))
    run("edit", "p")
    assert saved(vcfclick_home, "p")["filters"]["min_gq"] == 23

    edited["name"] = "other"
    monkeypatch.setattr("click.edit", lambda *a, **k: json.dumps(edited))
    assert "rename" in run("edit", "p", ok=False).output
    monkeypatch.setattr("click.edit", lambda *a, **k: None)
    assert "no changes" in run("edit", "p").output


def test_a_file_copied_in_under_another_name_is_the_preset_its_file_name_says(
    run, vcfclick_home, tmp_path, monkeypatch
):
    """A colleague's `q1.json` copied in as `lab-q.json` still declares `q1`. The
    file name is what the preset is called, so every writer must write that name:
    `set` must not leave the old one, `export` must not hand on a file that imports
    as `q1` (and replaces another preset with --force), and `edit` must be able to
    save without being told to rename."""
    d = vcfclick_home / "presets"
    run("save", "q1", "--min-gq", "21")
    shutil.copy(d / "q1.json", d / "lab-q.json")
    (d / "q1.json").unlink()
    mismatched = (d / "lab-q.json").read_text()
    assert json.loads(mismatched)["name"] == "q1"

    (row,) = store.list_presets("user")
    assert (row.name, row.doc["name"]) == ("lab-q", "lab-q")
    assert [p.code for p in row.problems] == ["name_mismatch"]  # still said, a warning
    assert list(row.doc)[:3] == ["schema", "schema_version", "name"]  # in its place

    run("set", "lab-q", "min_gq=22")
    assert saved(vcfclick_home, "lab-q")["name"] == "lab-q"

    (d / "lab-q.json").write_text(mismatched)  # declares q1 again
    exported = json.loads(run("export", "lab-q").output)
    assert exported["name"] == "lab-q"  # a plain import elsewhere gives lab-q back
    incoming = tmp_path / "in.json"
    incoming.write_text(json.dumps(exported))
    (item,) = exchange.read_incoming([str(incoming)], as_name=None, allow_unknown=False)
    assert item.name == "lab-q"

    shown: list[str] = []

    def editor(text, **kwargs):
        shown.append(json.loads(text)["name"])
        edited = json.loads(text)
        edited["filters"]["min_gq"] = 24
        return json.dumps(edited)

    monkeypatch.setattr("click.edit", editor)
    run("edit", "lab-q")  # not "the name cannot change here"
    assert shown == ["lab-q"]
    again = saved(vcfclick_home, "lab-q")
    assert (again["name"], again["filters"]["min_gq"]) == ("lab-q", 24)

    # no name at all, and a shared directory's file: the same rule
    bare = {"schema": ps.SCHEMA, "schema_version": 1, "filters": {"min_gq": 20}}
    (d / "bare.json").write_text(json.dumps(bare))
    assert store.find_preset("bare")[0].doc["name"] == "bare"
    shared = tmp_path / "shared"
    shared.mkdir()
    (shared / "team-q.json").write_text(mismatched)
    monkeypatch.setenv("VCFCLICK_PRESET_PATH", str(shared))
    (found,) = store.list_presets("shared")
    assert (found.name, found.doc["name"]) == ("team-q", "team-q")


# ───────────────────────── atomic writes, corrupt files ─────────────────────────


def test_a_write_that_dies_leaves_the_old_file_and_no_temp_file(
    vcfclick_home, monkeypatch
):
    path = vcfclick_home / "presets" / "p.json"
    store.write_text_atomic(path, "old\n")

    def boom(src, dst):
        raise OSError("killed between the temp file and the replace")

    monkeypatch.setattr(os, "replace", boom)
    with pytest.raises(OSError):
        store.write_text_atomic(path, "new\n")
    monkeypatch.undo()
    assert path.read_text() == "old\n"
    assert [p.name for p in path.parent.iterdir()] == ["p.json"]


def test_a_corrupt_file_is_listed_with_its_problem_and_breaks_nothing(
    run, vcfclick_home
):
    run("save", "good", "--min-gq", "21")
    d = vcfclick_home / "presets"
    (d / "broken.json").write_text("{ not json")
    (d / "Shouting.json").write_text(json.dumps(saved(vcfclick_home, "good")))
    (d / "rare-recessive.json").write_text(json.dumps(saved(vcfclick_home, "good")))
    (d / ".tmp-leftover.json").write_text("x")
    (d / "notes.txt").write_text("x")

    r = run("list")  # exits 0 whatever it finds
    assert "good" in r.output and "not valid JSON" in r.output
    rows = {p["name"]: p for p in json.loads(run("list", "--format", "json").stdout)["presets"]
            if p["source"] == "user"}  # fmt: skip
    assert set(rows) == {"good", "broken", "Shouting", "rare-recessive"}
    assert rows["good"]["available"] and not rows["broken"]["available"]
    assert rows["broken"]["problems"][0]["code"] == "invalid_json"
    assert rows["Shouting"]["problems"][0]["code"] == "bad_name"
    assert rows["rare-recessive"]["problems"][0]["code"] == "builtin_name"
    assert (d / "broken.json").exists()  # never deleted
    # ... and the built-in of that name is still the one that resolves
    assert store.load_named("rare-recessive").preset.source == "builtin"


def test_a_user_preset_shadows_a_shared_one_with_a_note(
    vcfclick_home, tmp_path, monkeypatch
):
    shared = tmp_path / "lab"
    shared.mkdir()
    doc = ps.validate_document(
        {"schema": ps.SCHEMA, "schema_version": 1, "filters": {"min_gq": 22}}
    ).doc
    (shared / "p.json").write_text(ps.to_text(doc))
    store.write_text_atomic(
        store.user_path("p"), ps.to_text({**doc, "filters": {"min_gq": 30}})
    )
    monkeypatch.setenv("VCFCLICK_PRESET_PATH", str(shared))
    found = store.load_named("p")
    assert found.preset.source == "user" and found.preset.filters == {"min_gq": 30}
    assert "using yours" in found.notes[0]
    listed = [p for p in store.list_presets() if p.name == "p"]
    assert [p.source for p in listed] == ["user", "shared"]
    assert [x.code for x in listed[1].problems] == ["shadowed"]


# ───────────────────────────── export and import ─────────────────────────────


def test_export_import_round_trip_records_a_basename_and_the_original_fingerprint(
    run, vcfclick_home, tmp_path
):
    run("save", "mine", "--from-preset", "rare-recessive", "--max-af", "0.004")
    out = tmp_path / "out"
    out.mkdir()
    run("export", "mine", "-o", str(out))
    exported = out / "mine.vcfclick-preset.json"
    sent = json.loads(exported.read_text())
    assert sent["meta"]["imported_from"] is None
    fingerprint = ps.fingerprint(sent)

    run("rm", "mine", "--yes")
    run("import", str(exported))
    got = saved(vcfclick_home, "mine")
    assert got["meta"]["imported_from"] == "mine.vcfclick-preset.json"  # not a path
    assert got["meta"]["original_fingerprint"] == fingerprint
    assert ps.fingerprint(got) == fingerprint
    assert str(tmp_path) not in json.dumps(got)

    assert "already exists" in run("import", str(exported), ok=False).output
    run("import", str(exported), "--force")
    run("import", str(exported), "--as", "copy")
    assert "copy.json" in files(vcfclick_home)


def test_import_refuses_a_built_in_name_and_is_all_or_nothing(
    run, vcfclick_home, tmp_path
):
    exported = tmp_path / "builtin.json"
    run("export", "rare-dominant", "-o", str(exported))
    r = run("import", str(exported), ok=False)
    assert "built-in" in r.output and "--as" in r.output
    run("import", str(exported), "--as", "dominant-copy")

    good = tmp_path / "good.json"
    bad = tmp_path / "bad.json"
    run("export", "rare-recessive", "-o", str(good))
    bundle = json.loads(good.read_text())
    bundle["name"] = "from-good"
    good.write_text(json.dumps(bundle))
    broken = {**bundle, "name": "from-bad", "filters": {"min_gq": 999}}
    bad.write_text(json.dumps(broken))
    before = files(vcfclick_home)
    r = run("import", str(good), str(bad), ok=False)
    assert "bad.json" in r.output and "min_gq" in r.output
    assert files(vcfclick_home) == before  # the good one was not written either

    # a name that is taken fails the whole import, before anything is written
    taken = {**bundle, "name": "dominant-copy"}
    (tmp_path / "taken.json").write_text(json.dumps(taken))
    assert (
        "already exists"
        in run("import", str(good), str(tmp_path / "taken.json"), ok=False).output
    )
    assert files(vcfclick_home) == before


def test_bundle_export_and_import(run, vcfclick_home, tmp_path):
    run("save", "a", "--min-gq", "21")
    run("save", "b", "--min-gq", "22")
    bundle = tmp_path / "pair.json"
    run("export", "a", "b", "--bundle", "-o", str(bundle))
    assert json.loads(bundle.read_text())["schema"] == ps.BUNDLE_SCHEMA
    assert "--bundle" in run("export", "a", "b", ok=False).output
    run("rm", "a", "--yes")
    run("rm", "b", "--yes")
    run("import", str(bundle))
    assert files(vcfclick_home) == ["a.json", "b.json"]
    assert (
        "several" in run("import", str(bundle), "--as", "x", "--force", ok=False).output
    )


def test_export_leaves_out_patient_specific_keys_unless_asked(
    vcfclick_home, monkeypatch
):
    from cli import filter_registry as fr
    from cli.preset_exchange import export_document

    monkeypatch.setattr(fr, "REGISTRY", dict(fr.REGISTRY))
    fr.register(fr.FilterKey(
        "hpo_terms", "hpo_terms", "--hpo", "category", [], "test", ("trio",),
        patient_specific=True,
    ))  # fmt: skip
    preset = store.Preset(
        "p", store.USER,
        {"schema": ps.SCHEMA, "schema_version": 1,
         "filters": {"min_gq": 20, "hpo_terms": ["HP:1"]}, "meta": {}},
    )  # fmt: skip
    assert export_document(preset)["filters"] == {"min_gq": 20}
    assert export_document(preset, with_patient_data=True)["filters"]["hpo_terms"] == [
        "HP:1"
    ]


def test_import_allow_unknown_keys_stores_it_and_it_is_unavailable_to_apply(
    run, vcfclick_home, tmp_path
):
    doc = {"schema": ps.SCHEMA, "schema_version": 1, "name": "future",
           "filters": {"min_gq": 20, "brand_new": 3},
           "meta": {"created_with": "9.9.9"}}  # fmt: skip
    f = tmp_path / "future.json"
    f.write_text(json.dumps(doc))
    r = run("import", str(f), ok=False)
    assert "did you mean" in r.output or "unknown filter key" in r.output
    assert "written by vcfclick 9.9.9" in r.output
    run("import", str(f), "--allow-unknown-keys")
    row = next(p for p in json.loads(run("list", "--format", "json").stdout)["presets"]
               if p["name"] == "future")  # fmt: skip
    assert row["available"] is False
    plan_err = pytest.raises(ps.PresetError)
    with plan_err as e:
        store.plan_apply(store.load_named("future").preset, "trio")
    assert "--preset-skip-unavailable" in str(e.value)
    plan = store.plan_apply(
        store.load_named("future").preset, "trio", skip_unavailable=True
    )
    assert plan.unavailable == ("brand_new",) and plan.applicable == {"min_gq": 20}


# ───────────────────────────── validate, from-result ─────────────────────────────


def test_validate_exits_1_on_an_error_and_names_the_nearest_key(run, tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"schema": ps.SCHEMA, "schema_version": 1,
                               "filters": {"min_gp": 20}}))  # fmt: skip
    r = run("validate", str(bad), ok=False)
    assert "did you mean min_gq" in r.output
    out = json.loads(run("validate", str(bad), "--format", "json", ok=False).stdout)
    assert out["valid"] is False and out["problems"][0]["code"] == "unknown_key"
    good = tmp_path / "good.json"
    good.write_text(json.dumps({"schema": ps.SCHEMA, "schema_version": 1,
                                "filters": {"min_gq": 20, "parent_alt_screen": True}}))  # fmt: skip
    assert "ok: fingerprint" in run("validate", str(good)).output
    r = run("validate", str(good), "--command", "segregate")
    assert "ok" in r.output
    assert "parent_alt_screen is not a filter of `segregate`" in r.output
    missing = run("validate", str(tmp_path / "nope.json"), ok=False)
    assert "cannot read" in missing.output


def _result_with(block) -> str:
    return json.dumps({"trio": {"proband": "C"}, "preset": block})


def test_save_from_result_rebuilds_the_filters_with_the_original_fingerprint(
    run, vcfclick_home, tmp_path
):
    p = store.builtin_preset("rare-recessive")
    block = {
        "name": "rare-recessive", "source": "builtin", "revision": 1,
        "fingerprint": p.fingerprint,
        "resolved": {
            k: {"value": v, "from": "preset"} for k, v in p.filters.items()
        } | {"pass_only": {"value": False, "from": "default"}},
    }  # fmt: skip
    result = tmp_path / "result.json"
    result.write_text(_result_with(block))
    run("save", "redo", "--from-result", str(result))
    got = saved(vcfclick_home, "redo")
    assert ps.fingerprint(got) == p.fingerprint  # the preset's own, not a lookalike
    assert got["filters"] == p.filters

    # a flag that overrode the preset is part of what ran; a default is not
    block["resolved"]["max_af"] = {"value": 0.005, "from": "flag"}
    result.write_text(_result_with(block))
    run("save", "redo2", "--from-result", str(result))
    assert saved(vcfclick_home, "redo2")["filters"]["max_af"] == 0.005


def test_save_from_result_says_when_a_built_in_has_changed_since(
    run, vcfclick_home, tmp_path
):
    block = {
        "name": "rare-recessive", "source": "builtin", "revision": 0,
        "fingerprint": "aaaaaaaaaaaa",
        "resolved": {"min_gq": {"value": 20, "from": "preset"}},
    }  # fmt: skip
    result = tmp_path / "r.json"
    result.write_text(_result_with(block))
    r = run("save", "redo", "--from-result", str(result))
    assert "rare-recessive rev 0 aaaaaaaaaaaa -> now rev 1 c19bb0a98ece" in r.output


def test_save_from_a_result_without_a_preset_block_says_what_to_do(run, tmp_path):
    result = tmp_path / "r.json"
    result.write_text(json.dumps({"trio": {}, "gates": {}}))
    r = run("save", "x", "--from-result", str(result), ok=False)
    assert "no preset block" in r.output and "same flags" in r.output
    docs = tmp_path / "many.json"
    docs.write_text(json.dumps({"trios": [{"trio": {}}, {"preset": {
        "name": "u", "source": "user", "resolved": {"min_gq": {"value": 21, "from": "flag"}}
    }}]}))  # fmt: skip
    run("save", "fromfirst", "--from-result", str(docs))
