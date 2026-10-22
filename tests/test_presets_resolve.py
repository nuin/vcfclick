"""Using a preset: precedence, pins, the name-only resolver, `--category` lists.

A dummy Click command with the real `PresetCommand`, `gate_options` and
`preset_option` stands in for `db trio`, so everything here is in-process and
opens no database. The precedence rule under test: a flag typed on the command
line (or in the environment) beats the preset, which beats the default; a
value that came from a preset is never "explicit".
"""

from __future__ import annotations

import json

import click
import pytest
from click.testing import CliRunner

from cli import filter_registry as fr
from cli import preset_schema as ps
from cli import presets as store
from cli.preset_options import (
    EXPLICIT_META,
    CategoryList,
    PresetCommand,
    gate_options,
    is_explicit,
    preset_option,
    preset_state,
    resolve_filters,
    scan_preset_args,
)
from cli.trio_analysis import CATEGORY_CHOICES


class LayeredCommand(PresetCommand):
    """A command with a layer between the defaults and the preset (the case
    file of the case command): min_gq 22, and max_ab counts as typed."""

    def default_layer(self, ctx, args):
        ctx.meta[EXPLICIT_META] = {"max_ab"}
        return {"min_gq": 22, "max_ab": 0.9}


def make(cls=PresetCommand):
    @click.command(name="trio", cls=cls)
    @click.option("--proband")
    @gate_options("trio")
    @preset_option()
    @click.pass_context
    def trio(ctx, **kw):
        state = preset_state(ctx)
        resolved, overridden = resolve_filters(ctx, "trio", state)
        click.echo(
            json.dumps(
                {
                    "values": {
                        k: list(v) if isinstance(v, tuple) else v for k, v in kw.items()
                    },
                    "source": {
                        k.param: ctx.get_parameter_source(k.param).name
                        for k in fr.keys_for("trio")
                    },
                    "explicit": {
                        k.param: is_explicit(ctx, k.param) for k in fr.keys_for("trio")
                    },
                    "resolved": {k: [r.value, r.source] for k, r in resolved.items()},
                    "overridden": overridden,
                    "skipped": list(state.plan.not_applicable) if state else None,
                    "pinned": state.located.pinned if state else None,
                    "source_kind": state.located.preset.source if state else None,
                }  # fmt: skip
            )
        )

    return trio


@pytest.fixture
def put(vcfclick_home, monkeypatch):
    monkeypatch.delenv("VCFCLICK_PRESET_PATH", raising=False)

    def write(name, filters, **extra):
        doc = ps.validate_document(
            {"schema": ps.SCHEMA, "schema_version": 1, "filters": filters, **extra},
            allow_unknown=True,
        ).doc
        store.write_text_atomic(
            store.user_path(name), ps.to_text({**doc, "name": name})
        )
        return ps.fingerprint(doc)

    return write


def invoke(*args, cls=PresetCommand, ok=True):
    r = CliRunner().invoke(make(cls), list(args))
    if ok:
        assert r.exit_code == 0, r.output
        return json.loads(r.stdout)
    assert r.exit_code != 0, r.output
    return r


# ───────────────────────────────── precedence ─────────────────────────────────


def test_no_preset_changes_nothing(vcfclick_home):
    out = invoke("--min-gq", "25")
    assert (
        out["source"]["min_gq"] == "COMMANDLINE"
        and out["source"]["max_af"] == "DEFAULT"
    )
    assert out["values"]["category"] == ["all"]
    assert out["pinned"] is None and out["skipped"] is None
    assert all(r[1] in ("flag", "default") for r in out["resolved"].values())


def test_a_preset_supplies_defaults_and_a_flag_wins(put):
    put("p", {"min_gq": 30, "max_af": 0.001, "category": ["recessive", "comphet"]})
    out = invoke("--preset", "p", "--max-af", "0.005")
    assert out["values"]["min_gq"] == 30 and out["source"]["min_gq"] == "DEFAULT_MAP"
    assert out["values"]["max_af"] == 0.005 and out["source"]["max_af"] == "COMMANDLINE"
    assert out["values"]["category"] == ["recessive", "comphet"]  # a tuple, as typed
    assert out["values"]["min_dp"] == 10 and out["source"]["min_dp"] == "DEFAULT"
    assert out["overridden"] == [{"key": "max_af", "preset": 0.001, "flag": 0.005}]
    # the registry's resolution says the same as Click's own sources
    got = {k: (v[0], v[1]) for k, v in out["resolved"].items()}
    assert got["min_gq"] == (30, "preset") and got["max_af"] == (0.005, "flag")
    assert got["min_dp"] == (10, "default")
    for key in fr.keys_for("trio"):
        value = out["values"][key.param]
        assert got[key.file_key][0] == value, key.file_key


def test_option_order_on_the_command_line_does_not_matter(put):
    put("p", {"min_gq": 30, "max_af": 0.001})
    first = invoke("--preset", "p", "--max-af", "0.005", "--min-gq", "31")
    second = invoke("--max-af", "0.005", "--min-gq", "31", "--preset", "p")
    third = invoke("--preset=p", "--max-af=0.005", "--min-gq=31")
    assert first == second == third
    assert first["values"]["min_gq"] == 31


def test_a_preset_value_is_never_explicit_a_typed_one_always_is(put):
    put("p", {"max_af": 0.001, "min_gq": 30})
    out = invoke("--preset", "p", "--min-gq", "30")  # equal to the preset's own
    assert out["explicit"]["max_af"] is False  # so the cohort-AF skip still applies
    assert out["explicit"]["min_gq"] is True
    assert out["source"]["min_gq"] == "COMMANDLINE"
    assert out["resolved"]["min_gq"] == [30, "flag"]
    assert out["overridden"] == [{"key": "min_gq", "preset": 30, "flag": 30}]


def test_a_boolean_flag_can_turn_a_preset_boolean_off(put):
    put("p", {"strict_missing": True, "pass_only": True})
    out = invoke("--preset", "p", "--no-strict-missing")
    assert (
        out["values"]["strict_missing"] is False and out["values"]["pass_only"] is True
    )
    assert out["overridden"] == [
        {"key": "strict_missing", "preset": True, "flag": False}
    ]


def test_preset_ignore_switches_a_key_off_whichever_side_it_is_typed_on(put):
    put("p", {"gnomad_max_af": 0.0001, "min_gq": 30})
    for args in (
        ("--preset", "p", "--preset-ignore", "gnomad_max_af"),
        ("--preset-ignore", "gnomad_max_af", "--preset", "p"),
    ):
        out = invoke(*args)
        assert out["values"]["gnomad_max_af"] is None
        assert out["values"]["min_gq"] == 30
        assert out["overridden"] == [
            {"key": "gnomad_max_af", "preset": 0.0001, "flag": "off"}
        ]
        assert out["resolved"]["gnomad_max_af"] == [None, "default"]


def test_a_typed_category_replaces_the_presets_list(put):
    put("p", {"category": ["recessive", "comphet"]})
    out = invoke("--preset", "p", "--category", "denovo")
    assert out["values"]["category"] == ["denovo"]
    assert out["overridden"][0]["preset"] == ["recessive", "comphet"]


def test_a_band_made_empty_by_a_flag_is_refused_and_a_flag_can_repair_it(put):
    put("p", {"min_ab": 0.8})
    assert "above max_ab" in invoke("--preset", "p", ok=False).output
    assert invoke("--preset", "p", "--max-ab", "0.9")["values"]["max_ab"] == 0.9


def test_a_layer_between_default_and_preset(put):
    put("p", {"min_gq": 30})
    out = invoke("--preset", "p", cls=LayeredCommand)
    assert out["values"]["min_gq"] == 30  # preset over the layer
    assert out["values"]["max_ab"] == 0.9 and out["explicit"]["max_ab"] is True
    assert invoke(cls=LayeredCommand)["values"]["min_gq"] == 22  # layer over default


def test_help_works_with_a_preset_that_does_not_exist():
    r = CliRunner().invoke(make(), ["--preset", "nope", "--help"])
    assert r.exit_code == 0 and "--preset-ignore" in r.output


# ─────────────────────── locating, pins, unavailable keys ───────────────────────


def test_scan_finds_the_four_options_and_stops_at_double_dash():
    scan = scan_preset_args(
        ["--min-gq", "3", "--preset=a", "--preset-ignore", "x", "--preset-ignore=y",
         "--preset-fingerprint", "abc", "--preset-skip-unavailable", "--", "--preset", "z"]
    )  # fmt: skip
    assert (scan.preset, scan.fingerprint, scan.ignore, scan.skip_unavailable) == (
        "a", "abc", ("x", "y"), True,
    )  # fmt: skip
    assert not scan_preset_args(["--min-gq", "3"]).any


def test_a_pin_must_match_and_the_message_says_what_changed(put):
    fp = put("p", {"min_gq": 30})
    assert invoke("--preset", f"p@{fp}")["pinned"] is True
    assert invoke("--preset", "p", "--preset-fingerprint", fp)["pinned"] is True
    out = invoke("--preset", "p@000000000000", ok=False).output
    assert f"fingerprint {fp}" in out and "pinned 000000000000" in out
    assert "not a fingerprint" in invoke("--preset", "p@xyz", ok=False).output
    assert (
        "two different"
        in invoke(
            "--preset", f"p@{fp}", "--preset-fingerprint", "0" * 12, ok=False
        ).output
    )


def test_a_built_ins_old_fingerprint_is_recognised(vcfclick_home, monkeypatch):
    from cli import builtin_presets as bp

    b = bp.BUILTINS["rare-recessive"]
    newer = bp.BuiltinPreset(
        **{**b.__dict__, "revision": 2, "filters": {**b.filters, "min_gq": 21},
           "history": (*b.history, {"revision": 2, "fingerprint": "f" * 12, "note": "x"})}
    )  # fmt: skip
    monkeypatch.setitem(bp.BUILTINS, "rare-recessive", newer)
    out = invoke("--preset", f"rare-recessive@{b.history[0]['fingerprint']}", ok=False)
    assert "revision 1 had that fingerprint; the current revision is 2" in out.output


def test_the_path_form_reads_a_file_and_bypasses_the_lookup(tmp_path, put):
    fp = put("p", {"min_gq": 30})
    f = tmp_path / "mine.json"
    f.write_text((store.user_path("p")).read_text())
    for value in (str(f), "./mine.json"):
        if value.startswith("."):
            continue
        out = invoke("--preset", value)
        assert out["source_kind"] == "file" and out["values"]["min_gq"] == 30
    assert invoke("--preset", str(f), "--preset-fingerprint", fp)["pinned"] is True
    assert (
        "no preset file"
        in invoke("--preset", str(tmp_path / "x.json"), ok=False).output
    )
    assert (
        str(tmp_path)
        not in invoke("--preset", str(tmp_path / "x.json"), ok=False).output
    )


def test_preset_companions_without_a_preset_are_a_usage_error():
    for opt in (("--preset-ignore", "min_gq"), ("--preset-skip-unavailable",)):
        r = invoke(*opt, ok=False)
        assert r.exit_code == 2 and "need --preset" in r.output


def test_ignoring_an_unknown_key_suggests_the_right_one(put):
    put("p", {"min_gq": 30})
    out = invoke("--preset", "p", "--preset-ignore", "gnomad_max_a", ok=False).output
    assert "Did you mean gnomad_max_af" in out


def test_unavailable_keys_are_an_error_unless_skipping_is_asked_for(put):
    put("future", {"min_gq": 30, "brand_new": 1})
    out = invoke("--preset", "future", ok=False).output
    assert "brand_new" in out and "--preset-skip-unavailable" in out
    ok = invoke("--preset", "future", "--preset-skip-unavailable")
    assert ok["values"]["min_gq"] == 30


def test_keys_the_command_does_not_take_are_skipped_not_applied():
    preset = store.Preset(
        "p", store.USER,
        {"schema": ps.SCHEMA, "schema_version": 1,
         "filters": {"parent_alt_screen": True, "min_gq": 30, "x_lab": 1}},
    )  # fmt: skip
    plan = store.plan_apply(preset, "segregate")
    assert plan.not_applicable == ("parent_alt_screen",)
    assert plan.applicable == {
        "min_gq": 30
    }  # x_ keys are kept in the file, never applied
    assert store.plan_apply(preset, "trio").not_applicable == ()


# ───────────────────────── the name-only resolver (web, MCP) ─────────────────────────


@pytest.mark.parametrize(
    "value",
    ["/etc/passwd", "../x", "a/b", "x.json", "C:\\x", "", "UPPER", ".hidden", "a b"],
)
def test_the_name_only_resolver_refuses_a_path_before_touching_the_file_system(
    value, monkeypatch, tmp_path
):
    def forbidden(*a, **k):
        raise AssertionError("the file system was touched")

    monkeypatch.setattr(store, "find_preset", forbidden)
    monkeypatch.setattr(store, "read_preset_file", forbidden)
    with pytest.raises(ps.PresetError):
        store.load_named(value)
    with pytest.raises(ps.PresetError):
        store.locate(value, allow_path=False)


def test_an_error_never_quotes_the_files_content(put, tmp_path):
    secret = tmp_path / "secret.json"
    secret.write_text("TOP-SECRET root:x:0:0")
    with pytest.raises(ps.PresetError) as e:
        store.locate(str(secret))
    assert "TOP-SECRET" not in str(e.value)
    put("p", {"min_gq": 1})
    store.user_path("p").write_text("TOP-SECRET")
    with pytest.raises(ps.PresetError) as e:
        store.load_named("p")
    assert "TOP-SECRET" not in str(e.value)


def test_load_named_takes_no_skip_option():
    """The resolver the web and MCP use cannot be told to skip unavailable keys:
    `skip_unavailable` belongs to the plan, which a front end calls without it."""
    import inspect

    assert "skip_unavailable" not in inspect.signature(store.load_named).parameters
    assert (
        inspect.signature(store.plan_apply).parameters["skip_unavailable"].default
        is False
    )


# ───────────────────────────────── --category ─────────────────────────────────


def test_category_list_converts_names_lists_and_all():
    t = CategoryList(CATEGORY_CHOICES)
    assert t.convert("denovo", None, None) == ("denovo",)
    assert t.convert("all", None, None) == ("all",)
    assert t.convert("recessive,comphet", None, None) == ("recessive", "comphet")
    assert t.convert(["recessive", "comphet"], None, None) == ("recessive", "comphet")
    assert t.convert("recessive, comphet,recessive", None, None) == (
        "recessive",
        "comphet",
    )


def test_category_list_keeps_the_choice_metavar_and_the_invalid_wording():
    plain = click.Choice(CATEGORY_CHOICES)
    mine = CategoryList(CATEGORY_CHOICES)
    param = click.Option(["--category"])
    assert mine.get_metavar(
        param, click.Context(click.Command("x"))
    ) == plain.get_metavar(param, click.Context(click.Command("x")))
    with pytest.raises(click.BadParameter) as old:
        plain.convert("nope", param, None)
    for bad in ("nope", "recessive,nope", "", ","):
        with pytest.raises(click.BadParameter) as new:
            mine.convert(bad, param, None)
        assert "is not one of" in new.value.message
    with pytest.raises(click.BadParameter) as single:
        mine.convert("nope", param, None)
    assert single.value.message == old.value.message
