"""The filter registry on its own: shape, value checks, precedence.

`cli/filter_registry.py` is pure (no Click, no SQL, no file system), so every
rule here is a unit test. That the Click options of `db trio` agree with the
registry is `tests/test_filter_registry_drift.py`.
"""

from __future__ import annotations

import pytest

from cli import filter_registry as fr


def test_registry_shape():
    assert len({k.param for k in fr.REGISTRY.values() if k.param}) == len(
        [k for k in fr.REGISTRY.values() if k.param]
    )
    for key in fr.REGISTRY.values():
        assert key.owner and key.commands
        assert not set(key.commands) - set(fr.COMMANDS)
    # a name cannot be both a filter and never presettable
    assert not set(fr.REGISTRY) & set(fr.NEVER_PRESET)
    assert all(fr.NEVER_PRESET.values()), "every refusal carries its reason"
    # the flags the registry says exist and the apps rely on keep their names
    assert fr.REGISTRY["require_bound"].flag == (
        "--require-bound-gates/--no-require-bound-gates"
    )
    assert fr.REGISTRY["parent_alt_screen"].commands == ("trio",)


def test_registering_twice_or_for_an_unknown_command_is_refused(monkeypatch):
    monkeypatch.setattr(fr, "REGISTRY", dict(fr.REGISTRY))
    with pytest.raises(ValueError, match="already registered"):
        fr.register(fr.REGISTRY["min_gq"])
    with pytest.raises(ValueError, match="unknown commands"):
        fr.register(fr.FilterKey("x", "x", "--x", "int", 1, "t", ("nope",), lo=0, hi=1))


def test_suggest_finds_a_misspelt_key_and_a_hyphenated_one():
    assert fr.suggest("min_gp")[0] == "min_gq"
    assert fr.suggest("max-af") == ["max_af"]
    assert fr.suggest("zzzzzz") == []


@pytest.mark.parametrize(
    "key,value,ok",
    [
        ("min_gq", 20, True),
        ("min_gq", 100, False),
        ("min_gq", True, False),  # Click's int type would take it
        ("min_gq", 1.5, False),
        ("min_gq", None, False),
        ("max_af", 0, True),
        ("max_af", 1, True),
        ("max_af", 1.0001, False),
        ("max_af", float("nan"), False),
        ("max_af", float("inf"), False),
        ("max_af", 10**400, False),  # does not fit a float: a value error, not a crash
        ("min_gq", 10**400, False),
        ("min_gq", -(10**400), False),
        ("max_af", "0.1", False),
        ("strict_missing", True, True),
        ("strict_missing", 1, False),
        ("category", ["recessive", "comphet"], True),
        ("category", "recessive,comphet", True),
        ("category", ["all"], True),
        ("category", ["recesive"], False),
        ("category", [], False),
        ("category", 3, False),
    ],
)
def test_check_value(key, value, ok):
    checked, why = fr.check_value(fr.REGISTRY[key], value)
    assert (why is None) is ok, why


def test_check_value_names_the_nearest_category():
    _, why = fr.check_value(fr.REGISTRY["category"], ["recesive"])
    assert "did you mean 'recessive'" in why


def test_precedence_flag_over_preset_over_default():
    preset = {"min_gq": 30, "max_af": 0.001, "category": ["denovo"]}
    resolved, overridden = fr.resolve_filters_from_values(
        "trio", preset, explicit={"max_af": 0.005, "min_dp": 12}
    )
    got = {k: (r.value, r.source) for k, r in resolved.items()}
    assert got["min_gq"] == (30, "preset")
    assert got["max_af"] == (0.005, "flag")  # a flag wins, even over a preset
    assert got["min_dp"] == (12, "flag")  # no preset value: still from=flag
    assert got["min_ab"] == (0.25, "default")
    assert got["category"] == (["denovo"], "preset")
    assert list(resolved) == [k.file_key for k in fr.keys_for("trio")]
    assert overridden == [{"key": "max_af", "preset": 0.001, "flag": 0.005}]


def test_a_flag_equal_to_the_preset_value_is_still_from_the_flag():
    resolved, overridden = fr.resolve_filters_from_values(
        "trio", {"min_gq": 30}, explicit={"min_gq": 30}
    )
    assert resolved["min_gq"].source == "flag"
    assert overridden == [{"key": "min_gq", "preset": 30, "flag": 30}]


def test_ignoring_a_key_switches_it_off_and_says_so():
    resolved, overridden = fr.resolve_filters_from_values(
        "trio", {"gnomad_max_af": 0.0001, "min_gq": 30}, {}, ignore=["gnomad_max_af"]
    )
    assert resolved["gnomad_max_af"].source == "default"
    assert resolved["gnomad_max_af"].value is None
    assert resolved["min_gq"].source == "preset"
    assert overridden == [{"key": "gnomad_max_af", "preset": 0.0001, "flag": "off"}]


def test_a_band_that_is_empty_after_resolving_is_reported():
    assert fr.cross_key_problems({"min_ab": 0.8, "max_ab": 0.2})
    assert not fr.cross_key_problems({"min_ab": 0.2, "max_ab": 0.8})
    assert not fr.cross_key_problems({})
