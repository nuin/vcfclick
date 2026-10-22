"""The built-in presets are pinned.

They are starting points, not validated clinical standards, and one name must
never silently mean different filters: every fingerprint is written down here,
so an edit to a built-in fails this file until the author bumps `revision`,
adds a `history` entry and writes the CHANGELOG line. A step that registers a
new inheritance model fails `test_every_registered_model_is_run_or_excluded`
until it has decided what each built-in does about it. In-process.
"""

from __future__ import annotations

import pytest

from cli import preset_schema as ps
from cli.builtin_presets import BUILTINS, STARTING_POINT
from cli.presets import builtin_preset
from cli.trio_output import MODELS

# name -> (revision, fingerprint). Change a number here only with the revision.
PINNED = {
    "broad-first-pass": (1, "cd3726cc3853"),
    "strict-denovo": (1, "30ffc72499a2"),
    "rare-recessive": (1, "c19bb0a98ece"),
    "rare-dominant": (1, "466b94313b5c"),
}


def test_the_built_ins_are_these_four_in_this_order():
    assert list(BUILTINS) == list(PINNED)


@pytest.mark.parametrize("name", PINNED)
def test_fingerprint_and_revision_are_pinned(name):
    revision, fingerprint = PINNED[name]
    preset = builtin_preset(name)
    assert preset.problems == ()
    assert preset.revision == revision
    assert preset.fingerprint == fingerprint
    # the history ends at the current revision with the current fingerprint
    last = BUILTINS[name].history[-1]
    assert (last["revision"], last["fingerprint"]) == (revision, fingerprint)
    assert [h["revision"] for h in BUILTINS[name].history] == list(
        range(1, revision + 1)
    )


@pytest.mark.parametrize("name", PINNED)
def test_every_built_in_is_a_valid_starting_point_with_its_limits_in_words(name):
    b = BUILTINS[name]
    assert ps.name_problem(name) is None
    assert b.status == STARTING_POINT
    assert b.covers.strip() and b.does_not_cover.strip() and b.description.strip()
    # common to all of them until a consequence filter exists
    assert "no consequence or impact filter" in b.does_not_cover.lower()
    # no genes, regions, HPO terms or a knob that does nothing
    assert set(b.filters) <= {
        "category", "min_gq", "min_dp", "min_ab", "max_ab", "max_af",
        "gnomad_max_af", "strict_missing", "require_bound", "pass_only",
        "parent_alt_screen",
    }  # fmt: skip


@pytest.mark.parametrize("name", PINNED)
def test_the_limits_in_words_say_what_the_frequency_gates_cannot_judge(name):
    """The text must not claim more than the filter does: a variant that gnomAD has
    no record for is KEPT, and the INFO/AF gate steps aside for a cohort's own AF."""
    b = BUILTINS[name]
    if "gnomad_max_af" in b.filters:
        assert "absent from the loaded gnomad slice" in b.does_not_cover.lower()
        assert "is kept as rare" in b.does_not_cover
    if "max_af" in b.filters:
        assert "skipped when INFO/AF is the cohort's own frequency" in b.does_not_cover


def test_strict_denovo_says_where_its_gates_do_not_reach():
    """Each of these is a place where the code is looser than the covers text reads."""
    b = BUILTINS["strict-denovo"]
    for limit in (
        "homozygous-alt proband call is not allele-balance gated",
        "fewer than 2 alt reads is trusted whatever the fraction",
    ):
        assert limit in b.does_not_cover, limit
    assert "or fewer than 2 alt reads" in b.covers
    assert 'unset or "." counts as PASS' in b.covers
    assert "heterozygous proband call" in b.covers


def test_every_registered_model_is_run_or_excluded():
    """A built-in whose category list omits a registered model must say why."""
    for name, b in BUILTINS.items():
        run = set(builtin_preset(name).categories)
        for model in MODELS:
            if "all" in run or model in run:
                assert model not in b.excludes, (name, model)
            else:
                assert b.excludes.get(model), f"{name} omits {model} without a reason"
        assert set(b.excludes) <= set(MODELS), f"{name} excludes a model that is gone"


def test_strict_denovo_insists_and_the_others_only_warn():
    strict = BUILTINS["strict-denovo"].filters
    assert strict["strict_missing"] and strict["require_bound"]
    assert strict["pass_only"] and strict["parent_alt_screen"]
    for name in ("broad-first-pass", "rare-recessive", "rare-dominant"):
        f = BUILTINS[name].filters
        assert not f.get("require_bound") and not f.get("strict_missing")


def test_the_user_visible_numbers_are_the_documented_ones():
    f = {n: BUILTINS[n].filters for n in BUILTINS}
    assert (f["broad-first-pass"]["min_gq"], f["broad-first-pass"]["min_dp"]) == (15, 8)
    assert f["strict-denovo"]["min_ab"] == 0.30 and f["strict-denovo"]["max_ab"] == 0.70
    assert f["strict-denovo"]["max_af"] == 0.001
    assert f["rare-recessive"]["max_af"] == f["rare-recessive"]["gnomad_max_af"] == 0.01
    assert f["rare-dominant"]["category"] == ["dominant", "denovo"]
    assert f["rare-recessive"]["category"] == ["recessive", "comphet"]


def test_a_built_in_title_names_the_parent_phenotype_limit():
    assert "parent phenotype not considered" in BUILTINS["rare-dominant"].title
