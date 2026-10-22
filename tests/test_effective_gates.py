"""The effective-gate report, without a database.

`build_report` is pure: every threshold (a gate binds at 90% of its rows, the
median-depth warning at 1.5 x min_dp) and every warning code is reached here
with synthetic counts. The DB-backed behaviour is in tests/test_presets_db.py.
"""

from __future__ import annotations

import pytest

from cli.effective_gates import (
    ENFORCEABLE_MIN,
    THIN_COVERAGE_FACTOR,
    SampleGateStats,
    VariantStats,
    build_report,
    median_from_histogram,
    refusal_message,
)
from cli.trio_analysis import Gates
from cli.warning_codes import CODES, warning

GATES = Gates(20, 10, 0.01, 0.25, 0.75)


def sample(role, rows=100, gq=None, dp=None, ad=None, depths=((40, 100),)):
    return SampleGateStats(
        role, role, rows,
        rows if gq is None else gq, rows if dp is None else dp,
        rows if ad is None else ad, depths,
    )  # fmt: skip


def trio(**by_role):
    return [by_role.get(r, sample(r)) for r in ("proband", "father", "mother")]


def report(stats=None, gates=GATES, **kw):
    kw.setdefault("variants", VariantStats(100, 100))
    kw.setdefault("gnomad_max_af", None)
    kw.setdefault("gnomad_loaded", False)
    return build_report(stats=stats or trio(), gates=gates, **kw)


def by_code(rep, code):
    return [w for w in rep.warnings if w["code"] == code]


# ───────────────────────────── a gate binds at 90% ─────────────────────────────


def test_a_gate_with_its_field_everywhere_binds_and_says_nothing():
    rep = report()
    assert rep.warnings == () and rep.errors == ()
    assert rep.bound == {"min_gq": True, "min_dp": True, "min_ab": True, "max_ab": True,
                         "max_af": True}  # fmt: skip
    assert rep.document["population_filter"] is True


def test_the_binding_threshold_is_ninety_percent_inclusive():
    assert ENFORCEABLE_MIN == 0.90
    at = report(trio(father=sample("father", rows=1000, gq=900)))
    below = report(trio(father=sample("father", rows=1000, gq=899)))
    assert at.bound["min_gq"] is True and not by_code(at, "gate_unbound")
    assert below.bound["min_gq"] is False
    (w,) = by_code(below, "gate_unbound")
    assert w["key"] == "min_gq" and w["level"] == "warning"
    assert (
        "min_gq=20 cannot be enforced: only 89.9% of father rows carry GQ"
        in w["message"]
    )


def test_a_gate_that_binds_on_some_rows_says_how_many_and_what_a_missing_one_does():
    stats = trio(father=sample("father", rows=1000, dp=980))
    (w,) = by_code(report(stats), "gate_partial")
    assert w["level"] == "info"
    assert "min_dp=10 enforced on 98% of rows" in w["message"]
    assert "counts as a pass" in w["message"]
    strict = report(stats, GATES._replace(strict_missing=True))
    assert "counts as a fail" in by_code(strict, "gate_partial")[0]["message"]
    # 99.96% is not shown as 100%
    nearly = trio(father=sample("father", rows=10000, dp=9996))
    assert "99.9%" in by_code(report(nearly), "gate_partial")[0]["message"]


def test_a_sample_with_no_stored_rows_says_nothing_but_no_rows_at_all_is_unbound():
    quiet = report(trio(father=sample("father", rows=0)))
    assert quiet.bound["min_gq"] is True
    none = report([sample(r, rows=0) for r in ("proband", "father", "mother")])
    assert not none.bound["min_gq"]
    assert "no stored genotype rows" in by_code(none, "gate_unbound")[0]["message"]


def test_allele_balance_reads_the_proband_unless_strict_then_the_parents_too():
    parents_without_ad = trio(
        father=sample("father", ad=0), mother=sample("mother", ad=0)
    )
    assert report(parents_without_ad).bound["min_ab"] is True
    strict = report(parents_without_ad, GATES._replace(strict_missing=True))
    assert strict.bound["min_ab"] is False
    off = report(trio(proband=sample("proband", ad=0)))
    assert off.bound["min_ab"] is False
    # a band that excludes nothing is not a gate
    none = report(
        trio(proband=sample("proband", ad=0)), GATES._replace(min_ab=0.0, max_ab=1.0)
    )
    assert "min_ab" not in none.bound and "max_ab" not in none.bound


def test_a_zero_threshold_is_not_a_requested_gate():
    rep = report(trio(proband=sample("proband", gq=0)), GATES._replace(min_gq=0))
    assert "min_gq" not in rep.bound and not by_code(rep, "gate_unbound")


def test_the_parent_alt_screen_needs_the_parents_allele_depths():
    gates = GATES._replace(parent_alt_screen=True)
    assert report(gates=gates).bound["parent_alt_screen"] is True
    rep = report(trio(father=sample("father", ad=0)), gates)
    assert rep.bound["parent_alt_screen"] is False
    (w,) = [w for w in by_code(rep, "gate_unbound") if w["key"] == "parent_alt_screen"]
    assert "every parent 0/0 would be trusted" in w["message"]


# ───────────────────────────── the frequency filters ─────────────────────────────


def test_info_af_must_be_on_most_variants_and_not_the_cohorts_own():
    assert report(variants=VariantStats(100, 90)).bound["max_af"] is True
    rep = report(variants=VariantStats(100, 89))
    assert rep.bound["max_af"] is False
    assert (
        "only 89% of variants carry INFO/AF"
        in by_code(rep, "gate_unbound")[0]["message"]
    )
    skipped = report(gates=GATES._replace(af_gate="skipped_cohort_derived"))
    assert skipped.bound["max_af"] is False
    assert not [w for w in by_code(skipped, "gate_unbound") if w["key"] == "max_af"]


def test_no_population_filter_is_said_with_the_reason():
    cohort = report(gates=GATES._replace(af_gate="skipped_cohort_derived"))
    (w,) = by_code(cohort, "no_population_frequency_filter")
    assert w["message"].startswith("NO population frequency filter applied")
    assert (
        "cohort's own frequency" in w["message"]
        and "no gnomAD data is loaded" in w["message"]
    )
    loaded = report(
        gates=GATES._replace(af_gate="skipped_cohort_derived"), gnomad_loaded=True
    )
    assert (
        "--gnomad-max-af is not set"
        in by_code(loaded, "no_population_frequency_filter")[0]["message"]
    )
    absent = report(variants=VariantStats(100, 0))
    assert (
        "missing from most variants"
        in by_code(absent, "no_population_frequency_filter")[0]["message"]
    )
    # gnomAD loaded and asked for: that IS a population filter
    ok = report(
        gates=GATES._replace(af_gate="skipped_cohort_derived"),
        gnomad_max_af=0.01,
        gnomad_loaded=True,
    )
    assert not by_code(ok, "no_population_frequency_filter")
    assert (
        ok.document["population_filter"] is True and ok.bound["gnomad_max_af"] is True
    )


def test_gnomad_asked_for_and_not_loaded():
    rep = report(gnomad_max_af=0.01)
    (w,) = by_code(rep, "gnomad_not_loaded")
    assert "gnomad_max_af=0.01 cannot be enforced" in w["message"]
    assert rep.bound["gnomad_max_af"] is False
    assert rep.document["gnomad"] == {
        "loaded": False, "bound": False, "field": "popmax (grpmax else AF)",
        "version": None, "candidates_without_record": None,
    }  # fmt: skip


# ───────────────────────────── the other warnings ─────────────────────────────


def test_thin_coverage_is_one_and_a_half_times_min_dp():
    assert THIN_COVERAGE_FACTOR == 1.5
    thin = trio(proband=sample("proband", depths=((31, 100),)))
    (w,) = by_code(
        report(thin, GATES._replace(min_dp=25), assay="exome"), "thin_coverage"
    )
    assert "median proband DP is 31 but min_dp is 25" in w["message"]
    assert "too thin for the preset" in w["message"] and "exome" in w["message"]
    just_enough = trio(proband=sample("proband", depths=((38, 100),)))
    assert not by_code(report(just_enough, GATES._replace(min_dp=25)), "thin_coverage")
    plain = report(thin, GATES._replace(min_dp=25), subject="gates")
    assert "too thin for the gates" in by_code(plain, "thin_coverage")[0]["message"]


def test_skipped_keys_that_widen_the_search_are_warnings_others_are_not():
    skipped = [
        {
            "key": "pass_only",
            "reason": "not_applicable",
            "widening": True,
            "message": "m1",
        },
        {"key": "x", "reason": "not_applicable", "widening": False, "message": "m2"},
    ]
    rep = report(skipped=skipped)
    assert [w["message"] for w in by_code(rep, "skipped_widening")] == ["m1"]


def test_a_duo_asking_for_de_novo_is_told_it_is_weaker():
    assert not by_code(report(kind="trio"), "duo_denovo")
    assert by_code(report(kind="duo", categories=["denovo"]), "duo_denovo")
    assert by_code(report(kind="duo", categories=["all"]), "duo_denovo")
    assert not by_code(report(kind="duo", categories=["recessive"]), "duo_denovo")


def test_a_presets_max_af_that_was_skipped_is_explained_as_info():
    rep = report(
        gates=GATES._replace(af_gate="skipped_cohort_derived"), preset_max_af=0.001
    )
    (w,) = by_code(rep, "preset_max_af_not_applied")
    assert (
        w["level"] == "info"
        and "max_af 0.001 from the preset was not applied" in w["message"]
    )
    assert not by_code(report(preset_max_af=0.001), "preset_max_af_not_applied")


# ───────────────────────────── refusing ─────────────────────────────


def broken():
    return trio(
        proband=sample("proband", gq=0, dp=0, ad=0),
        father=sample("father", gq=0, dp=0, ad=0),
        mother=sample("mother", gq=0, dp=0, ad=0),
    )


def test_require_bound_makes_unbound_gates_and_a_missing_population_filter_errors():
    gates = GATES._replace(af_gate="skipped_cohort_derived")
    rep = report(broken(), gates, enforce_bound=True)
    assert {e["code"] for e in rep.errors} == {
        "gate_unbound",
        "no_population_frequency_filter",
    }
    assert {e["key"] for e in rep.errors if e["code"] == "gate_unbound"} == {
        "min_gq", "min_dp", "min_ab",
    }  # fmt: skip
    assert report(broken(), gates).errors == ()  # warnings only, by default
    gn = report(gnomad_max_af=0.01, enforce_bound=True)
    assert [e["code"] for e in gn.errors] == ["gnomad_not_loaded"]


def test_strict_missing_alone_refuses_only_for_the_gq_dp_ad_gates():
    gates = GATES._replace(af_gate="skipped_cohort_derived", strict_missing=True)
    rep = report(broken(), gates, enforce_format=True)
    assert {e["key"] for e in rep.errors} == {"min_gq", "min_dp", "min_ab"}
    assert "no_population_frequency_filter" not in {e["code"] for e in rep.errors}


def test_the_refusal_names_each_problem_once_and_says_how_to_run_anyway():
    rep = report(
        broken(), GATES._replace(af_gate="skipped_cohort_derived"),
        gnomad_max_af=0.01, enforce_bound=True,
    )  # fmt: skip
    text = refusal_message(rep.errors, subject="preset 'strict-denovo'")
    assert text.startswith("preset 'strict-denovo' needs every gate to bind")
    assert text.count("Remedy:") == 2  # the population remedy is said once, not twice
    assert "annotations load-gnomad" in text and "broad-first-pass" in text
    assert "--no-require-bound-gates" in text


def test_warning_codes_are_a_closed_list():
    assert set(CODES) >= {
        "gate_unbound", "gate_partial", "thin_coverage", "gnomad_not_loaded",
        "no_population_frequency_filter", "skipped_widening", "duo_denovo",
        "preset_max_af_not_applied",
    }  # fmt: skip
    assert warning("gate_unbound", "min_gq", "m") == {
        "code": "gate_unbound", "level": "warning", "key": "min_gq", "message": "m",
    }  # fmt: skip
    with pytest.raises(ValueError):
        warning("made_up", None, "m")


def test_the_document_has_exactly_the_documented_keys():
    doc = report().document
    assert list(doc) == [
        "rows", "gq", "dp", "ad", "median_dp", "info_af", "filter_column",
        "gnomad", "population_filter", "warnings",
    ]  # fmt: skip
    assert list(doc["gq"]) == ["proband", "father", "mother", "bound"]
    assert doc["rows"] == {"proband": 100, "father": 100, "mother": 100}
    assert doc["median_dp"] == {"proband": 40, "father": 40, "mother": 40}
    assert doc["info_af"] == {"state": "applied", "bound": True}


def test_median_from_a_histogram():
    assert median_from_histogram([]) is None
    assert median_from_histogram([(30, 5)]) == 30
    assert median_from_histogram([(10, 1), (20, 1), (30, 1)]) == 20
    assert median_from_histogram([(10, 1), (20, 1)]) == 15  # even: the two middles
    assert median_from_histogram([(40, 3), (10, 1), (50, 1)]) == 40
    assert median_from_histogram([(1, 1000), (100, 1)]) == 1
