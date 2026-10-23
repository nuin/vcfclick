"""The effective-gate report, without a database.

`build_report` is pure: every threshold (a gate binds at 90% of its rows, the
median-depth warning at 1.5 x min_dp) and every warning code is reached here
with synthetic counts. The DB-backed behaviour is in tests/test_presets_db.py.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from cli import trio_gate_output as tgo
from cli.effective_gates import (
    ENFORCEABLE_MIN,
    GNOMAD_MIN_HIT_RATE,
    GNOMAD_PROBE_SIZE,
    STRICT_KEYS,
    THIN_COVERAGE_FACTOR,
    CandidateStats,
    GenotypeClassStats,
    GnomadProbe,
    SampleGateStats,
    VariantStats,
    build_report,
    candidate_models,
    candidate_presence_sql,
    collect_candidate_stats,
    gnomad_record_warning,
    judged_gate_keys,
    median_from_histogram,
    refusal_message,
    sample_cohort_keys,
    unbound_strict_keys,
)
from cli.trio_analysis import MODEL_PREDICATES, Gates, Trio, trio_sql
from cli.warning_codes import CODES, warning

GATES = Gates(20, 10, 0.01, 0.25, 0.75)


def sample(role, rows=100, gq=None, dp=None, ad=None, depths=((40, 100),)):
    return SampleGateStats(
        role, role, rows,
        rows if gq is None else gq, rows if dp is None else dp,
        rows if ad is None else ad, depths,
    )  # fmt: skip


def split(role, **classes):
    """A sample whose stored rows are split by genotype class: `ref=(rows, gq, dp,
    ad)` is the number of its 0/0 rows and how many of them carry each field (`het`
    is 0/1, `alt` is 1/1). The whole-sample totals are the sums."""
    parts = tuple(GenotypeClassStats(c, *counts) for c, counts in classes.items())
    return SampleGateStats(
        role, role, *(sum(p[i] for p in classes.values()) for i in range(4)),
        ((40, sum(p[0] for p in classes.values())),), parts,
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


def test_pass_only_binds_only_where_the_data_has_a_filter_label_to_act_on():
    """FILTER is NULL for both PASS and '.', so on a VCF that labels nothing the
    gate removes nothing: it must not read as a gate that worked."""
    gates = GATES._replace(pass_only=True)
    labelled = report(gates=gates, variants=VariantStats(100, 100, non_pass=3))
    assert labelled.bound["pass_only"] is True and labelled.warnings == ()
    bare = report(
        gates=gates, variants=VariantStats(100, 100, non_pass=0), enforce_bound=True
    )
    assert bare.bound["pass_only"] is False
    (w,) = by_code(bare, "pass_only_no_labels")
    assert (w["level"], w["key"]) == ("warning", "pass_only")
    assert w["message"] == (
        "pass_only removes nothing: no variant in this database carries a FILTER "
        "label other than PASS or unset"
    )
    # a warning and never a refusal: a VCF that was hard-filtered before it was
    # ingested is all PASS, and is not a reason to stop a strict preset
    assert bare.errors == ()
    # not counted (a stand-in, or an older caller): nothing is claimed either way
    unknown = report(gates=gates, variants=VariantStats(100, 100))
    assert unknown.bound["pass_only"] is True and unknown.warnings == ()
    # not asked for: not judged
    off = report(variants=VariantStats(100, 100, non_pass=0))
    assert "pass_only" not in off.bound and not by_code(off, "pass_only_no_labels")
    # no filter column at all: the gate cannot be enforced, and a run that insists refuses
    nocol = report(gates=gates, filter_present=False, enforce_bound=True)
    assert nocol.bound["pass_only"] is False
    assert [e["key"] for e in nocol.errors] == ["pass_only"]
    assert nocol.errors[0]["code"] == "gate_unbound"


def test_the_filter_column_object_counts_the_labelled_variants():
    doc = report(variants=VariantStats(100, 100, non_pass=4)).document["filter_column"]
    assert doc == {"present": True, "non_pass_variants": 4}
    assert report().document["filter_column"] == {
        "present": True, "non_pass_variants": None,
    }  # fmt: skip


def test_a_report_limited_to_pass_only_says_only_that():
    gates = GATES._replace(af_gate="skipped_cohort_derived", pass_only=True)
    scoped = report(
        broken(), gates, variants=VariantStats(100, 100, non_pass=0), only={"pass_only"}
    )
    assert [w["code"] for w in scoped.warnings] == ["pass_only_no_labels"]
    assert scoped.bound == {"pass_only": False}
    assert scoped.document["gq"]["bound"] is None


# ─────────────────── a gate is judged on the rows it reads ───────────────────


def mismatched(field):
    """A trio whose parents' 0/0 rows lack `field` (10 rows) while their 0/1 rows
    carry it (200 rows): over all rows the share is 95.2%, over the 0/0 rows a de
    novo gate reads it is 0%. The proband's rows are complete."""
    have = {"gq": 0, "dp": 1, "ad": 2}[field]
    ref = [10, 10, 10, 10]
    ref[1 + have] = 0
    parent = dict(ref=tuple(ref), het=(200, 200, 200, 200))
    return [
        split("proband", het=(10, 10, 10, 10), alt=(5, 5, 5, 5), ref=(200, 200, 200, 200)),
        split("father", **parent),
        split("mother", **parent),
    ]  # fmt: skip


def test_the_whole_sample_share_hides_a_field_missing_where_the_gate_reads():
    """Over every stored row of the parents the share is 95.2% and the old report
    called the gate bound; at the sites a de novo call is made the parents are 0/0
    and none of those rows carries the field."""
    stats = mismatched("ad")
    assert stats[1].fraction("ad") > ENFORCEABLE_MIN  # the number that misled
    gates = GATES._replace(parent_alt_screen=True)
    rep = report(stats, gates, categories=["denovo"], enforce_bound=True)
    assert rep.bound["parent_alt_screen"] is False
    (w,) = [w for w in by_code(rep, "gate_unbound") if w["key"] == "parent_alt_screen"]
    assert "only 0% of father/mother 0/0 rows carry AD" in w["message"]
    assert [e["key"] for e in rep.errors] == ["parent_alt_screen"]
    # the other models never read it, and the document keeps the whole-sample facts
    assert report(stats, gates, categories=["recessive"]).bound["parent_alt_screen"]
    assert report(stats, gates, categories=["all"]).bound["parent_alt_screen"] is False
    # AD is read by two gates: the document's ad.bound is bound only if both are
    assert rep.bound["min_ab"] is True
    assert rep.document["ad"]["bound"] is False
    assert report(stats, gates, categories=["recessive"]).document["ad"]["bound"]
    assert rep.document["ad"]["father"] == round(200 / 210, 4)
    # allele balance (strict: the parents' het calls too) is judged on het rows,
    # which have their depths
    strict = GATES._replace(strict_missing=True)
    denovo = report(stats, strict, categories=["denovo"])
    assert denovo.bound["min_ab"] is True
    # a de novo candidate's parents are 0/0: their het rows are not asked for, so
    # the 95.2% over all their rows is not quoted as the gate's share either
    assert not [w for w in denovo.warnings if w["key"] == "min_ab"]
    assert report(stats, strict, categories=["recessive"]).bound["min_ab"] is True


def test_gq_and_dp_are_judged_on_the_class_each_model_reads():
    for field, key, label in (("gq", "min_gq", "GQ"), ("dp", "min_dp", "DP")):
        stats = mismatched(field)
        denovo = report(stats, categories=["denovo"])
        assert denovo.bound[key] is False, key
        (w,) = [w for w in by_code(denovo, "gate_unbound") if w["key"] == key]
        assert f"only 0% of father/mother 0/0 rows carry {label}" in w["message"]
        # recessive reads the parents' 0/1 rows, dominant both: the 0/1 rows are
        # complete, and dominant's share is over the pooled classes
        assert report(stats, categories=["recessive"]).bound[key] is True
        assert report(stats, categories=["dominant"]).bound[key] is True
        # a run of all models is only as bound as its weakest
        assert report(stats, categories=["all"]).bound[key] is False
        assert report(stats, categories=["recessive", "denovo"]).bound[key] is False
        assert report(stats, categories=["recessive", "dominant"]).bound[key] is True


def test_the_proband_is_judged_on_its_non_reference_rows():
    """A caller whose hom-ref calls carry AD and whose het calls do not: the
    allele-balance gate reads the het rows only."""
    stats = [
        split("proband", ref=(200, 200, 200, 200), het=(10, 10, 10, 0)),
        split("father", ref=(10, 10, 10, 10)),
        split("mother", ref=(10, 10, 10, 10)),
    ]
    assert stats[0].fraction("ad") > ENFORCEABLE_MIN
    rep = report(stats, categories=["denovo"])
    assert rep.bound["min_ab"] is False
    (w,) = [w for w in by_code(rep, "gate_unbound") if w["key"] == "min_ab"]
    assert "only 0% of proband 0/1 rows carry AD" in w["message"]
    # a hom-alt proband is not allele-balance gated, so recessive reads no proband row
    # and the gate is judged as it always was, on the whole sample
    assert report(stats, categories=["recessive"]).bound["min_ab"] is True


def test_a_class_without_rows_falls_back_to_the_whole_sample():
    """No stored 0/0 rows (a database ingested without --keep-reference, or a
    parent with only het calls): nothing to read, so the whole sample decides."""
    only_het = dict(het=(100, 100, 100, 100))
    stats = [split("proband", het=(100, 100, 100, 100)), split("father", **only_het),
             split("mother", **only_het)]  # fmt: skip
    assert report(stats, categories=["denovo"]).bound["min_gq"] is True
    poor = dict(het=(100, 50, 100, 100))
    low = [split("proband", het=(100, 100, 100, 100)), split("father", **poor),
           split("mother", **poor)]  # fmt: skip
    rep = report(low, categories=["denovo"])
    assert rep.bound["min_gq"] is False
    (w,) = by_code(rep, "gate_unbound")
    assert "only 50% of father/mother rows carry GQ" in w["message"]
    # a sample with no class split at all (a stand-in) is judged whole, as before
    assert (
        report(trio(father=sample("father", rows=1000, gq=899))).bound["min_gq"]
        is False
    )


def test_only_the_models_the_run_selected_decide():
    stats = mismatched("gq")
    assert report(stats, categories=["recessive"]).warnings == ()
    assert report(stats, categories=["denovo"]).bound["min_gq"] is False
    # a model name this version does not know reads everything, never nothing
    assert report(stats, categories=["xlinked"]).bound["min_gq"] is True


# ─────────────────── a gate is judged on the candidate rows ───────────────────


def at_sites(model, sites, **roles):
    """What the rows at a model's candidate sites carry: `proband=dict(het=(rows,
    gq, dp, ad))` is the number of the proband's 0/1 rows at those sites and how
    many carry each field."""
    return CandidateStats(
        model,
        sites,
        {
            role: tuple(GenotypeClassStats(c, *n) for c, n in classes.items())
            for role, classes in roles.items()
        },
    )


def shaped(field, *lacking):
    """A trio whose `lacking` roles carry `field` on the inherited calls (95.2%
    of the class the gate reads) and on none of their rows at the de novo
    candidate sites. The proband has 10 de novo hets and 200 inherited hom-alts;
    each parent has 10 0/0 rows at the de novo sites and 200 0/0 rows beside
    them. Returns the whole-sample stats and the candidate stats of `denovo`."""
    have = {"gq": 0, "dp": 1, "ad": 2}[field]

    def counts(role, rows, inherited):
        """(rows, gq, dp, ad) of `rows` candidate rows and `inherited` others."""
        out = [rows + inherited] * 4
        if role in lacking:
            out[1 + have] = inherited
        return tuple(out)

    def at_denovo(role, rows):
        out = [rows] * 4
        if role in lacking:
            out[1 + have] = 0
        return tuple(out)

    whole = [
        split("proband", het=at_denovo("proband", 10), alt=(200, 200, 200, 200)),
        split("father", ref=counts("father", 10, 200)),
        split("mother", ref=counts("mother", 10, 200)),
    ]  # fmt: skip
    denovo = at_sites(
        "denovo", 10, proband=dict(het=at_denovo("proband", 10)),
        father=dict(ref=at_denovo("father", 10)),
        mother=dict(ref=at_denovo("mother", 10)),
    )  # fmt: skip
    return whole, {"denovo": denovo}


def test_a_field_present_on_the_class_but_absent_on_the_candidates_does_not_bind():
    """The reproduced case: 10 de novo hets without GQ and 200 inherited hom-alts
    with it. The proband's non-reference rows carry GQ on 95.2% of them, the old
    verdict ("bound"), and a strict_missing run then failed all 10 candidates and
    printed `denovo candidates: 0`."""
    whole, candidates = shaped("gq", "proband")
    assert whole[0].share("gq", ("het", "alt")) == 200 / 210  # the number that misled
    gates = GATES._replace(strict_missing=True)
    old = report(whole, gates, categories=["denovo"], enforce_format=True)
    assert old.bound["min_gq"] is True and old.errors == ()  # without the candidates
    rep = report(
        whole, gates, categories=["denovo"], candidates=candidates, enforce_format=True
    )
    assert rep.bound["min_gq"] is False
    (w,) = [w for w in by_code(rep, "gate_unbound") if w["key"] == "min_gq"]
    assert (
        "min_gq=20 cannot be enforced: only 0% of proband 0/1 or 1/1 rows at denovo "
        "candidate sites carry GQ" in w["message"]
    )
    assert "strict_missing counts a call without GQ as a fail" in w["message"]
    assert [e["key"] for e in rep.errors] == ["min_gq"]
    # the remedy names the field
    text = refusal_message(rep.errors, subject="preset 'strict-denovo'")
    assert "Remedy: GQ is missing where the gates read it" in text
    # the document keeps the whole-sample facts and the verdict of the candidates
    assert rep.document["gq"]["proband"] == round(200 / 210, 4)
    assert rep.document["gq"]["bound"] is False
    assert rep.document["dp"]["bound"] is True
    # a warn-only run says it too, and does not refuse
    warn = report(whole, GATES, categories=["denovo"], candidates=candidates)
    assert by_code(warn, "gate_unbound")[0]["level"] == "warning"
    assert warn.errors == ()


def test_the_parents_reference_rows_are_judged_at_the_candidate_sites_too():
    whole, candidates = shaped("ad", "father", "mother")
    assert whole[1].share("ad", ("ref",)) == 200 / 210
    gates = GATES._replace(parent_alt_screen=True)
    rep = report(whole, gates, categories=["denovo"], candidates=candidates)
    assert rep.bound["parent_alt_screen"] is False
    (w,) = [w for w in by_code(rep, "gate_unbound") if w["key"] == "parent_alt_screen"]
    assert (
        "only 0% of father/mother 0/0 rows at denovo candidate sites carry AD"
        in (w["message"])
    )
    assert rep.bound["min_ab"] is True  # the proband's hets at the sites have it
    # without the candidates it reads as bound
    assert report(whole, gates, categories=["denovo"]).bound["parent_alt_screen"]


def test_the_weakest_model_decides_and_each_is_judged_on_its_own_candidates():
    whole, candidates = shaped("gq", "proband")
    candidates["recessive"] = at_sites(
        "recessive", 200, proband=dict(alt=(200, 200, 200, 200)),
        father=dict(het=(200, 200, 200, 200)), mother=dict(het=(200, 200, 200, 200)),
    )  # fmt: skip
    assert not report(whole, categories=["recessive"], candidates=candidates).warnings
    both = report(whole, categories=["recessive", "denovo"], candidates=candidates)
    assert both.bound["min_gq"] is False
    (w,) = by_code(both, "gate_unbound")
    assert "rows at denovo candidate sites" in w["message"]
    assert "recessive" not in w["message"]


def test_only_the_classes_the_gate_reads_count_at_the_candidate_sites():
    """Allele balance reads the proband's het rows: a hom-alt candidate row without
    AD is not gated, so it does not make the gate unbound, and the het ones do."""
    candidates = {
        "denovo": at_sites(
            "denovo",
            15,
            proband=dict(het=(10, 10, 10, 10), alt=(5, 5, 5, 0)),
            father=dict(ref=(15, 15, 15, 15)),
            mother=dict(ref=(15, 15, 15, 15)),
        )  # fmt: skip
    }
    whole = [split("proband", het=(10, 10, 10, 10), alt=(5, 5, 5, 0)),
             split("father", ref=(15, 15, 15, 15)),
             split("mother", ref=(15, 15, 15, 15))]  # fmt: skip
    assert report(whole, categories=["denovo"], candidates=candidates).bound["min_ab"]
    lacking = {
        "denovo": at_sites(
            "denovo",
            15,
            proband=dict(het=(10, 10, 10, 4), alt=(5, 5, 5, 5)),
            father=dict(ref=(15, 15, 15, 15)),
            mother=dict(ref=(15, 15, 15, 15)),
        )  # fmt: skip
    }
    rep = report(whole, categories=["denovo"], candidates=lacking)
    assert rep.bound["min_ab"] is False
    assert (
        "only 40% of proband 0/1 rows at denovo candidate sites carry AD"
        in (by_code(rep, "gate_unbound")[0]["message"])
    )


def test_a_model_without_candidate_rows_is_judged_on_the_class_share():
    """No candidate site (the model finds nothing, or nothing was counted for it):
    the stored rows of the class are what there is to go by, as before."""
    whole, _ = shaped("gq", "proband")
    none = {"denovo": at_sites("denovo", 0)}
    for candidates in (none, {}, None):
        rep = report(whole, categories=["denovo"], candidates=candidates)
        assert rep.bound["min_gq"] is True, candidates  # the whole-sample class share
    poor = [split("proband", het=(100, 100, 100, 100)),
            split("father", ref=(100, 50, 100, 100)),
            split("mother", ref=(100, 50, 100, 100))]  # fmt: skip
    rep = report(poor, categories=["denovo"], candidates=none)
    assert (
        "only 50% of father/mother 0/0 rows carry GQ"
        in (by_code(rep, "gate_unbound")[0]["message"])
    )
    # sites exist but none of the classes the gate reads: the same fallback
    elsewhere = {"denovo": at_sites("denovo", 3, proband=dict(alt=(3, 0, 0, 0)))}
    assert report(poor, categories=["denovo"], candidates=elsewhere).bound["min_ab"]


def test_a_partial_gate_says_where_its_share_was_taken():
    candidates = {
        "denovo": at_sites(
            "denovo",
            100,
            proband=dict(het=(100, 100, 98, 100)),
            father=dict(ref=(100, 100, 100, 100)),
            mother=dict(ref=(100, 100, 100, 100)),
        )  # fmt: skip
    }
    whole = [split("proband", het=(100, 100, 98, 100)),
             split("father", ref=(100, 100, 100, 100)),
             split("mother", ref=(100, 100, 100, 100))]  # fmt: skip
    rep = report(whole, categories=["denovo"], candidates=candidates)
    (w,) = by_code(rep, "gate_partial")
    assert "min_dp=10 enforced on 98% of rows at denovo candidate sites" in w["message"]


def test_candidate_rows_are_tallied_from_the_grouped_count():
    """[proband gt, father gt, mother gt, rows, then gq, dp, ad for each of the
    three samples]: two groups of the de novo join."""
    from cli.effective_gates import _candidate_stats

    rows = [
        [1, 0, 0, 7, 7, 6, 5, 7, 7, 7, 7, 7, 0],
        [2, 0, 0, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3],
    ]
    c = _candidate_stats("denovo", rows)
    assert c.sites == 10 and c.model == "denovo"
    assert c.share("proband", "gq", ("het",)) == 1.0
    assert c.share("proband", "gq", ("alt",)) == 1 / 3
    assert c.share("proband", "gq", ("het", "alt")) == 8 / 10
    assert c.share("proband", "dp", ("het", "alt")) == 9 / 10
    assert c.share("mother", "ad", ("ref",)) == 3 / 10
    assert c.share("father", "ad", ("ref",)) == 1.0
    assert c.share("father", "gq", ("het",)) is None  # no such row
    assert c.share("nobody", "gq", ("ref",)) is None


class FakeSession:
    def __init__(self):
        self.sql = []

    def query(self, sql, fmt):
        self.sql.append(sql)
        return self

    def bytes(self):
        return b'{"data": [[1, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]]}'


def test_one_count_per_distinct_predicate_and_it_reads_the_models_own_join():
    trio_ = Trio("i1", "CHILD", "FATHER", "MOTHER")
    sess = FakeSession()
    got = collect_candidate_stats(
        sess, trio_, ["denovo", "dominant", "comphet", "xlinked"]
    )
    assert set(got) == {"denovo", "dominant", "comphet"}  # xlinked has no predicate
    assert len(sess.sql) == 2  # comphet shares dominant's candidates
    assert got["comphet"].model == "comphet" and got["comphet"].sites == 4
    sql = candidate_presence_sql("denovo", trio_)
    assert MODEL_PREDICATES["denovo"] in sql and "GROUP BY g.gt, f.gt, m.gt" in sql
    assert "ingest_id = 'i1'" in sql and "g.sample_id = 'CHILD'" in sql
    # the very predicate and parent joins `db trio` runs, and no quality gate
    model_sql = trio_sql("denovo", trio_, Gates(0, 0, 1.0, 0.0, 1.0), count_only=True)
    assert MODEL_PREDICATES["denovo"] in model_sql
    join = "INNER JOIN genotypes f ON"
    assert join in sql and join in model_sql
    assert "gq >=" not in sql and "info_AF" not in sql
    # a trio without both parents has no candidate sites to count
    assert (
        collect_candidate_stats(FakeSession(), Trio("i1", "C", "F", None), ["denovo"])
        == {}
    )


def test_only_the_models_that_read_a_judged_gate_have_their_candidates_counted():
    only_gq = GATES._replace(min_dp=0, min_ab=0.0, max_ab=1.0)
    assert judged_gate_keys(only_gq) == ("min_gq",)
    assert judged_gate_keys(GATES) == ("min_gq", "min_dp", "min_ab")
    screen = GATES._replace(parent_alt_screen=True)
    assert judged_gate_keys(screen, only={"parent_alt_screen"}) == (
        "parent_alt_screen",
    )
    assert candidate_models(screen, ["all"], {"parent_alt_screen"}) == ("denovo",)
    assert candidate_models(only_gq, ["recessive"]) == ("recessive",)
    assert candidate_models(GATES, ["all"]) == (
        "denovo",
        "recessive",
        "dominant",
        "comphet",
    )
    nothing = GATES._replace(min_gq=0, min_dp=0, min_ab=0.0, max_ab=1.0)
    assert candidate_models(nothing, ["all"]) == ()


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
        "version": None, "candidates_checked": None,
        "candidates_without_record": None, "cohort_variants_sampled": None,
        "cohort_variants_with_record": None,
    }  # fmt: skip


def test_gnomad_bound_is_reported_only_for_a_gate_that_was_asked_for():
    """`bound` says the filter was asked for and probed. A gnomAD that is loaded
    but not asked for (the preset's key switched off, or no `--gnomad-max-af`) was
    never looked at, so it is neither bound nor unbound."""
    for loaded in (True, False):
        doc = report(gnomad_loaded=loaded).document["gnomad"]
        assert doc["loaded"] is loaded and doc["bound"] is None, loaded
        assert "gnomad_max_af" not in report(gnomad_loaded=loaded).bound
    asked = report(gnomad_max_af=0.01, gnomad_loaded=True).document["gnomad"]
    assert asked["bound"] is True
    assert report(gnomad_max_af=0.01).document["gnomad"]["bound"] is False


def cohort(**per_chrom):
    """VariantStats for a cohort with this many variants on each chromosome."""
    chroms = tuple(sorted(per_chrom.items(), key=lambda t: (-t[1], t[0])))
    n = sum(per_chrom.values())
    return VariantStats(n, n, chroms)


def gnomad_report(variants, chroms, **kw):
    return report(
        variants=variants,
        gates=GATES._replace(af_gate="skipped_cohort_derived"),
        gnomad_max_af=0.001,
        gnomad_loaded=True,
        gnomad_chroms=chroms,
        **kw,
    )


def test_gnomad_binds_where_it_has_records_for_ninety_percent_of_the_cohort():
    v = cohort(chr1=90, chr2=10)
    both = gnomad_report(v, {"chr1", "chr2"})
    at = gnomad_report(v, {"chr1"})  # 90%: inclusive
    for rep in (both, at):
        assert rep.bound["gnomad_max_af"] is True and rep.warnings == ()
        assert rep.document["gnomad"]["bound"] is True
    # not looked at (None): any loaded gnomAD is taken to cover the cohort
    assert gnomad_report(v, None).bound["gnomad_max_af"] is True
    # no variants to judge by: not a reason to refuse
    assert gnomad_report(VariantStats(0, 0), {"chr9"}).bound["gnomad_max_af"] is True


def test_a_gnomad_with_no_record_on_the_cohorts_chromosomes_binds_nothing():
    rep = gnomad_report(cohort(chr1=90, chr2=10), {"9"}, enforce_bound=True)
    (w,) = by_code(rep, "gnomad_no_coverage")
    assert w["level"] == "warning" and w["key"] == "gnomad_max_af"
    assert "no record on any chromosome this cohort has variants on" in w["message"]
    assert "(cohort: chr1, chr2; gnomAD: 9)" in w["message"]
    assert "chr1 against 1" in w["message"]
    assert rep.bound["gnomad_max_af"] is False
    assert rep.document["gnomad"]["loaded"] is True
    assert rep.document["gnomad"]["bound"] is False
    # so it is no population filter either, and the reason says gnomAD does not reach
    assert rep.document["population_filter"] is False
    (none,) = by_code(rep, "no_population_frequency_filter")
    assert "does not cover this cohort's variants" in none["message"]
    assert {e["code"] for e in rep.errors} == {
        "gnomad_no_coverage",
        "no_population_frequency_filter",
    }
    text = refusal_message(rep.errors, subject="preset 'strict-denovo'")
    assert text.count("Remedy:") == 1  # one remedy for both, said once
    assert "this cohort's chromosomes and genome build" in text
    assert gnomad_report(cohort(chr1=90, chr2=10), {"9"}).errors == ()  # warn only


def test_a_gnomad_on_few_of_the_cohorts_chromosomes_says_how_few():
    rep = gnomad_report(cohort(chr1=90, chr2=10), {"chr2"})
    (w,) = by_code(rep, "gnomad_no_coverage")
    assert "chromosomes that hold only 10% of this cohort's variants" in w["message"]
    assert rep.bound["gnomad_max_af"] is False
    many = gnomad_report(cohort(**{f"chr{i}": 100 - i for i in range(1, 9)}), {"x"})
    assert (
        "chr1, chr2, chr3, chr4, ..."
        in by_code(many, "gnomad_no_coverage")[0]["message"]
    )


def test_gnomad_not_asked_for_is_not_judged_for_coverage():
    rep = report(variants=cohort(chr1=10), gnomad_loaded=True, gnomad_chroms={"9"})
    assert not by_code(rep, "gnomad_no_coverage") and "gnomad_max_af" not in rep.bound


def probed(found, sampled=1000, chroms=None, variants=None, **kw):
    """A gnomAD report whose pre-run probe found `found` of `sampled` variants."""
    return gnomad_report(
        variants or cohort(chr1=100),
        chroms if chroms is not None else {"chr1"},
        gnomad_probe=GnomadProbe(sampled, found),
        **kw,
    )


def test_gnomad_binds_only_when_it_has_a_record_for_a_tenth_of_the_cohorts_variants():
    """Rows on the right chromosome are not enough: a one-row load, a slice of one
    gene or another build judges almost nothing there."""
    assert GNOMAD_MIN_HIT_RATE == 0.10
    at = probed(100)  # 10%: inclusive
    assert at.bound["gnomad_max_af"] is True and at.warnings == ()
    low = probed(99, enforce_bound=True)
    (w,) = by_code(low, "gnomad_no_coverage")
    assert w["level"] == "warning" and w["key"] == "gnomad_max_af"
    assert (
        "a record for only 9.9% of a sample of 1000 of this cohort's variants"
        in (w["message"])
    )
    assert "at least 10% are needed" in w["message"] and "genome build" in w["message"]
    assert low.bound["gnomad_max_af"] is False
    assert low.document["gnomad"]["bound"] is False
    assert low.document["population_filter"] is False  # so it is no filter at all
    assert {e["code"] for e in low.errors} == {
        "gnomad_no_coverage",
        "no_population_frequency_filter",
    }
    assert probed(99).errors == ()  # warn only, unless the run insists
    # the one-row store: a record, on the right chromosome, for almost nothing
    one_row = probed(1, sampled=5000)
    assert one_row.bound["gnomad_max_af"] is False
    # the counts are in the document, for the app to show
    assert one_row.document["gnomad"]["cohort_variants_sampled"] == 5000
    assert one_row.document["gnomad"]["cohort_variants_with_record"] == 1


def test_a_probe_that_has_nothing_to_judge_by_is_not_a_reason_to_refuse():
    for probe in (None, GnomadProbe(0, 0)):
        rep = gnomad_report(cohort(chr1=100), {"chr1"}, gnomad_probe=probe)
        assert rep.bound["gnomad_max_af"] is True and rep.warnings == ()


def test_the_chromosomes_are_still_judged_first_and_say_the_most():
    rep = probed(0, chroms={"9"})
    (w,) = by_code(rep, "gnomad_no_coverage")
    assert "no record on any chromosome this cohort has variants on" in w["message"]
    assert "chr1 against 1" in w["message"]


def test_the_probe_samples_every_nth_position_of_a_big_cohort_and_all_of_a_small_one():
    class Sess:
        def __init__(self):
            self.sql = []

        def query(self, sql, fmt):
            self.sql.append(sql)
            return self

        def bytes(self):
            return b'{"data": [["chr1", 5, "A", "G"]]}'

    small, big = Sess(), Sess()
    assert sample_cohort_keys(small, "i1", GNOMAD_PROBE_SIZE) == [("chr1", 5, "A", "G")]
    assert "pos %" not in small.sql[0] and "ingest_id = 'i1'" in small.sql[0]
    sample_cohort_keys(big, "i1", 3 * GNOMAD_PROBE_SIZE + 1)
    assert "pos % 4 = 0" in big.sql[0]  # ceil(15001 / 5000)
    assert f"LIMIT {2 * GNOMAD_PROBE_SIZE}" in big.sql[0]


def test_a_probe_rate_is_found_over_sampled():
    assert GnomadProbe(200, 50).rate == 0.25
    assert GnomadProbe(0, 0).rate is None


# ───────────────────── after the run: a filter that judged nothing ─────────────────────


def test_a_gnomad_filter_that_found_no_record_for_any_candidate_is_a_warning():
    w = gnomad_record_warning(
        checked=15, without_record=15, max_af=0.0001, info_af_bound=False
    )
    assert (w["code"], w["level"], w["key"]) == (
        "gnomad_no_record", "warning", "gnomad_max_af",
    )  # fmt: skip
    assert "judged none of the 15 candidates" in w["message"]
    assert "all are kept as rare" in w["message"]
    assert "INFO/AF is not applied either" in w["message"]
    assert "genome build and chromosome naming" in w["message"]
    covered = gnomad_record_warning(
        checked=15, without_record=15, max_af=0.0001, info_af_bound=True
    )
    assert "INFO/AF" not in covered["message"]
    one = gnomad_record_warning(
        checked=1, without_record=1, max_af=0.0001, info_af_bound=True
    )
    assert (
        "none of the 1 candidate:" in one["message"] and "it is kept" in one["message"]
    )
    # when the pre-run probe says gnomAD does judge this cohort, the candidates are
    # probably novel, and the advice is not "check the build"
    reach = gnomad_record_warning(
        checked=15, without_record=15, max_af=0.0001, info_af_bound=True,
        cohort_rate=0.62,
    )  # fmt: skip
    assert "a record for 62% of a sample of the cohort's variants" in reach["message"]
    assert "these are probably novel" in reach["message"]
    assert "genome build" not in reach["message"]
    for checked, without in ((15, 14), (15, 0), (0, 0)):
        assert (
            gnomad_record_warning(
                checked=checked, without_record=without, max_af=0.1, info_af_bound=False
            )
            is None
        )


def stub_analysis(
    seen,
    absent,
    *,
    loaded=True,
    af_bound=False,
    max_af=0.001,
    probe=None,
    reported=True,
):
    """What the gnomAD hooks read from a finished run."""
    effective = report(
        gates=GATES._replace(af_gate="skipped_cohort_derived"),
        gnomad_max_af=max_af,
        gnomad_loaded=loaded,
        gnomad_probe=probe,
    ).document
    effective["gnomad"]["loaded"] = loaded
    effective["info_af"]["bound"] = af_bound
    stats = {
        "effective_block": effective,
        "gnomad_seen": set(range(seen)),
        "gnomad_absent": set(range(absent)),
        "gnomad_reported": reported,
    }
    return SimpleNamespace(gnomad_max_af=max_af, stats=stats)


def test_the_text_footer_says_what_the_gnomad_filter_could_not_judge():
    (all_lost,) = tgo.gnomad_footer(stub_analysis(15, 15))
    assert all_lost.startswith("  WARNING: gnomad_max_af=0.001 judged none of the 15")
    assert tgo.gnomad_footer(stub_analysis(15, 3)) == [
        "  gnomAD: no record for 3 of 15 candidates (kept as rare)"
    ]
    assert tgo.gnomad_footer(stub_analysis(15, 0)) == []
    assert tgo.gnomad_footer(stub_analysis(0, 0)) == []
    # not loaded: its own warning says so, the footer does not repeat it
    assert tgo.gnomad_footer(stub_analysis(15, 15, loaded=False)) == []
    assert tgo.gnomad_footer(stub_analysis(15, 15, max_af=None)) == []
    bare = SimpleNamespace(gnomad_max_af=0.001, stats={})
    assert tgo.gnomad_footer(bare) == []
    # a run whose report is not about gnomAD (it asked for one opt-in gate only)
    quiet = stub_analysis(15, 15, reported=False)
    assert tgo.gnomad_footer(quiet) == []


def test_the_footer_says_how_much_of_the_cohort_gnomad_can_judge():
    (line,) = tgo.gnomad_footer(stub_analysis(9, 8, probe=GnomadProbe(1000, 620)))
    assert line == (
        "  gnomAD: no record for 8 of 9 candidates (kept as rare); it has a "
        "record for 62% of a sample of the cohort's variants"
    )
    (all_lost,) = tgo.gnomad_footer(stub_analysis(9, 9, probe=GnomadProbe(1000, 620)))
    assert "probably novel" in all_lost


def test_the_json_effective_object_carries_the_counts_and_the_late_warning():
    analysis = stub_analysis(15, 15)
    doc = tgo._preset_document({}, SimpleNamespace(analysis=analysis))["effective"]
    assert doc["gnomad"]["candidates_checked"] == 15
    assert doc["gnomad"]["candidates_without_record"] == 15
    assert [w["code"] for w in doc["warnings"]][-1] == "gnomad_no_record"
    assert doc["population_filter"] is False
    # the stored block is not changed by a render
    assert analysis.stats["effective_block"]["gnomad"]["candidates_checked"] is None
    partial = stub_analysis(15, 2)
    doc = tgo._preset_document({}, SimpleNamespace(analysis=partial))["effective"]
    assert doc["gnomad"]["candidates_without_record"] == 2
    assert "gnomad_no_record" not in [w["code"] for w in doc["warnings"]]
    # nothing loaded: the filter ran over an empty store, and a count over a store
    # that held nothing (every candidate "without a record") says nothing
    empty = stub_analysis(9, 9, loaded=False)
    doc = tgo._preset_document({}, SimpleNamespace(analysis=empty))["effective"]
    assert doc["gnomad"]["loaded"] is False
    assert doc["gnomad"]["candidates_checked"] is None
    assert doc["gnomad"]["candidates_without_record"] is None
    # no --gnomad-max-af: the placeholders stay None
    off = stub_analysis(0, 0, max_af=None)
    doc = tgo._preset_document({}, SimpleNamespace(analysis=off))["effective"]
    assert doc["gnomad"]["candidates_checked"] is None
    # a run whose report is not about gnomAD: no counts, no late warning
    quiet = stub_analysis(15, 15, reported=False)
    doc = tgo._preset_document({}, SimpleNamespace(analysis=quiet))["effective"]
    assert doc["gnomad"]["candidates_checked"] is None
    assert "gnomad_no_record" not in [w["code"] for w in doc["warnings"]]


# ─────────────────────── a selected model that cannot run ───────────────────────


NO_REF = ("denovo", "needs_keep_reference")
NO_GENES = ("comphet", "needs_genes")


def test_a_selected_model_that_cannot_run_is_a_warning_and_a_refusal_when_insisted():
    """Its count is 0 because of the database and not the genome: that is never
    a quiet result. It is a warning, and an error for a run that insists."""
    rep = report(blocked_models=[NO_REF, NO_GENES])
    assert [w["message"] for w in by_code(rep, "model_blocked")] == [
        "denovo cannot run: this database has no stored hom-reference calls; "
        "re-ingest with --keep-reference",
        "comphet cannot run: no gene annotations are loaded "
        "(`vcfclick annotations load`)",
    ]
    for w in by_code(rep, "model_blocked"):
        assert (w["level"], w["key"]) == ("warning", "category")
    assert rep.errors == () and rep.document["warnings"] == list(rep.warnings)
    insisting = report(blocked_models=[NO_REF, NO_GENES], enforce_bound=True)
    assert [e["code"] for e in insisting.errors] == ["model_blocked"] * 2
    text = refusal_message(insisting.errors, subject="preset 'strict-denovo'")
    assert text.startswith("preset 'strict-denovo' needs every gate to bind and")
    assert "denovo cannot run" in text and text.count("Remedy:") == 1
    assert "--keep-reference" in text and "--no-require-bound-gates" in text
    # strict_missing alone is about the format fields, not about the models
    assert report(blocked_models=[NO_REF], enforce_format=True).errors == ()
    # no model blocked: nothing said
    assert not by_code(report(), "model_blocked")


def test_a_report_limited_to_one_gate_does_not_speak_for_the_models():
    rep = report(
        gates=GATES._replace(parent_alt_screen=True),
        blocked_models=[NO_REF],
        only={"parent_alt_screen"},
    )
    assert not by_code(rep, "model_blocked")


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
    # ... and the parent-allele screen, which reads the same field
    screen = report(
        broken(), gates._replace(parent_alt_screen=True), enforce_format=True
    )
    assert {e["key"] for e in screen.errors} == set(STRICT_KEYS)


def test_strict_missing_is_judged_only_for_the_gates_that_are_asked_for():
    """A field nobody asked about is not held against a call, so it cannot make
    a gate unbound (and the run refuse) either."""
    nothing = GATES._replace(
        min_gq=0, min_dp=0, min_ab=0.0, max_ab=1.0, strict_missing=True
    )
    rep = report(broken(), nothing, enforce_format=True)
    assert rep.errors == () and not by_code(rep, "gate_unbound")
    assert unbound_strict_keys(rep) == ()
    only_gq = report(broken(), nothing._replace(min_gq=20), enforce_format=True)
    assert [e["key"] for e in only_gq.errors] == ["min_gq"]


def test_an_unbound_gate_says_what_strict_missing_does_to_it():
    strict = GATES._replace(strict_missing=True)
    rep = report(broken(), strict)
    for key, label in (("min_gq", "GQ"), ("min_dp", "DP"), ("min_ab", "AD")):
        (w,) = [w for w in by_code(rep, "gate_unbound") if w["key"] == key]
        assert f"strict_missing counts a call without {label} as a fail" in w["message"]
        assert "removes every candidate it cannot judge" in w["message"]
        assert "passes" not in w["message"]
    # without strict_missing, the message is the tolerant one it always was
    plain = report(broken())
    assert not any("strict_missing" in w["message"] for w in plain.warnings)
    # exempted for a run that goes ahead anyway: said the other way round
    gone = strict._replace(strict_exempt=("min_gq", "min_dp", "min_ab"))
    rep = report(broken(), gone)
    for key, label in (("min_gq", "GQ"), ("min_dp", "DP"), ("min_ab", "AD")):
        (w,) = [w for w in by_code(rep, "gate_unbound") if w["key"] == key]
        assert "strict_missing is not applied to it for this run" in w["message"]
        assert f"a call without {label} passes" in w["message"]
    assert (
        "counts as a pass"
        in by_code(
            report(trio(father=sample("father", rows=1000, dp=980)), gone),
            "gate_partial",
        )[0]["message"]
    )


def test_the_parent_alt_screen_says_what_strict_missing_does_to_it():
    gates = GATES._replace(parent_alt_screen=True, strict_missing=True)
    stats = trio(father=sample("father", ad=0))
    (w,) = [w for w in by_code(report(stats, gates), "gate_unbound")
            if w["key"] == "parent_alt_screen"]  # fmt: skip
    assert "no de novo candidate can pass" in w["message"]
    gone = gates._replace(strict_exempt=("parent_alt_screen",))
    (w,) = [w for w in by_code(report(stats, gone), "gate_unbound")
            if w["key"] == "parent_alt_screen"]  # fmt: skip
    assert "every parent 0/0 would be trusted" in w["message"]
    assert "strict_missing is not applied to it for this run" in w["message"]


def test_the_strict_gates_a_go_ahead_run_exempts_are_the_unbound_ones():
    gates = GATES._replace(strict_missing=True, parent_alt_screen=True)
    assert unbound_strict_keys(report(broken(), gates)) == STRICT_KEYS
    # a gate whose field is there is not exempted
    mixed = trio(
        proband=sample("proband", gq=0), father=sample("father", gq=0),
        mother=sample("mother", gq=0),
    )  # fmt: skip
    assert unbound_strict_keys(report(mixed, gates)) == ("min_gq",)


def test_a_report_limited_to_one_gate_says_only_that_gate():
    """`--parent-alt-screen` typed alone is told about the screen, and not about
    the default GQ / DP gates or the missing population filter it did not ask for."""
    gates = GATES._replace(af_gate="skipped_cohort_derived", parent_alt_screen=True)
    everything = report(broken(), gates)
    assert {w["key"] for w in everything.warnings} > {"parent_alt_screen", "min_gq"}
    scoped = report(broken(), gates, only={"parent_alt_screen"})
    assert [w["key"] for w in scoped.warnings] == ["parent_alt_screen"]
    assert "every parent 0/0 would be trusted" in scoped.warnings[0]["message"]
    assert scoped.document["warnings"] == list(scoped.warnings)
    assert scoped.bound["parent_alt_screen"] is False
    fine = report(trio(), gates, only={"parent_alt_screen"})
    assert fine.warnings == () and fine.bound["parent_alt_screen"] is True


def test_a_report_limited_to_one_gate_claims_nothing_about_the_gates_it_did_not_examine():
    """A typed `--parent-alt-screen` does not probe gnomAD, so a gnomAD that judges
    nothing must not read as bound, nor as a population filter in effect. Every
    `bound` outside the scope, and `population_filter`, is null (not examined),
    and the screen's own verdict is the one `bound` that is given."""
    gates = GATES._replace(af_gate="skipped_cohort_derived", parent_alt_screen=True)
    asked = {"gnomad_max_af": 0.001, "gnomad_loaded": True}
    everything = report(trio(), gates, **asked).document
    assert everything["gnomad"]["bound"] is True  # an unscoped run takes it as probed
    assert everything["population_filter"] is True and everything["gq"]["bound"]

    for stats, screen_binds, share in ((trio(), True, 1.0), (broken(), False, 0.0)):
        scoped = report(stats, gates, only={"parent_alt_screen"}, **asked)
        doc = scoped.document
        assert doc["gnomad"]["bound"] is None and doc["gnomad"]["loaded"] is True
        assert doc["gnomad"]["cohort_variants_sampled"] is None
        assert doc["population_filter"] is None
        assert doc["info_af"] == {"state": "skipped_cohort_derived", "bound": None}
        assert doc["gq"]["bound"] is None and doc["dp"]["bound"] is None
        assert doc["ad"]["bound"] is screen_binds  # the one gate that was examined
        assert scoped.bound == {"parent_alt_screen": screen_binds}
        assert doc["rows"] == {"proband": 100, "father": 100, "mother": 100}
        assert doc["gq"]["proband"] == share  # the facts about the data stay
    # no gate in scope: every verdict is null, none is invented
    none = report(trio(), GATES, only=frozenset(), **asked).document
    assert none["ad"]["bound"] is None and none["population_filter"] is None


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
        "no_population_frequency_filter", "skipped_widening",
        "preset_max_af_not_applied", "gnomad_no_coverage", "gnomad_no_record",
        "pass_only_no_labels", "model_blocked",
    }  # fmt: skip
    # a code is registered by the change that first emits it: none for a run that
    # `db trio` refuses (a duo's de novo answer comes with duo analysis)
    assert "duo_denovo" not in CODES
    assert warning("gate_unbound", "min_gq", "m") == {
        "code": "gate_unbound", "level": "warning", "key": "min_gq", "message": "m",
    }  # fmt: skip
    with pytest.raises(ValueError):
        warning("made_up", None, "m")


def test_the_document_has_exactly_the_documented_keys():
    doc = report().document
    assert list(doc["gnomad"]) == [
        "loaded", "bound", "field", "version", "candidates_checked",
        "candidates_without_record", "cohort_variants_sampled",
        "cohort_variants_with_record",
    ]  # fmt: skip
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
