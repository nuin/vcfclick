"""Sibling segregation: the pure rules (no database server, no CLI).

`cli/segregation_verdict.py` decides every status from member states, and
`cli/segregation_sql.py` prunes rows with a predicate built from the same state
sets. The exhaustive test at the bottom evaluates that predicate in an
in-process DuckDB over every combination of states and requires it to agree
with the verdict, so the two cannot drift.
"""

from __future__ import annotations

import itertools

import pytest

from cli import segregation_verdict as sv
from cli.segregation_sql import ChildInfo, Sibship, SqlOptions, predicate
from cli.segregation_verdict import (
    ABSENT_CALLED,
    ABSENT_UNKNOWN,
    HET,
    HOM_ALT,
    HOM_REF,
    LOW_QUALITY,
    NO_CALL,
    REF_ALT_READS,
    REF_PARTIAL,
    Child,
    Policy,
)

FA, MO = "FA", "MO"


def kids(*spec):
    """kids(("S1","affected",HET), ("S3","unaffected",HOM_REF))"""
    return [Child(*s) for s in spec]


def judge(model, f, m, children, *, policy=None, pc=False, **kw):
    return sv.judge(
        model,
        father_sample=FA,
        mother_sample=MO,
        father=f,
        mother=m,
        children=children,
        partial_site=pc,
        policy=policy or Policy(),
        **kw,
    )


A, U = "affected", "unaffected"


# ───────────────────────────── the status ladder ─────────────────────────────


def test_dominant_consistent_when_every_child_is_positively_resolved():
    j = judge(
        "dominant", HET, HOM_REF,
        kids(("S1", A, HET), ("S2", A, HET), ("S3", U, HOM_REF), ("S4", U, HOM_REF)),
    )  # fmt: skip
    assert j.status == sv.CONSISTENT
    assert j.carrier_parent == FA
    assert (j.n_affected_carrying, j.n_unaffected_carrying, j.n_unaffected_ref) == (
        2,
        0,
        2,
    )
    assert j.n_informative_meioses == 4 and j.p_random == 0.0625


def test_a_no_call_sibling_makes_it_compatible_not_consistent():
    j = judge(
        "dominant", HET, HOM_REF,
        kids(("S1", A, HET), ("S3", U, NO_CALL)),
    )  # fmt: skip
    assert j.status == sv.COMPATIBLE and j.unresolved == ("S3",)


@pytest.mark.parametrize("state", [LOW_QUALITY, ABSENT_UNKNOWN, REF_PARTIAL])
def test_every_unresolved_unaffected_state_is_compatible(state):
    j = judge("dominant", HET, HOM_REF, kids(("S1", A, HET), ("S3", U, state)))
    assert j.status == sv.COMPATIBLE


def test_an_unaffected_carrier_excludes_unless_tolerated():
    children = kids(("S1", A, HET), ("S3", U, HET), ("S4", U, HOM_REF))
    assert judge("dominant", HET, HOM_REF, children) is None
    j = judge(
        "dominant", HET, HOM_REF, children, policy=Policy(max_unaffected_carriers=1)
    )
    assert j.status == sv.REDUCED_PENETRANCE
    assert j.contradicting == ("S3",)


def test_reduced_penetrance_is_never_consistent_whatever_else_is_resolved():
    children = kids(
        ("S1", A, HET), ("S2", A, HET), ("S3", U, HOM_ALT), ("S4", U, HOM_REF)
    )
    j = judge(
        "dominant", HET, HOM_REF, children, policy=Policy(max_unaffected_carriers=1)
    )
    assert j.status == sv.REDUCED_PENETRANCE and j.status != sv.CONSISTENT


def test_a_het_unaffected_at_a_partial_call_site_is_unresolved_not_a_carrier():
    children = kids(("S1", A, HET), ("S3", U, HET))
    assert judge("dominant", HET, HOM_REF, children) is None
    j = judge("dominant", HET, HOM_REF, children, pc=True)
    assert j.status == sv.COMPATIBLE and j.unresolved == ("S3",)


def test_a_partial_site_parent_proof_is_only_compatible():
    j = judge(
        "dominant", HET, REF_PARTIAL, kids(("S1", A, HET), ("S3", U, HOM_REF)), pc=True
    )
    assert j.status == sv.COMPATIBLE


def test_absent_called_is_inferred_and_never_consistent():
    children = kids(("S1", A, HET), ("S3", U, ABSENT_CALLED))
    strict = judge("dominant", HET, HOM_REF, children)
    assert strict.status == sv.COMPATIBLE  # state 3 proves nothing under strict
    called = judge(
        "dominant", HET, ABSENT_CALLED, children, policy=Policy(absent_called=True)
    )
    assert called.status == sv.CONSISTENT_INFERRED
    assert called.inferred_reference and set(called.inferred) == {MO, "S3"}


def test_denovo_parents_must_be_stored_reference_even_under_absent_called():
    children = kids(("S1", A, HET), ("S3", U, HOM_REF))
    called = Policy(absent_called=True)
    assert judge("denovo", ABSENT_CALLED, HOM_REF, children, policy=called) is None
    assert (
        judge("denovo", HOM_REF, HOM_REF, children, policy=called).status
        == sv.CONSISTENT
    )
    # an unaffected sibling's reference may be inferred, the parents' may not
    j = judge(
        "denovo",
        HOM_REF,
        HOM_REF,
        kids(("S1", A, HET), ("S3", U, ABSENT_CALLED)),
        policy=called,
    )
    assert j.status == sv.CONSISTENT_INFERRED


def test_denovo_reference_with_alt_reads_is_a_contradiction_for_a_sibling():
    assert (
        judge(
            "denovo", HOM_REF, HOM_REF, kids(("S1", A, HET), ("S3", U, REF_ALT_READS))
        )
        is None
    )
    assert (
        judge(
            "denovo", REF_ALT_READS, HOM_REF, kids(("S1", A, HET), ("S3", U, HOM_REF))
        )
        is None
    )


def test_denovo_mosaic_suspect_from_parental_alt_reads_or_two_affected_groups():
    one = kids(("S1", A, HET), ("S3", U, HOM_REF))
    assert not judge(
        "denovo", HOM_REF, HOM_REF, one, parent_alt_reads=(0, 0)
    ).parental_mosaic_suspect
    assert judge(
        "denovo", HOM_REF, HOM_REF, one, parent_alt_reads=(1, 0)
    ).parental_mosaic_suspect
    two = kids(("S1", A, HET), ("S2", A, HET), ("S3", U, HOM_REF))
    assert judge("denovo", HOM_REF, HOM_REF, two).parental_mosaic_suspect
    # an MZ pair is one group: no recurrence evidence
    mz = Policy(groups={"S1": 0, "S2": 0})
    assert not judge("denovo", HOM_REF, HOM_REF, two, policy=mz).parental_mosaic_suspect


def test_recessive_needs_a_positive_observation_of_each_unaffected_child():
    base = [("S1", A, HOM_ALT), ("S2", A, HOM_ALT)]
    ok = judge("recessive", HET, HET, kids(*base, ("S3", U, HOM_REF), ("S4", U, HET)))
    assert ok.status == sv.CONSISTENT
    assert ok.n_informative_meioses == 4 and ok.p_random == round(0.25**2 * 0.75**2, 10)
    assert (
        judge("recessive", HET, HET, kids(*base, ("S3", U, NO_CALL))).status
        == sv.COMPATIBLE
    )
    assert judge("recessive", HET, HET, kids(*base, ("S3", U, HOM_ALT))) is None
    # absent_called satisfies "not hom-alt" with no inference label
    j = judge("recessive", HET, HET, kids(*base, ("S3", U, ABSENT_CALLED)))
    assert j.status == sv.CONSISTENT and not j.inferred_reference


def test_recessive_parents_must_both_be_het_and_affected_hom_alt():
    children = kids(("S1", A, HOM_ALT), ("S3", U, HOM_REF))
    assert judge("recessive", HET, HOM_REF, children) is None
    assert (
        judge("recessive", HET, HET, kids(("S1", A, HET), ("S3", U, HOM_REF))) is None
    )


def test_at_least_one_affected_child_must_be_strictly_satisfied():
    all_open = kids(("S1", A, NO_CALL), ("S2", A, LOW_QUALITY), ("S3", U, HOM_REF))
    assert judge("dominant", HET, HOM_REF, all_open) is None
    one = kids(("S1", A, NO_CALL), ("S2", A, HET), ("S3", U, HOM_REF))
    assert judge("dominant", HET, HOM_REF, one).status == sv.COMPATIBLE


def test_the_result_does_not_depend_on_child_order():
    children = kids(
        ("S1", A, NO_CALL), ("S2", A, HET), ("S3", U, HOM_REF), ("S4", U, NO_CALL)
    )
    forward = judge("dominant", HET, HOM_REF, children)
    backward = judge("dominant", HET, HOM_REF, list(reversed(children)))
    assert forward.status == backward.status
    assert set(forward.unresolved) == set(backward.unresolved)


def test_a_hom_alt_carrier_parent_is_allowed_and_uninformative():
    j = judge("dominant", HOM_ALT, HOM_REF, kids(("S1", A, HET), ("S3", U, HOM_REF)))
    assert j.status == sv.CONSISTENT and j.carrier_parent_hom_alt
    assert j.n_informative_meioses == 0 and j.p_random == 1.0


def test_parent_phenotype_use_downgrades_an_unaffected_carrier_parent():
    children = kids(("S1", A, HET), ("S3", U, HOM_REF))
    pheno = {FA: "unaffected", MO: None}
    ignore = judge("dominant", HET, HOM_REF, children, parent_phenotype=pheno)
    assert ignore.status == sv.CONSISTENT
    use = judge(
        "dominant", HET, HOM_REF, children, parent_phenotype=pheno,
        policy=Policy(parent_phenotype_use=True),
    )  # fmt: skip
    assert use.status == sv.REDUCED_PENETRANCE and use.contradicting == (FA,)
    affected_parent = judge(
        "dominant", HET, HOM_REF, children, parent_phenotype={FA: "affected"},
        policy=Policy(parent_phenotype_use=True),
    )  # fmt: skip
    assert affected_parent.status == sv.CONSISTENT


def test_mz_pair_counts_as_one_informative_meiosis():
    children = kids(("S1", A, HET), ("S2", A, HET), ("S3", U, HOM_REF))
    plain = judge("dominant", HET, HOM_REF, children)
    mz = judge(
        "dominant", HET, HOM_REF, children, policy=Policy(groups={"S1": 0, "S2": 0})
    )
    assert plain.n_informative_meioses == 3 and plain.p_random == 0.125
    assert mz.n_informative_meioses == 2 and mz.p_random == 0.25


def test_denovo_has_no_informativeness_number():
    j = judge("denovo", HOM_REF, HOM_REF, kids(("S1", A, HET), ("S3", U, HOM_REF)))
    assert j.n_informative_meioses is None and j.p_random is None


def test_judge_refuses_comphet():
    with pytest.raises(ValueError):
        judge("comphet", HET, HOM_REF, kids(("S1", A, HET)))


# ───────────────────────────── compound het pairs ─────────────────────────────


def pair(
    p_states,
    m_states,
    *,
    pheno=(A, A, U, U),
    partial=(False, False),
    parents=(HET, HOM_REF, HOM_REF, HET),
    policy=None,
):
    names = ["S1", "S2", "S3", "S4"]
    mk = lambda states: [Child(n, ph, st) for n, ph, st in zip(names, pheno, states)]  # noqa: E731
    return sv.judge_pair(
        paternal=mk(p_states),
        maternal=mk(m_states),
        paternal_partial=partial[0],
        maternal_partial=partial[1],
        parent_states=parents,
        policy=policy or Policy(),
    )


def test_pair_consistent_when_each_unaffected_child_is_reference_at_one_hit():
    j = pair([HET, HET, HOM_REF, HET], [HET, HET, HET, HOM_REF])
    assert j.status == sv.CONSISTENT
    assert j.n_informative_meioses == 4


def test_pair_excluded_when_an_unaffected_child_carries_both():
    assert pair([HET, HET, HET, HOM_REF], [HET, HET, HET, HOM_REF]) is None


def test_pair_excluded_when_an_affected_child_misses_a_hit():
    assert pair([HET, HOM_REF, HOM_REF, HOM_REF], [HET, HET, HOM_REF, HOM_REF]) is None
    assert pair([HET, HOM_ALT, HOM_REF, HOM_REF], [HET, HET, HOM_REF, HOM_REF]) is None


def test_pair_compatible_when_a_child_is_unresolved():
    assert (
        pair([HET, HET, HOM_REF, NO_CALL], [HET, HET, HOM_REF, HET]).status
        == sv.COMPATIBLE
    )
    assert (
        pair([HET, NO_CALL, HOM_REF, HOM_REF], [HET, HET, HOM_REF, HOM_REF]).status
        == sv.COMPATIBLE
    )


def test_pair_hets_at_a_partial_site_are_unresolved_not_both():
    j = pair([HET, HET, HET, HOM_REF], [HET, HET, HET, HOM_REF], partial=(True, False))
    assert j is not None and j.status == sv.COMPATIBLE


def test_pair_inferred_reference_lowers_the_status():
    j = pair(
        [HET, HET, ABSENT_CALLED, HET], [HET, HET, HET, HOM_REF],
        policy=Policy(absent_called=True),
    )  # fmt: skip
    assert j.status == sv.CONSISTENT_INFERRED


def test_pair_needs_one_affected_child_carrying_both():
    assert (
        pair([NO_CALL, NO_CALL, HOM_REF, HOM_REF], [HET, HET, HOM_REF, HOM_REF]) is None
    )


# ───────────────── the SQL predicate and the verdict agree ─────────────────


SIB = Sibship(
    ingest_id="i1",
    family_id="f",
    father=FA,
    mother=MO,
    children=(ChildInfo("A", A), ChildInfo("U1", U), ChildInfo("U2", U)),
)
STATES = range(9)


def _expected(model, f, m, a, u1, u2, pc, *, absent_called, n):
    """Does the verdict list the row (ignoring --parent-phenotype)?"""
    if model == "comphet":
        return (
            sv.parents_fit("comphet", f, m, absent_called)
            and sv.classify_affected("comphet", a) != sv.CONTRADICTS
            and sv.classify_affected("comphet", a) == sv.SATISFIED
        )
    j = sv.judge(
        model,
        father_sample=FA,
        mother_sample=MO,
        father=f,
        mother=m,
        children=[Child("A", A, a), Child("U1", U, u1), Child("U2", U, u2)],
        partial_site=bool(pc),
        policy=Policy(absent_called=absent_called, max_unaffected_carriers=n),
    )
    return j is not None


@pytest.mark.parametrize("model", sv.MODELS)
@pytest.mark.parametrize("absent_called", [False, True])
@pytest.mark.parametrize("n", [0, 1])
def test_the_sql_predicate_equals_the_verdict_over_every_state(model, absent_called, n):
    duckdb = pytest.importorskip("duckdb")
    opts = SqlOptions(
        n_samples=7, absent_called=absent_called, max_unaffected_carriers=n
    )
    sql_pred = predicate(model, SIB, opts)
    con = duckdb.connect()
    con.execute(
        "CREATE TABLE t AS SELECT * FROM (VALUES "
        + ", ".join(
            f"({f},{m},{a},{u1},{u2},{pc})"
            for f, m, a, u1, u2, pc in itertools.product(
                STATES, STATES, STATES, STATES, STATES, (0, 1)
            )
        )
        + ") v(st_f, st_m, st_0, st_1, st_2, pc)"
    )
    got = {
        tuple(r)
        for r in con.execute(
            f"SELECT st_f, st_m, st_0, st_1, st_2, pc FROM t WHERE {sql_pred}"
        ).fetchall()
    }
    # the predicate prunes, the verdict decides: every listed row must pass the
    # predicate, and (for the row models) nothing else may
    for f, m, a, u1, u2, pc in itertools.product(
        STATES, STATES, STATES, STATES, STATES, (0, 1)
    ):
        want = _expected(model, f, m, a, u1, u2, pc, absent_called=absent_called, n=n)
        assert ((f, m, a, u1, u2, pc) in got) == want, (model, f, m, a, u1, u2, pc)
