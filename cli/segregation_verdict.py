"""Sibling segregation: member states, predicates, status ladder.

Pure (no SQL, no storage), so every rule is a unit test over state tuples.
`cli/segregation_sql.py` builds its WHERE clauses from the same state sets
below, and `judge` / `judge_pair` re-derive each candidate from the states the
SQL returned, so the two cannot read different rules.

A member's STATE at a candidate site says what the stored data proves about
that sample:

  0 hom_ref           stored gt=0 row that passes the gates (--keep-reference)
  1 het, 2 hom_alt    stored row that passes the gates
  3 absent_called     no row anywhere and the site's call accounting is
                      complete: called, no ALT allele, but no GQ/DP evidence
  4 no_call           a missing_genotypes row
  5 low_quality       a stored row that fails GQ/DP or allele balance
  6 absent_unknown    no row anywhere and the accounting is not complete
  7 ref_partial_site  stored gt=0 at a site where some sample has a partial
                      call ('0/.', './1'): the 0 cannot be trusted
  8 ref_alt_reads     (de novo only) stored gt=0 whose alt reads say otherwise

`consistent` means every constrained child is positively satisfied on stored,
gated evidence. A child the data cannot resolve makes the candidate
`compatible`. `reduced_penetrance` is an unaffected child that carries the
variant (only with --max-unaffected-carriers N or --parent-phenotype use) and
is never labelled consistent. `consistent_inferred` used `--absent called`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

HOM_REF, HET, HOM_ALT, ABSENT_CALLED, NO_CALL = 0, 1, 2, 3, 4
LOW_QUALITY, ABSENT_UNKNOWN, REF_PARTIAL, REF_ALT_READS = 5, 6, 7, 8

STATE_NAMES: tuple[str, ...] = (
    "hom_ref",
    "het",
    "hom_alt",
    "absent_called",
    "no_call",
    "low_quality",
    "absent_unknown",
    "ref_partial_site",
    "ref_alt_reads",
)

CONSISTENT = "consistent"
CONSISTENT_INFERRED = "consistent_inferred"
COMPATIBLE = "compatible"
REDUCED_PENETRANCE = "reduced_penetrance"
STATUSES: tuple[str, ...] = (
    CONSISTENT,
    CONSISTENT_INFERRED,
    COMPATIBLE,
    REDUCED_PENETRANCE,
)

MODELS: tuple[str, ...] = ("denovo", "recessive", "dominant", "comphet")
ROW_MODELS: tuple[str, ...] = ("denovo", "recessive", "dominant")

SATISFIED, UNRESOLVED, CONTRADICTS = "satisfied", "unresolved", "contradicts"

# States that leave an affected child's claim open: nothing is proven either way.
_OPEN = frozenset({NO_CALL, LOW_QUALITY, ABSENT_UNKNOWN, REF_PARTIAL, REF_ALT_READS})

# What an affected child may be in (anything else contradicts the model) and
# what positively satisfies it. A comphet half is the dominant pattern.
AFFECTED_OK: dict[str, frozenset[int]] = {
    "denovo": frozenset({HET, HOM_ALT}) | _OPEN,
    "recessive": frozenset({HOM_ALT}) | _OPEN,
    "dominant": frozenset({HET}) | _OPEN,
    "comphet": frozenset({HET}) | _OPEN,
}
AFFECTED_STRICT: dict[str, frozenset[int]] = {
    "denovo": frozenset({HET, HOM_ALT}),
    "recessive": frozenset({HOM_ALT}),
    "dominant": frozenset({HET}),
    "comphet": frozenset({HET}),
}
CARRIES = frozenset({HET, HOM_ALT})


def reference_states(absent_called: bool) -> frozenset[int]:
    """States that prove a PARENT or unaffected sibling does not carry the
    allele: a stored 0/0, one at a partial-call site (the candidate is then
    only `compatible`) and, under --absent called, `absent_called`."""
    base = {HOM_REF, REF_PARTIAL}
    if absent_called:
        base.add(ABSENT_CALLED)
    return frozenset(base)


def unaffected_contradicts(model: str, state: int, partial_site: bool) -> bool:
    """Does an unaffected child's state contradict the model? Dominant, de
    novo and a comphet half: it carries (a het at a partial-call site is
    unresolved, not a carrier; for de novo a 0/0 with alt reads is one too).
    Recessive: it is hom-alt."""
    if model == "recessive":
        return state == HOM_ALT
    if state == HOM_ALT or (state == HET and not partial_site):
        return True
    return model == "denovo" and state == REF_ALT_READS


def classify_affected(model: str, state: int) -> str:
    if state in AFFECTED_STRICT[model]:
        return SATISFIED
    return UNRESOLVED if state in AFFECTED_OK[model] else CONTRADICTS


def classify_unaffected(
    model: str, state: int, partial_site: bool, absent_called: bool
) -> str:
    if unaffected_contradicts(model, state, partial_site):
        return CONTRADICTS
    if model == "recessive":
        # not hom-alt needs no reference inference, but still a positive
        # observation: a no-call or unknown stays unresolved
        return SATISFIED if state in (HOM_REF, HET, ABSENT_CALLED) else UNRESOLVED
    if state == HOM_REF or (state == ABSENT_CALLED and absent_called):
        return SATISFIED
    return UNRESOLVED


def parents_fit(model: str, father: int, mother: int, absent_called: bool) -> bool:
    """The parental pattern of a model. A parent may be unresolved only by
    state 7 (reference not proven at a partial-call site)."""
    if model == "denovo":
        return father in (HOM_REF, REF_PARTIAL) and mother in (HOM_REF, REF_PARTIAL)
    if model == "recessive":
        return father == HET and mother == HET
    ref = reference_states(absent_called)
    return (father in CARRIES and mother in ref) or (
        father in ref and mother in CARRIES
    )


@dataclass(frozen=True)
class Policy:
    """The options a verdict depends on."""

    absent_called: bool = False
    max_unaffected_carriers: int = 0
    parent_phenotype_use: bool = False
    # sample -> group id: members of a declared MZ pair share one id
    groups: Mapping[str, int] = field(default_factory=dict)

    def group(self, sample: str) -> object:
        return self.groups.get(sample, sample)


@dataclass(frozen=True)
class Child:
    sample: str
    phenotype: str | None  # affected | unaffected | None
    state: int


@dataclass(frozen=True)
class Judgement:
    status: str
    unresolved: tuple[str, ...]
    contradicting: tuple[str, ...]
    inferred: tuple[str, ...]  # members whose reference proof was state 3
    inferred_reference: bool
    n_affected_carrying: int
    n_unaffected_carrying: int
    n_unaffected_ref: int
    n_informative_meioses: int | None
    p_random: float | None
    parental_mosaic_suspect: bool
    carrier_parent: str | None  # sample id of the one parent that carries
    carrier_parent_hom_alt: bool = False


def _p_dominant(k: int) -> float:
    return round(0.5**k, 10)


def _p_recessive(affected: int, unaffected: int) -> float:
    return round(0.25**affected * 0.75**unaffected, 10)


def _independent(policy: Policy, samples: Sequence[str]) -> int:
    return len({policy.group(s) for s in samples})


def judge(
    model: str,
    *,
    father_sample: str,
    mother_sample: str,
    father: int,
    mother: int,
    children: Sequence[Child],
    partial_site: bool,
    policy: Policy,
    parent_phenotype: Mapping[str, str | None] | None = None,
    parent_alt_reads: tuple[int, int] | None = None,
) -> Judgement | None:
    """The verdict on one site for one of denovo, recessive, dominant: None
    when it is not a candidate, else the status and the evidence behind it."""
    if model not in ROW_MODELS:
        raise ValueError(f"judge() takes a row model, not {model!r}")
    R = policy.absent_called
    if not parents_fit(model, father, mother, R):
        return None
    affected = [c for c in children if c.phenotype == "affected"]
    unaffected = [c for c in children if c.phenotype == "unaffected"]

    verdicts = {c.sample: classify_affected(model, c.state) for c in affected}
    if CONTRADICTS in verdicts.values() or SATISFIED not in verdicts.values():
        return None
    uverdicts = {
        c.sample: classify_unaffected(model, c.state, partial_site, R)
        for c in unaffected
    }
    contradicting = [s for s, v in uverdicts.items() if v == CONTRADICTS]
    if len(contradicting) > policy.max_unaffected_carriers:
        return None

    unresolved = [s for s, v in {**verdicts, **uverdicts}.items() if v == UNRESOLVED]
    parent_partial = REF_PARTIAL in (father, mother)

    # who proved a reference by `absent_called` (inference, not a stored row)
    inferred: list[str] = []
    if model != "recessive":
        for sample, state in ((father_sample, father), (mother_sample, mother)):
            if state == ABSENT_CALLED:
                inferred.append(sample)
    for c in unaffected:
        if c.state == ABSENT_CALLED and uverdicts[c.sample] == SATISFIED:
            inferred.append(c.sample)
    # recessive: state 3 on an unaffected sibling proves "not hom-alt" by the
    # accounting itself, so it is listed as evidence but never lowers the status
    inferred_reference = bool(inferred) and model != "recessive"

    carrier_parent = None
    carrier_hom_alt = False
    if model == "dominant":
        if father in CARRIES:
            carrier_parent, carrier_hom_alt = father_sample, father == HOM_ALT
        else:
            carrier_parent, carrier_hom_alt = mother_sample, mother == HOM_ALT

    contradicting_all = list(contradicting)
    if (
        policy.parent_phenotype_use
        and model == "dominant"
        and carrier_parent is not None
        and (parent_phenotype or {}).get(carrier_parent) == "unaffected"
    ):
        contradicting_all.append(carrier_parent)

    if contradicting_all:
        status = REDUCED_PENETRANCE
    elif unresolved or parent_partial:
        status = COMPATIBLE
    elif inferred_reference:
        status = CONSISTENT_INFERRED
    else:
        status = CONSISTENT

    resolved_affected = [s for s, v in verdicts.items() if v == SATISFIED]
    resolved_unaffected = [s for s, v in uverdicts.items() if v == SATISFIED]
    n_inf: int | None
    p: float | None
    if model == "denovo":
        n_inf, p = None, None
    elif model == "dominant":
        k = _independent(policy, resolved_affected + resolved_unaffected)
        if carrier_hom_alt:
            k = 0  # a hom-alt carrier parent transmits to every child: uninformative
        n_inf, p = k, _p_dominant(k)
    else:
        a = _independent(policy, resolved_affected)
        u = _independent(policy, resolved_unaffected)
        n_inf, p = a + u, _p_recessive(a, u)

    mosaic = False
    if model == "denovo":
        carrying = [c.sample for c in affected if c.state in CARRIES]
        mosaic = (
            any((parent_alt_reads or (0, 0))) or _independent(policy, carrying) >= 2
        )

    return Judgement(
        status=status,
        unresolved=tuple(unresolved),
        contradicting=tuple(contradicting_all),
        inferred=tuple(inferred),
        inferred_reference=inferred_reference,
        n_affected_carrying=sum(1 for c in affected if c.state in CARRIES),
        n_unaffected_carrying=sum(1 for c in unaffected if c.state in CARRIES),
        n_unaffected_ref=sum(
            1
            for c in unaffected
            if c.state == HOM_REF or (c.state == ABSENT_CALLED and R)
        ),
        n_informative_meioses=n_inf,
        p_random=p,
        parental_mosaic_suspect=mosaic,
        carrier_parent=carrier_parent,
        carrier_parent_hom_alt=carrier_hom_alt,
    )


@dataclass(frozen=True)
class PairJudgement:
    status: str
    unresolved: tuple[str, ...]
    contradicting: tuple[str, ...]
    n_informative_meioses: int
    p_random: float


def _carries_for_sure(state: int, partial_site: bool) -> bool:
    return state == HOM_ALT or (state == HET and not partial_site)


def judge_pair(
    *,
    paternal: Sequence[Child],
    maternal: Sequence[Child],
    paternal_partial: bool,
    maternal_partial: bool,
    parent_states: Sequence[int],
    policy: Policy,
) -> PairJudgement | None:
    """The verdict on a (paternal hit, maternal hit) pair in one gene: both
    hits het in every affected child (an affected child unresolved at either
    makes the pair compatible, at least one carries both), and no unaffected
    child carries both. `paternal` / `maternal` hold each child's state at the
    two sites; `parent_states` the parents' states at both (a state 3 or 7
    there lowers the status). --max-unaffected-carriers is ignored."""
    R = policy.absent_called
    by_p = {c.sample: c for c in paternal}
    by_m = {c.sample: c for c in maternal}
    unresolved: list[str] = []
    contradicting: list[str] = []
    both = 0
    resolved_aff: list[str] = []
    resolved_unaff: list[str] = []
    inferred = ABSENT_CALLED in parent_states
    for sample, p in by_p.items():
        m = by_m[sample]
        if p.phenotype == "affected":
            a, b = (
                classify_affected("comphet", p.state),
                classify_affected("comphet", m.state),
            )
            if CONTRADICTS in (a, b):
                return None
            if a == b == SATISFIED:
                both += 1
                resolved_aff.append(sample)
            else:
                unresolved.append(sample)
        elif p.phenotype == "unaffected":
            if _carries_for_sure(p.state, paternal_partial) and _carries_for_sure(
                m.state, maternal_partial
            ):
                contradicting.append(sample)
                continue
            ref = {HOM_REF} | ({ABSENT_CALLED} if R else set())
            proof = [s for s in (p.state, m.state) if s in ref]
            if proof:
                resolved_unaff.append(sample)
                if all(s == ABSENT_CALLED for s in proof):
                    inferred = True
            else:
                unresolved.append(sample)
    if both == 0 or contradicting:
        return None
    if unresolved or REF_PARTIAL in parent_states:
        status = COMPATIBLE
    elif inferred:
        status = CONSISTENT_INFERRED
    else:
        status = CONSISTENT
    a = _independent(policy, resolved_aff)
    u = _independent(policy, resolved_unaff)
    return PairJudgement(status, tuple(unresolved), (), a + u, _p_recessive(a, u))
