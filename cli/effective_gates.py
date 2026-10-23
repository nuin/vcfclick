"""The effective-gate report: does each gate actually bind in THIS data?

The quality gates are NULL-tolerant (a call with no GQ passes), a VCF's INFO/AF
may be the cohort's own frequency, and gnomAD may not be loaded. A preset
named "strict" can therefore be a no-op, and used to be silent about it. This
report says, per gate, whether the field it reads exists in the data, and turns
that into warnings (printed on stdout, carried in the JSON `effective`
object) and, on request, into a refusal.

A gate BINDS when its field is present on at least ENFORCEABLE_MIN of the rows
it reads:

  * GQ, DP and AD on the CANDIDATE ROWS of each selected model: the rows of the
    proband and of the parents at the sites where the model's genotype pattern
    holds (`MODEL_PREDICATES`: for de novo the sites where the proband carries
    and both parents are 0/0), counted over the very join `db trio` runs, before
    any quality gate. The gate is judged on the genotype classes of those rows
    that it reads (`gate_reads`: allele balance reads the proband's het rows, the
    parent-allele screen the parents' 0/0 rows, GQ and DP every row), and the
    weakest model decides. A caller can fill a field on its inherited calls and
    leave it empty on the de novo ones, and neither a share over all of a sample's
    rows nor one over a genotype class says anything about the rows a candidate
    rests on. Only a model with no candidate rows at all (or none of the classes
    the gate reads) is judged on the sample's stored rows of those classes
    (`collect_sample_stats`), and a role with none of those on its whole sample;
  * the INFO/AF gate on the variants (and its probe not saying the AF is the
    cohort's own);
  * gnomAD when it is loaded, has records on the chromosomes that hold
    ENFORCEABLE_MIN of the cohort's variants (a slice on other chromosomes, or
    under another naming, judges nothing) AND has a record for at least
    GNOMAD_MIN_HIT_RATE of a sample of the cohort's own variants (a sparse test
    load, a slice of one gene, another genome build judge almost nothing, whatever
    chromosomes they have a row on).

`--require-bound-gates` makes a gate that cannot bind, a selected model that
cannot run, and the absence of any population-frequency filter, an error;
`--strict-missing` makes that error for the GQ / DP / AD gates (and the
parent-allele screen) alone, because it turns a missing field into a fail and a
field missing everywhere would then fail everything. `--no-require-bound-gates`
typed on the command line overrides both for one run, and for that run
strict_missing is not applied to the gates that cannot bind (`Gates.strict_exempt`):
they are NULL-tolerant, as their warnings say, rather than failing every call.

The FILTER gate (`pass_only`) binds when some variant of the cohort carries a
FILTER label other than PASS or unset: `variants.filter` is NULL for both PASS and
'.', so on a VCF that labels nothing the gate removes nothing. That is a warning
(`pass_only_no_labels`) and never a refusal, because a VCF that was hard-filtered
before it was ingested is all PASS and the gate is right to remove nothing.

A strict_missing gate is judged only when it is asked for (min_gq > 0, min_dp
> 0, an allele-balance band narrower than 0..1, the parent-allele screen on):
an unasked-for field is not held against a call, so it cannot make a gate
"unbound" either.

The gnomAD probe is a fact about the COHORT, known before the run, so a
run that insists (`require_bound`, `strict-denovo`) can refuse on it. What is
left for after the run is a fact about the CANDIDATES: a de novo candidate is
novel by nature, so one that gnomAD has no record for is the expected case and
not a reason to refuse; `gnomad_record_warning` says it when it is all of them.

`collect_sample_stats` (one pass over the three samples' rows),
`collect_candidate_stats` (one grouped count per selected model over its
candidate join) and the variant and gnomAD probes are the queries; `build_report`
is pure so every threshold is a unit test. The numbers below are the author's
choices, named so a reviewer can change one place.
"""

from __future__ import annotations

import json
import math
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import NamedTuple

from cli.trio_analysis import (
    AF_GATE_SKIPPED,
    CATEGORIES,
    MODEL_PREDICATES,
    PASS_ONLY_SQL,
    STRICT_GATE_KEYS,
    Gates,
    Trio,
    parent_sample_joins,
)
from cli.warning_codes import warning

# A gate binds when its field is present on at least this share of the rows.
ENFORCEABLE_MIN = 0.90
# The proband's median depth below this many times min_dp: the gate is close to
# the assay's own depth, so the preset may be too strict for it.
THIN_COVERAGE_FACTOR = 1.5
# gnomAD binds only when it has a record for at least this share of a sample of
# the cohort's own variants. A real load has one for most of a cohort's common
# variants; the floor is low on purpose (a cohort's private variants are
# legitimately absent), so only a load that judges almost nothing fails it.
GNOMAD_MIN_HIT_RATE = 0.10
# how many of the cohort's variants the probe looks up (one Arrow join)
GNOMAD_PROBE_SIZE = 5000

GNOMAD_FIELD = "popmax (grpmax else AF)"

# which warning codes a refusal is made of
_BOUND_CODES = (
    "gate_unbound",
    "model_blocked",
    "gnomad_not_loaded",
    "gnomad_no_coverage",
    "no_population_frequency_filter",
)
# the gates strict_missing turns into "fail when the field is missing": the
# GQ / DP / AD gates and the parent-allele screen
STRICT_KEYS = STRICT_GATE_KEYS
# why a model cannot run (the `blocked_reason` values a run can have), said as
# the end of "<model> cannot run: ..."
_BLOCKED_WHY = {
    "needs_keep_reference": "this database has no stored hom-reference calls; "
    "re-ingest with --keep-reference",
    "needs_genes": "no gene annotations are loaded (`vcfclick annotations load`)",
}
_NAMES_SHOWN = 4  # chromosome names listed in a message


# The genotype classes of a stored row. A gate reads different rows in different
# models (a parent's 0/0 for de novo, a parent's 0/1 for recessive, ...), and a
# caller can fill GQ, DP and AD on one class and not another, so a share over all
# of a sample's rows says little about the rows a gate reads. (Nor does a share
# over one class of them: see `CandidateStats`.)
REF, HET, ALT = "ref", "het", "alt"
_ANY = (REF, HET, ALT)
_LABEL = {REF: "0/0", HET: "0/1", ALT: "1/1"}


def genotype_class(gt: int) -> str:
    """`genotypes.gt` is the alt-allele count: 0 hom-ref, 1 het, 2 hom-alt (a
    haploid ALT call is stored as 2)."""
    return HET if gt == 1 else ALT if gt >= 2 else REF


@dataclass(frozen=True)
class GenotypeClassStats:
    """What one sample's stored rows of one genotype class carry."""

    cls: str  # ref | het | alt
    rows: int
    gq: int
    dp: int
    ad: int


@dataclass(frozen=True)
class SampleGateStats:
    """What one sample's stored rows carry."""

    sample: str
    role: str  # proband | father | mother
    rows: int
    gq: int  # rows with a GQ
    dp: int  # rows with a DP
    ad: int  # rows with both ad_ref and ad_alt
    depths: tuple[tuple[int, int], ...] = ()  # (dp, rows) histogram
    # the same counts by genotype class; empty for a stand-in with no split
    classes: tuple[GenotypeClassStats, ...] = ()

    def fraction(self, name: str) -> float | None:
        """Share of this sample's rows that carry `name` (gq | dp | ad); None
        for a sample with no stored rows, which says nothing either way."""
        return getattr(self, name) / self.rows if self.rows else None

    def share(self, name: str, classes: Collection[str]) -> float | None:
        """Share of this sample's rows of these genotype `classes` that carry
        `name`; None when it has none of those rows (or no class split)."""
        picked = [c for c in self.classes if c.cls in classes]
        rows = sum(c.rows for c in picked)
        return sum(getattr(c, name) for c in picked) / rows if rows else None

    @property
    def median_dp(self) -> float | None:
        return median_from_histogram(self.depths)


def median_from_histogram(pairs: Sequence[tuple[int, int]]) -> float | None:
    """The median of values given as (value, count) pairs. (The engines differ
    in their quantile syntax, so the histogram comes back and the median is
    taken here.)"""
    total = sum(n for _, n in pairs)
    if not total:
        return None
    ordered = sorted(pairs)

    def at(rank: int) -> int:  # the value at 0-based rank
        seen = 0
        for value, n in ordered:
            seen += n
            if rank < seen:
                return value
        return ordered[-1][0]

    lo, hi = at((total - 1) // 2), at(total // 2)
    return (lo + hi) / 2


def _num(v: float | None) -> float | int | None:
    if v is None:
        return None
    return int(v) if float(v).is_integer() else v


def collect_sample_stats(
    sess, ingest_id: str, samples: Sequence[tuple[str, str]]
) -> list[SampleGateStats]:
    """For each (sample, role): row count, rows carrying GQ / DP / AD (in all and
    by genotype class), and the DP histogram. One statement over the samples'
    stored rows, however many samples there are."""
    from storage import sql_quote_str

    listed = ", ".join(sql_quote_str(s) for s, _ in samples)
    sql = (
        "SELECT sample_id, gt, dp, count(*), "
        "sum(CASE WHEN gq IS NOT NULL THEN 1 ELSE 0 END), "
        "sum(CASE WHEN ad_ref IS NOT NULL AND ad_alt IS NOT NULL THEN 1 ELSE 0 END) "
        f"FROM genotypes WHERE ingest_id = {sql_quote_str(ingest_id)} "
        f"AND sample_id IN ({listed}) GROUP BY sample_id, gt, dp"
    )
    data = json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]
    acc = {s: [0, 0, 0, 0, {}] for s, _ in samples}  # rows, gq, dp, ad, hist
    by_class: dict[str, dict[str, list[int]]] = {s: {} for s, _ in samples}
    for sample, gt, dp, n, n_gq, n_ad in data:
        a = acc[sample]
        n, n_gq, n_ad = int(n), int(n_gq or 0), int(n_ad or 0)
        n_dp = n if dp is not None else 0
        c = by_class[sample].setdefault(genotype_class(int(gt)), [0, 0, 0, 0])
        for tally in (a, c):
            tally[0] += n
            tally[1] += n_gq
            tally[3] += n_ad
        a[2] += n_dp
        c[2] += n_dp
        if dp is not None:
            a[4][int(dp)] = a[4].get(int(dp), 0) + n
    return [
        SampleGateStats(
            s,
            role,
            a[0],
            a[1],
            a[2],
            a[3],
            tuple(sorted(a[4].items())),
            tuple(
                GenotypeClassStats(k, *by_class[s][k]) for k in _ANY if k in by_class[s]
            ),
        )
        for (s, role) in samples
        for a in (acc[s],)
    ]


@dataclass(frozen=True)
class CandidateStats:
    """What the rows a model reads at its candidate sites carry.

    A candidate site is one where the model's genotype pattern holds
    (`MODEL_PREDICATES`: for de novo the proband carries and both parents are
    0/0). `roles` gives, for the proband and each parent, the rows at those
    sites by genotype class with the number that carry GQ, DP and AD. These are
    the rows the quality gates are applied to, so a field's presence on them is
    what decides whether a gate can act; the same field on the sample's other
    rows (the inherited calls of a de novo search) decides nothing."""

    model: str
    sites: int  # candidate sites: the model's predicate only, no quality gate
    roles: Mapping[str, tuple[GenotypeClassStats, ...]]

    def share(self, role: str, name: str, classes: Collection[str]) -> float | None:
        """Share of this role's candidate rows of these genotype `classes` that
        carry `name` (gq | dp | ad); None when there are none."""
        picked = [c for c in self.roles.get(role, ()) if c.cls in classes]
        rows = sum(c.rows for c in picked)
        return sum(getattr(c, name) for c in picked) / rows if rows else None


def candidate_presence_sql(model: str, trio: Trio) -> str:
    """One grouped count over the join `trio_sql` runs for `model` (proband `g`,
    parents `f` and `m`, the model's genotype predicate and no quality or rarity
    gate): per combination of the three genotypes, the number of rows and, for
    each of the three samples, how many carry GQ, DP and both allele depths."""
    from storage import count_expr, sql_quote_str

    cols = []
    for alias in ("g", "f", "m"):
        cols += [
            f"sum(CASE WHEN {alias}.gq IS NOT NULL THEN 1 ELSE 0 END)",
            f"sum(CASE WHEN {alias}.dp IS NOT NULL THEN 1 ELSE 0 END)",
            f"sum(CASE WHEN {alias}.ad_ref IS NOT NULL AND {alias}.ad_alt IS NOT NULL "
            "THEN 1 ELSE 0 END)",
        ]
    return (
        f"SELECT g.gt, f.gt, m.gt, {count_expr()}, {', '.join(cols)} "
        f"FROM genotypes g {parent_sample_joins(trio.father, trio.mother)} "
        f"WHERE g.ingest_id = {sql_quote_str(trio.ingest_id)} "
        f"AND g.sample_id = {sql_quote_str(trio.proband)} "
        f"AND {MODEL_PREDICATES[model]} GROUP BY g.gt, f.gt, m.gt"
    )


def _candidate_stats(model: str, data: Sequence[Sequence]) -> CandidateStats:
    """A `CandidateStats` from the rows of `candidate_presence_sql`."""
    tallies: dict[str, dict[str, list[int]]] = {role: {} for role in _ROLES}
    sites = 0
    for row in data:
        n = int(row[3])
        sites += n
        for i, role in enumerate(_ROLES):  # proband g, father f, mother m
            tally = tallies[role].setdefault(genotype_class(int(row[i])), [0, 0, 0, 0])
            tally[0] += n
            for j in range(3):  # gq, dp, ad
                tally[1 + j] += int(row[4 + 3 * i + j] or 0)
    return CandidateStats(
        model,
        sites,
        {
            role: tuple(
                GenotypeClassStats(k, *tallies[role][k])
                for k in _ANY
                if k in tallies[role]
            )
            for role in _ROLES
        },
    )


def collect_candidate_stats(
    sess, trio: Trio, models: Sequence[str]
) -> dict[str, CandidateStats]:
    """For each of `models` that has a genotype predicate, what the rows at its
    candidate sites carry. One grouped count per distinct predicate (compound het
    reads the dominant pattern and shares its count). A model without one, or a
    trio without both parents, is left out, and its gates are judged on the
    sample's stored rows instead."""
    out: dict[str, CandidateStats] = {}
    if trio.father is None or trio.mother is None:
        return out
    by_predicate: dict[str, CandidateStats] = {}
    for model in models:
        predicate = MODEL_PREDICATES.get(model)
        if predicate is None:
            continue
        if predicate not in by_predicate:
            raw = sess.query(candidate_presence_sql(model, trio), "JSONCompact")
            data = json.loads(raw.bytes().decode())["data"]
            by_predicate[predicate] = _candidate_stats(model, data)
        out[model] = replace(by_predicate[predicate], model=model)
    return out


@dataclass(frozen=True)
class VariantStats:
    """How many of the cohort's variants carry INFO/AF, and a FILTER label."""

    n: int
    with_af: int
    # (chrom, variants on it), most variants first: where gnomAD must have
    # records to judge anything
    chroms: tuple[tuple[str, int], ...] = ()
    # variants whose FILTER is a label other than PASS, unset or '.' (None: not
    # counted). variants.filter is NULL for both PASS and '.', so a gate that keeps
    # PASS has something to remove only where this is above 0.
    non_pass: int | None = None

    @property
    def af_fraction(self) -> float | None:
        return self.with_af / self.n if self.n else None


def collect_variant_stats(
    sess, ingest_id: str, *, filter_present: bool = True
) -> VariantStats:
    """One pass over the cohort's variants: how many, how many carry INFO/AF,
    how many carry a FILTER label other than PASS, and how many sit on each
    chromosome."""
    from storage import count_expr, sql_quote_str

    labelled = (
        f"sum(CASE WHEN NOT {PASS_ONLY_SQL} THEN 1 ELSE 0 END)"
        if filter_present
        else "0"
    )
    sql = (
        f"SELECT chrom, {count_expr()}, "
        "sum(CASE WHEN info_AF IS NOT NULL THEN 1 ELSE 0 END), "
        f"{labelled} "
        f"FROM variants v WHERE ingest_id = {sql_quote_str(ingest_id)} GROUP BY chrom"
    )
    data = json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]
    per = sorted(((c, int(n or 0)) for c, n, _, _ in data), key=lambda t: (-t[1], t[0]))
    return VariantStats(
        sum(n for _, n in per),
        sum(int(w or 0) for _, _, w, _ in data),
        tuple(per),
        sum(int(x or 0) for _, _, _, x in data),
    )


@dataclass(frozen=True)
class GnomadProbe:
    """How many of a sample of the cohort's own variants the loaded gnomAD data
    can judge (has a record with a frequency for)."""

    sampled: int
    found: int

    @property
    def rate(self) -> float | None:
        return self.found / self.sampled if self.sampled else None


def sample_cohort_keys(sess, ingest_id: str, n: int) -> list[tuple[str, int, str, str]]:
    """About GNOMAD_PROBE_SIZE of the cohort's `n` variants as (chrom, pos, ref,
    alt), spread over the whole genome: every stride-th position (the same on
    both engines, and the same sample on every run). A cohort of at most
    GNOMAD_PROBE_SIZE variants is taken whole."""
    from storage import sql_quote_str

    stride = max(1, -(-n // GNOMAD_PROBE_SIZE))
    where = f"ingest_id = {sql_quote_str(ingest_id)}"
    if stride > 1:
        where += f" AND pos % {stride} = 0"
    sql = (
        f"SELECT chrom, pos, ref, alt FROM variants WHERE {where} "
        f"LIMIT {2 * GNOMAD_PROBE_SIZE}"
    )
    data = json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]
    return [(c, int(p), r, a) for c, p, r, a in data]


def collect_gnomad_probe(sess, ingest_id: str, variants: VariantStats) -> GnomadProbe:
    """Look the cohort's sampled variants up in the loaded gnomAD data, exactly as
    the `--gnomad-max-af` filter will (an exact chrom / pos / ref / alt join), and
    count the ones it can judge. A load that holds one row, a slice of one gene or
    another genome build finds almost none, whatever chromosomes it has rows on."""
    from annotations import gnomad_af_many

    keys = sample_cohort_keys(sess, ingest_id, variants.n) if variants.n else []
    unique = list(dict.fromkeys(keys))
    found = gnomad_af_many(unique)
    judged = sum(1 for g in found.values() if g.popmax is not None)
    return GnomadProbe(len(unique), judged)


# ───────────────────────────────── the report ─────────────────────────────────


@dataclass(frozen=True)
class EffectiveReport:
    document: dict  # the JSON `effective` object (warnings included)
    warnings: tuple[dict, ...]
    bound: dict[str, bool]  # file key -> does the gate bind in this data
    errors: tuple[dict, ...] = ()  # the warnings a refusal is made of


def pct(fraction: float) -> str:
    """A share as a percent, never rounded up to 100 when some rows lack it."""
    value = math.floor(fraction * 1000) / 10
    return f"{value:g}%"


def _weakest(
    stats: Sequence[SampleGateStats], name: str, roles: Sequence[str] | None = None
) -> tuple[float | None, list[str]]:
    """The lowest share of rows carrying `name` among the samples that have
    rows (limited to `roles`), and the roles below ENFORCEABLE_MIN."""
    shares = [
        (s.fraction(name), s.role)
        for s in stats
        if (roles is None or s.role in roles) and s.rows
    ]
    if not shares:
        return None, []
    lowest = min(f for f, _ in shares)
    below = [r for f, r in shares if f < ENFORCEABLE_MIN]
    return lowest, below


# The rows each model reads at its candidate sites, by role: the join in
# cli/trio_analysis.trio_sql (proband g, parents f / m) as genotype classes. A
# de novo call needs both parents' stored 0/0; a recessive one needs both 0/1;
# dominant and compound het need one parent with the allele and the other 0/0.
# These are the classes of `MODEL_PREDICATES`, which `collect_candidate_stats`
# counts the rows of.
_SITE_ROWS: dict[str, dict[str, tuple[str, ...]]] = {
    "denovo": {"proband": (HET, ALT), "father": (REF,), "mother": (REF,)},
    "recessive": {"proband": (ALT,), "father": (HET,), "mother": (HET,)},
    "dominant": {"proband": (HET,), "father": _ANY, "mother": _ANY},
    "comphet": {"proband": (HET,), "father": _ANY, "mother": _ANY},
}
_ROLES = ("proband", "father", "mother")

# one selected model's reading of a gate: (model, {role: genotype classes})
Read = tuple[str, dict[str, tuple[str, ...]]]


def gate_reads(
    key: str, categories: Sequence[str], *, strict: bool = False
) -> list[Read]:
    """For gate `key`, one (model, {role: genotype classes}) for each selected
    model that reads it: the rows its WHERE clause looks at on a candidate (not
    the sample's rows at large).

      min_gq, min_dp      every row the model joins
      min_ab              het rows only: the proband's (a hom-alt or hom-ref call
                          is not gated), and, under strict_missing, a parent's
                          het call must carry its depths too
      parent_alt_screen   the parents' 0/0 rows, in the de novo model only

    A model this version does not know reads every row (never none)."""
    models = CATEGORIES if "all" in categories else tuple(categories)
    out: list[Read] = []
    for model in models:
        rows = _SITE_ROWS.get(model) or dict.fromkeys(_ROLES, _ANY)
        if key == "parent_alt_screen":
            read = (
                {r: (REF,) for r in ("father", "mother")} if model == "denovo" else {}
            )
        elif key == "min_ab":
            read = {
                role: (HET,)
                for role, classes in rows.items()
                if HET in classes and (strict or role == "proband")
            }
        else:
            read = dict(rows)
        if read:
            out.append((model, read))
    return out


def judged_gate_keys(
    gates: Gates, only: Collection[str] | None = None
) -> tuple[str, ...]:
    """The gate keys a report judges on the rows of a model's candidates: the GQ /
    DP / AD gates and the parent-allele screen that this run asks for (and, in a
    report limited to some gates, that are among them)."""
    asked = [
        key
        for key, on in (
            ("min_gq", gates.min_gq > 0),
            ("min_dp", gates.min_dp > 0),
            ("min_ab", gates.ab_requested),
            ("parent_alt_screen", gates.parent_alt_screen),
        )
        if on
    ]
    return tuple(k for k in asked if only is None or k in only)


def candidate_models(
    gates: Gates, categories: Sequence[str], only: Collection[str] | None = None
) -> tuple[str, ...]:
    """The selected models whose candidate rows a report of `gates` reads: the
    ones that read a gate it judges. (strict_missing makes a parent's het rows
    count for allele balance; the models that read them are included whether or
    not the gate turns out to be exempt, which only costs a count.)"""
    return tuple(
        dict.fromkeys(
            model
            for key in judged_gate_keys(gates, only)
            for model, _ in gate_reads(key, categories, strict=gates.strict_missing)
        )
    )


def strict_gates_in_force(gates: Gates) -> tuple[str, ...]:
    """The strict_missing gates that fail a call for a missing field in this run:
    asked for, strict_missing on, and not exempted for want of the field."""
    return tuple(k for k in judged_gate_keys(gates) if gates.strict_for(k))


def strict_check_models(gates: Gates, categories: Sequence[str]) -> tuple[str, ...]:
    """The selected models that a strict_missing gate in force can drop a candidate
    from: the ones that read such a gate. A model that shares another's candidates
    (compound het reads the dominant pattern) is left out when that one is
    selected: the same variants would be counted and said twice."""
    models = dict.fromkeys(
        model
        for key in strict_gates_in_force(gates)
        for model, _ in gate_reads(key, categories, strict=True)
    )
    by_candidates = {MODEL_PREDICATES.get(m, m): m for m in reversed(list(models))}
    return tuple(m for m in models if by_candidates[MODEL_PREDICATES.get(m, m)] == m)


def strict_missing_warning(model: str, dropped: int) -> dict | None:
    """After the run: None, or the warning that strict_missing dropped `dropped`
    candidates of `model` that pass every other gate (the ones a run without
    strict_missing lists beyond this one), only because a field a gate needs is
    missing. A candidate is a variant, for compound het too."""
    if dropped <= 0:
        return None
    noun = "candidate variants" if model == "comphet" else "candidates"
    return warning(
        "strict_missing_dropped",
        "strict_missing",
        f"strict_missing dropped {dropped} {model} {noun} that pass every other "
        "gate, only because a GQ, DP or allele depth they are judged on is "
        "missing from their rows; --no-strict-missing lists them",
    )


def _who(roles: Sequence[str], classes: Collection[str] | None) -> str:
    names = "/".join(roles)
    if classes is None or set(classes) == set(_ANY):
        return names
    return f"{names} {' or '.join(_LABEL[c] for c in _ANY if c in classes)}"


class Judged(NamedTuple):
    """The verdict of `_weakest_read` on one field."""

    lowest: float | None  # the lowest share among the rows the gate reads
    below: str  # the rows under ENFORCEABLE_MIN ("father/mother 0/0 rows at denovo
    # candidate sites"), empty when none is
    where: str  # where the lowest share was taken (" at denovo candidate sites"),
    # empty for rows at large


def _at(model: str | None) -> str:
    return f" at {model} candidate sites" if model else ""


def _weakest_read(
    stats: Sequence[SampleGateStats],
    name: str,
    reads: Sequence[Read],
    candidates: Mapping[str, CandidateStats] | None = None,
) -> Judged:
    """The lowest share of rows carrying `name` among the rows the gate reads, and
    a phrase for the roles and rows below ENFORCEABLE_MIN.

    A role's share is taken, for each selected model that reads the gate, over the
    model's CANDIDATE rows of the genotype classes it reads (`candidates`), and the
    weakest model decides. A model with no candidate rows of those classes (none
    at all, or no candidate counted for it) is judged on the sample's stored rows of
    the classes; a role a model reads but that has none of those on its whole
    sample. A role no model reads is not judged (a de novo candidate's parents are
    0/0, so their het rows' depths are not asked for). `reads` is never empty: a
    gate that no selected model reads is not applicable (`gate_not_applicable`)
    and is not judged at all."""
    present = {s.role: s for s in stats}
    # role -> (share, the model whose candidates it was taken over | None, classes)
    worst: dict[str, tuple[float, str | None, tuple[str, ...] | None]] = {}
    read_roles: set[str] = set()
    for model, read in reads:
        at = candidates.get(model) if candidates else None
        for role, classes in read.items():
            read_roles.add(role)
            share, source = None, None
            if at is not None and at.sites:
                share, source = at.share(role, name, classes), model
            if share is None:
                s = present.get(role)
                share, source = (
                    (s.share(name, classes) if s is not None else None),
                    None,
                )
            if share is not None and (role not in worst or share < worst[role][0]):
                worst[role] = (share, source, classes)
    for s in stats:
        if s.role in read_roles and s.role not in worst and s.rows:
            worst[s.role] = (s.fraction(name), None, None)
    if not worst:
        return Judged(None, "", "")
    roles = [r for r in _ROLES if r in worst]
    first = min(roles, key=lambda r: worst[r][0])
    groups: dict[tuple[str | None, tuple[str, ...] | None], list[str]] = {}
    for role in roles:
        share, source, classes = worst[role]
        if share < ENFORCEABLE_MIN:
            groups.setdefault((source, classes), []).append(role)
    below = ", ".join(
        f"{_who(rs, classes)} rows{_at(source)}"
        for (source, classes), rs in groups.items()
    )
    return Judged(worst[first][0], below, _at(worst[first][1]))


def build_report(
    *,
    stats: Sequence[SampleGateStats],
    variants: VariantStats,
    gates: Gates,
    gnomad_max_af: float | None,
    gnomad_loaded: bool,
    gnomad_chroms: Collection[str] | None = None,
    gnomad_probe: GnomadProbe | None = None,
    filter_present: bool = True,
    categories: Sequence[str] = ("all",),
    candidates: Mapping[str, CandidateStats] | None = None,
    blocked_models: Sequence[tuple[str, str]] = (),
    assay: str | None = None,
    subject: str = "preset",
    skipped: Sequence[dict] = (),
    preset_max_af: float | None = None,
    enforce_bound: bool = False,
    enforce_format: bool = False,
    only: Collection[str] | None = None,
) -> EffectiveReport:
    """Whether each requested gate binds, and what to say about it.

    `gnomad_chroms` are the chromosomes the loaded gnomAD data has records on
    and `gnomad_probe` how much of a sample of the cohort's variants it can
    judge (None for either: not looked at, so any loaded gnomAD is taken to
    cover the cohort). `candidates` are what the rows at each selected model's
    candidate sites carry (`collect_candidate_stats`): the GQ / DP / AD gates are
    judged on them, and on the sample's stored rows only for a model that has none
    (None: no candidate was counted, as in a unit test with synthetic counts).
    `only` limits the report to these gate keys: a run that
    asked for one opt-in gate, and nothing else, is told about that gate only.
    The warnings are those about these keys, and in the document every `bound`
    (and `population_filter`) of a gate outside them is null: it was not
    examined, which is not the same as bound.
    `blocked_models` are the (model, reason) pairs of the models the run selected
    that cannot run on this database (reasons from `BLOCKED_REASONS`): an empty
    model is not a clean result, so each is a `model_blocked` warning, and an
    error for a run that insists on its gates.
    `skipped` are the preset's skipped keys (from the echo); `preset_max_af` is
    the max_af a preset supplied, set only so a skip of it can be explained;
    `enforce_bound` / `enforce_format` choose which warnings are also errors
    (see the module docstring). `subject` names what the user chose, in the
    thin-coverage sentence."""
    warnings: list[dict] = []
    bound: dict[str, bool] = {}
    selected = CATEGORIES if "all" in categories else tuple(categories)

    def not_applicable(key: str, shown: str, why: str) -> None:
        """A requested gate that no selected model reads: it removes nothing, so
        it is not bound. A warning and not a refusal (not in `_BOUND_CODES`):
        nothing about the data is wrong, the run's models just do not use it."""
        bound[key] = False
        many = len(selected) > 1
        warnings.append(
            warning(
                "gate_not_applicable",
                key,
                f"{shown} is not applied by the model{'s' if many else ''} this run "
                f"selects ({', '.join(selected)}): {why}, so it removes nothing",
            )
        )

    def strict_note(key: str, label: str) -> str:
        """What the gate does with a call that lacks the field, said when that
        is not the NULL-tolerant default."""
        if gates.strict_for(key):
            return (
                f"; strict_missing counts a call without {label} as a fail, so "
                "this gate removes every candidate it cannot judge"
            )
        if gates.strict_missing:  # asked for, and switched off for this run
            return (
                f"; strict_missing is not applied to it for this run, so a call "
                f"without {label} passes"
            )
        return ""

    def check_format(key: str, field_name: str, label: str, shown: str, why: str):
        """One GQ / DP / AD gate: bound, partial or unbound, judged on the rows
        the selected models read it from (their candidate rows), or not applicable
        when none of them reads it (`why` says what the gate judges)."""
        reads = gate_reads(key, categories, strict=gates.strict_for(key))
        if not reads:
            return not_applicable(key, shown, why)
        judged = _weakest_read(stats, field_name, reads, candidates)
        lowest = judged.lowest
        if lowest is None:
            bound[key] = False
            warnings.append(
                warning(
                    "gate_unbound",
                    key,
                    f"{shown} cannot be enforced: no stored genotype rows to judge",
                )
            )
        elif lowest < ENFORCEABLE_MIN:
            bound[key] = False
            warnings.append(
                warning(
                    "gate_unbound",
                    key,
                    f"{shown} cannot be enforced: only {pct(lowest)} of "
                    f"{judged.below} carry {label}" + strict_note(key, label),
                )
            )
        else:
            bound[key] = True
            if lowest < 1.0:
                outcome = "fail" if gates.strict_for(key) else "pass"
                warnings.append(
                    warning(
                        "gate_partial",
                        key,
                        f"{shown} enforced on {pct(lowest)} of rows{judged.where} "
                        f"(a call with no {label} counts as a {outcome})",
                    )
                )

    if gates.min_gq > 0:
        check_format(
            "min_gq", "gq", "GQ", f"min_gq={gates.min_gq}", "no model reads a call"
        )
    if gates.min_dp > 0:
        check_format(
            "min_dp", "dp", "DP", f"min_dp={gates.min_dp}", "no model reads a call"
        )
    if gates.ab_requested:
        check_format(
            "min_ab",
            "ad",
            "AD",
            f"min_ab={gates.min_ab:g}/max_ab={gates.max_ab:g}",
            "the band judges heterozygous calls and they read none (a homozygous-alt "
            "proband call is not allele-balance gated)",
        )
        bound["max_ab"] = bound["min_ab"]

    screen_reads = gate_reads("parent_alt_screen", categories)
    if gates.parent_alt_screen and not screen_reads:
        not_applicable(
            "parent_alt_screen",
            "parent_alt_screen",
            "the screen is read by the denovo model only",
        )
    elif gates.parent_alt_screen:
        judged = _weakest_read(stats, "ad", screen_reads, candidates)
        lowest = judged.lowest
        ok = lowest is not None and lowest >= ENFORCEABLE_MIN
        bound["parent_alt_screen"] = ok
        if not ok:
            why = (
                "no stored rows for the parents"
                if lowest is None
                else f"only {pct(lowest)} of {judged.below} carry AD"
            )
            if gates.strict_for("parent_alt_screen"):
                outcome = (
                    "; strict_missing counts a parent 0/0 without allele depths "
                    "as a fail, so no de novo candidate can pass"
                )
            else:
                outcome = ", so every parent 0/0 would be trusted" + (
                    "; strict_missing is not applied to it for this run"
                    if gates.strict_missing
                    else ""
                )
            warnings.append(
                warning(
                    "gate_unbound",
                    "parent_alt_screen",
                    f"parent_alt_screen cannot be enforced: {why}{outcome}",
                )
            )
    if gates.pass_only:
        # the gate keeps PASS and unset, and the column holds NULL for both, so it
        # has something to remove only where some variant carries another label.
        # A warning and not a refusal (_BOUND_CODES): a VCF that was hard-filtered
        # before it was ingested is all PASS, and the gate is right to remove nothing.
        if not filter_present:
            bound["pass_only"] = False
            warnings.append(
                warning(
                    "gate_unbound",
                    "pass_only",
                    "pass_only cannot be enforced: the variants table has no "
                    "filter column",
                )
            )
        elif variants.non_pass == 0:
            bound["pass_only"] = False
            warnings.append(
                warning(
                    "pass_only_no_labels",
                    "pass_only",
                    "pass_only removes nothing: no variant in this database "
                    "carries a FILTER label other than PASS or unset",
                )
            )
        else:
            bound["pass_only"] = True

    # the population-frequency filters
    af_run = gates.af_gate != AF_GATE_SKIPPED and gates.max_af < 1
    af_share = variants.af_fraction
    af_ok = af_run and af_share is not None and af_share >= ENFORCEABLE_MIN
    if af_run:
        bound["max_af"] = af_ok
        if not af_ok:
            have = (
                "no variants"
                if af_share is None
                else f"only {pct(af_share)} of variants carry INFO/AF"
            )
            warnings.append(
                warning(
                    "gate_unbound",
                    "max_af",
                    f"max_af={gates.max_af:g} cannot be enforced: {have}",
                )
            )
    elif gates.af_gate == AF_GATE_SKIPPED:
        bound["max_af"] = False
    gnomad_bound = gnomad_loaded
    if gnomad_max_af is not None:
        if not gnomad_loaded:
            warnings.append(
                warning(
                    "gnomad_not_loaded",
                    "gnomad_max_af",
                    f"gnomad_max_af={gnomad_max_af:g} cannot be enforced: no "
                    "gnomAD data is loaded (`vcfclick annotations load-gnomad`)",
                )
            )
        else:
            why_not = _gnomad_unbound(
                gnomad_max_af, variants, gnomad_chroms, gnomad_probe
            )
            if why_not:
                gnomad_bound = False
                warnings.append(warning("gnomad_no_coverage", "gnomad_max_af", why_not))
        bound["gnomad_max_af"] = gnomad_bound
    population = af_ok or (gnomad_max_af is not None and gnomad_bound)
    if not population:
        warnings.append(
            warning(
                "no_population_frequency_filter",
                None,
                "NO population frequency filter applied: "
                + _why_no_population(
                    gates, af_run, af_share, gnomad_max_af, gnomad_loaded, gnomad_bound
                ),
            )
        )

    proband = next((s for s in stats if s.role == "proband"), None)
    median = proband.median_dp if proband else None
    if (
        gates.min_dp > 0
        and bound.get("min_dp")
        and median is not None
        and median < THIN_COVERAGE_FACTOR * gates.min_dp
    ):
        hint = f" (the {subject} says it was written for {assay} data)" if assay else ""
        warnings.append(
            warning(
                "thin_coverage",
                "min_dp",
                f"median proband DP is {_num(median):g} but min_dp is "
                f"{gates.min_dp}; this assay may be too thin for the "
                f"{subject}{hint}",
            )
        )
    for model, reason in blocked_models:
        warnings.append(
            warning(
                "model_blocked",
                "category",
                f"{model} cannot run: {_BLOCKED_WHY.get(reason, reason)}",
            )
        )
    for item in skipped:
        if item.get("widening"):
            warnings.append(warning("skipped_widening", item["key"], item["message"]))
    if preset_max_af is not None and gates.af_gate == AF_GATE_SKIPPED:
        warnings.append(
            warning(
                "preset_max_af_not_applied",
                "max_af",
                f"max_af {preset_max_af:g} from the preset was not applied: "
                "INFO/AF is this cohort's own frequency (gates.af_gate = "
                "skipped_cohort_derived); pass --max-af to force it, or use "
                "--gnomad-max-af",
            )
        )

    if only is not None:
        warnings = [w for w in warnings if w["key"] in only]
        # what was not examined is not claimed either: no `bound` for a gate
        # outside the scope (the document says null for it)
        bound = {k: v for k, v in bound.items() if k in only}
    errors = tuple(
        w
        for w in warnings
        if (enforce_bound and w["code"] in _BOUND_CODES)
        or (enforce_format and w["code"] == "gate_unbound" and w["key"] in STRICT_KEYS)
    )
    document = _document(
        stats,
        variants,
        gates,
        gnomad_loaded,
        gnomad_bound,
        filter_present,
        bound,
        population,
        warnings,
        gnomad_probe,
        only,
    )
    return EffectiveReport(document, tuple(warnings), bound, errors)


def unbound_strict_keys(report: EffectiveReport) -> tuple[str, ...]:
    """The strict_missing gates whose field the data lacks: the keys a run that
    was told to go ahead anyway exempts from strict_missing (`Gates.strict_exempt`)."""
    return tuple(
        dict.fromkeys(
            w["key"]
            for w in report.warnings
            if w["code"] == "gate_unbound" and w["key"] in STRICT_KEYS
        )
    )


def _gnomad_unbound(
    gnomad_max_af: float,
    variants: VariantStats,
    gnomad_chroms: Collection[str] | None,
    probe: GnomadProbe | None,
) -> str | None:
    """Why the loaded gnomAD data cannot judge this cohort, or None when it can
    (or when nothing was looked at). The chromosomes are tested first: they say
    the most when the naming differs."""
    if gnomad_chroms is not None:
        share, cohort = _gnomad_coverage(variants, gnomad_chroms)
        if share is not None and share < ENFORCEABLE_MIN:
            return _no_coverage_message(gnomad_max_af, share, cohort, gnomad_chroms)
    rate = probe.rate if probe is not None else None
    if rate is not None and rate < GNOMAD_MIN_HIT_RATE:
        return (
            f"gnomad_max_af={gnomad_max_af:g} cannot be enforced: the loaded gnomAD "
            f"data has a record for only {pct(rate)} of a sample of "
            f"{probe.sampled} of this cohort's variants (at least "
            f"{pct(GNOMAD_MIN_HIT_RATE)} are needed), so most are kept as rare "
            "without being judged. Check that its genome build and chromosome "
            "naming match the cohort's and that it covers more than a slice of one "
            "region"
        )
    return None


def _gnomad_coverage(
    variants: VariantStats, gnomad_chroms: Collection[str]
) -> tuple[float | None, list[str]]:
    """(the share of the cohort's variants on a chromosome gnomAD has records
    on, the cohort's chromosomes most variants first); None for the share when
    the cohort has no variants to judge it by."""
    if not variants.n or not variants.chroms:
        return None, []
    covered = sum(n for chrom, n in variants.chroms if chrom in gnomad_chroms)
    return covered / variants.n, [c for c, _ in variants.chroms]


def _names(names: Sequence[str]) -> str:
    shown = ", ".join(n[:20] for n in list(names)[:_NAMES_SHOWN])
    return shown + (", ..." if len(names) > _NAMES_SHOWN else "")


def _no_coverage_message(
    gnomad_max_af: float,
    share: float,
    cohort_chroms: Sequence[str],
    gnomad_chroms: Collection[str],
) -> str:
    head = f"gnomad_max_af={gnomad_max_af:g} cannot be enforced: "
    if share == 0:
        return (
            head + "the loaded gnomAD data has no record on any chromosome this "
            f"cohort has variants on (cohort: {_names(cohort_chroms)}; gnomAD: "
            f"{_names(sorted(gnomad_chroms))}). Check the chromosome naming "
            "(chr1 against 1) and the genome build"
        )
    return (
        head + f"the loaded gnomAD data has records on chromosomes that hold only "
        f"{pct(share)} of this cohort's variants, so the rest are kept as rare "
        "without being judged"
    )


def gnomad_record_warning(
    *,
    checked: int,
    without_record: int,
    max_af: float,
    info_af_bound: bool,
    cohort_rate: float | None = None,
) -> dict | None:
    """After the gnomAD filter has run: None, or the warning that it found a
    record for none of the `checked` candidates. `cohort_rate` is the share of
    the cohort's sampled variants gnomAD could judge (the pre-run probe). Only when
    it is known and reaches GNOMAD_MIN_HIT_RATE does the filter work on this cohort
    and the candidates are most likely novel variants, which is what a de novo
    candidate is. When it is not known, or the probe already said gnomAD cannot
    judge the cohort (under GNOMAD_MIN_HIT_RATE: a slice that does not reach it,
    another build or chromosome naming), the candidates are not "probably novel":
    nothing was looked up for most of the cohort either, and the advice is to check
    the build and the naming."""
    if not checked or without_record < checked:
        return None
    one = checked == 1
    tail = (
        ""
        if info_af_bound
        else ", and INFO/AF is not applied either, so no population frequency "
        "filter took effect"
    )
    if cohort_rate is None or cohort_rate < GNOMAD_MIN_HIT_RATE:
        advice = (
            ". Check that its genome build and chromosome naming match the cohort's"
        )
    else:
        advice = (
            f". It does hold a record for {pct(cohort_rate)} of a sample of the "
            f"cohort's variants, so {'this one is' if one else 'these are'} "
            "probably novel"
        )
    return warning(
        "gnomad_no_record",
        "gnomad_max_af",
        f"gnomad_max_af={max_af:g} judged none of the {checked} "
        f"candidate{'' if one else 's'}: the loaded gnomAD data has no record for "
        f"{'it' if one else 'any of them'}, so {'it is' if one else 'all are'} "
        f"kept as rare{tail}{advice}",
    )


def _why_no_population(
    gates: Gates,
    af_run: bool,
    af_share: float | None,
    gnomad_max_af: float | None,
    gnomad_loaded: bool,
    gnomad_bound: bool,
) -> str:
    if gates.af_gate == AF_GATE_SKIPPED:
        af = "INFO/AF is this cohort's own frequency, so the --max-af gate was skipped"
    elif not af_run:
        af = "the --max-af gate is off"
    else:
        af = "INFO/AF is missing from most variants"
    if gnomad_max_af is None:
        gn = "--gnomad-max-af is not set" + (
            "" if gnomad_loaded else " and no gnomAD data is loaded"
        )
    elif gnomad_loaded and not gnomad_bound:
        gn = "the loaded gnomAD data does not cover this cohort's variants"
    else:
        gn = "no gnomAD data is loaded"
    return f"{af}, and {gn}"


def _document(
    stats: Sequence[SampleGateStats],
    variants: VariantStats,
    gates: Gates,
    gnomad_loaded: bool,
    gnomad_bound: bool,
    filter_present: bool,
    bound: dict[str, bool],
    population: bool,
    warnings: list[dict],
    probe: GnomadProbe | None = None,
    only: Collection[str] | None = None,
) -> dict:
    """The `effective` object. In a report limited to some gates (`only`) the
    shares of rows, the medians and the facts about the store stay, but every
    `bound`, and `population_filter`, is reported only for what the scope looked
    at and is null for the rest: a gate nobody examined (gnomAD, say, which a
    scoped run does not probe) is neither bound nor unbound. The same holds for
    `gnomad.bound` in a full report: it is given only when `--gnomad-max-af` was
    asked for (so the probe ran), and is null for a loaded store nobody consulted."""

    def examined(*keys: str) -> bool:
        return only is None or any(k in only for k in keys)

    def per_sample(name: str, key: str | None, *scoped: str) -> dict:
        """The share of rows per sample (the whole sample's: a fact about the
        data) and `bound`: the verdict of the gate(s) that read this field,
        `scoped` being the extra ones (AD is read by both min_ab and the
        parent-allele screen, and is bound only if every one of them is). `bound`
        holds only the gates the report judged."""
        out: dict = {s.role: _round(s.fraction(name)) for s in stats}
        judged = [k for k in (key, *scoped) if k and k in bound]
        if judged:
            out["bound"] = all(bound[k] for k in judged)
        elif only is not None:
            out["bound"] = None  # nobody examined it
        else:
            # not requested: bound in the sense of "the data has it"
            lowest, _ = _weakest(stats, name)
            out["bound"] = lowest is not None and lowest >= ENFORCEABLE_MIN
        return out

    return {
        "rows": {s.role: s.rows for s in stats},
        "gq": per_sample("gq", "min_gq"),
        "dp": per_sample("dp", "min_dp"),
        "ad": per_sample("ad", "min_ab", "parent_alt_screen"),
        "median_dp": {s.role: _num(s.median_dp) for s in stats},
        "info_af": {
            "state": gates.af_gate,
            "bound": bound.get("max_af", False) if examined("max_af") else None,
        },
        "filter_column": {
            "present": filter_present,
            "non_pass_variants": variants.non_pass,
        },
        "gnomad": {
            "loaded": gnomad_loaded,
            # asked for and probed, or null: a loaded store nobody asked to
            # consult (the key switched off, no --gnomad-max-af, a scoped run) was
            # not looked at, which is neither bound nor unbound
            "bound": gnomad_bound if "gnomad_max_af" in bound else None,
            "field": GNOMAD_FIELD,
            "version": None,
            "candidates_checked": None,
            "candidates_without_record": None,
            "cohort_variants_sampled": probe.sampled if probe else None,
            "cohort_variants_with_record": probe.found if probe else None,
        },
        "population_filter": (
            population if examined("max_af", "gnomad_max_af") else None
        ),
        "warnings": list(warnings),
    }


def _round(v: float | None) -> float | None:
    return None if v is None else round(v, 4)


# ───────────────────────────────── refusal ─────────────────────────────────

_POPULATION_REMEDY = (
    "load gnomAD for this cohort's chromosomes and genome build with "
    "`vcfclick annotations load-gnomad` and set --gnomad-max-af, or use "
    "broad-first-pass"
)
_REMEDY = {
    "gate_unbound": (
        "the VCF lacks the field: use flags or a preset that does not need it "
        "(for example broad-first-pass), or ingest a VCF that carries it"
    ),
    "gnomad_not_loaded": _POPULATION_REMEDY,
    "gnomad_no_coverage": _POPULATION_REMEDY,
    "no_population_frequency_filter": _POPULATION_REMEDY,
    "model_blocked": (
        "re-ingest the cohort with --keep-reference (a de novo, dominant or "
        "compound-het call needs the parents' stored 0/0), load gene annotations "
        "with `vcfclick annotations load` for compound het, or leave the model "
        "out with --category"
    ),
}
# the field each GQ / DP / AD gate reads, as the remedy names it
_GATE_FIELD = {
    "min_gq": "GQ",
    "min_dp": "DP",
    "min_ab": "AD",
    "max_ab": "AD",
    "parent_alt_screen": "AD",
}


def _unbound_remedy(errors: Sequence[dict]) -> str:
    """The remedy for gates whose field is missing, naming the fields: GQ, DP or
    AD absent where the gates read them (on the rows at the selected models'
    candidate sites, which is not the same as absent from the whole VCF)."""
    fields = list(
        dict.fromkeys(
            _GATE_FIELD[e["key"]]
            for e in errors
            if e["code"] == "gate_unbound" and e["key"] in _GATE_FIELD
        )
    )
    if not fields:
        return _REMEDY["gate_unbound"]
    many = len(fields) > 1
    names = f"{', '.join(fields[:-1])} and {fields[-1]}" if many else fields[0]
    it = "them" if many else "it"
    return (
        f"{names} {'are' if many else 'is'} missing where the gates read {it} (on "
        f"the rows at the candidate sites of the selected models): use flags or a "
        f"preset that does not need {it} (for example broad-first-pass), or ingest "
        f"a VCF that carries {it}"
    )


def refusal_message(errors: Sequence[dict], *, subject: str) -> str:
    """The text of the error that stops a run whose gates cannot bind."""
    needs = "every gate to bind"
    if any(e["code"] == "model_blocked" for e in errors):
        needs += " and every model it selects to run"
    lines = [f"{subject} needs {needs}, and these cannot in this data:"]
    lines += [f"  - {e['message']}" for e in errors]
    remedies = list(
        dict.fromkeys(
            _unbound_remedy(errors)
            if e["code"] == "gate_unbound"
            else _REMEDY[e["code"]]
            for e in errors
        )
    )
    lines += [f"Remedy: {r}." for r in remedies]
    lines.append(
        "To run anyway, pass --no-require-bound-gates: the warnings stay in the "
        "output, and a call whose field is missing passes the gate instead of "
        "failing it."
    )
    return "\n".join(lines)
