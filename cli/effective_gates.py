"""The effective-gate report: does each gate actually bind in THIS data?

The quality gates are NULL-tolerant (a call with no GQ passes), a VCF's INFO/AF
may be the cohort's own frequency, and gnomAD may not be loaded. A preset
named "strict" can therefore be a no-op, and used to be silent about it. This
report says, per gate, whether the field it reads exists in the data, and turns
that into warnings (printed on stdout, carried in the JSON `effective`
object) and, on request, into a refusal.

A gate BINDS when its field is present on at least ENFORCEABLE_MIN of the rows
it reads: GQ, DP and AD on the stored rows of the proband and both parents, the
INFO/AF gate on the variants (and its probe not saying the AF is the cohort's
own), gnomAD when it is loaded, has records on the chromosomes that hold
ENFORCEABLE_MIN of the cohort's variants (a slice on other chromosomes, or
under another naming, judges nothing) AND has a record for at least
GNOMAD_MIN_HIT_RATE of a sample of the cohort's own variants (a sparse test
load, a slice of one gene, another genome build judge almost nothing, whatever
chromosomes they have a row on). `--require-bound-gates` makes a gate
that cannot bind, and the absence of any population-frequency filter, an error;
`--strict-missing` makes that error for the GQ / DP / AD gates (and the
parent-allele screen) alone, because it turns a missing field into a fail and a
field missing everywhere would then fail everything. `--no-require-bound-gates`
typed on the command line overrides both for one run, and for that run
strict_missing is not applied to the gates that cannot bind (`Gates.strict_exempt`):
they are NULL-tolerant, as their warnings say, rather than failing every call.

A strict_missing gate is judged only when it is asked for (min_gq > 0, min_dp
> 0, an allele-balance band narrower than 0..1, the parent-allele screen on):
an unasked-for field is not held against a call, so it cannot make a gate
"unbound" either.

The gnomAD probe is a fact about the COHORT, known before the run, so a
run that insists (`require_bound`, `strict-denovo`) can refuse on it. What is
left for after the run is a fact about the CANDIDATES: a de novo candidate is
novel by nature, so one that gnomAD has no record for is the expected case and
not a reason to refuse; `gnomad_record_warning` says it when it is all of them.

`collect_sample_stats` is the only query (one pass over the three samples'
rows); `build_report` is pure so every threshold is a unit test. The numbers
below are the author's choices, named so a reviewer can change one place.
"""

from __future__ import annotations

import json
import math
from collections.abc import Collection, Sequence
from dataclasses import dataclass

from cli.trio_analysis import AF_GATE_SKIPPED, Gates
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
    "gnomad_not_loaded",
    "gnomad_no_coverage",
    "no_population_frequency_filter",
)
# the gates strict_missing turns into "fail when the field is missing": the
# GQ / DP / AD gates and the parent-allele screen
STRICT_KEYS = ("min_gq", "min_dp", "min_ab", "parent_alt_screen")
_NAMES_SHOWN = 4  # chromosome names listed in a message


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

    def fraction(self, name: str) -> float | None:
        """Share of this sample's rows that carry `name` (gq | dp | ad); None
        for a sample with no stored rows, which says nothing either way."""
        return getattr(self, name) / self.rows if self.rows else None

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
    """For each (sample, role): row count, rows carrying GQ / DP / AD, and the
    DP histogram. One statement over the samples' stored rows, however many
    samples there are."""
    from storage import sql_quote_str

    listed = ", ".join(sql_quote_str(s) for s, _ in samples)
    sql = (
        "SELECT sample_id, dp, count(*), "
        "sum(CASE WHEN gq IS NOT NULL THEN 1 ELSE 0 END), "
        "sum(CASE WHEN ad_ref IS NOT NULL AND ad_alt IS NOT NULL THEN 1 ELSE 0 END) "
        f"FROM genotypes WHERE ingest_id = {sql_quote_str(ingest_id)} "
        f"AND sample_id IN ({listed}) GROUP BY sample_id, dp"
    )
    data = json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]
    acc = {s: [0, 0, 0, 0, {}] for s, _ in samples}  # rows, gq, dp, ad, hist
    for sample, dp, n, n_gq, n_ad in data:
        a = acc[sample]
        n, n_gq, n_ad = int(n), int(n_gq or 0), int(n_ad or 0)
        a[0] += n
        a[1] += n_gq
        a[3] += n_ad
        if dp is not None:
            a[2] += n
            a[4][int(dp)] = a[4].get(int(dp), 0) + n
    return [
        SampleGateStats(s, role, a[0], a[1], a[2], a[3], tuple(sorted(a[4].items())))
        for (s, role) in samples
        for a in (acc[s],)
    ]


@dataclass(frozen=True)
class VariantStats:
    """How many of the cohort's variants carry INFO/AF."""

    n: int
    with_af: int
    # (chrom, variants on it), most variants first: where gnomAD must have
    # records to judge anything
    chroms: tuple[tuple[str, int], ...] = ()

    @property
    def af_fraction(self) -> float | None:
        return self.with_af / self.n if self.n else None


def collect_variant_stats(sess, ingest_id: str) -> VariantStats:
    """One pass over the cohort's variants: how many, how many carry INFO/AF,
    and how many sit on each chromosome."""
    from storage import count_expr, sql_quote_str

    sql = (
        f"SELECT chrom, {count_expr()}, "
        "sum(CASE WHEN info_AF IS NOT NULL THEN 1 ELSE 0 END) "
        f"FROM variants WHERE ingest_id = {sql_quote_str(ingest_id)} GROUP BY chrom"
    )
    data = json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]
    per = sorted(((c, int(n or 0)) for c, n, _ in data), key=lambda t: (-t[1], t[0]))
    return VariantStats(
        sum(n for _, n in per), sum(int(w or 0) for _, _, w in data), tuple(per)
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
    cover the cohort). `only` limits the report to these gate keys: a run that
    asked for one opt-in gate, and nothing else, is told about that gate only.
    The warnings are those about these keys, and in the document every `bound`
    (and `population_filter`) of a gate outside them is null: it was not
    examined, which is not the same as bound.
    `skipped` are the preset's skipped keys (from the echo); `preset_max_af` is
    the max_af a preset supplied, set only so a skip of it can be explained;
    `enforce_bound` / `enforce_format` choose which warnings are also errors
    (see the module docstring). `subject` names what the user chose, in the
    thin-coverage sentence."""
    warnings: list[dict] = []
    bound: dict[str, bool] = {}
    all_roles = [s.role for s in stats]
    parents = [r for r in all_roles if r in ("father", "mother")]

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

    def check_format(key: str, field_name: str, label: str, shown: str, roles):
        """One GQ / DP / AD gate: bound, partial or unbound."""
        lowest, below = _weakest(stats, field_name, roles)
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
                    f"{'/'.join(below)} rows carry {label}" + strict_note(key, label),
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
                        f"{shown} enforced on {pct(lowest)} of rows (a call "
                        f"with no {label} counts as a {outcome})",
                    )
                )

    ab_roles = None if gates.strict_for("min_ab") else ["proband"]
    if gates.min_gq > 0:
        check_format("min_gq", "gq", "GQ", f"min_gq={gates.min_gq}", None)
    if gates.min_dp > 0:
        check_format("min_dp", "dp", "DP", f"min_dp={gates.min_dp}", None)
    if gates.ab_requested:
        check_format(
            "min_ab",
            "ad",
            "AD",
            f"min_ab={gates.min_ab:g}/max_ab={gates.max_ab:g}",
            ab_roles,
        )
        bound["max_ab"] = bound["min_ab"]

    if gates.parent_alt_screen:
        lowest, below = _weakest(stats, "ad", parents)
        ok = lowest is not None and lowest >= ENFORCEABLE_MIN
        bound["parent_alt_screen"] = ok
        if not ok:
            why = (
                "no stored rows for the parents"
                if lowest is None
                else f"only {pct(lowest)} of {'/'.join(below)} rows carry AD"
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
        bound["pass_only"] = filter_present

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
    the cohort's sampled variants gnomAD could judge (the pre-run probe): when
    it is known the filter does work on this cohort and the candidates are most
    likely novel variants, which is what a de novo candidate is; when it is not,
    a slice that does not reach them or another build is the suspect."""
    if not checked or without_record < checked:
        return None
    one = checked == 1
    tail = (
        ""
        if info_af_bound
        else ", and INFO/AF is not applied either, so no population frequency "
        "filter took effect"
    )
    if cohort_rate is None:
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
    scoped run does not probe) is neither bound nor unbound."""

    def examined(*keys: str) -> bool:
        return only is None or any(k in only for k in keys)

    def per_sample(name: str, key: str | None, *scoped: str) -> dict:
        out: dict = {s.role: _round(s.fraction(name)) for s in stats}
        if only is not None:
            # `scoped` are the extra gates whose verdict this field decides
            hit = next((k for k in (key, *scoped) if k and k in only), None)
            out["bound"] = bound.get(hit) if hit else None
            return out
        out["bound"] = bound.get(key) if key else None
        if out["bound"] is None:
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
        "filter_column": {"present": filter_present},
        "gnomad": {
            "loaded": gnomad_loaded,
            "bound": gnomad_bound if examined("gnomad_max_af") else None,
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
}


def refusal_message(errors: Sequence[dict], *, subject: str) -> str:
    """The text of the error that stops a run whose gates cannot bind."""
    lines = [f"{subject} needs every gate to bind, and these cannot in this data:"]
    lines += [f"  - {e['message']}" for e in errors]
    remedies = list(dict.fromkeys(_REMEDY[e["code"]] for e in errors))
    lines += [f"Remedy: {r}." for r in remedies]
    lines.append(
        "To run anyway, pass --no-require-bound-gates: the warnings stay in the "
        "output, and a call whose field is missing passes the gate instead of "
        "failing it."
    )
    return "\n".join(lines)
