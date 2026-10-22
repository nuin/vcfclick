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
own), gnomAD when any is loaded. `--require-bound-gates` makes a gate that
cannot bind, and the absence of any population-frequency filter, an error;
`--strict-missing` makes that error for the GQ / DP / AD gates alone, because
it turns a missing field into a fail and a field missing everywhere would then
fail everything. `--no-require-bound-gates` typed on the command line overrides
both for one run.

`collect_sample_stats` is the only query (one pass over the three samples'
rows); `build_report` is pure so every threshold is a unit test. The numbers
below are the author's choices, named so a reviewer can change one place.
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import dataclass

from cli.trio_analysis import AF_GATE_SKIPPED, Gates
from cli.warning_codes import warning

# A gate binds when its field is present on at least this share of the rows.
ENFORCEABLE_MIN = 0.90
# The proband's median depth below this many times min_dp: the gate is close to
# the assay's own depth, so the preset may be too strict for it.
THIN_COVERAGE_FACTOR = 1.5

GNOMAD_FIELD = "popmax (grpmax else AF)"

# which warning codes a refusal is made of
_BOUND_CODES = ("gate_unbound", "gnomad_not_loaded", "no_population_frequency_filter")
_FORMAT_KEYS = ("min_gq", "min_dp", "min_ab")  # the GQ / DP / AD gates


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

    @property
    def af_fraction(self) -> float | None:
        return self.with_af / self.n if self.n else None


def collect_variant_stats(sess, ingest_id: str) -> VariantStats:
    from storage import count_expr, sql_quote_str

    sql = (
        f"SELECT {count_expr()}, "
        "sum(CASE WHEN info_AF IS NOT NULL THEN 1 ELSE 0 END) "
        f"FROM variants WHERE ingest_id = {sql_quote_str(ingest_id)}"
    )
    n, with_af = json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"][0]
    return VariantStats(int(n or 0), int(with_af or 0))


# ───────────────────────────────── the report ─────────────────────────────────


@dataclass(frozen=True)
class EffectiveReport:
    document: dict  # the JSON `effective` object (warnings included)
    warnings: tuple[dict, ...]
    bound: dict[str, bool]  # file key -> does the gate bind in this data
    errors: tuple[dict, ...] = ()  # the warnings a refusal is made of


def _pct(fraction: float) -> str:
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
    filter_present: bool = True,
    kind: str = "trio",
    categories: Sequence[str] = ("all",),
    assay: str | None = None,
    subject: str = "preset",
    skipped: Sequence[dict] = (),
    preset_max_af: float | None = None,
    enforce_bound: bool = False,
    enforce_format: bool = False,
) -> EffectiveReport:
    """Whether each requested gate binds, and what to say about it.

    `skipped` are the preset's skipped keys (from the echo); `preset_max_af` is
    the max_af a preset supplied, set only so a skip of it can be explained;
    `enforce_bound` / `enforce_format` choose which warnings are also errors
    (see the module docstring). `subject` names what the user chose, in the
    thin-coverage sentence."""
    warnings: list[dict] = []
    bound: dict[str, bool] = {}
    all_roles = [s.role for s in stats]
    parents = [r for r in all_roles if r in ("father", "mother")]

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
                    f"{shown} cannot be enforced: only {_pct(lowest)} of "
                    f"{'/'.join(below)} rows carry {label}",
                )
            )
        else:
            bound[key] = True
            if lowest < 1.0:
                outcome = "fail" if gates.strict_missing else "pass"
                warnings.append(
                    warning(
                        "gate_partial",
                        key,
                        f"{shown} enforced on {_pct(lowest)} of rows (a call "
                        f"with no {label} counts as a {outcome})",
                    )
                )

    ab_requested = gates.min_ab > 0 or gates.max_ab < 1
    ab_roles = None if gates.strict_missing else ["proband"]
    if gates.min_gq > 0:
        check_format("min_gq", "gq", "GQ", f"min_gq={gates.min_gq}", None)
    if gates.min_dp > 0:
        check_format("min_dp", "dp", "DP", f"min_dp={gates.min_dp}", None)
    if ab_requested:
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
                else f"only {_pct(lowest)} of {'/'.join(below)} rows carry AD"
            )
            warnings.append(
                warning(
                    "gate_unbound",
                    "parent_alt_screen",
                    f"parent_alt_screen cannot be enforced: {why}, so every "
                    "parent 0/0 would be trusted",
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
                else f"only {_pct(af_share)} of variants carry INFO/AF"
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
    if gnomad_max_af is not None:
        bound["gnomad_max_af"] = gnomad_loaded
        if not gnomad_loaded:
            warnings.append(
                warning(
                    "gnomad_not_loaded",
                    "gnomad_max_af",
                    f"gnomad_max_af={gnomad_max_af:g} cannot be enforced: no "
                    "gnomAD data is loaded (`vcfclick annotations load-gnomad`)",
                )
            )
    population = af_ok or (gnomad_max_af is not None and gnomad_loaded)
    if not population:
        warnings.append(
            warning(
                "no_population_frequency_filter",
                None,
                "NO population frequency filter applied: "
                + _why_no_population(
                    gates, af_run, af_share, gnomad_max_af, gnomad_loaded
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
    if kind == "duo" and any(c in ("all", "denovo") for c in categories):
        warnings.append(
            warning(
                "duo_denovo",
                "category",
                "with one parent the de novo answer is possible_denovo: weaker "
                "evidence than a trio's",
            )
        )
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

    errors = tuple(
        w
        for w in warnings
        if (enforce_bound and w["code"] in _BOUND_CODES)
        or (enforce_format and w["code"] == "gate_unbound" and w["key"] in _FORMAT_KEYS)
    )
    document = _document(
        stats,
        variants,
        gates,
        gnomad_loaded,
        filter_present,
        bound,
        population,
        warnings,
    )
    return EffectiveReport(document, tuple(warnings), bound, errors)


def _why_no_population(
    gates: Gates,
    af_run: bool,
    af_share: float | None,
    gnomad_max_af: float | None,
    gnomad_loaded: bool,
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
    else:
        gn = "no gnomAD data is loaded"
    return f"{af}, and {gn}"


def _document(
    stats: Sequence[SampleGateStats],
    variants: VariantStats,
    gates: Gates,
    gnomad_loaded: bool,
    filter_present: bool,
    bound: dict[str, bool],
    population: bool,
    warnings: list[dict],
) -> dict:
    def per_sample(name: str, key: str | None) -> dict:
        out: dict = {s.role: _round(s.fraction(name)) for s in stats}
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
        "ad": per_sample("ad", "min_ab"),
        "median_dp": {s.role: _num(s.median_dp) for s in stats},
        "info_af": {
            "state": gates.af_gate,
            "bound": bound.get("max_af", False),
        },
        "filter_column": {"present": filter_present},
        "gnomad": {
            "loaded": gnomad_loaded,
            "bound": gnomad_loaded,
            "field": GNOMAD_FIELD,
            "version": None,
            "candidates_without_record": None,
        },
        "population_filter": population,
        "warnings": list(warnings),
    }


def _round(v: float | None) -> float | None:
    return None if v is None else round(v, 4)


# ───────────────────────────────── refusal ─────────────────────────────────

_POPULATION_REMEDY = (
    "load gnomAD with `vcfclick annotations load-gnomad` and set "
    "--gnomad-max-af, or use broad-first-pass"
)
_REMEDY = {
    "gate_unbound": (
        "the VCF lacks the field: use flags or a preset that does not need it "
        "(for example broad-first-pass), or ingest a VCF that carries it"
    ),
    "gnomad_not_loaded": _POPULATION_REMEDY,
    "no_population_frequency_filter": _POPULATION_REMEDY,
}


def refusal_message(errors: Sequence[dict], *, subject: str) -> str:
    """The text of the error that stops a run whose gates cannot bind."""
    lines = [f"{subject} needs every gate to bind, and these cannot in this data:"]
    lines += [f"  - {e['message']}" for e in errors]
    remedies = list(dict.fromkeys(_REMEDY[e["code"]] for e in errors))
    lines += [f"Remedy: {r}." for r in remedies]
    lines.append(
        "To run anyway, pass --no-require-bound-gates (the warnings stay in the "
        "output)."
    )
    return "\n".join(lines)
