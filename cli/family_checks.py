"""Checks every family analysis shares: does a sample have stored reference
calls, and is the VCF's INFO/AF the cohort's own?

Both are about what the stored data can and cannot prove, so they live
once here and each model asks them instead of inventing its own.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import dataclass


def _rows(sess, sql: str) -> list:
    return json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]


# ─────────────────────────── reference calls ───────────────────────────


def samples_with_reference(
    sess, ingest_id: str, samples: Iterable[str]
) -> dict[str, bool]:
    """For each sample, does it have at least one stored hom-reference
    (gt = 0) row in this ingestion?

    `genotypes` is sparse: unless the cohort was ingested with
    --keep-reference there are no gt = 0 rows, and a sample absent at a site
    is 0/0 or a no-call. A sample with no stored reference row can therefore
    never be proven 0/0 at any site.

    One statement for all the samples, however many (a 32-sample case is one
    query, not 32). A sample with no stored reference row makes a per-sample
    `LIMIT 1` scan all of that sample's rows before it can answer, so a
    sparse store paid one scan per sample. `DISTINCT ... LIMIT n` asks for at
    most the `n` samples listed, so an engine that can stops reading as soon
    as each has shown a reference row. The result keeps the order the
    samples were first listed in.
    """
    from storage import sql_quote_str

    wanted = list(dict.fromkeys(samples))
    if not wanted:
        return {}
    listed = ", ".join(sql_quote_str(s) for s in wanted)
    rows = _rows(
        sess,
        "SELECT DISTINCT sample_id FROM genotypes WHERE ingest_id = "
        f"{sql_quote_str(ingest_id)} AND sample_id IN ({listed}) AND gt = 0 "
        f"LIMIT {len(wanted)}",
    )
    have = {r[0] for r in rows}
    return {s: s in have for s in wanted}


def family_has_reference(sess, ingest_id: str, samples: Iterable[str]) -> bool:
    """True when every listed sample has stored reference calls (and at
    least one sample was listed). This is the per-family test the models
    that prove a parent is 0/0 use; the database-wide `_has_reference_rows`
    keeps driving only the four classic models' `blocked` flag."""
    have = samples_with_reference(sess, ingest_id, samples)
    return bool(have) and all(have.values())


# ───────────────────── cohort-derived allele frequency ─────────────────────

# A VCF's INFO/AF is meant to be a population frequency. Family VCFs from
# GATK-style joint callers instead carry the AF of the callset itself, which
# for a trio makes every private het 1/(2N) = 0.167: a --max-af 0.01 rarity
# gate then removes every candidate, silently. The probe recognises that
# case from the data: at how many sites does INFO/AF equal ac_called /
# an_called (the AF of this cohort, computed at ingest)?
#
# "Equal" has to be relative. The cohort's own AF of a rare site is tiny
# (1/(2N): 0.005 at N = 100, 0.0017 at N = 300), and a real population
# frequency of 0.0001 is still a different number from it. A fixed absolute
# tolerance of 0.005 would call every rare population AF "cohort-like" from
# about N = 100 on, and the rarity gate would then be skipped for a large
# cohort that carries real, common variants. So the tolerance scales with
# the value it is compared to; only the rounding a writer applies to
# printed values is allowed on top: a three-decimal AF (0.167 for 1/6) is
# off by up to 0.0005, but only values from 0.01 up are printed that way.
PROBE_MIN_SITES = 20  # fewer sites with AF and called counts: cannot judge
PROBE_FRACTION = 0.9  # at least this share of sites equals the cohort's AF
PROBE_REL_TOLERANCE = 0.01  # |INFO/AF - ac/an| <= this share of ac/an ...
PROBE_ROUNDING = 0.0005  # ... plus this, when ac/an >= PROBE_ROUNDING_FROM
PROBE_ROUNDING_FROM = 0.01  # below this, writers print significant digits


@dataclass(frozen=True)
class CohortAfProbe:
    """What the probe saw. `available` is False when it could not run."""

    available: bool
    n_sites: int = 0  # sites with INFO/AF and called counts
    n_cohort_like: int = 0  # of those, INFO/AF equals ac/an (see PROBE_* above)
    n_singleton_like: int = 0  # cohort-like sites with ac_called = 1 (AF = 1/an)
    reason: str = ""  # why it is unavailable

    @property
    def fraction(self) -> float:
        return self.n_cohort_like / self.n_sites if self.n_sites else 0.0

    @property
    def suspected(self) -> bool:
        return self.available and self.fraction >= PROBE_FRACTION


def cohort_like_sql() -> str:
    """The predicate "this row's INFO/AF is the cohort's own frequency",
    over the `variants` columns info_AF, ac_called and an_called (see the
    PROBE_* notes above for the tolerance)."""
    from storage.sql import to_double

    ratio = f"({to_double('ac_called')} / {to_double('an_called')})"
    tolerance = (
        f"({PROBE_REL_TOLERANCE!r} * {ratio} + CASE WHEN {ratio} >= "
        f"{PROBE_ROUNDING_FROM!r} THEN {PROBE_ROUNDING!r} ELSE 0.0 END)"
    )
    return f"abs(info_AF - {ratio}) <= {tolerance}"


def cohort_af_probe(sess, ingest_id: str) -> CohortAfProbe:
    """Probe whether INFO/AF in `ingest_id` is the cohort's own frequency.

    Unavailable (never an error) when the store predates the called-count
    columns, or fewer than PROBE_MIN_SITES sites carry both INFO/AF and
    called counts. Reads `variants` only; no genotype scan.
    """
    from storage import count_expr, sql_quote_str, table_columns

    cols = set(table_columns("variants"))
    if not {"ac_called", "an_called"} <= cols:
        return CohortAfProbe(False, reason="no called counts (older database)")

    like = cohort_like_sql()
    sql = (
        f"SELECT {count_expr()}, "
        f"sum(CASE WHEN {like} THEN 1 ELSE 0 END), "
        f"sum(CASE WHEN {like} AND ac_called = 1 THEN 1 ELSE 0 END) "
        f"FROM variants WHERE ingest_id = {sql_quote_str(ingest_id)} "
        "AND info_AF IS NOT NULL AND ac_called IS NOT NULL "
        "AND an_called IS NOT NULL AND an_called > 0"
    )
    n, n_like, n_single = (int(v or 0) for v in _rows(sess, sql)[0])
    if n < PROBE_MIN_SITES:
        return CohortAfProbe(
            False, n_sites=n, reason=f"only {n} sites carry INFO/AF and called counts"
        )
    return CohortAfProbe(True, n, n_like, n_single)


def resolve_af_gate(probe: CohortAfProbe, *, explicit: bool) -> str:
    """The `gates.af_gate` value. `explicit` is True when the user gave
    --max-af: an explicit threshold is always honoured, so a cohort-derived
    AF only changes the report, not the gate."""
    if not probe.available:
        return "unavailable"
    if not probe.suspected:
        return "applied"
    return "cohort_derived_suspected" if explicit else "skipped_cohort_derived"


def af_gate_message(probe: CohortAfProbe, state: str, max_af: float) -> str | None:
    """One sentence for a suspected cohort-derived AF (None for the other
    states): what the probe saw, what the gate did, what to do about it."""
    if state not in ("cohort_derived_suspected", "skipped_cohort_derived"):
        return None
    seen = (
        f"INFO/AF equals this cohort's own allele frequency at "
        f"{probe.n_cohort_like} of {probe.n_sites} sites "
        f"({probe.fraction:.0%}), so it is not a population frequency"
    )
    if state == "skipped_cohort_derived":
        return (
            f"{seen}: a private het in a small cohort has a high "
            f"AF, so the --max-af {max_af:g} gate was skipped. Pass --max-af "
            "to apply it anyway, or --gnomad-max-af for a population "
            "frequency."
        )
    return (
        f"{seen}; --max-af {max_af:g} is applied as requested and may "
        "remove real candidates (a private het in a small cohort has a high "
        "AF). Use --gnomad-max-af for a population frequency."
    )


def trio_af_gate(
    sess, ingest_id: str, max_af: float, *, explicit: bool
) -> tuple[str, str | None]:
    """Run the probe and resolve it: the `gates.af_gate` value for a trio
    analysis of `ingest_id`, and the sentence to show the user (None when
    there is nothing to say). `explicit` is True when the caller passed
    --max-af (`db trio`) or max_af (the web API); the CLI and the web UI
    both resolve the gate here so they cannot drift."""
    probe = cohort_af_probe(sess, ingest_id)
    state = resolve_af_gate(probe, explicit=explicit)
    return state, af_gate_message(probe, state, max_af)
