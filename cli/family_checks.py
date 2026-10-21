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
    never be proven 0/0 at any site. One cheap `LIMIT 1` per sample, so the
    cost does not depend on the size of the cohort.
    """
    from storage import sql_quote_str

    out: dict[str, bool] = {}
    for sample in dict.fromkeys(samples):
        rows = _rows(
            sess,
            "SELECT 1 FROM genotypes WHERE ingest_id = "
            f"{sql_quote_str(ingest_id)} AND sample_id = {sql_quote_str(sample)} "
            "AND gt = 0 LIMIT 1",
        )
        out[sample] = bool(rows)
    return out


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
PROBE_MIN_SITES = 20  # fewer sites with AF and called counts: cannot judge
PROBE_FRACTION = 0.9  # at least this share of sites equals the cohort's AF
PROBE_TOLERANCE = 0.005  # |INFO/AF - ac_called/an_called| within this


@dataclass(frozen=True)
class CohortAfProbe:
    """What the probe saw. `available` is False when it could not run."""

    available: bool
    n_sites: int = 0  # sites with INFO/AF and called counts
    n_cohort_like: int = 0  # of those, INFO/AF within tolerance of ac/an
    n_singleton_like: int = 0  # cohort-like sites with ac_called = 1 (AF = 1/an)
    reason: str = ""  # why it is unavailable

    @property
    def fraction(self) -> float:
        return self.n_cohort_like / self.n_sites if self.n_sites else 0.0

    @property
    def suspected(self) -> bool:
        return self.available and self.fraction >= PROBE_FRACTION


def cohort_af_probe(sess, ingest_id: str) -> CohortAfProbe:
    """Probe whether INFO/AF in `ingest_id` is the cohort's own frequency.

    Unavailable (never an error) when the store predates the called-count
    columns, or fewer than PROBE_MIN_SITES sites carry both INFO/AF and
    called counts. Reads `variants` only; no genotype scan.
    """
    from storage import count_expr, sql_quote_str, table_columns
    from storage.sql import to_double

    cols = set(table_columns("variants"))
    if not {"ac_called", "an_called"} <= cols:
        return CohortAfProbe(False, reason="no called counts (older database)")

    ratio = f"{to_double('ac_called')} / {to_double('an_called')}"
    like = f"abs(info_AF - ({ratio})) <= {PROBE_TOLERANCE!r}"
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


def af_gate_note(probe: CohortAfProbe, state: str, max_af: float) -> str | None:
    """The stderr / `note` line for a suspected cohort-derived AF, or None."""
    if state not in ("cohort_derived_suspected", "skipped_cohort_derived"):
        return None
    seen = (
        f"INFO/AF equals this cohort's own allele frequency at "
        f"{probe.n_cohort_like} of {probe.n_sites} sites "
        f"({probe.fraction:.0%}), so it is not a population frequency"
    )
    if state == "skipped_cohort_derived":
        return (
            f"note: {seen}: a private het in a small cohort has a high "
            f"AF, so the --max-af {max_af:g} gate was skipped. Pass --max-af "
            "to apply it anyway, or --gnomad-max-af for a population "
            "frequency."
        )
    return (
        f"note: {seen}; --max-af {max_af:g} is applied as requested and may "
        "remove real candidates (a private het in a small cohort has a high "
        "AF). Use --gnomad-max-af for a population frequency."
    )
