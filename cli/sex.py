"""Sex from chrX, and the sex each member of a trio is analysed as.

Two pieces of evidence, both read from non-reference calls:

  * the chrX het fraction: hets / non-reference calls. A male is
    hemizygous on non-PAR X, so he has (almost) no hets there. The
    thresholds are the ones `db qc` has always used.
  * X transmission, for a child: at non-PAR sites where the father is
    hemizygous ALT (stored gt 2) and the mother has a stored 0/0, a
    daughter must be het (she has her father's X) and a son must be 0/0
    (his X is his mother's). It is immune to long runs of homozygosity,
    which can make a daughter's het fraction look male.

`resolve_trio_sexes` combines them with the pedigree. No import of cli.main
(cli.db_qc imports this module and is itself imported by cli.main).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

from cli.warning_codes import warning

# Below this chrX het fraction a sample reads as male (hemizygous: few
# hets); at or above the upper bound, female; between, ambiguous.
MALE_MAX_HETFRAC = 0.10
FEMALE_MIN_HETFRAC = 0.20
MIN_CHRX_CALLS = 20  # too few chrX calls to infer sex

# X transmission (a child): informative sites needed, and the share of them
# at which the child is het that reads female / male. Provisional values,
# calibrated only on the synthetic fixture.
MIN_TRANSMISSION_SITES = 5
TRANSMISSION_FEMALE_MIN = 0.75
TRANSMISSION_MALE_MAX = 0.25

SEXES = ("male", "female")


def infer_sex(chrx_het: int, chrx_calls: int) -> tuple[str, float | None]:
    """(male | female | ambiguous | unknown, het fraction) from the chrX het
    and non-reference call counts."""
    if chrx_calls < MIN_CHRX_CALLS:
        return "unknown", None
    frac = chrx_het / chrx_calls
    if frac <= MALE_MAX_HETFRAC:
        return "male", frac
    if frac >= FEMALE_MIN_HETFRAC:
        return "female", frac
    return "ambiguous", frac


def transmission_call(
    informative: int, child_het: int
) -> tuple[str | None, float | None]:
    """(male | female | None, share of informative sites where the child is
    het). None when there are too few sites, or the share is in between."""
    if informative < MIN_TRANSMISSION_SITES:
        return None, None
    frac = child_het / informative
    if frac >= TRANSMISSION_FEMALE_MIN:
        return "female", frac
    if frac <= TRANSMISSION_MALE_MAX:
        return "male", frac
    return None, frac


# ─────────────────────────────── SQL ───────────────────────────────


def inference_sql(ingest_id: str, samples, build: str) -> str:
    """Per sample: (sample_id, chrX hets, chrX non-reference calls) on sure
    non-PAR X of `build`."""
    from cli.xchrom import npx_sql
    from storage.sql import sql_quote_str, str_list

    return (
        "SELECT sample_id, sum(CASE WHEN gt = 1 THEN 1 ELSE 0 END), "
        "sum(CASE WHEN gt > 0 THEN 1 ELSE 0 END) FROM genotypes "
        f"WHERE ingest_id = {sql_quote_str(ingest_id)} "
        f"AND sample_id IN ({str_list(samples)}) AND gt > 0 AND {npx_sql(build)} "
        "GROUP BY sample_id"
    )


def transmission_sql(
    ingest_id: str, child: str, father: str, mother: str, build: str
) -> str:
    """(informative sites, child het there): non-PAR sites with the father
    stored hemizygous ALT (gt 2), the mother stored 0/0 and the child stored
    0/0 or het."""
    from cli.xchrom import npx_sql
    from storage.sql import genotype_pivot_sql

    pivot = genotype_pivot_sql(
        ingest_id, {"p": child, "f": father, "m": mother}, where=npx_sql(build)
    )
    informative = "x.f_gt = 2 AND x.m_gt = 0 AND x.p_gt IN (0, 1)"
    return (
        f"SELECT sum(CASE WHEN {informative} THEN 1 ELSE 0 END), "
        f"sum(CASE WHEN {informative} AND x.p_gt = 1 THEN 1 ELSE 0 END) "
        f"FROM ({pivot}) x"
    )


def _rows(session, sql: str) -> list:
    return json.loads(session.query(sql, "JSONCompact").bytes().decode())["data"]


def chrx_counts(session, ingest_id: str, samples, build: str) -> dict:
    """{sample: (hets, non-reference calls)} on non-PAR X; a sample with no
    such call is absent."""
    return {
        s: (int(het or 0), int(calls or 0))
        for s, het, calls in _rows(session, inference_sql(ingest_id, samples, build))
    }


def transmission(
    session, ingest_id: str, child: str, father: str, mother: str, build: str
) -> tuple[int, int]:
    data = _rows(session, transmission_sql(ingest_id, child, father, mother, build))
    if not data:
        return 0, 0
    informative, het = data[0]
    return int(informative or 0), int(het or 0)


# ───────────────────────────── resolution ─────────────────────────────


@dataclass(frozen=True)
class SexEvidence:
    """What is known about one member's sex."""

    declared: str | None  # pedigree.sex male/female, else None
    inferred: str  # infer_sex: male | female | ambiguous | unknown
    het_frac: float | None
    # a child only: the X transmission call and its share
    transmitted: str | None = None
    transmitted_frac: float | None = None


@dataclass
class TrioSexes:
    proband_sex: str | None = None
    proband_sex_source: str = "unknown"  # pedigree | inferred | unknown
    father_sex: str | None = None
    father_sex_source: str = "unknown"  # pedigree | inferred | role | unknown
    mother_sex: str | None = None
    mother_sex_source: str = "unknown"
    warnings: list[dict] = field(default_factory=list)
    blocked_reason: str | None = None
    # one line for the text output when a parent's sex blocks the model
    conflict: str | None = None


def _frac(f: float | None) -> str:
    return "n/a" if f is None else f"{f:.2f}"


def _parent(
    out: TrioSexes, role: str, sample: str, ev: SexEvidence
) -> tuple[str | None, str]:
    """(sex, source) of a parent, recording a conflict with the role."""
    role_sex = "male" if role == "father" else "female"
    confident = ev.inferred if ev.inferred in SEXES else None
    if ev.declared in SEXES:
        sex, source = ev.declared, "pedigree"
        if confident and confident != ev.declared:
            out.warnings.append(
                warning(
                    "sex_conflict",
                    None,
                    f"{sample}: pedigree={ev.declared} but chrX het fraction="
                    f"{_frac(ev.het_frac)} (possible swap)",
                )
            )
    elif confident:
        sex, source = confident, "inferred"
    else:
        sex, source = role_sex, "role"
    if ev.inferred == "ambiguous":
        out.warnings.append(
            warning(
                "sex_ambiguous",
                None,
                f"{sample}: chrX het fraction={_frac(ev.het_frac)} is between the "
                "male and female thresholds (possible XXY, Turner or low depth)",
            )
        )
    wrong = [
        f"{label} {s}"
        for label, s in (("declared", ev.declared), ("chrX says", confident))
        if s in SEXES and s != role_sex
    ]
    if wrong and out.blocked_reason is None:
        out.blocked_reason = "parent_sex_conflict"
        out.conflict = f"{sample} is the {role}, " + ", ".join(wrong)
    return sex, source


def resolve_trio_sexes(
    proband: str,
    father: str,
    mother: str,
    evidence: dict[str, SexEvidence],
) -> TrioSexes:
    """The sex each member is analysed as, in priority order: the pedigree,
    a confident chrX inference (a child's checked against X transmission),
    for a parent the role, else unknown. A parent whose declared or
    confidently inferred sex contradicts the role blocks the model
    (`parent_sex_conflict`): the father's hemizygosity is what the patterns
    rest on. An unresolved proband is `needs_sex`."""
    out = TrioSexes()
    out.father_sex, out.father_sex_source = _parent(
        out, "father", father, evidence[father]
    )
    out.mother_sex, out.mother_sex_source = _parent(
        out, "mother", mother, evidence[mother]
    )
    ev = evidence[proband]
    confident = ev.inferred if ev.inferred in SEXES else None
    if ev.declared in SEXES:
        out.proband_sex, out.proband_sex_source = ev.declared, "pedigree"
        if confident and confident != ev.declared:
            out.warnings.append(
                warning(
                    "sex_conflict",
                    None,
                    f"{proband}: pedigree={ev.declared} but chrX het fraction="
                    f"{_frac(ev.het_frac)} (possible swap)",
                )
            )
    elif confident and ev.transmitted not in (None, confident):
        out.warnings.append(
            warning(
                "sex_ambiguous",
                None,
                f"{proband}: het fraction says {confident}, X transmission says "
                f"{ev.transmitted} (long homozygosity?)",
            )
        )
    elif confident:
        out.proband_sex, out.proband_sex_source = confident, "inferred"
    elif ev.inferred == "ambiguous":
        out.warnings.append(
            warning(
                "sex_ambiguous",
                None,
                f"{proband}: chrX het fraction={_frac(ev.het_frac)} is between the "
                "male and female thresholds (possible XXY, Turner or low depth)",
            )
        )
    if out.proband_sex is None and out.blocked_reason is None:
        out.blocked_reason = "needs_sex"
    return out
