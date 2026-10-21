"""Inheritance queries and family resolution for trio analysis."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, NamedTuple

import click


# The inheritance models, in the order a trio document lists them. This is
# the one definition: the CLI's --category choices, the web UI's category
# list and the MCP briefing all derive from it (tests/test_trio_categories.py
# ties them together).
ROW_CATEGORIES: tuple[str, ...] = (
    "denovo",
    "recessive",
    "dominant",
)  # one row per variant
GENE_CATEGORIES: tuple[str, ...] = ("comphet",)  # reported per gene
CATEGORIES: tuple[str, ...] = ROW_CATEGORIES + GENE_CATEGORIES
CATEGORY_CHOICES: tuple[str, ...] = CATEGORIES + ("all",)

# The default for `db trio --max-af` (and the web API's max_af).
DEFAULT_MAX_AF = 0.01

# Why a model reports `blocked_reason` instead of results: the closed
# vocabulary of the JSON `blocked_reason` key (null when not blocked).
# Defined once here; each desktop app mirrors it in one reason -> text map.
# The four classic models use the first two; the rest are reserved for the
# models that land later in the family program.
BLOCKED_REASONS: tuple[str, ...] = (
    "needs_keep_reference",  # no stored hom-ref rows to prove a 0/0 parent
    "needs_genes",  # gene annotations not loaded (comphet)
    "needs_phase",
    "needs_parents",
    "needs_sex",
    "needs_ad",
    "needs_parent_ad",
    "no_autosomes",
    "no_x_data",
    "parent_sex_conflict",
    "build_conflict",
)


# What happened to the --max-af gate, reported as `gates.af_gate`:
#   applied                   the gate ran; the cohort-AF probe found no sign
#                             the VCF's INFO/AF is the cohort's own
#   cohort_derived_suspected  INFO/AF looks cohort-derived but the gate ran
#                             because --max-af was given explicitly
#   skipped_cohort_derived    INFO/AF looks cohort-derived and --max-af was
#                             not given, so the gate was not applied
#   unavailable               the probe could not run (called counts missing
#                             or too few sites to judge); the gate ran
AF_GATE_VALUES: tuple[str, ...] = (
    "applied",
    "cohort_derived_suspected",
    "skipped_cohort_derived",
    "unavailable",
)
AF_GATE_SKIPPED = "skipped_cohort_derived"


def check_blocked_reason(reason: str | None) -> str | None:
    """Return `reason` if it is null or in the closed enum, else raise."""
    if reason is not None and reason not in BLOCKED_REASONS:
        raise ValueError(f"unknown blocked_reason {reason!r}")
    return reason


class Trio(NamedTuple):
    """A resolved family: the cohort, the proband and the declared parents.

    A parent is None when the pedigree declares none (a duo has one parent
    and only the analyses written for it accept that). The classic models
    need both and refuse a trio without them."""

    ingest_id: str
    proband: str
    father: str | None
    mother: str | None

    @property
    def parents(self) -> tuple[str, ...]:
        """The declared parents, father first."""
        return tuple(p for p in (self.father, self.mother) if p is not None)

    @property
    def kind(self) -> str:
        """`trio` (both parents), `duo` (one) or `none`."""
        return {2: "trio", 1: "duo"}.get(len(self.parents), "none")


class Gates(NamedTuple):
    """The tunable quality / rarity thresholds, passed together.

    Only ever gains trailing defaulted fields, so positional construction
    (cli/db_trio.py, vcfclick_web/app.py) keeps working."""

    min_gq: int
    min_dp: int
    max_af: float
    min_ab: float
    max_ab: float
    # A parent's stored 0/0 row is trusted as "parent does not carry it" only
    # when its allele balance stays under this (see parent_ref_ok). The four
    # classic models do not read it.
    max_parent_alt_frac: float = 0.05
    # See AF_GATE_VALUES. Only `skipped_cohort_derived` changes the SQL.
    af_gate: str = "applied"


def _sole_ingest_id(name: str) -> str | None:
    """Return the db's ingest_id if it has exactly one ingestion, else
    None. Lets `db ped` infer the target when unambiguous."""
    from storage import get_session

    sess = get_session(name)
    raw = (
        sess.query("SELECT DISTINCT ingest_id FROM ingestions FORMAT TabSeparated")
        .bytes()
        .decode()
    )
    ids = [s for s in raw.splitlines() if s.strip()]
    return ids[0] if len(ids) == 1 else None


def _quality_gate(alias: str, min_gq: int, min_dp: int) -> str:
    return (
        f"({alias}.gq IS NULL OR {alias}.gq >= {min_gq}) "
        f"AND ({alias}.dp IS NULL OR {alias}.dp >= {min_dp})"
    )


def _het_ab_gate(alias: str, min_ab: float, max_ab: float) -> str:
    # Allele balance for a het call: ad_alt / (ad_ref + ad_alt) should sit
    # near 0.5. Applies only to het rows (gt=1) with AD present.
    denom = f"({alias}.ad_ref + {alias}.ad_alt)"
    frac = f"({alias}.ad_alt * 1.0 / {denom})"
    return (
        f"({alias}.gt != 1 OR {alias}.ad_ref IS NULL OR {alias}.ad_alt IS NULL "
        f"OR {denom} = 0 OR {frac} BETWEEN {min_ab} AND {max_ab})"
    )


def _parent_joins(father: str | None, mother: str | None) -> str:
    from storage.sql import sample_join

    if father is None or mother is None:
        raise ValueError(
            "this model needs both parents declared (father="
            f"{father!r}, mother={mother!r})"
        )
    return (
        f"{sample_join('f', father)} {sample_join('m', mother)} "
        "LEFT JOIN variants v ON v.ingest_id = g.ingest_id AND v.chrom = g.chrom "
        "AND v.pos = g.pos AND v.ref = g.ref AND v.alt = g.alt"
    )


def parent_ref_ok(alias: str, gates: Gates) -> str:
    """SQL boolean: a parent's STORED 0/0 row at this site can be trusted as
    "the parent does not carry the allele" (alias = the parent's genotypes
    row). It is not trusted when the reads say otherwise, since a caller can
    emit 0/0 for a site with real alt support.

        ad_alt IS NULL                    no allele depths recorded: trusted
        ad_alt < 2                        zero or one alt read: noise, trusted
        ad_ref IS NOT NULL AND
          alt / (ref + alt) < frac        otherwise the alt fraction must stay
                                          strictly below gates.max_parent_alt_frac

    A NULL ad_ref with two or more alt reads is NOT trusted: the fraction
    cannot be computed, and two alt reads are the least evidence we accept as
    contradicting a 0/0. The expression is never NULL (a NULL would silently
    drop the row), so it can be negated or combined freely. The four classic
    models do not use it; the models that prove a parent's absence do.
    """
    frac = gates.max_parent_alt_frac
    if not 0.0 <= frac <= 1.0:
        raise ValueError(f"max_parent_alt_frac must be in [0, 1], got {frac!r}")
    from storage.sql import allele_fraction_sql

    ratio = allele_fraction_sql(f"{alias}.ad_alt", f"{alias}.ad_ref")
    return (
        f"({alias}.ad_alt IS NULL OR {alias}.ad_alt < 2 OR "
        f"({alias}.ad_ref IS NOT NULL AND {ratio} < {frac!r}))"
    )


def _where(trio: Trio, gates: Gates, model: str) -> str:
    """Shared WHERE clause: the genotype model plus the quality / rarity
    gates, scoped to this proband and ingest."""
    from storage import sql_quote_str
    from storage.sql import af_gate_sql

    parts = [
        model,
        _quality_gate("g", gates.min_gq, gates.min_dp),
        _quality_gate("f", gates.min_gq, gates.min_dp),
        _quality_gate("m", gates.min_gq, gates.min_dp),
        _het_ab_gate("g", gates.min_ab, gates.max_ab),
    ]
    if gates.af_gate != AF_GATE_SKIPPED:
        parts.append(af_gate_sql("v.info_AF", gates.max_af))
    predicate = " AND ".join(parts)
    return (
        f"g.ingest_id = {sql_quote_str(trio.ingest_id)} "
        f"AND g.sample_id = {sql_quote_str(trio.proband)} AND {predicate}"
    )


def trio_sql(category: str, trio: Trio, gates: Gates, *, count_only: bool) -> str:
    """Build the SQL for one inheritance model. The proband row is `g`;
    parents are joined as `f`/`m`; `v` brings population AF."""
    from storage import count_expr

    if category == "denovo":
        model = "g.gt > 0 AND f.gt = 0 AND m.gt = 0"
    elif category == "recessive":
        model = "g.gt = 2 AND f.gt = 1 AND m.gt = 1"
    elif category == "dominant":
        model = "g.gt = 1 AND ((f.gt > 0 AND m.gt = 0) OR (f.gt = 0 AND m.gt > 0))"
    else:
        raise ValueError(f"unknown category {category!r}")

    joins = _parent_joins(trio.father, trio.mother)
    where = _where(trio, gates, model)
    if count_only:
        return f"SELECT {count_expr()} FROM genotypes g {joins} WHERE {where}"
    return (
        "SELECT g.chrom, g.pos, g.ref, g.alt, g.gt AS proband_gt, "
        "f.gt AS father_gt, m.gt AS mother_gt, v.info_AF AS af "
        f"FROM genotypes g {joins} WHERE {where} ORDER BY g.chrom, g.pos"
    )


def _comphet_sql(trio: Trio, gates: Gates) -> str:
    """Candidate variants for compound-het: each is a rare proband het
    inherited from exactly one parent (the dominant pattern), tagged with
    parent-of-origin. The gene grouping happens in Python — genes live in
    the annotation store, which can't be SQL-joined to the cohort."""
    model = "g.gt = 1 AND ((f.gt > 0 AND m.gt = 0) OR (f.gt = 0 AND m.gt > 0))"
    where = _where(trio, gates, model)
    origin = (
        "CASE WHEN f.gt > 0 AND m.gt = 0 THEN 'paternal' "
        "WHEN m.gt > 0 AND f.gt = 0 THEN 'maternal' END AS origin"
    )
    return (
        f"SELECT g.chrom, g.pos, g.ref, g.alt, v.info_AF AS af, {origin} "
        f"FROM genotypes g {_parent_joins(trio.father, trio.mother)} "
        f"WHERE {where} ORDER BY g.chrom, g.pos"
    )


def _comphet_genes(candidate_rows: list) -> dict:
    """Group origin-tagged candidate variants by gene, keeping only genes
    that carry BOTH a paternal and a maternal het (trans configuration →
    both gene copies hit). A variant overlapping several genes counts for
    each. Returns {gene_symbol: {"paternal": [...], "maternal": [...]}}."""
    from annotations import gene_at

    genes: dict = {}
    for chrom, pos, ref, alt, af, origin in candidate_rows:
        if origin not in ("paternal", "maternal"):
            continue
        for gr in gene_at(chrom, int(pos)):
            entry = genes.setdefault(gr.gene_symbol, {"paternal": [], "maternal": []})
            entry[origin].append((chrom, pos, ref, alt, af))
    return {sym: e for sym, e in genes.items() if e["paternal"] and e["maternal"]}


def _resolve_parents(sess, ingest_id: str, proband: str) -> tuple[str, str]:
    from storage import sql_quote_str

    raw = (
        sess.query(
            f"SELECT father_id, mother_id FROM pedigree WHERE ingest_id = "
            f"{sql_quote_str(ingest_id)} AND sample_id = {sql_quote_str(proband)} "
            "FORMAT JSONCompact"
        )
        .bytes()
        .decode()
    )
    data = json.loads(raw)["data"]
    if not data:
        raise click.ClickException(
            f"no pedigree entry for proband {proband!r}. Load a PED first: "
            f"vcfclick db ped <name> <file.ped>"
        )
    father, mother = data[0]
    if father in ("0", "", None) or mother in ("0", "", None):
        raise click.ClickException(
            f"proband {proband!r} is missing a parent in the pedigree "
            f"(father={father!r}, mother={mother!r}); trio analysis needs both."
        )
    return father, mother


def _children(sess, ingest_id: str) -> list[tuple[str, str, str]]:
    """Every (child, father, mother) in the pedigree with both parents
    declared, in sample-id order: the probands `db trio` runs without
    --proband."""
    from storage import sql_quote_str

    raw = (
        sess.query(
            "SELECT sample_id, father_id, mother_id FROM pedigree "
            f"WHERE ingest_id = {sql_quote_str(ingest_id)} "
            "AND father_id IS NOT NULL AND father_id NOT IN ('0', '') "
            "AND mother_id IS NOT NULL AND mother_id NOT IN ('0', '') "
            "ORDER BY sample_id FORMAT JSONCompact"
        )
        .bytes()
        .decode()
    )
    return [(c, f, m) for c, f, m in json.loads(raw)["data"]]


def _has_reference_rows(sess) -> bool:
    from storage import count_expr

    raw = (
        sess.query(
            f"SELECT {count_expr()} FROM genotypes WHERE gt = 0 FORMAT JSONCompact"
        )
        .bytes()
        .decode()
    )
    return int(json.loads(raw)["data"][0][0]) > 0


def _gnomad_keep(chrom: str, pos, ref: str, alt: str, max_af: float) -> bool:
    """Keep a candidate whose gnomAD popmax AF is <= max_af, or that is
    absent from the loaded gnomAD slice — absence is treated as rare
    (the slice may simply not cover the locus), never as AF 0."""
    from annotations import gnomad_af

    g = gnomad_af(chrom, int(pos), ref, alt)
    return g is None or g.popmax is None or g.popmax <= max_af


@dataclass(frozen=True)
class TrioAnalysis:
    """Query a resolved trio using the same gates for every output format."""

    session: Any
    trio: Trio
    gates: Gates
    has_reference: bool
    gnomad_max_af: float | None

    def _gnomad(self, rows: list) -> list:
        if self.gnomad_max_af is None:
            return rows
        return [row for row in rows if _gnomad_keep(*row[:4], self.gnomad_max_af)]

    def detail_rows(self, category: str) -> list:
        sql = trio_sql(category, self.trio, self.gates, count_only=False)
        rows = json.loads(self.session.query(sql, "JSONCompact").bytes().decode())[
            "data"
        ]
        return self._gnomad(rows)

    def count(self, category: str) -> int:
        if self.gnomad_max_af is not None:
            return len(self.detail_rows(category))
        sql = trio_sql(category, self.trio, self.gates, count_only=True)
        return json.loads(self.session.query(sql, "JSONCompact").bytes().decode())[
            "data"
        ][0][0]

    def comphet_genes(self) -> dict:
        sql = _comphet_sql(self.trio, self.gates)
        rows = json.loads(self.session.query(sql, "JSONCompact").bytes().decode())[
            "data"
        ]
        return _comphet_genes(self._gnomad(rows))
