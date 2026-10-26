"""Inheritance queries and family resolution for trio analysis."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, NamedTuple

import click

from cli.filter_registry import DEFAULT_MAX_AF as DEFAULT_MAX_AF
from cli.filter_registry import DEFAULTS

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
# Sex- and region-aware models with their own row shape (cli/xlinked.py).
X_CATEGORIES: tuple[str, ...] = ("xlinked",)
CATEGORIES: tuple[str, ...] = ROW_CATEGORIES + GENE_CATEGORIES + X_CATEGORIES
CATEGORY_CHOICES: tuple[str, ...] = CATEGORIES + ("all",)

# DEFAULT_MAX_AF, the default for `db trio --max-af` (and the web API's
# max_af), is declared in cli/filter_registry.py with every other default.

# Why a model reports `blocked_reason` instead of results: the closed
# vocabulary of the JSON `blocked_reason` key (null when not blocked).
# Defined once here; each desktop app mirrors it in one reason -> text map.
# The four classic models use the first two; xlinked also uses needs_sex,
# no_x_data, parent_sex_conflict and build_conflict; the rest are reserved
# for the family models that are not implemented yet (phasing, ...).
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
    # The proband's family in the pedigree (None when not looked up). It
    # only labels the report; no query reads it.
    family_id: str | None = None

    @property
    def parents(self) -> tuple[str, ...]:
        """The declared parents, father first."""
        return tuple(p for p in (self.father, self.mother) if p is not None)

    @property
    def kind(self) -> str:
        """`trio` (both parents), `duo` (one) or `none`."""
        return {2: "trio", 1: "duo"}.get(len(self.parents), "none")


# The gates strict_missing turns into "fail when the field is missing", by file
# key: the GQ / DP / AD gates and the parent-allele screen.
STRICT_GATE_KEYS: tuple[str, ...] = ("min_gq", "min_dp", "min_ab", "parent_alt_screen")


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
    # when its allele balance stays under this (see parent_ref_ok). The
    # classic models read it only for the opt-in `parent_alt_screen`.
    max_parent_alt_frac: float = DEFAULTS["max_parent_alt_frac"]
    # See AF_GATE_VALUES. Only `skipped_cohort_derived` changes the SQL.
    af_gate: str = "applied"
    # Opt-in gates (cli/filter_registry.py). All off: the SQL is the text it
    # always was. strict_missing turns a NULL GQ / DP / AD into a fail,
    # pass_only keeps VCF FILTER = PASS (or unset), parent_alt_screen makes
    # the de novo model distrust a parent's 0/0 that has alt reads.
    strict_missing: bool = False
    pass_only: bool = False
    parent_alt_screen: bool = False
    # The gates strict_missing is NOT applied to in this run, by file key
    # (STRICT_GATE_KEYS). Set only by a run told to go
    # ahead anyway (`--no-require-bound-gates`) for a gate whose field the data
    # lacks: strict there would fail every call, and the run would report "no
    # candidates" for a reason that has nothing to do with the genome.
    strict_exempt: tuple[str, ...] = ()
    # The X-linked model (cli/xlinked.py). min_male_x_dp None: half min_dp.
    min_male_x_dp: int | None = None
    max_af_x: float = DEFAULTS["max_af_x"]
    min_hemi_ab: float = DEFAULTS["min_hemi_ab"]
    max_parent_alt_reads: int = DEFAULTS["max_parent_alt_reads"]

    def strict_for(self, key: str) -> bool:
        """Does strict_missing apply to the gate `key`?"""
        return self.strict_missing and key not in self.strict_exempt

    @property
    def ab_requested(self) -> bool:
        """Is the allele-balance band a gate at all? A band of 0..1 excludes
        nothing."""
        return self.min_ab > 0 or self.max_ab < 1


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


def _quality_gate(
    alias: str,
    min_gq: int,
    min_dp: int,
    *,
    strict_gq: bool = False,
    strict_dp: bool = False,
) -> str:
    """GQ and DP gate. Tolerant by default: a call with no GQ (or DP) passes.
    `strict_gq` / `strict_dp` (strict_missing) count a NULL as a fail, but only
    for a gate that is asked for (a minimum above 0): a field nobody asked
    about is not held against the call."""

    def one(column: str, minimum: int, strict: bool) -> str:
        if strict and minimum > 0:
            return f"({alias}.{column} IS NOT NULL AND {alias}.{column} >= {minimum})"
        return f"({alias}.{column} IS NULL OR {alias}.{column} >= {minimum})"

    return f"{one('gq', min_gq, strict_gq)} AND {one('dp', min_dp, strict_dp)}"


def _het_ab_gate(alias: str, min_ab: float, max_ab: float) -> str:
    # Allele balance for a het call: ad_alt / (ad_ref + ad_alt) should sit
    # near 0.5. Applies only to het rows (gt=1) with AD present.
    denom = f"({alias}.ad_ref + {alias}.ad_alt)"
    frac = f"({alias}.ad_alt * 1.0 / {denom})"
    return (
        f"({alias}.gt != 1 OR {alias}.ad_ref IS NULL OR {alias}.ad_alt IS NULL "
        f"OR {denom} = 0 OR {frac} BETWEEN {min_ab} AND {max_ab})"
    )


def _het_ad_present(alias: str) -> str:
    """strict_missing with an allele-balance band: a het call must carry allele
    depths. Only the proband is held to the band; this asks that no het call,
    the proband's or a parent's, has its depths missing."""
    return (
        f"({alias}.gt != 1 OR "
        f"({alias}.ad_ref IS NOT NULL AND {alias}.ad_alt IS NOT NULL))"
    )


# variants.filter holds NULL for both PASS and '.' (cyvcf2 reports neither),
# and the label otherwise. The literals cover a store written another way.
PASS_ONLY_SQL = "(v.filter IS NULL OR v.filter IN ('PASS', '.', ''))"


def parent_sample_joins(father: str | None, mother: str | None) -> str:
    """The parents' calls at the proband's sites: `f` and `m` joined to the
    proband's row `g`. Inner joins, so a site where a parent has no stored row
    is not a candidate."""
    from storage.sql import sample_join

    if father is None or mother is None:
        raise ValueError(
            "this model needs both parents declared (father="
            f"{father!r}, mother={mother!r})"
        )
    return f"{sample_join('f', father)} {sample_join('m', mother)}"


def parent_joins(father: str | None, mother: str | None) -> str:
    """`parent_sample_joins` and the site's variant `v` (population AF, FILTER)."""
    return (
        f"{parent_sample_joins(father, mother)} "
        "LEFT JOIN variants v ON v.ingest_id = g.ingest_id AND v.chrom = g.chrom "
        "AND v.pos = g.pos AND v.ref = g.ref AND v.alt = g.alt"
    )


# The genotype predicate of each model: the sites it reads (proband `g`,
# parents `f` and `m`), before any quality or rarity gate. One definition, used
# by the models' SQL below and by the effective-gate report, which counts what
# the rows at these sites carry (cli/effective_gates.collect_candidate_stats),
# so the two can never read different sites. Compound het takes the dominant
# pattern (a het inherited from exactly one parent) and groups it by gene.
MODEL_PREDICATES: dict[str, str] = {
    "denovo": "g.gt > 0 AND f.gt = 0 AND m.gt = 0",
    "recessive": "g.gt = 2 AND f.gt = 1 AND m.gt = 1",
    "dominant": "g.gt = 1 AND ((f.gt > 0 AND m.gt = 0) OR (f.gt = 0 AND m.gt > 0))",
}
MODEL_PREDICATES["comphet"] = MODEL_PREDICATES["dominant"]


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


def _where_parts(
    trio: Trio,
    gates: Gates,
    model: str,
    *,
    parent_alt_screen: bool = False,
    x_guard: str | None = None,
) -> list[tuple[str, str]]:
    """The predicates of the shared WHERE clause, in application order, each
    with its name: model, x_guard, quality_proband, quality_father,
    quality_mother, allele_balance, population_af, pass_only,
    parent_alt_screen. `x_guard` (cli/xlinked.x_guard_sql) drops the
    genotypes that cannot exist on non-PAR X. `_where`
    joins them; a diagnostic that wants to count what each removes (the
    funnel) walks them. A part is absent when its gate is off. (`trio` is
    unused today; the scope predicate of a later step needs it.)"""
    from storage.sql import af_gate_sql

    def quality(alias: str) -> str:
        return _quality_gate(
            alias,
            gates.min_gq,
            gates.min_dp,
            strict_gq=gates.strict_for("min_gq"),
            strict_dp=gates.strict_for("min_dp"),
        )

    allele_balance = _het_ab_gate("g", gates.min_ab, gates.max_ab)
    if gates.ab_requested and gates.strict_for("min_ab"):
        allele_balance += " AND " + " AND ".join(
            _het_ad_present(alias) for alias in ("g", "f", "m")
        )
    parts = [("model", model)]
    if x_guard:
        parts.append(("x_guard", x_guard))
    parts += [
        ("quality_proband", quality("g")),
        ("quality_father", quality("f")),
        ("quality_mother", quality("m")),
        ("allele_balance", allele_balance),
    ]
    if gates.af_gate != AF_GATE_SKIPPED:
        parts.append(("population_af", af_gate_sql("v.info_AF", gates.max_af)))
    if gates.pass_only:
        parts.append(("pass_only", PASS_ONLY_SQL))
    if parent_alt_screen:
        screen = f"{parent_ref_ok('f', gates)} AND {parent_ref_ok('m', gates)}"
        if gates.strict_for("parent_alt_screen"):
            # the screen reads the parents' allele depths: a 0/0 without them
            # cannot be screened, and strict_missing counts that as a fail
            present = " AND ".join(
                f"({alias}.ad_ref IS NOT NULL AND {alias}.ad_alt IS NOT NULL)"
                for alias in ("f", "m")
            )
            screen = f"{present} AND {screen}"
        parts.append(("parent_alt_screen", screen))
    return parts


# The parts of `_where_parts` that decide WHICH SITES a model lists, not whether
# a call is trusted: the rarity and FILTER gates. (GQ, DP, allele balance and the
# parent-allele screen are the quality gates; the gnomAD filter is not SQL at all,
# it runs on the rows this query returns.)
SITE_GATE_PARTS: tuple[str, ...] = ("population_af", "pass_only")


def site_gate_sql(
    trio: Trio, gates: Gates, model: str, *, filter_present: bool = True
) -> list[str]:
    """The SQL predicates (on the variant `v`) that remove a site before any
    quality gate is looked at: INFO/AF above `max_af` and, with `pass_only`, a
    FILTER label other than PASS. The effective-gate report counts what the rows at
    a model's candidate sites carry, and the sites these two remove are not
    candidates: a field missing there says nothing about the rows a clinician is
    shown. `filter_present` False (a variants table without the column) leaves out
    the FILTER gate, which the report says it cannot enforce."""
    parts = dict(_where_parts(trio, gates, model))
    return [
        parts[name]
        for name in SITE_GATE_PARTS
        if name in parts and (filter_present or name != "pass_only")
    ]


def _where(
    trio: Trio,
    gates: Gates,
    model: str,
    *,
    parent_alt_screen: bool = False,
    x_guard: str | None = None,
) -> str:
    """Shared WHERE clause: the genotype model plus the quality / rarity
    gates, scoped to this proband and ingest."""
    from storage import sql_quote_str

    parts = _where_parts(
        trio, gates, model, parent_alt_screen=parent_alt_screen, x_guard=x_guard
    )
    predicate = " AND ".join(sql for _, sql in parts)
    return (
        f"g.ingest_id = {sql_quote_str(trio.ingest_id)} "
        f"AND g.sample_id = {sql_quote_str(trio.proband)} AND {predicate}"
    )


def trio_sql(
    category: str,
    trio: Trio,
    gates: Gates,
    *,
    count_only: bool,
    x_guard: str | None = None,
) -> str:
    """Build the SQL for one inheritance model. The proband row is `g`;
    parents are joined as `f`/`m`; `v` brings population AF. `x_guard` is
    the X guard predicate (None: the text it always was)."""
    from storage import count_expr

    if category not in ROW_CATEGORIES:
        raise ValueError(f"unknown category {category!r}")
    model = MODEL_PREDICATES[category]

    joins = parent_joins(trio.father, trio.mother)
    # parent_alt_screen is about a parent's ABSENCE, so only the model that
    # claims a variant is new in the child reads it.
    where = _where(
        trio,
        gates,
        model,
        parent_alt_screen=gates.parent_alt_screen and category == "denovo",
        x_guard=x_guard,
    )
    if count_only:
        return f"SELECT {count_expr()} FROM genotypes g {joins} WHERE {where}"
    return (
        "SELECT g.chrom, g.pos, g.ref, g.alt, g.gt AS proband_gt, "
        "f.gt AS father_gt, m.gt AS mother_gt, v.info_AF AS af "
        f"FROM genotypes g {joins} WHERE {where} ORDER BY g.chrom, g.pos"
    )


def strict_dropped_sql(
    category: str,
    trio: Trio,
    gates: Gates,
    *,
    count_only: bool,
    x_guard: str | None = None,
) -> str:
    """The candidates of `category` that pass every gate of a run WITHOUT
    strict_missing and fail only because strict_missing counts a missing GQ, DP or
    allele depth as a fail: the rows of the tolerant query that the strict gates
    reject. (The strict gates are a subset of the tolerant ones: a NULL fails
    instead of passing, a value is judged the same, and each strict predicate is a
    plain boolean, never NULL, so the NOT below loses no row.) A candidate here is
    a variant, for compound het too."""
    from storage import count_expr

    model = MODEL_PREDICATES[category]
    screen = gates.parent_alt_screen and category == "denovo"
    tolerant = gates._replace(strict_exempt=STRICT_GATE_KEYS)
    changed = [
        sql
        for name, sql in _where_parts(
            trio, gates, model, parent_alt_screen=screen, x_guard=x_guard
        )
        if name.startswith("quality") or name in ("allele_balance", "parent_alt_screen")
    ]
    select = count_expr() if count_only else "g.chrom, g.pos, g.ref, g.alt"
    return (
        f"SELECT {select} FROM genotypes g "
        f"{parent_joins(trio.father, trio.mother)} "
        f"WHERE {_where(trio, tolerant, model, parent_alt_screen=screen, x_guard=x_guard)} "
        f"AND NOT ({' AND '.join(changed)})"
    )


def _comphet_sql(
    trio: Trio,
    gates: Gates,
    *,
    x_guard: str | None = None,
    x_only: str | None = None,
) -> str:
    """Candidate variants for compound-het: each is a rare proband het
    inherited from exactly one parent (the dominant pattern), tagged with
    parent-of-origin. The gene grouping happens in Python — genes live in
    the annotation store, which can't be SQL-joined to the cohort.

    `x_only`, an NPX predicate on `g` (cli/xchrom.npx_sql), makes this the
    X model's female compound het: NPX sites only, a paternal het only from
    a hemizygous-ALT father (cli/xlinked.comphet_x_model)."""
    model = MODEL_PREDICATES["comphet"]
    if x_only is not None:
        from cli.xlinked import comphet_x_model

        model = f"{comphet_x_model()} AND {x_only}"
    where = _where(trio, gates, model, x_guard=x_guard)
    origin = (
        "CASE WHEN f.gt > 0 AND m.gt = 0 THEN 'paternal' "
        "WHEN m.gt > 0 AND f.gt = 0 THEN 'maternal' END AS origin"
    )
    return (
        f"SELECT g.chrom, g.pos, g.ref, g.alt, v.info_AF AS af, {origin} "
        f"FROM genotypes g {parent_joins(trio.father, trio.mother)} "
        f"WHERE {where} ORDER BY g.chrom, g.pos"
    )


def _comphet_genes(candidate_rows: list) -> dict:
    """Group origin-tagged candidate variants by gene, keeping only genes
    that carry BOTH a paternal and a maternal het (trans configuration →
    both gene copies hit). A variant overlapping several genes counts for
    each. Returns {gene_symbol: {"paternal": [...], "maternal": [...]}}."""
    from annotations.db import GeneIndex

    genes: dict = {}
    index = None  # read the gene table once, at the first candidate
    for chrom, pos, ref, alt, af, origin in candidate_rows:
        if origin not in ("paternal", "maternal"):
            continue
        if index is None:
            index = GeneIndex()
        for gr in index.at(chrom, int(pos)):
            entry = genes.setdefault(gr.gene_symbol, {"paternal": [], "maternal": []})
            entry[origin].append((chrom, pos, ref, alt, af))
    return {sym: e for sym, e in genes.items() if e["paternal"] and e["maternal"]}


def _resolve_parents(sess, ingest_id: str, proband: str) -> tuple[str, str]:
    """The proband's declared parents. Deliberately its own narrow query
    (not the Pedigree graph): it keeps resolving a proband whatever else is
    wrong with the pedigree, such as a cycle that loaded years ago."""
    from cli.pedigree_load import NO_PEDIGREE, NoPedigree
    from storage import sql_quote_str, table_exists

    if not table_exists("pedigree"):
        raise NoPedigree(NO_PEDIGREE)
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


def _gnomad_filter(rows: list, max_af: float, stats: dict | None = None) -> list:
    """The rows (chrom, pos, ref, alt, ...) whose gnomAD popmax AF is
    <= max_af, or that are absent from the loaded gnomAD slice: absence is
    treated as rare (the slice may simply not cover the locus), never as
    AF 0. One batched lookup for all rows, not one per row: with the AF
    gate skipped on a family VCF, a trio has thousands of candidates.

    `stats`, when given, collects every variant the filter looked at in
    `stats["gnomad_seen"]` and those it kept for want of a record in
    `stats["gnomad_absent"]` (sets of (chrom, pos, ref, alt)), so a report can
    say how many candidates the filter could not judge, out of how many."""
    from annotations import gnomad_af_many

    seen = None if stats is None else stats.setdefault("gnomad_seen", set())
    absent = None if stats is None else stats.setdefault("gnomad_absent", set())
    if seen is not None:
        seen.update((r[0], int(r[1]), r[2], r[3]) for r in rows)
    found = gnomad_af_many((row[0], row[1], row[2], row[3]) for row in rows)
    if not found:
        # nothing loaded, or none of the candidates in it
        if absent is not None:
            absent.update((r[0], int(r[1]), r[2], r[3]) for r in rows)
        return list(rows)
    keep = []
    for row in rows:
        key = (row[0], int(row[1]), row[2], row[3])
        g = found.get(key)
        if g is None and absent is not None:
            absent.add(key)
        if g is None or g.popmax is None or g.popmax <= max_af:
            keep.append(row)
    return keep


@dataclass(frozen=True)
class TrioAnalysis:
    """Query a resolved trio using the same gates for every output format."""

    session: Any
    trio: Trio
    gates: Gates
    has_reference: bool
    gnomad_max_af: float | None
    # The proband's siblings as (sample, kind) pairs, kind being `full`,
    # `half_paternal` or `half_maternal`. Reported, never queried.
    siblings: tuple[tuple[str, str], ...] = ()
    # Facts about this run that a report reads back: the preset echo and the
    # effective-gate report `db trio` computed up front, the variants the
    # gnomAD filter kept for want of a record. A deliberate exception to the
    # frozen dataclass (a mutable dict, not part of equality), written by the
    # command and by the filter that has the fact, read by the document and
    # text decorators.
    stats: dict = field(default_factory=dict, compare=False)
    # The chrX handling (cli/xchrom.py, cli/xlinked.py): the --build option
    # (auto | GRCh37 | GRCh38) and whether the X guard is on.
    build: str = "auto"
    x_guard: bool = False
    # Answers computed once per analysis (the X context, its rows), private.
    _cache: dict = field(default_factory=dict, compare=False, repr=False)

    def _gnomad(self, rows: list) -> list:
        if self.gnomad_max_af is None:
            return rows
        return _gnomad_filter(rows, self.gnomad_max_af, self.stats)

    def _query(self, sql: str) -> list:
        return json.loads(self.session.query(sql, "JSONCompact").bytes().decode())[
            "data"
        ]

    # ─────────────────────────── chrX ───────────────────────────

    def x_context(self):
        """The X model's context (cli/xlinked.XContext): build, sexes, why it
        is blocked. Computed once; a trio without both parents has none."""
        if "x_context" not in self._cache:
            from cli.xlinked import XContext, build_context

            if self.trio.father is None or self.trio.mother is None:
                ctx = XContext(has_x=False, blocked_reason="needs_parents")
            else:
                ctx = build_context(self.session, self.trio, self.build)
            self._cache["x_context"] = ctx
        return self._cache["x_context"]

    def x_guard_sql(self) -> str | None:
        """The X guard predicate for the classic models, or None (guard off,
        no chrX data): a het father on sure non-PAR X, and a het proband
        there when the proband is male."""
        if not self.x_guard:
            return None
        ctx = self.x_context()
        if not ctx.has_x:
            return None
        from cli.xlinked import x_guard_sql

        sexes = ctx.sexes
        return x_guard_sql(
            ctx.region_build,
            father=sexes.father_sex != "female",
            proband_male=sexes.proband_sex == "male",
        )

    def x_rows(self) -> list:
        """The X-linked candidates (cli/xlinked.XRow), gnomAD-filtered, empty
        when the context blocks the model."""
        if "x_rows" not in self._cache:
            from cli.xlinked import decode_rows, xlinked_sql

            ctx = self.x_context()
            rows: list = []
            if ctx.has_x and ctx.blocked_reason is None:
                sql = xlinked_sql(
                    self.trio, self.gates, ctx.proband_sex, ctx.region_build
                )
                rows = self._gnomad(decode_rows(self._query(sql), ctx.affected))
            self._cache["x_rows"] = rows
        return self._cache["x_rows"]

    def x_excluded(self) -> dict:
        """`context.excluded` of the X model (cli/xlinked.excluded_counts)."""
        if "x_excluded" not in self._cache:
            from cli.xlinked import excluded_counts

            self._cache["x_excluded"] = excluded_counts(
                self.session, self.trio, self.x_context()
            )
        return self._cache["x_excluded"]

    def x_genes(self) -> dict:
        """Female compound het on non-PAR X: {gene: {paternal, maternal}}."""
        if "x_genes" not in self._cache:
            from cli.xchrom import npx_sql

            ctx = self.x_context()
            genes: dict = {}
            if (
                ctx.has_x
                and ctx.blocked_reason is None
                and ctx.proband_sex == "female"
                and self.has_reference
            ):
                where = npx_sql(ctx.region_build, chrom="g.chrom", pos="g.pos")
                sql = _comphet_sql(self.trio, self.x_gates(), x_only=where)
                genes = _comphet_genes(self._gnomad(self._query(sql)))
            self._cache["x_genes"] = genes
        return self._cache["x_genes"]

    def x_gates(self) -> Gates:
        """The gates with the X rarity cap: max_af = min(max_af, max_af_x)."""
        from cli.xlinked import x_max_af

        return self.gates._replace(max_af=x_max_af(self.gates))

    def x_guard_excluded(self, category: str) -> int | None:
        """How many candidates (genes for comphet) the X guard removed from a
        classic model. 0 with no chrX data; None when the guard is off or the
        gnomAD filter is set (the count would need every row)."""
        if not self.x_guard:
            return None
        guard = self.x_guard_sql()
        if guard is None:
            return 0
        if self.gnomad_max_af is not None:
            return None
        if category == "comphet":
            plain = _comphet_genes(self._query(_comphet_sql(self.trio, self.gates)))
            return len(plain) - len(self.comphet_genes())
        sql = trio_sql(category, self.trio, self.gates, count_only=True)
        return int(self._query(sql)[0][0]) - int(self.count(category))

    # ─────────────────────── classic models ───────────────────────

    def detail_rows(self, category: str) -> list:
        sql = trio_sql(
            category,
            self.trio,
            self.gates,
            count_only=False,
            x_guard=self.x_guard_sql(),
        )
        return self._gnomad(self._query(sql))

    def count(self, category: str) -> int:
        if self.gnomad_max_af is not None:
            return len(self.detail_rows(category))
        sql = trio_sql(
            category, self.trio, self.gates, count_only=True, x_guard=self.x_guard_sql()
        )
        return self._query(sql)[0][0]

    def dropped_by_strict(self, category: str) -> int:
        """How many candidates of this model pass every gate (the gnomAD filter
        too) and are dropped only because strict_missing counts a missing GQ, DP
        or allele depth as a fail: what a run without strict_missing would list
        beyond this one. Variants, for compound het as for the others. 0 when no
        strict gate is in force, and for a model without a classic genotype
        predicate (xlinked)."""
        gates = self.gates
        if category not in MODEL_PREDICATES:
            return 0
        if not any(gates.strict_for(k) for k in STRICT_GATE_KEYS):
            return 0
        guard = self.x_guard_sql()
        if self.gnomad_max_af is None:
            sql = strict_dropped_sql(
                category, self.trio, gates, count_only=True, x_guard=guard
            )
            return int(self._query(sql)[0][0])
        sql = strict_dropped_sql(
            category, self.trio, gates, count_only=False, x_guard=guard
        )
        return len(_gnomad_filter(self._query(sql), self.gnomad_max_af))

    def comphet_genes(self) -> dict:
        sql = _comphet_sql(self.trio, self.gates, x_guard=self.x_guard_sql())
        return _comphet_genes(self._gnomad(self._query(sql)))
