"""The SQL of `db segregate`: one query per sibship per model.

Candidate sites are the UNION of every affected child's alt rows, ungated, so
the result does not depend on child order or on one anchor child's quality;
the gates are applied per member in its state (cli/segregation_verdict.py
documents the nine states). The query pivots the members into one row per
site with conditional aggregation (`max(CASE WHEN sample_id = ... )`, the
idiom of cli/db_qc.py) and keeps the rows whose states satisfy the model's
pruning predicate, built from the same state sets the Python verdict uses.

Portable between chDB and DuckDB: gt+1 and 0/1 flags instead of NULLs (a
LEFT JOIN default is 0, not NULL, on chDB), counts written additively
(`n_called + n_miss < n_samples`, never `n_samples - n_called`, so UInt32
cannot wrap), every outer column explicitly aliased (the chDB analyzer
raises UNKNOWN_IDENTIFIER on an ORDER BY over an unaliased inner column).
"""

from __future__ import annotations

from dataclasses import dataclass

from cli import segregation_verdict as sv
from cli.trio_analysis import (
    PASS_ONLY_SQL,
    Gates,
    _het_ab_gate,
    _het_ad_present,
    _quality_gate,
    parent_ref_ok,
)


@dataclass(frozen=True)
class ChildInfo:
    sample: str
    affected: str | None  # affected | unaffected | None
    sex: str | None = None


@dataclass(frozen=True)
class Sibship:
    """Children with the same declared father and mother in one family."""

    ingest_id: str
    family_id: str | None
    father: str
    mother: str
    children: tuple[ChildInfo, ...]
    # the parents' declared phenotypes (affected | unaffected | None), read
    # only by --parent-phenotype use
    father_affected: str | None = None
    mother_affected: str | None = None

    @property
    def affected(self) -> tuple[ChildInfo, ...]:
        return tuple(c for c in self.children if c.affected == "affected")

    @property
    def unaffected(self) -> tuple[ChildInfo, ...]:
        return tuple(c for c in self.children if c.affected == "unaffected")

    @property
    def unknown(self) -> tuple[ChildInfo, ...]:
        return tuple(
            c for c in self.children if c.affected not in ("affected", "unaffected")
        )

    @property
    def members(self) -> tuple[str, ...]:
        return (self.father, self.mother, *(c.sample for c in self.children))


@dataclass(frozen=True)
class Caps:
    """What the database can account for. A database ingested before
    `missing_genotypes` or `variants.n_called` existed gets the reduced
    builder: state 3 and 7 never arise, every absent member is state 6."""

    have_miss: bool = True
    have_ncalled: bool = True
    have_filter: bool = True


@dataclass(frozen=True)
class SqlOptions:
    n_samples: int
    absent_called: bool = False
    max_unaffected_carriers: int = 0
    partial_strict: bool = True


def _q(s: str) -> str:
    from storage import sql_quote_str

    return sql_quote_str(s)


def _ids(samples) -> str:
    return ", ".join(_q(s) for s in samples)


def _in(col: str, states) -> str:
    return f"{col} IN ({', '.join(str(s) for s in sorted(states))})"


def _key(a: str, b: str) -> str:
    return " AND ".join(f"{a}.{c} = {b}.{c}" for c in ("chrom", "pos", "ref", "alt"))


def _scope_part(model: str) -> str:
    # compound-het inclusions act after the gene grouping
    return "exclude" if model == "comphet" else "all"


def _homalt_floor(alias: str, floor: float) -> str:
    from storage.sql import allele_fraction_sql, to_int64

    alt, ref = f"{alias}.ad_alt", f"{alias}.ad_ref"
    zero = f"({to_int64(alt)} + {to_int64(ref)}) = 0"
    frac = allele_fraction_sql(alt, ref)
    return (
        f"({alias}.gt != 2 OR {ref} IS NULL OR {alt} IS NULL OR {zero} "
        f"OR {frac} >= {float(floor)!r})"
    )


def cand_sql(sib: Sibship, model: str, scope=None) -> str:
    """The candidate sites: every site where some affected child has an alt
    row, off the sex and mitochondrial contigs (haploid calls are stored as
    hom-alt, so no model here can read them) and inside the scope."""
    from cli.scope import scope_sql
    from storage.contigs import NON_AUTOSOMAL
    from storage.sql import chrom_in_sql

    where = [
        f"a.ingest_id = {_q(sib.ingest_id)}",
        f"a.sample_id IN ({_ids(c.sample for c in sib.affected)})",
        "a.gt > 0",
        f"NOT {chrom_in_sql('a.chrom', sorted(NON_AUTOSOMAL))}",
    ]
    scoped = scope_sql("a", scope, part=_scope_part(model))
    if scoped:
        where.append(scoped)
    return (
        "SELECT DISTINCT a.chrom AS chrom, a.pos AS pos, a.ref AS ref, a.alt AS alt "
        f"FROM genotypes a WHERE {' AND '.join(where)}"
    )


def skipped_sites_sql(sib: Sibship) -> str:
    """How many candidate sites lie on contigs the models do not read."""
    from storage import count_expr
    from storage.contigs import NON_AUTOSOMAL
    from storage.sql import chrom_in_sql

    return (
        f"SELECT {count_expr()} FROM (SELECT DISTINCT a.chrom, a.pos, a.ref, a.alt "
        f"FROM genotypes a WHERE a.ingest_id = {_q(sib.ingest_id)} "
        f"AND a.sample_id IN ({_ids(c.sample for c in sib.affected)}) AND a.gt > 0 "
        f"AND {chrom_in_sql('a.chrom', sorted(NON_AUTOSOMAL))}) s"
    )


def partial_sites_sql(sib: Sibship, n_samples: int, scope=None) -> str:
    """How many candidate sites have a partial call (cohort-wide): a sample
    with one allele called and one missing. n_called counts the fully called
    samples and missing_genotypes the fully missing ones, so the rest are
    partial: n_called + n_missing < n_samples."""
    from storage import count_expr

    return (
        f"SELECT {count_expr()} FROM ({cand_sql(sib, 'dominant', scope)}) c "
        "INNER JOIN variants v ON "
        f"v.ingest_id = {_q(sib.ingest_id)} AND {_key('v', 'c')} "
        "LEFT JOIN (SELECT m.chrom AS chrom, m.pos AS pos, m.ref AS ref, "
        "m.alt AS alt, count(DISTINCT m.sample_id) AS n_miss "
        f"FROM missing_genotypes m WHERE m.ingest_id = {_q(sib.ingest_id)} "
        "GROUP BY m.chrom, m.pos, m.ref, m.alt) x ON "
        f"{_key('x', 'c')} "
        f"WHERE v.n_called IS NOT NULL AND v.n_called + coalesce(x.n_miss, 0) < {int(n_samples)}"
    )


def _flags(sib: Sibship, gates: Gates, opts: SqlOptions, model: str) -> tuple[str, str]:
    """(ok, ra) expressions on the genotype row `g`: does it pass the gates for
    its role, and (de novo) is it a 0/0 whose alt reads say otherwise."""
    quality = _quality_gate(
        "g",
        gates.min_gq,
        gates.min_dp,
        strict_gq=gates.strict_for("min_gq"),
        strict_dp=gates.strict_for("min_dp"),
    )
    kids = _ids(c.sample for c in sib.children)
    aff = _ids(c.sample for c in sib.affected)
    ok = f"{quality} AND (g.sample_id NOT IN ({kids}) OR {_het_ab_gate('g', gates.min_ab, gates.max_ab)})"
    ok += f" AND (g.sample_id NOT IN ({aff}) OR {_homalt_floor('g', gates.min_homalt_ab)})"
    if gates.ab_requested and gates.strict_for("min_ab"):
        ok += f" AND {_het_ad_present('g')}"
    ra = "0 = 1"
    if model == "denovo":
        ra = f"g.gt = 0 AND NOT {parent_ref_ok('g', gates)}"
    return ok, ra


def member_columns(sib: Sibship) -> list[tuple[str, str]]:
    """(column suffix, sample id) in output order: father, mother, children."""
    return [("f", sib.father), ("m", sib.mother)] + [
        (str(i), c.sample) for i, c in enumerate(sib.children)
    ]


def predicate(model: str, sib: Sibship, opts: SqlOptions) -> str:
    """The pruning predicate over the state columns: the parental pattern, no
    affected child contradicting, at least one affected child strictly
    satisfied, and (not compound het) unaffected contradictions within
    --max-unaffected-carriers. The Python verdict re-derives each row."""
    R = opts.absent_called
    ref = sv.reference_states(R)
    f, m = "st_f", "st_m"
    if model == "denovo":
        parents = f"{_in(f, (0, 7))} AND {_in(m, (0, 7))}"
    elif model == "recessive":
        parents = f"{f} = 1 AND {m} = 1"
    else:
        parents = (
            f"(({_in(f, sv.CARRIES)} AND {_in(m, ref)}) "
            f"OR ({_in(f, ref)} AND {_in(m, sv.CARRIES)}))"
        )
    parts = [parents]
    idx = {c.sample: str(i) for i, c in enumerate(sib.children)}
    aff = [f"st_{idx[c.sample]}" for c in sib.affected]
    parts += [_in(col, sv.AFFECTED_OK[model]) for col in aff]
    parts.append(
        "(" + " OR ".join(_in(col, sv.AFFECTED_STRICT[model]) for col in aff) + ")"
    )
    if model != "comphet" and sib.unaffected:
        contra = []
        for c in sib.unaffected:
            col = f"st_{idx[c.sample]}"
            if model == "recessive":
                expr = f"{col} = 2"
            else:
                expr = f"{col} = 2 OR ({col} = 1 AND pc = 0)"
                if model == "denovo":
                    expr += f" OR {col} = 8"
            contra.append(expr)
        if opts.max_unaffected_carriers <= 0:
            parts += [f"NOT ({e})" for e in contra]
        else:
            total = " + ".join(f"(CASE WHEN {e} THEN 1 ELSE 0 END)" for e in contra)
            parts.append(f"({total}) <= {int(opts.max_unaffected_carriers)}")
    return " AND ".join(parts)


def build_sql(
    model: str,
    sib: Sibship,
    gates: Gates,
    opts: SqlOptions,
    caps: Caps = Caps(),
    scope=None,
) -> str:
    """One query: a row per candidate site with the state of each member."""
    if model not in sv.MODELS:
        raise ValueError(f"unknown model {model!r}")
    members = member_columns(sib)
    ok, ra = _flags(sib, gates, opts, model)
    cand = cand_sql(sib, model, scope)
    ing = _q(sib.ingest_id)
    from storage.sql import to_int64

    fam = (
        "SELECT g.chrom AS chrom, g.pos AS pos, g.ref AS ref, g.alt AS alt, "
        "g.sample_id AS sample_id, g.gt AS gt, "
        f"coalesce({to_int64('g.ad_alt')}, 0) AS aa, "
        f"CASE WHEN {ok} THEN 1 ELSE 0 END AS ok, "
        f"CASE WHEN {ra} THEN 1 ELSE 0 END AS ra "
        f"FROM genotypes g INNER JOIN cand c ON {_key('c', 'g')} "
        f"WHERE g.ingest_id = {ing} AND g.sample_id IN ({_ids(s for _, s in members)})"
    )
    site_cols = []
    for sfx, sample in members:
        s = _q(sample)
        site_cols.append(
            f"max(CASE WHEN sample_id = {s} THEN gt + 1 ELSE 0 END) AS gp_{sfx}"
        )
        site_cols.append(
            f"max(CASE WHEN sample_id = {s} THEN ok ELSE 0 END) AS ok_{sfx}"
        )
        site_cols.append(
            f"max(CASE WHEN sample_id = {s} THEN ra ELSE 0 END) AS ra_{sfx}"
        )
        if sfx in ("f", "m"):
            site_cols.append(
                f"max(CASE WHEN sample_id = {s} THEN aa ELSE 0 END) AS aa_{sfx}"
            )
    site = (
        "SELECT chrom, pos, ref, alt, "
        + ", ".join(site_cols)
        + " FROM fam GROUP BY chrom, pos, ref, alt"
    )
    ctes = [f"cand AS ({cand})", f"fam AS ({fam})", f"site AS ({site})"]
    n = int(opts.n_samples)
    if caps.have_miss:
        miss_cols = ", ".join(
            f"max(CASE WHEN m.sample_id = {_q(sample)} THEN 1 ELSE 0 END) AS nc_{sfx}"
            for sfx, sample in members
        )
        ctes.append(
            "miss AS (SELECT m.chrom AS chrom, m.pos AS pos, m.ref AS ref, "
            f"m.alt AS alt, count(DISTINCT m.sample_id) AS n_miss, {miss_cols} "
            f"FROM missing_genotypes m INNER JOIN cand c ON {_key('c', 'm')} "
            f"WHERE m.ingest_id = {ing} GROUP BY m.chrom, m.pos, m.ref, m.alt)"
        )
        miss_join = f"LEFT JOIN miss x ON {_key('x', 's')}"
        n_miss = "coalesce(x.n_miss, 0)"
    else:
        miss_join, n_miss = "", "0"
    if caps.have_ncalled:
        complete = f"(v.n_called IS NOT NULL AND v.n_called + {n_miss} = {n})"
        partial = (
            f"(CASE WHEN v.n_called IS NOT NULL AND v.n_called + {n_miss} < {n} "
            "THEN 1 ELSE 0 END)"
        )
    else:
        complete, partial = "(1 = 0)", "0"
    if not opts.partial_strict:
        partial = "0"

    def state(sfx: str) -> str:
        nc = f"coalesce(x.nc_{sfx}, 0) = 1" if caps.have_miss else "1 = 0"
        return (
            f"CASE WHEN s.gp_{sfx} = 0 THEN (CASE WHEN {nc} THEN 4 "
            f"WHEN {complete} THEN 3 ELSE 6 END) "
            f"WHEN s.ok_{sfx} = 0 THEN 5 "
            f"WHEN s.gp_{sfx} = 1 THEN (CASE WHEN s.ra_{sfx} = 1 THEN 8 "
            f"WHEN {partial} = 1 THEN 7 ELSE 0 END) "
            f"WHEN s.gp_{sfx} = 2 THEN 1 ELSE 2 END"
        )

    inner_cols = [
        "s.chrom AS chrom",
        "s.pos AS pos",
        "s.ref AS ref",
        "s.alt AS alt",
        "v.info_AF AS af",
        f"{partial} AS pc",
    ]
    inner_cols += [f"{state(sfx)} AS st_{sfx}" for sfx, _ in members]
    outer_cols = ["chrom", "pos", "ref", "alt", "af", "pc"]
    outer_cols += [f"st_{sfx}" for sfx, _ in members]
    inner_cols += [f"s.gp_{sfx} AS gp_{sfx}" for sfx, _ in members]
    outer_cols += [f"gp_{sfx}" for sfx, _ in members]
    if model == "denovo":
        inner_cols += ["s.aa_f AS aa_f", "s.aa_m AS aa_m"]
        outer_cols += ["aa_f", "aa_m"]
    site_where = ""
    if gates.pass_only and caps.have_filter:
        site_where = f" WHERE {PASS_ONLY_SQL}"
    inner = (
        f"SELECT {', '.join(inner_cols)} FROM site s {miss_join} "
        f"LEFT JOIN variants v ON v.ingest_id = {ing} AND {_key('v', 's')}" + site_where
    )
    return (
        f"WITH {', '.join(ctes)} "
        f"SELECT {', '.join(outer_cols)} FROM ({inner}) t "
        f"WHERE {predicate(model, sib, opts)} ORDER BY chrom, pos, ref, alt"
    )
