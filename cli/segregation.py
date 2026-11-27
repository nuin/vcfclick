"""Sibling segregation: sibship resolution and the per-model analysis.

`db segregate` reads one sibship (children with the same declared father and
mother) and reports, per model, the variants the affected children carry and
the unaffected children do not contradict, with every member's state at the
site (cli/segregation_verdict.py) and a status ladder: consistent >
consistent_inferred > compatible > reduced_penetrance. The SQL is in
cli/segregation_sql.py, the text and JSON in cli/segregation_output.py, the
command in cli/db_segregate.py; the web endpoint and the MCP tool call
`segregation_document` directly.

Definitions:

  * sibship: children of one family with the same two declared parents, both
    parents and the children in the data. Needs >= 2 children and >= 1
    affected child, else the question is a trio's (`db trio`).
  * a child's phenotype is affected, unaffected or unknown; unknown children
    are shown and never constrain a model. Half siblings are separate sibships.
  * a sibship has at most MAX_CHILDREN children (the SQL carries one column
    set per member).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import click

from cli import segregation_verdict as sv
from cli.segregation_sql import (
    Caps,
    ChildInfo,
    Sibship,
    SqlOptions,
    build_sql,
    partial_sites_sql,
    skipped_sites_sql,
)
from cli.trio_analysis import AF_GATE_SKIPPED, Gates, _gnomad_filter
from storage import query_rows
from storage.pedigree_graph import Pedigree
from storage.sql import AF_EPSILON

MAX_CHILDREN = 12
MAX_PAIRS_PER_GENE = 500
ABSENT_CHOICES = ("strict", "called")
# The models of this command: the four of `db trio` that read autosomes. The
# X-linked model is not here: sex contigs are skipped (see segregation_sql).
MODELS = sv.MODELS
NO_ACCOUNTING_NOTE = (
    "this database has no complete missing-call accounting (no "
    "missing_genotypes table or variants.n_called); re-ingest to enable "
    "--absent called, and the partial-call check"
)


@dataclass(frozen=True)
class SegOptions:
    """What a run asks for, resolved (flag > preset > default)."""

    models: tuple[str, ...]
    gates: Gates
    gnomad_max_af: float | None = None
    absent: str = "strict"
    partial_calls: str = "strict"
    allow_unresolved: bool = False
    max_unaffected_carriers: int = 0
    parent_phenotype: str = "ignore"
    mz_pairs: tuple[tuple[str, str], ...] = ()
    limit: int = 1000

    @property
    def absent_called(self) -> bool:
        return self.absent == "called"

    def policy(self, sib: Sibship) -> sv.Policy:
        return sv.Policy(
            absent_called=self.absent_called,
            max_unaffected_carriers=self.max_unaffected_carriers,
            parent_phenotype_use=self.parent_phenotype == "use",
            groups=mz_groups(self.mz_pairs, sib),
        )


def select_models(category: Sequence[str] | str) -> tuple[str, ...]:
    """The models `--category` names, in document order. `all` anywhere is all
    of them; `xlinked` is not a model of this command."""
    tokens = (category,) if isinstance(category, str) else tuple(category)
    if "xlinked" in tokens:
        raise click.UsageError(
            "db segregate has no xlinked model: sex and mitochondrial contigs "
            "are skipped (a haploid male call is stored as hom-alt). Use "
            "`db trio --category xlinked` for a proband."
        )
    if "all" in tokens:
        return MODELS
    unknown = [t for t in tokens if t not in MODELS]
    if unknown:
        raise click.UsageError(
            f"unknown category {unknown[0]!r}; choose from {', '.join(MODELS)}, all"
        )
    return tuple(m for m in MODELS if m in tokens)


def mz_groups(pairs: Sequence[tuple[str, str]], sib: Sibship) -> dict[str, int]:
    """sample -> group id for the declared MZ pairs inside this sibship (a
    pair and a third member joined by another pair share one group)."""
    names = {c.sample for c in sib.children}
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in pairs:
        if a in names and b in names:
            parent[find(a)] = find(b)
    roots = {s: find(s) for s in parent}
    ids = {r: i for i, r in enumerate(sorted(set(roots.values())))}
    return {s: ids[r] for s, r in roots.items()}


# ─────────────────────────── sibship resolution ───────────────────────────


def _label(f: str, m: str, kids: Sequence[ChildInfo]) -> str:
    return f"{f} x {m}: {', '.join(c.sample for c in kids)}"


def find_sibships(
    graph: Pedigree,
    *,
    family: str | None = None,
    father: str | None = None,
    mother: str | None = None,
) -> tuple[list[Sibship], list[str]]:
    """The sibships to analyse and the notes for stderr. Without a selector,
    every eligible sibship; with one, exactly one must match (a ClickException
    lists the candidates of an ambiguous one)."""
    groups: dict[tuple[str, str, str], list[ChildInfo]] = {}
    for sid in graph.members():
        p = graph.person(sid)
        if p is None or not (p.father and p.mother) or not graph.is_genotyped(sid):
            continue
        groups.setdefault((p.family_id, p.father, p.mother), []).append(
            ChildInfo(sid, p.affected, p.sex)
        )

    def problem(key, kids) -> str | None:
        fam, f, m = key
        gone = [x for x in (f, m) if not graph.is_genotyped(x)]
        if gone:
            return f"parent {gone[0]!r} is not a sample in the data"
        if len(kids) < 2:
            return (
                f"one child ({kids[0].sample}): a trio, not a sibship; use "
                f"`vcfclick db trio --proband {kids[0].sample}`"
            )
        if not any(c.affected == "affected" for c in kids):
            return "no affected child (phenotype 2); set one with `vcfclick db ped`"
        if len(kids) > MAX_CHILDREN:
            return f"{len(kids)} children; at most {MAX_CHILDREN} are supported"
        return None

    selecting = any(x is not None for x in (family, father, mother))
    keys = sorted(groups)
    if selecting:
        keys = [
            k
            for k in keys
            if (family is None or k[0] == family)
            and (father is None or k[1] == father)
            and (mother is None or k[2] == mother)
        ]
        if not keys:
            raise click.ClickException(
                "no sibship matches the selection "
                f"(family={family}, father={father}, mother={mother}). "
                "Families with children: "
                + (", ".join(sorted({k[0] for k in groups})) or "(none)")
            )
        if len(keys) > 1:
            raise click.ClickException(
                "the selection matches several sibships (half siblings); add "
                "--father and --mother to pick one:\n  "
                + "\n  ".join(_label(k[1], k[2], groups[k]) for k in keys)
            )
        why = problem(keys[0], groups[keys[0]])
        if why:
            raise click.ClickException(f"family {keys[0][0]!r}: {why}")
    sibships: list[Sibship] = []
    notes: list[str] = []
    for key in keys:
        why = problem(key, groups[key])
        if why:
            notes.append(f"skipping family {key[0]!r} ({key[1]} x {key[2]}): {why}.")
            continue
        fam, f, m = key
        sibships.append(
            Sibship(
                ingest_id="",
                family_id=fam,
                father=f,
                mother=m,
                children=tuple(groups[key]),
                father_affected=graph.person(f).affected,
                mother_affected=graph.person(m).affected,
            )
        )
    if not sibships:
        raise click.ClickException(
            "no sibship with >=2 children, both parents declared and >=1 "
            "affected child. Load one with `vcfclick db ped`."
        )
    return sibships, notes


def check_mz_pairs(
    pairs: Sequence[tuple[str, str]], sibships: Sequence[Sibship]
) -> None:
    """Each declared MZ pair must be two children of one sibship with the
    same phenotype."""
    for a, b in pairs:
        for sib in sibships:
            kids = {c.sample: c for c in sib.children}
            if a in kids and b in kids:
                if kids[a].affected != kids[b].affected:
                    raise click.ClickException(
                        f"--mz-pair {a},{b}: the two must share a phenotype "
                        f"({a} is {kids[a].affected or 'unknown'}, {b} is "
                        f"{kids[b].affected or 'unknown'})"
                    )
                break
        else:
            raise click.ClickException(
                f"--mz-pair {a},{b}: not two children of one sibship in this run"
            )


def check_unaffected_carriers(sib: Sibship, n: int) -> None:
    """--max-unaffected-carriers N must leave the unaffected filter in force:
    0 <= N <= (unaffected children) - 1."""
    if n < 0:
        raise click.UsageError("--max-unaffected-carriers must be >= 0")
    if n > 0 and n > len(sib.unaffected) - 1:
        raise click.ClickException(
            f"--max-unaffected-carriers {n}: family {sib.family_id!r} has "
            f"{len(sib.unaffected)} unaffected children, so the filter would "
            f"exclude nothing; use at most {max(len(sib.unaffected) - 1, 0)}."
        )


# ─────────────────────────── what the data can say ───────────────────────────


def probe_caps() -> Caps:
    """Which accounting tables and columns this database has (the active one).
    A database ingested before `missing_genotypes` or `variants.n_called`
    existed gets the reduced builder."""
    from storage import table_columns, table_exists

    cols = set(table_columns("variants"))
    return Caps(
        have_miss=table_exists("missing_genotypes"),
        have_ncalled="n_called" in cols,
        have_filter="filter" in cols,
    )


def ingest_n_samples(sess, ingest_id: str) -> int:
    from storage import sql_quote_str

    q = sql_quote_str(ingest_id)
    rows = query_rows(
        sess, f"SELECT n_samples FROM ingestions WHERE ingest_id = {q} LIMIT 1"
    )
    if rows and rows[0][0] is not None:
        return int(rows[0][0])
    rows = query_rows(
        sess,
        f"SELECT count(DISTINCT sample_id) FROM samples WHERE ingest_id = {q}",
    )
    return int(rows[0][0])


# ─────────────────────────── sites and candidates ───────────────────────────


@dataclass(frozen=True)
class Site:
    chrom: str
    pos: int
    ref: str
    alt: str
    af: float | None
    pc: bool
    states: tuple[int, ...]  # father, mother, children in sibship order
    gts: tuple[int | None, ...]  # the stored gt of each member, None: no row
    parent_alt: tuple[int, int] | None = None

    @property
    def key(self) -> tuple[str, int, str, str]:
        return (self.chrom, self.pos, self.ref, self.alt)


def decode_rows(
    rows: Sequence[Sequence], n_members: int, *, denovo: bool
) -> list[Site]:
    out = []
    for r in rows:
        states = tuple(int(x) for x in r[6 : 6 + n_members])
        gps = [int(x) for x in r[6 + n_members : 6 + 2 * n_members]]
        tail = r[6 + 2 * n_members :]
        out.append(
            Site(
                r[0],
                int(r[1]),
                r[2],
                r[3],
                None if r[4] is None else float(r[4]),
                bool(int(r[5])),
                states,
                tuple(None if g == 0 else g - 1 for g in gps),
                (int(tail[0]), int(tail[1])) if denovo and len(tail) >= 2 else None,
            )
        )
    return out


def af_passes(af: float | None, gates: Gates) -> bool:
    """The rarity gate on INFO/AF: a missing AF passes; skipped entirely when
    the cohort-AF probe says INFO/AF is the cohort's own frequency."""
    if gates.af_gate == AF_GATE_SKIPPED or af is None:
        return True
    return af <= gates.max_af + AF_EPSILON


def children_of(sib: Sibship, site: Site) -> list[sv.Child]:
    return [
        sv.Child(c.sample, c.affected, site.states[2 + i])
        for i, c in enumerate(sib.children)
    ]


def state_name(state: int) -> str:
    return sv.STATE_NAMES[state]


_STORED = frozenset(
    {sv.HOM_REF, sv.HET, sv.HOM_ALT, sv.LOW_QUALITY, sv.REF_PARTIAL, sv.REF_ALT_READS}
)


def call_entry(sample: str, state: int, gt: int | None, inferred: bool) -> dict:
    """One member's `calls[]` object. gt is the stored value for a stored row
    and null otherwise; evidence says whether the state rests on a stored row
    or on the missing-call accounting."""
    return {
        "sample": sample,
        "gt": gt if state in _STORED else None,
        "state": state_name(state),
        "evidence": "stored"
        if state in _STORED
        else ("inferred" if state == sv.ABSENT_CALLED and inferred else None),
    }


def calls_for(sib: Sibship, site: Site, inferred: Sequence[str]) -> list[dict]:
    return [
        call_entry(c.sample, site.states[2 + i], site.gts[2 + i], c.sample in inferred)
        for i, c in enumerate(sib.children)
    ]


def by_status_zero() -> dict[str, int]:
    return dict.fromkeys(sv.STATUSES, 0)


@dataclass
class RunContext:
    """What every model of one sibship shares."""

    sess: Any
    sib: Sibship
    opts: SegOptions
    has_reference: bool
    caps: Caps
    n_samples: int
    scope: Any = None
    stats: dict = field(default_factory=dict)  # facts the effective report reads
    scope_stats: dict = field(default_factory=dict)

    @property
    def sql_options(self) -> SqlOptions:
        return SqlOptions(
            n_samples=self.n_samples,
            absent_called=self.opts.absent_called,
            max_unaffected_carriers=self.opts.max_unaffected_carriers,
            partial_strict=self.opts.partial_calls == "strict",
        )


def blocked_reason(model: str, ctx: RunContext, genes_loaded: bool) -> str | None:
    """Why a model cannot run on this database, or None. De novo needs stored
    parental 0/0 under every policy (a de novo claim needs positive reference
    evidence); dominant and comphet need it unless --absent called."""
    needs_ref = model == "denovo" or (
        model in ("dominant", "comphet") and not ctx.opts.absent_called
    )
    if needs_ref and not ctx.has_reference:
        return "needs_keep_reference"
    if model == "comphet" and not genes_loaded:
        return "needs_genes"
    return None


def _scoped(rows: list[Site], ctx: RunContext) -> list[Site]:
    if ctx.scope is None:
        return rows
    return [r for r in rows if ctx.scope.contains(r.chrom, r.pos)]


def _query_sites(model: str, ctx: RunContext) -> list[Site]:
    sql = build_sql(
        model, ctx.sib, ctx.opts.gates, ctx.sql_options, ctx.caps, scope=ctx.scope
    )
    rows = query_rows(ctx.sess, sql)
    return decode_rows(rows, len(ctx.sib.members), denovo=model == "denovo")


def _gnomad(items: list, ctx: RunContext, key) -> list:
    """Keep the items whose variant passes the gnomAD popmax filter (absent
    from the loaded slice: kept as rare)."""
    if ctx.opts.gnomad_max_af is None or not items:
        return items
    keys = [(*key(it)[:4], i) for i, it in enumerate(items)]
    kept = _gnomad_filter(keys, ctx.opts.gnomad_max_af, ctx.stats)
    keep = {k[4] for k in kept}
    return [it for i, it in enumerate(items) if i in keep]


def row_candidate(model: str, ctx: RunContext, site: Site, j: sv.Judgement) -> dict:
    """One candidate of a row model, before the annotation columns."""
    sib = ctx.sib
    out: dict = {
        "chrom": site.chrom,
        "pos": site.pos,
        "ref": site.ref,
        "alt": site.alt,
        "af": None if site.af is None else _num(site.af),
    }
    carrier = None
    if j.carrier_parent is not None:
        carrier = {
            "sample": j.carrier_parent,
            "affected": sib.father_affected
            if j.carrier_parent == sib.father
            else sib.mother_affected,
        }
    rest = {
        "status": j.status,
        "father_gt": site.gts[0] if site.states[0] in _STORED else None,
        "mother_gt": site.gts[1] if site.states[1] in _STORED else None,
        "father_state": state_name(site.states[0]),
        "mother_state": state_name(site.states[1]),
        "calls": calls_for(sib, site, j.inferred),
        "n_affected_carrying": j.n_affected_carrying,
        "n_unaffected_carrying": j.n_unaffected_carrying,
        "n_unaffected_ref": j.n_unaffected_ref,
        "unresolved": list(j.unresolved),
        "contradicting": list(j.contradicting),
        "inferred_reference": j.inferred_reference,
        "reduced_penetrance": j.status == sv.REDUCED_PENETRANCE,
        "partial_call_site": site.pc,
        "parental_mosaic_suspect": j.parental_mosaic_suspect,
        "parent_alt_reads": None
        if site.parent_alt is None
        else {"father": site.parent_alt[0], "mother": site.parent_alt[1]},
        "carrier_parent": carrier,
        "n_informative_meioses": j.n_informative_meioses,
        "p_random": j.p_random,
    }
    return {"base": out, "rest": rest}


def _num(v: float) -> float | int:
    f = float(v)
    return int(f) if f.is_integer() else f


def listed(status: str, opts: SegOptions) -> bool:
    return status != sv.COMPATIBLE or opts.allow_unresolved


def run_row_model(model: str, ctx: RunContext, render) -> dict:
    """The document object of denovo, recessive or dominant."""
    from cli.trio_output import decorate_rows

    opts, sib = ctx.opts, ctx.sib
    sites = _scoped(_query_sites(model, ctx), ctx)
    n_af = sum(1 for s in sites if not af_passes(s.af, opts.gates))
    sites = [s for s in sites if af_passes(s.af, opts.gates)]
    policy = opts.policy(sib)
    parent_pheno = {sib.father: sib.father_affected, sib.mother: sib.mother_affected}
    judged = []
    for s in sites:
        j = sv.judge(
            model,
            father_sample=sib.father,
            mother_sample=sib.mother,
            father=s.states[0],
            mother=s.states[1],
            children=children_of(sib, s),
            partial_site=s.pc,
            policy=policy,
            parent_phenotype=parent_pheno,
            parent_alt_reads=s.parent_alt,
        )
        if j is not None:
            judged.append((s, j))
    judged = _gnomad(judged, ctx, key=lambda it: it[0].key)
    by_status = by_status_zero()
    for _, j in judged:
        by_status[j.status] += 1
    shown = [(s, j) for s, j in judged if listed(j.status, opts)]
    page = shown[: opts.limit]
    built = [row_candidate(model, ctx, s, j) for s, j in page]
    outs = decorate_rows([b["base"] for b in built], [s.key for s, _ in page], render)
    candidates = [{**o, **b["rest"]} for o, b in zip(outs, built, strict=True)]
    return {
        "count": len(shown),
        "n_compatible": by_status[sv.COMPATIBLE],
        "by_status": by_status,
        "n_af_filtered": n_af,
        "blocked": False,
        "blocked_reason": None,
        "truncated": len(shown) > opts.limit,
        "candidates": candidates,
    }


def blocked_result(model: str, reason: str) -> dict:
    doc = {
        "count": 0,
        "n_compatible": 0,
        "by_status": by_status_zero(),
        "n_af_filtered": 0,
        "blocked": True,
        "blocked_reason": reason,
        "truncated": False,
    }
    if model == "comphet":
        doc["needs_genes"] = reason == "needs_genes"
        doc["genes"] = []
    else:
        doc["candidates"] = []
    return doc


# ───────────────────────────── compound het ─────────────────────────────


def _pair_id(s: Site) -> str:
    return f"{s.chrom}:{s.pos}:{s.ref}:{s.alt}"


def run_comphet(ctx: RunContext, render) -> dict:
    """Per gene, pairs (paternal hit, maternal hit): rare hets with one parent
    carrying and the other proven reference, both het in every affected child,
    no unaffected child carrying both. Python groups by gene (genes live in
    the annotation store), as `db trio` does."""
    from annotations.db import GeneIndex
    from cli.scope import genes_in_scope

    opts, sib = ctx.opts, ctx.sib
    policy = opts.policy(sib)
    sites = _scoped_exclusions(_query_sites("comphet", ctx), ctx)
    n_af = sum(1 for s in sites if not af_passes(s.af, opts.gates))
    sites = [s for s in sites if af_passes(s.af, opts.gates)]
    sites = [
        s
        for s in sites
        if sv.parents_fit("comphet", s.states[0], s.states[1], opts.absent_called)
    ]
    sites = _gnomad(sites, ctx, key=lambda s: s.key)

    index = GeneIndex(render.store)
    genes: dict[str, dict[str, list[Site]]] = {}
    for s in sites:
        origin = "paternal" if s.states[0] in sv.CARRIES else "maternal"
        for g in index.at(s.chrom, s.pos):
            entry = genes.setdefault(g.gene_symbol, {"paternal": [], "maternal": []})
            entry[origin].append(s)
    genes = {g: e for g, e in genes.items() if e["paternal"] and e["maternal"]}
    if ctx.scope is not None:
        wrapped = {
            g: {o: [(s.chrom, s.pos, s) for s in e[o]] for o in e}
            for g, e in genes.items()
        }
        kept, dropped = genes_in_scope(wrapped, ctx.scope)
        ctx.scope_stats["scope_genes_dropped"] = {"comphet": dropped}
        genes = {g: {o: [t[2] for t in e[o]] for o in e} for g, e in kept.items()}
        if ctx.scope.flags_variants:
            ctx.scope_stats["scope_outside"] = {
                "comphet": sum(
                    1
                    for e in genes.values()
                    for o in e
                    for s in e[o]
                    if not ctx.scope.contains(s.chrom, s.pos)
                )
            }

    by_status = by_status_zero()
    entries = []
    truncated = False
    n_genes = 0
    for sym in sorted(genes):
        entry = _gene_entry(sym, genes[sym], ctx, policy, by_status)
        if entry is None:
            continue
        n_genes += 1
        truncated = truncated or entry.pop("_truncated")
        entries.append(entry)
    shown = entries[: opts.limit]
    from cli.trio_output import decorate_rows

    # annotate every reported variant once (genes, gnomAD, ClinVar)
    variants = [v for e in shown for o in ("paternal", "maternal") for v in e[o]]
    decorate_rows(
        variants, [(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in variants], render
    )
    return {
        "count": n_genes,
        "n_compatible": by_status[sv.COMPATIBLE],
        "by_status": by_status,
        "n_af_filtered": n_af,
        "blocked": False,
        "blocked_reason": None,
        "needs_genes": False,
        "truncated": truncated or len(entries) > opts.limit,
        "genes": shown,
    }


def _scoped_exclusions(rows: list[Site], ctx: RunContext) -> list[Site]:
    """Comphet halves: only the scope's exclusions act before the grouping."""
    scope = ctx.scope
    if scope is None or not scope.exclude:
        return rows
    return [r for r in rows if not scope.excluded(r.chrom, r.pos)]


def _half_dict(sib: Sibship, s: Site, ctx: RunContext) -> dict:
    d = {
        "chrom": s.chrom,
        "pos": s.pos,
        "ref": s.ref,
        "alt": s.alt,
        "af": None if s.af is None else _num(s.af),
        "calls": calls_for(sib, s, ()),
    }
    scope = ctx.scope
    if scope is not None and scope.flags_variants:
        d["in_scope"] = scope.contains(s.chrom, s.pos)
    return d


def _gene_entry(
    sym: str,
    sides: dict[str, list[Site]],
    ctx: RunContext,
    policy: sv.Policy,
    by_status: dict[str, int],
) -> dict | None:
    sib, opts = ctx.sib, ctx.opts
    pairs = []
    used: dict[str, dict[tuple, Site]] = {"paternal": {}, "maternal": {}}
    truncated = False
    n = 0
    for p in sides["paternal"]:
        for m in sides["maternal"]:
            if (p.chrom, p.pos, p.ref, p.alt) == (m.chrom, m.pos, m.ref, m.alt):
                continue
            n += 1
            if n > MAX_PAIRS_PER_GENE:
                truncated = True
                break
            kids_p, kids_m = children_of(sib, p), children_of(sib, m)
            j = sv.judge_pair(
                paternal=kids_p,
                maternal=kids_m,
                paternal_partial=p.pc,
                maternal_partial=m.pc,
                parent_states=[*p.states[:2], *m.states[:2]],
                policy=policy,
            )
            if j is None:
                continue
            by_status[j.status] += 1
            if not listed(j.status, opts):
                continue
            used["paternal"][p.key] = p
            used["maternal"][m.key] = m
            pairs.append(
                {
                    "paternal": _pair_id(p),
                    "maternal": _pair_id(m),
                    "status": j.status,
                    "contradicting": list(j.contradicting),
                    "unresolved": list(j.unresolved),
                    "n_informative_meioses": j.n_informative_meioses,
                    "p_random": j.p_random,
                }
            )
        if truncated:
            break
    if not pairs:
        return None
    return {
        "gene": sym,
        "paternal": [_half_dict(sib, s, ctx) for s in used["paternal"].values()],
        "maternal": [_half_dict(sib, s, ctx) for s in used["maternal"].values()],
        "pairs": pairs,
        "_truncated": truncated,
    }


# ───────────────────────────── the document ─────────────────────────────


def gates_json(opts: SegOptions) -> dict:
    g = opts.gates
    out = {
        "min_gq": g.min_gq,
        "min_dp": g.min_dp,
        "max_af": g.max_af,
        "af_gate": g.af_gate,
        "min_ab": g.min_ab,
        "max_ab": g.max_ab,
        "min_homalt_ab": g.min_homalt_ab,
        "max_parent_alt_frac": g.max_parent_alt_frac,
        "gnomad_max_af": opts.gnomad_max_af,
    }
    if g.strict_missing:
        out["strict_missing"] = True
        if g.strict_exempt:
            out["strict_exempt"] = list(g.strict_exempt)
    if g.pass_only:
        out["pass_only"] = True
    return out


def family_json(sib: Sibship, opts: SegOptions) -> dict:
    mz = [list(p) for p in opts.mz_pairs if p[0] in {c.sample for c in sib.children}]
    return {
        "ingest_id": sib.ingest_id,
        "family_id": sib.family_id,
        "father": sib.father,
        "mother": sib.mother,
        "children": [
            {"sample": c.sample, "affected": c.affected, "sex": c.sex}
            for c in sib.children
        ],
        "mz_pairs": mz,
    }


def segregation_document(ctx: RunContext, render, notes: list[str]) -> dict:
    """The JSON document of one sibship: every selected model, run."""
    sib, opts = ctx.sib, ctx.opts
    genes_loaded = bool(render.have["genes"])
    counts = {
        "children": len(sib.children),
        "affected": len(sib.affected),
        "unaffected": len(sib.unaffected),
        "unknown": len(sib.unknown),
        "skipped_sites": int(query_rows(ctx.sess, skipped_sites_sql(sib))[0][0]),
        "partial_call_sites": 0,
    }
    if ctx.caps.have_miss and ctx.caps.have_ncalled and opts.partial_calls == "strict":
        counts["partial_call_sites"] = int(
            query_rows(ctx.sess, partial_sites_sql(sib, ctx.n_samples, ctx.scope))[0][0]
        )
    notes = list(notes)
    if counts["skipped_sites"]:
        notes.append(
            f"{counts['skipped_sites']} candidate sites on sex or mitochondrial "
            "contigs were skipped (haploid calls are stored as hom-alt, so no "
            "model here can read them)."
        )
    if not sib.unaffected:
        notes.append(
            "no unaffected child: candidates are filtered by being shared by "
            "every affected child only."
        )
    if sib.unknown:
        notes.append(
            "children with unknown phenotype are shown and never constrain a model: "
            + ", ".join(c.sample for c in sib.unknown)
            + "."
        )
    if not (ctx.caps.have_miss and ctx.caps.have_ncalled):
        notes.append(NO_ACCOUNTING_NOTE + ".")
    if opts.max_unaffected_carriers and "comphet" in opts.models:
        notes.append(
            "--max-unaffected-carriers is ignored for comphet: a pair carried "
            "by an unaffected child is excluded."
        )
    if opts.absent_called:
        notes.append(
            "--absent called: a member with no row at a fully accounted site "
            "counts as reference (status consistent_inferred); de novo still "
            "needs stored parental 0/0 rows."
        )
    models: dict[str, dict] = {}
    for model in opts.models:
        reason = blocked_reason(model, ctx, genes_loaded)
        if reason is not None:
            models[model] = blocked_result(model, reason)
        elif model == "comphet":
            models[model] = run_comphet(ctx, render)
        else:
            models[model] = run_row_model(model, ctx, render)
    return {
        "family": family_json(sib, opts),
        "counts": counts,
        "keep_reference": ctx.has_reference,
        "absent": opts.absent,
        "partial_calls": opts.partial_calls,
        "parent_phenotype": opts.parent_phenotype,
        "gates": gates_json(opts),
        "options": {
            "allow_unresolved": opts.allow_unresolved,
            "max_unaffected_carriers": opts.max_unaffected_carriers,
            "limit": opts.limit,
        },
        "annotations": render.have,
        "notes": notes,
        "models": models,
    }


def with_ingest(sib: Sibship, ingest_id: str) -> Sibship:
    return replace(sib, ingest_id=ingest_id)


def scope_builder(trio: Any, gates: Gates, scope: Any) -> list[tuple[str, str, str]]:
    """`cli.scope.SCOPE_AWARE_BUILDERS` entry: every statement a model runs
    under a scope, for a stand-in sibship of the given trio's parents."""
    sib = Sibship(
        ingest_id=trio.ingest_id,
        family_id="fam",
        father=trio.father or "F",
        mother=trio.mother or "M",
        children=(
            ChildInfo(trio.proband, "affected"),
            ChildInfo(f"{trio.proband}_sib", "unaffected"),
        ),
    )
    opts = SqlOptions(n_samples=3)
    out = [
        (
            build_sql(model, sib, gates, opts, Caps(), scope=scope),
            "exclude" if model == "comphet" else "all",
            "a",
        )
        for model in MODELS
    ]
    out.append((partial_sites_sql(sib, 3, scope), "all", "a"))
    return out


def _register_scope_builder() -> None:
    from cli.scope import register_scope_builder

    register_scope_builder("segregate", scope_builder)


_register_scope_builder()

__all__ = [
    "SegOptions",
    "find_sibships",
    "segregation_document",
]
