"""`vcfclick db segregate` — sibling segregation.

For a sibship (children with the same declared father and mother) it reports,
per model, the variants carried by the affected children and not contradicted
by the unaffected ones, with every member's state at the site and a status
ladder (consistent > consistent_inferred > compatible > reduced_penetrance).

    vcfclick db ped fam1 fam1.ped
    vcfclick db segregate fam1 --category dominant

Models (all autosomal, slivar-style gates as `db trio`):

  * denovo     an affected child carries, both parents are stored 0/0, no
               unaffected child carries.
  * dominant   every affected child het, exactly one parent carries, the other
               and every unaffected child proven reference.
  * recessive  every affected child hom-alt, both parents het, no unaffected
               child hom-alt.
  * comphet    per gene, a paternal and a maternal het in every affected child,
               no unaffected child carrying both.

A model a child's data cannot decide is `compatible`, listed with
--allow-unresolved. See docs/TRIO.md#sibling-segregation.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from types import SimpleNamespace

import click

from cli.family_checks import trio_af_gate
from cli.main import _set_db, db
from cli.options import command_options
from cli.preset_options import (
    PresetCommand,
    gate_options,
    is_explicit,
    preset_option,
    preset_state,
    resolve_filters,
)
from cli.trio_analysis import (
    Gates,
    Trio,
    TrioAnalysis,
    _has_reference_rows,
    _sole_ingest_id,
)
from cli.xchrom import BUILD_CHOICES


@dataclass(frozen=True)
class SegregateOptions:
    name: str
    family: str | None
    father: str | None
    mother: str | None
    category: tuple[str, ...]
    min_gq: int
    min_dp: int
    max_af: float
    min_ab: float
    max_ab: float
    gnomad_max_af: float | None
    strict_missing: bool
    pass_only: bool
    max_parent_alt_frac: float
    require_bound: bool
    min_homalt_ab: float
    partial_calls: str
    parent_phenotype: str
    allow_unresolved: bool
    genes: tuple[str, ...]
    gene_file: str | None
    exclude_genes: tuple[str, ...]
    exclude_gene_file: str | None
    gene_flank: int
    region: tuple[str, ...]
    region_file: str | None
    exclude_region: tuple[str, ...]
    allow_unresolved_genes: bool
    build: str
    absent: str
    max_unaffected_carriers: int
    mz_pair: tuple[str, ...]
    out_format: str
    limit: int
    preset: str | None
    preset_fingerprint: str | None
    preset_ignore: tuple[str, ...]
    preset_skip_unavailable: bool

    @property
    def gates(self) -> Gates:
        return Gates(
            self.min_gq,
            self.min_dp,
            self.max_af,
            self.min_ab,
            self.max_ab,
            max_parent_alt_frac=self.max_parent_alt_frac,
            strict_missing=self.strict_missing,
            pass_only=self.pass_only,
            min_homalt_ab=self.min_homalt_ab,
        )


def parse_mz_pairs(values) -> tuple[tuple[str, str], ...]:
    """`--mz-pair A,B` (repeatable) as pairs of two different sample ids."""
    out = []
    for v in values:
        parts = [p.strip() for p in str(v).split(",")]
        if len(parts) != 2 or not all(parts) or parts[0] == parts[1]:
            raise click.UsageError(
                f"--mz-pair {v!r}: give two different sample ids, A,B"
            )
        out.append((parts[0], parts[1]))
    return tuple(out)


@db.command(name="segregate", cls=PresetCommand)
@click.argument("name")
@click.option(
    "--family",
    default=None,
    help="Only the sibship(s) of this family_id. Without --family, --father or "
    "--mother, every eligible sibship is analysed (json: a `families` list).",
)
@click.option(
    "--father",
    default=None,
    help="The sibship whose declared father is this sample (with --family or "
    "--mother to pick one half-sibship).",
)
@click.option(
    "--mother",
    default=None,
    help="The sibship whose declared mother is this sample.",
)
@gate_options("segregate")
@click.option(
    "--allow-unresolved-genes",
    is_flag=True,
    help="Run even when tokens of --genes / --exclude-genes resolve to no gene "
    "(by default a list of 20 or fewer may miss none, a longer one at most 5%); "
    "they are left out of the scope and reported.",
)
@click.option(
    "--build",
    type=click.Choice(BUILD_CHOICES),
    default="auto",
    show_default=True,
    help="Assembly of the database, for --genes / --region (auto: from chrX "
    "positions when the data prove it).",
)
@click.option(
    "--absent",
    type=click.Choice(["strict", "called"]),
    default="strict",
    show_default=True,
    help="strict: only a stored 0/0 proves a parent or an unaffected child "
    "does not carry the allele. called: a member with no row at a site whose "
    "calls are all accounted for (n_called + missing = n_samples) counts as "
    "reference; the candidate is then consistent_inferred, never consistent. "
    "De novo still needs stored parental 0/0 rows.",
)
@click.option(
    "--max-unaffected-carriers",
    type=int,
    default=0,
    show_default=True,
    help="Reduced penetrance (non-standard): tolerate up to N unaffected "
    "children that contradict the model (carry the variant; hom-alt for "
    "recessive). Such candidates are status reduced_penetrance, never "
    "consistent. 0 <= N <= unaffected children - 1. Ignored for comphet.",
)
@click.option(
    "--mz-pair",
    "mz_pair",
    multiple=True,
    metavar="A,B",
    help="Two children who are monozygotic twins or duplicate samples "
    "(repeatable): they count as one informative meiosis and one affected "
    "group. Both must share a phenotype. Not detected automatically.",
)
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
    help="json: one document per sibship (`families` list unless a selector "
    "names one sibship).",
)
@click.option(
    "--limit",
    type=int,
    default=1000,
    show_default=True,
    help="json: at most this many candidates per model.",
)
@preset_option()
@command_options(SegregateOptions)
def db_segregate(options: SegregateOptions) -> None:
    """Report, per inheritance model, the variants the affected children of a
    sibship carry and the unaffected children do not contradict, with every
    member's state and a status (consistent, consistent_inferred, compatible,
    reduced_penetrance)."""
    from cli import segregation as seg
    from cli.pedigree_load import load_pedigree_graph
    from cli.trio_gate_report import attach_gate_reports, needs_report
    from cli.trio_output import open_render
    from storage import db_path, get_session

    if not db_path(options.name).exists():
        raise click.ClickException(f"db {options.name!r} does not exist.")
    if options.limit < 1:
        raise click.UsageError("--limit must be at least 1")
    cat_ctx = click.get_current_context()
    cat_typed = is_explicit(cat_ctx, "category")
    models = seg.select_models(options.category, typed=cat_typed)
    if not cat_typed and "xlinked" in options.category:
        click.echo(
            "note: the preset's xlinked category is not a segregation model; skipped.",
            err=True,
        )
    mz_pairs = parse_mz_pairs(options.mz_pair)
    _set_db(options.name)
    sess = get_session(options.name)
    ingest_id = _sole_ingest_id(options.name)
    if ingest_id is None:
        raise click.ClickException(
            "database has multiple ingestions; segregation analysis assumes one "
            "joint-called cohort. Re-ingest the family as a single VCF "
            "(see `vcfclick merge`)."
        )
    graph = load_pedigree_graph(sess, ingest_id)
    sibships, notes = seg.find_sibships(
        graph, family=options.family, father=options.father, mother=options.mother
    )
    sibships = [seg.with_ingest(s, ingest_id) for s in sibships]
    seg.check_mz_pairs(mz_pairs, sibships)
    for sib in sibships:
        seg.check_unaffected_carriers(sib, options.max_unaffected_carriers)
    for note in notes:
        click.echo(f"note: {note}", err=True)

    ctx = click.get_current_context()
    af_gate, af_message = trio_af_gate(
        sess, ingest_id, options.max_af, explicit=is_explicit(ctx, "max_af")
    )
    run_notes: list[str] = []
    if af_message:
        click.echo(f"note: {af_message}", err=True)
        run_notes.append(af_message)
    gates = options.gates._replace(af_gate=af_gate)
    has_ref = _has_reference_rows(sess)
    scope = _build_scope(ctx, options, sess, ingest_id)
    caps = seg.probe_caps()
    n_samples = seg.ingest_n_samples(sess, ingest_id)
    opts = seg.SegOptions(
        models=models,
        gates=gates,
        gnomad_max_af=options.gnomad_max_af,
        absent=options.absent,
        partial_calls=options.partial_calls,
        allow_unresolved=options.allow_unresolved,
        max_unaffected_carriers=options.max_unaffected_carriers,
        parent_phenotype=options.parent_phenotype,
        mz_pairs=mz_pairs,
        limit=options.limit,
    )

    # The effective-gate report (does each gate bind in this data?) and the
    # preset echo: computed for a stand-in trio of each sibship (its first
    # affected child and the parents), because the fields a gate reads are a
    # property of the data, not of the sibling. Only for a run that asks.
    state = preset_state(ctx)
    proxies: list[TrioAnalysis] = []
    if needs_report(
        state, options.require_bound, options.strict_missing, False, options.pass_only
    ):
        resolved, overridden = (
            resolve_filters(ctx, "segregate", state) if state else (None, [])
        )
        report_options = SimpleNamespace(
            require_bound=options.require_bound,
            strict_missing=options.strict_missing,
            parent_alt_screen=False,
            pass_only=options.pass_only,
            category=models,
        )
        proxies = attach_gate_reports(
            [_proxy(sess, sib, gates, has_ref, opts, scope) for sib in sibships],
            ctx,
            options=report_options,
            state=state,
            resolved=resolved,
            overridden=overridden,
        )
        if scope is not None:
            from cli.db_trio import _attach_scope

            for proxy in proxies:
                _attach_scope(proxy, scope)

    docs: list[dict] = []
    headers: list[list[str]] = []
    with open_render(None, options.limit) as render:
        for i, sib in enumerate(sibships):
            run_opts = replace(opts, gates=proxies[i].gates) if proxies else opts
            rctx = seg.RunContext(
                sess=sess,
                sib=sib,
                opts=run_opts,
                has_reference=has_ref,
                caps=caps,
                n_samples=n_samples,
                scope=scope,
                stats=proxies[i].stats if proxies else {},
            )
            doc = seg.segregation_document(rctx, render, run_notes)
            if proxies:
                stats = proxies[i].stats
                if "preset_block" in stats:
                    doc["preset"] = stats["preset_block"]
                if "effective_block" in stats:
                    doc["effective"] = stats["effective_block"]
                headers.append(list(stats.get("text_header", ())))
            elif scope is not None:
                from cli.scope import text_lines

                headers.append(text_lines(scope))
            else:
                headers.append([])
            if scope is not None:
                from cli.scope import scope_document

                doc["scope"] = scope_document(scope, rctx.scope_stats)
            docs.append(doc)
    for doc in docs:
        for note in doc["notes"]:
            if note not in run_notes:
                click.echo(f"note: {note}", err=True)
    selecting = any((options.family, options.father, options.mother))
    if options.out_format == "json":
        result = docs[0] if selecting else {"families": docs}
        click.echo(json.dumps(result, indent=2))
        return
    from cli.segregation_output import text_report

    for i, (doc, head) in enumerate(zip(docs, headers, strict=True)):
        if i:
            click.echo()
        text_report(doc, options.category, head)


def _proxy(sess, sib, gates, has_ref, opts, scope) -> TrioAnalysis:
    first = sib.affected[0].sample
    return TrioAnalysis(
        session=sess,
        trio=Trio(sib.ingest_id, first, sib.father, sib.mother, sib.family_id),
        gates=gates,
        has_reference=has_ref,
        gnomad_max_af=opts.gnomad_max_af,
        scope=scope,
    )


def _build_scope(ctx: click.Context, options: SegregateOptions, sess, ingest_id: str):
    """The run's cli.scope.Scope (genes, regions) as `db trio` builds it."""
    from cli.db_trio import _build_scope as build

    return build(ctx, options, sess, ingest_id, "segregate")
