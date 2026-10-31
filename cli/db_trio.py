"""`vcfclick db trio` — family-based analysis.

A PED file (or `db ped` flags, see cli/db_ped.py) declares parents among
already-ingested samples. Load it, then `db trio` resolves a proband's
parents and reports candidate variants under each Mendelian inheritance
model.

    vcfclick db ped  fam1 fam1.ped
    vcfclick db trio fam1 --proband CHILD

Inheritance models, with slivar-style genotype quality (GQ, depth,
allele balance) and population-AF rarity:

  * de novo     proband carries; BOTH parents provably hom-ref (gt=0).
                Requires a `--keep-reference` ingest — without stored
                parent 0/0 rows, "parent absent" is a no-call, not
                reference, so de novo is undecidable.
  * recessive   proband hom-alt; both parents heterozygous carriers.
                Works on a normal sparse ingest.
  * dominant    proband het; exactly one parent carries, the other
                provably hom-ref. Also needs --keep-reference (to prove
                the non-carrier parent is reference, not no-call).
  * comphet     two rare proband hets in the SAME gene, one inherited
                from each parent (trans → both gene copies hit). Needs
                gene annotations loaded (`vcfclick annotations load`) and
                --keep-reference (to prove each non-carrier parent is
                hom-ref). Reported per gene, not per variant.
  * xlinked     sex- and build-aware patterns on non-PAR chrX (a son read
                with his mother only, a daughter with a hemizygous father);
                cli/xlinked.py. The X guard drops a het father and a het
                male proband there from the four models above
                (--no-x-guard turns it off).
"""

from __future__ import annotations

import json
from dataclasses import dataclass

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
from cli.trio_analysis import DEFAULT_MAX_AF as DEFAULT_MAX_AF
from cli.trio_analysis import Gates as Gates
from cli.trio_analysis import Trio as Trio
from cli.trio_analysis import (
    TrioAnalysis,
    _has_reference_rows,
    _resolve_parents,
    _sole_ingest_id,
)
from cli.trio_analysis import _comphet_genes as _comphet_genes
from cli.trio_analysis import trio_sql as trio_sql
from cli.trio_gate_output import show_parent_reads
from cli.trio_gate_report import attach_gate_reports, needs_report
from cli.trio_output import _trio_json, omitted_models, trio_text
from cli.trio_select import make_pick, select_trios
from cli.xchrom import BUILD_CHOICES


@dataclass(frozen=True)
class TrioOptions:
    name: str
    proband: str | None
    family: str | None
    all_trios: bool
    category: tuple[str, ...]
    min_gq: int
    min_dp: int
    max_af: float
    min_ab: float
    max_ab: float
    gnomad_max_af: float | None
    strict_missing: bool
    pass_only: bool
    parent_alt_screen: bool
    max_parent_alt_frac: float
    require_bound: bool
    min_male_x_dp: int | None
    max_af_x: float
    min_hemi_ab: float
    max_parent_alt_reads: int
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
    x_guard: bool
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
            parent_alt_screen=self.parent_alt_screen,
            min_male_x_dp=self.min_male_x_dp,
            max_af_x=self.max_af_x,
            min_hemi_ab=self.min_hemi_ab,
            max_parent_alt_reads=self.max_parent_alt_reads,
        )


@db.command(name="trio", cls=PresetCommand)
@click.argument("name")
@click.option(
    "--proband",
    default=None,
    help="Sample id of the affected child. Omitted: every proband in the "
    "pedigree, one report each. A proband is a child with both parents "
    "declared who is not marked unaffected and who either has no children "
    "of their own or is affected (`db families` lists them).",
)
@click.option(
    "--family",
    default=None,
    help="Analyse every proband of this family (see `db families`). Not "
    "with --proband.",
)
@click.option(
    "--all-trios",
    is_flag=True,
    help="Without --proband: analyse every child with both parents "
    "declared, including children marked unaffected and parents who have "
    "children of their own (the rule before probands were defined).",
)
@gate_options("trio")
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
    help="Assembly of chrX, for its non-PAR region (xlinked and the X guard). "
    "auto: GRCh38 when a chrX position exceeds the GRCh37 chrX length "
    "(155,270,560), else unknown: the two bands where GRCh37 and GRCh38 "
    "disagree are left out and counted. CHM13/T2T is not supported.",
)
@click.option(
    "--x-guard/--no-x-guard",
    "x_guard",
    default=True,
    show_default=True,
    help="On sure non-PAR chrX, drop a het father, and a het male proband, "
    "from denovo, recessive, dominant and comphet: a male has one X, so "
    "those calls are artefacts. --no-x-guard restores the unguarded models.",
)
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
    help="json: every model's count and candidates (with gene, gnomAD and "
    "ClinVar when loaded) in one document; without --proband, a `trios` "
    "list with one document per proband, and the `families` they are in.",
)
@click.option(
    "--limit",
    type=int,
    default=1000,
    show_default=True,
    help="json: at most this many candidates per model.",
)
@preset_option()
@command_options(TrioOptions)
def db_trio(options: TrioOptions) -> None:
    """Report candidate variants under Mendelian inheritance models for
    a trio, with genotype quality gates and an AF rarity filter. Without
    --proband, every proband in the pedigree is analysed in turn (--family
    limits that to one family, --all-trios widens it to every child)."""
    from cli.pedigree_load import load_pedigree_graph
    from storage import db_path, get_session

    if not db_path(options.name).exists():
        raise click.ClickException(f"db {options.name!r} does not exist.")
    if options.proband is not None and (options.family or options.all_trios):
        raise click.UsageError(
            "--family and --all-trios choose which children to analyse; "
            "they cannot be combined with --proband."
        )
    _set_db(options.name)
    sess = get_session(options.name)

    ingest_id = _sole_ingest_id(options.name)
    if ingest_id is None:
        raise click.ClickException(
            "database has multiple ingestions; trio analysis assumes one "
            "joint-called cohort. Re-ingest the trio as a single VCF "
            "(see `vcfclick merge`)."
        )
    if options.proband is not None:
        # The proband's own row decides the parents, as it always has; the
        # graph (read tolerantly) only adds the family and the siblings.
        father, mother = _resolve_parents(sess, ingest_id, options.proband)
        graph = load_pedigree_graph(sess, ingest_id)
        picks = [make_pick(graph, options.proband, father, mother)]
        notes: list[str] = []
    else:
        graph = load_pedigree_graph(sess, ingest_id)
        selection = select_trios(
            graph, family=options.family, all_trios=options.all_trios
        )
        picks, notes = selection.picks, selection.notes
    for note in notes:
        click.echo(f"note: {note}", err=True)
    has_ref = _has_reference_rows(sess)
    # Is INFO/AF this cohort's own frequency? Then --max-af would remove
    # every private candidate: say so, and skip the gate unless --max-af was
    # typed on the command line (an explicit threshold is always honoured).
    # A value a preset supplied is not typed.
    ctx = click.get_current_context()
    af_gate, af_message = trio_af_gate(
        sess, ingest_id, options.max_af, explicit=is_explicit(ctx, "max_af")
    )
    if af_message:
        click.echo(f"note: {af_message}", err=True)
    gates = options.gates._replace(af_gate=af_gate)
    scope = _build_scope(ctx, options, sess, ingest_id)
    analyses = [
        TrioAnalysis(
            session=sess,
            trio=Trio(ingest_id, p.child, p.father, p.mother, p.family_id),
            gates=gates,
            has_reference=has_ref,
            gnomad_max_af=options.gnomad_max_af,
            siblings=p.siblings,
            build=options.build,
            x_guard=options.x_guard,
            scope=scope,
        )
        for p in picks
    ]
    for analysis in analyses:
        for model in omitted_models(analysis, options.category):
            click.echo(
                f"note: the {model} model does not apply to proband "
                f"{analysis.trio.proband}; left out",
                err=True,
            )
    state = preset_state(ctx)
    if state or options.parent_alt_screen:
        for analysis in analyses:
            show_parent_reads(analysis)
    if needs_report(
        state,
        options.require_bound,
        options.strict_missing,
        options.parent_alt_screen,
        options.pass_only,
    ):
        resolved, overridden = (
            resolve_filters(ctx, "trio", state) if state else (None, [])
        )
        analyses = attach_gate_reports(
            analyses,
            ctx,
            options=options,
            state=state,
            resolved=resolved,
            overridden=overridden,
        )
    if scope is not None:
        for analysis in analyses:
            _attach_scope(analysis, scope)
    if options.out_format == "json":
        docs = [_trio_json(a, options.category, options.limit) for a in analyses]
        if options.proband is not None:
            result = docs[0]
        else:
            result = {
                "trios": docs,
                "families": list(dict.fromkeys(p.family_id for p in picks)),
            }
        click.echo(json.dumps(result, indent=2))
        return
    several = len({p.family_id for p in picks}) > 1
    for i, (pick, analysis) in enumerate(zip(picks, analyses, strict=True)):
        if i:
            click.echo()
        if several:
            click.echo(f"family: {pick.family_id}")
        trio_text(analysis, options.category)


def _build_scope(ctx: click.Context, options: TrioOptions, sess, ingest_id: str):
    """The run's cli.scope.Scope from the resolved scope values (flag > preset >
    default), or None. Gene lists are resolved here, once for every trio; a
    refusal (unresolved genes, no gene table, a GRCh37 database) is a
    ClickException."""
    from cli import filter_registry as fr
    from cli.scope import SCOPE_KEYS, ScopeError, ScopeRequest, build_scope

    state = preset_state(ctx)
    resolved, _ = resolve_filters(ctx, "trio", state)
    values = {k: resolved[k].value for k in SCOPE_KEYS if k in resolved}
    if state is not None:
        # the descriptive companions follow the list they describe: a typed
        # --genes replaces the preset's list, and its label with it
        filters = state.located.preset.filters
        for meta, of in (("genes_label", "genes"), ("regions_build", "regions")):
            if meta in filters and resolved[of].source == fr.FROM_PRESET:
                values[meta] = filters[meta]
    try:
        request = ScopeRequest.from_values(
            values, allow_unresolved=options.allow_unresolved_genes
        )
        if not request.any:
            return None
        needs_build = request.genes or request.exclude_genes or request.regions_build
        build = _database_build(options, sess, ingest_id) if needs_build else None
        scope = build_scope(request, build=build)
    except ScopeError as e:
        raise click.ClickException(str(e)) from None
    for note in scope.notes:
        click.echo(f"note: {note}", err=True)
    return scope


def _database_build(options: TrioOptions, sess, ingest_id: str) -> str:
    """GRCh37 / GRCh38 when --build says so or the chrX data proves it, else
    `unknown` (cli/xchrom.resolve_build)."""
    from cli.xchrom import BUILDS, resolve_build
    from cli.xlinked import _max_x_pos

    if options.build in BUILDS:
        return options.build
    build, _, _ = resolve_build(options.build, _max_x_pos(sess, ingest_id))
    return build


def _attach_scope(analysis: TrioAnalysis, scope) -> None:
    """The scope's text lines under the run's header, and its warnings in the
    `effective` object when the run has one."""
    from cli.scope import text_lines

    stats = analysis.stats
    stats["text_header"] = [*stats.get("text_header", ()), *text_lines(scope)]
    effective = stats.get("effective_block")
    if effective is not None:
        have = {w["message"] for w in effective["warnings"]}
        effective["warnings"].extend(
            w for w in scope.warnings if w["message"] not in have
        )
