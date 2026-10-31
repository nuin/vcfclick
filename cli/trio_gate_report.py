"""`db trio` with a preset or a strictness flag: compute, enforce, leave the echo.

For a run that uses a preset, `--strict-missing`, `--require-bound-gates` or
`--parent-alt-screen`, `attach_gate_reports` asks the data whether each gate binds
(`cli.effective_gates`), stops the run when it was told to insist and a gate
cannot bind, and leaves the answers in each analysis' `stats` for the output
hooks of `cli.trio_gate_output` (the `preset` and `effective` JSON objects, the
text header). A plain run never gets here, so its output and its queries are
unchanged.

The GQ / DP / AD gates are judged on the rows each selected model reads at its
candidate sites, counted here with one grouped query per model
(`collect_candidate_stats`, only for a model that reads a gate the report judges
and that can run), not on the share over all of a sample's rows.

A model the run selected that cannot run on this database (no stored hom-reference
calls, no gene table) is a `model_blocked` warning, and an error for a run that
insists: a clean "0 candidates" from a model that could not run is the quiet
result this report exists to prevent.

`--parent-alt-screen` typed alone is an opt-in gate that exists to catch false
de novos, and on a VCF without allele depths it screens nothing: every parent
0/0 is trusted. `--pass-only` typed alone keeps FILTER = PASS, and on a VCF whose
FILTER column carries no label it removes nothing. Those runs are reported too,
but only about the gates that were typed (the others were not asked for), so
they add their warning and the `effective` object.

A run told to go ahead anyway (`--no-require-bound-gates`) does not fail every
call on a field the data lacks: for each strict_missing gate that cannot bind,
the analysis is re-made with that gate exempt from strict_missing
(`Gates.strict_exempt`), and the report says so. Without that, "run anyway"
would mean "no candidates", whatever the genome holds.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace
from functools import partial

import click

from cli import filter_registry as fr
from cli.effective_gates import (
    build_report,
    candidate_models,
    collect_candidate_stats,
    collect_gnomad_probe,
    collect_sample_stats,
    collect_variant_stats,
    refusal_message,
    strict_check_models,
    unbound_strict_keys,
)
from cli.preset_echo import header_lines, preset_block, skipped_items
from cli.preset_options import PresetState, is_explicit
from cli.trio_analysis import TrioAnalysis
from cli.trio_gate_output import run_footer
from cli.trio_output import (
    MODELS,
    _annotation_status,
    _applicable,
    model_blockers,
    selected_categories,
)

# info-level codes the user does not need to see twice
_SILENT_NOTES = ("preset_max_af_not_applied",)


def needs_report(
    state: PresetState | None,
    require_bound: bool,
    strict: bool,
    parent_alt_screen: bool = False,
    pass_only: bool = False,
) -> bool:
    return (
        state is not None or require_bound or strict or parent_alt_screen or pass_only
    )


def report_scope(
    state: PresetState | None,
    require_bound: bool,
    strict: bool,
    parent_alt_screen: bool = False,
    pass_only: bool = False,
) -> frozenset[str] | None:
    """The gate keys a run's warnings are limited to: None (all of them) for a
    preset, `--strict-missing` or `--require-bound-gates`; only the opt-in gates
    that were typed (the parent-allele screen, the FILTER gate) for a run that
    asked for nothing else."""
    if state is not None or require_bound or strict:
        return None
    return frozenset(
        key
        for key, on in (
            ("parent_alt_screen", parent_alt_screen),
            ("pass_only", pass_only),
        )
        if on
    )


def _gnomad_chromosomes() -> frozenset[str] | None:
    """The chromosomes the loaded gnomAD data has records on; None when the
    store cannot be read (the coverage is then not judged)."""
    try:
        from annotations.db import get_connection

        conn = get_connection()
    except Exception:  # noqa: BLE001 - no annotation store at all
        return None
    try:
        rows = conn.execute("SELECT DISTINCT chrom FROM gnomad_af").fetchall()
    except Exception:  # noqa: BLE001 - table absent in an older store
        return None
    finally:
        conn.close()
    return frozenset(r[0] for r in rows)


def attach_gate_reports(
    analyses: Sequence[TrioAnalysis],
    ctx: click.Context,
    *,
    options,
    state: PresetState | None,
    resolved: dict[str, fr.Resolved] | None,
    overridden: list[dict],
) -> list[TrioAnalysis]:
    """Fill `analysis.stats` for every trio, or raise a ClickException naming
    the gates that cannot bind when the run insists on it. Returns the analyses
    to run: the ones given, except that a trio whose strict_missing gates could
    not bind under `--no-require-bound-gates` comes back with those gates
    exempt."""
    first = analyses[0]
    sess, ingest_id = first.session, first.trio.ingest_id
    from storage import table_columns

    filter_present = "filter" in set(table_columns("variants"))
    variants = collect_variant_stats(sess, ingest_id, filter_present=filter_present)
    status = _annotation_status()
    gnomad_loaded = bool(status["gnomad"])
    genes_loaded = bool(status["genes"])
    wants_gnomad = first.gnomad_max_af is not None
    # --no-require-bound-gates typed on the command line overrides both kinds
    # of refusal for this run
    overridden_off = is_explicit(ctx, "require_bound") and not options.require_bound
    enforce_bound = bool(options.require_bound)
    enforce_format = bool(options.strict_missing) and not overridden_off
    only = report_scope(
        state,
        enforce_bound,
        bool(options.strict_missing),
        bool(options.parent_alt_screen),
        bool(options.pass_only),
    )
    # gnomAD is looked at (its chromosomes, and a sample of the cohort's own
    # variants looked up in it) only when the report is about it
    looked_up = wants_gnomad and gnomad_loaded and only is None
    gnomad_chroms = _gnomad_chromosomes() if looked_up else None
    gnomad_probe = (
        collect_gnomad_probe(sess, ingest_id, variants) if looked_up else None
    )
    subject = f"preset {state.located.preset.name!r}" if state else "this run"
    skipped = []
    preset_max_af = None
    notes: list[str] = list(state.notes) if state else []
    if state:
        skipped = skipped_items(
            state.plan, state.command, skip_unavailable=state.skip_unavailable
        )
        r = resolved["max_af"] if resolved else None
        if r is not None and r.source == fr.FROM_PRESET:
            preset_max_af = r.value

    # one pass over the rows of every sample of every trio (siblings share
    # parents), then each trio takes its own three with their roles
    unique = dict.fromkeys(
        s for a in analyses for s in (a.trio.proband, a.trio.father, a.trio.mother) if s
    )
    counted = {
        c.sample: c
        for c in collect_sample_stats(
            sess, ingest_id, [(s, "") for s in unique], scope=first.scope
        )
    }

    def blocked_models(analysis: TrioAnalysis) -> list[tuple[str, str]]:
        """(model, reason) for each selected model that cannot run on this
        database: its count would be 0 for a reason that is not the genome."""
        return [
            (spec.name, reason)
            for spec in _applicable(analysis, options.category)
            for reason in model_blockers(analysis, spec, genes_loaded)
        ]

    def candidate_rows(analysis: TrioAnalysis) -> dict:
        """What the rows at each runnable selected model's candidate sites carry:
        the rows the GQ / DP / AD gates are judged on. Counted only for the
        models that read a gate this report judges, and not for one that cannot
        run (its candidate sites do not exist). The sites the rarity (INFO/AF) and
        FILTER gates remove are not counted: nobody is shown them."""
        blocked = {model for model, _ in blocked_models(analysis)}
        wanted = [
            m
            for m in candidate_models(
                analysis.gates, selected_categories(options.category), only
            )
            if m not in blocked
        ]
        return collect_candidate_stats(
            sess,
            analysis.trio,
            wanted,
            analysis.gates,
            filter_present=filter_present,
            x_guard=analysis.x_guard_sql(),
            scope=analysis.scope,
        )

    def strict_check(analysis: TrioAnalysis) -> tuple[str, ...]:
        """The runnable selected models a strict_missing gate in force can drop a
        candidate from: after the run, the candidates it dropped for a missing
        field alone are counted (`trio_gate_output.strict_missing_lost`)."""
        blocked = {model for model, _ in blocked_models(analysis)}
        return tuple(
            m
            for m in strict_check_models(
                analysis.gates, selected_categories(options.category)
            )
            if m not in blocked
        )

    def make_report(analysis: TrioAnalysis, gates, candidates):
        trio = analysis.trio
        roles = (
            (trio.proband, "proband"),
            (trio.father, "father"),
            (trio.mother, "mother"),
        )
        return build_report(
            stats=[replace(counted[s], role=role) for s, role in roles if s],
            variants=variants,
            gates=gates,
            gnomad_max_af=analysis.gnomad_max_af,
            gnomad_loaded=gnomad_loaded,
            gnomad_chroms=gnomad_chroms,
            gnomad_probe=gnomad_probe,
            filter_present=filter_present,
            categories=selected_categories(options.category),
            candidates=candidates,
            blocked_models=blocked_models(analysis),
            assay=state.located.preset.assay if state else None,
            subject="preset" if state else "gates",
            skipped=skipped,
            preset_max_af=preset_max_af,
            enforce_bound=enforce_bound,
            enforce_format=enforce_format,
            only=only,
        )

    printed: set[str] = set()
    refusals: list[str] = []
    ready: list[TrioAnalysis] = []
    for analysis in analyses:
        trio = analysis.trio
        candidates = candidate_rows(analysis)
        report = make_report(analysis, analysis.gates, candidates)
        if overridden_off and analysis.gates.strict_missing:
            exempt = unbound_strict_keys(report)
            if exempt:
                gates = analysis.gates._replace(strict_exempt=exempt)
                analysis = replace(analysis, gates=gates)
                report = make_report(analysis, gates, candidates)
        ready.append(analysis)
        strict_models = strict_check(analysis) if only is None else ()
        if strict_models:
            analysis.stats["strict_models"] = strict_models
        if wants_gnomad and only is None:
            analysis.stats["gnomad_reported"] = True
        if strict_models or (wants_gnomad and only is None):
            analysis.stats["text_footer"] = partial(run_footer, analysis)
        if report.errors:
            who = f" (proband {trio.proband})" if len(analyses) > 1 else ""
            refusals.append(refusal_message(report.errors, subject=subject + who))
        analysis.stats["effective_block"] = report.document
        warnings = report.warnings
        for w in warnings:
            if w["level"] == "info" and w["code"] not in _SILENT_NOTES:
                if w["message"] not in printed:
                    printed.add(w["message"])
                    click.echo(f"note: {w['message']}", err=True)
        if state:
            wanted = selected_categories(options.category)
            not_run = [m for m in MODELS if "all" not in wanted and m not in wanted]
            block = preset_block(
                state.located,
                state.plan,
                resolved,
                overridden,
                command=state.command,
                skip_unavailable=state.skip_unavailable,
                bound=report.bound,
                notes=notes
                + [w["message"] for w in warnings if w["code"] in _SILENT_NOTES],
                models_not_run=not_run,
            )
            analysis.stats["preset_block"] = block
            analysis.stats["text_header"] = header_lines(
                block, resolved, analysis.gates.af_gate, warnings
            )
        else:
            analysis.stats["text_header"] = [
                f"  WARNING: {w['message']}"
                for w in warnings
                if w["level"] == "warning"
            ]
    if refusals:
        raise click.ClickException("\n\n".join(refusals))
    for note in state.notes if state else ():
        click.echo(f"note: {note}", err=True)
    return ready
