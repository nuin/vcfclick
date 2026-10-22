"""`db trio` with a preset or a strictness flag: compute, enforce, leave the echo.

For a run that uses a preset, `--strict-missing` or `--require-bound-gates`,
`attach_gate_reports` asks the data whether each gate binds
(`cli.effective_gates`), stops the run when it was told to insist and a gate
cannot bind, and leaves the answers in each analysis' `stats` for the output
hooks of `cli.trio_gate_output` (the `preset` and `effective` JSON objects, the
text header). A plain run never gets here, so its output and its queries are
unchanged.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace

import click

from cli import filter_registry as fr
from cli.effective_gates import (
    build_report,
    collect_sample_stats,
    collect_variant_stats,
    refusal_message,
)
from cli.preset_echo import header_lines, preset_block, skipped_items
from cli.preset_options import PresetState, is_explicit
from cli.trio_analysis import TrioAnalysis
from cli.trio_output import MODELS, _annotation_status, selected_categories

# info-level codes the user does not need to see twice
_SILENT_NOTES = ("preset_max_af_not_applied",)


def needs_report(state: PresetState | None, require_bound: bool, strict: bool) -> bool:
    return state is not None or require_bound or strict


def attach_gate_reports(
    analyses: Sequence[TrioAnalysis],
    ctx: click.Context,
    *,
    options,
    state: PresetState | None,
    resolved: dict[str, fr.Resolved] | None,
    overridden: list[dict],
) -> None:
    """Fill `analysis.stats` for every trio, or raise a ClickException naming
    the gates that cannot bind when the run insists on it."""
    first = analyses[0]
    sess, ingest_id = first.session, first.trio.ingest_id
    variants = collect_variant_stats(sess, ingest_id)
    gnomad_loaded = bool(_annotation_status()["gnomad"])
    from storage import table_columns

    filter_present = "filter" in set(table_columns("variants"))

    # --no-require-bound-gates typed on the command line overrides both kinds
    # of refusal for this run
    overridden_off = is_explicit(ctx, "require_bound") and not options.require_bound
    enforce_bound = bool(options.require_bound)
    enforce_format = bool(options.strict_missing) and not overridden_off

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
        for c in collect_sample_stats(sess, ingest_id, [(s, "") for s in unique])
    }

    printed: set[str] = set()
    refusals: list[str] = []
    for analysis in analyses:
        trio = analysis.trio
        roles = (
            (trio.proband, "proband"),
            (trio.father, "father"),
            (trio.mother, "mother"),
        )
        report = build_report(
            stats=[replace(counted[s], role=role) for s, role in roles if s],
            variants=variants,
            gates=analysis.gates,
            gnomad_max_af=analysis.gnomad_max_af,
            gnomad_loaded=gnomad_loaded,
            filter_present=filter_present,
            kind=trio.kind,
            categories=selected_categories(options.category),
            assay=state.located.preset.assay if state else None,
            subject="preset" if state else "gates",
            skipped=skipped,
            preset_max_af=preset_max_af,
            enforce_bound=enforce_bound,
            enforce_format=enforce_format,
        )
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
