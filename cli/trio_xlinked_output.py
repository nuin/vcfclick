"""Text and JSON of the X-linked model (`models.xlinked`, `--category xlinked`).

Registered by cli/trio_output.py after the four classic models. Its rows are
`cli.xlinked.XRow`s (one more column than a classic row), so they never go
through the 8-column candidate decorators: they get the per-variant columns
(genes, gnomAD popmax, ClinVar) and then the X keys (mode, pattern,
hemizygous, father_flag, parent_affected).

`models.xlinked` is present for every trio. With no chrX data it is
`blocked: true, blocked_reason: no_x_data`, and the `all` text summary has
no xlinked line.
"""

from __future__ import annotations

import click

from cli.trio_output import (
    ModelSpec,
    Render,
    _comphet_entry,
    _num,
    decorate_rows,
    outside_note,
)
from cli.xlinked import (
    APPLICABLE,
    MODES,
    NEEDS_REFERENCE,
    father_flag,
    hemizygous_roles,
    parent_affected,
)

# What a blocked X model says on its `all` line and in its listing.
_BLOCKED_TEXT = {
    "needs_sex": "needs proband sex: set it in the PED, pass db ped --sex, or "
    "ingest more chrX calls",
    "needs_parents": "needs both parents",
}


def _context_json(analysis) -> dict:
    ctx = analysis.x_context()
    sexes = ctx.sexes
    return {
        "build": ctx.build,
        "build_source": ctx.build_source,
        "proband_sex": sexes.proband_sex,
        "proband_sex_source": sexes.proband_sex_source,
        "father_sex": sexes.father_sex,
        "father_sex_source": sexes.father_sex_source,
        "mother_sex": sexes.mother_sex,
        "mother_sex_source": sexes.mother_sex_source,
        "warnings": list(ctx.warnings),
        "excluded": analysis.x_excluded(),
    }


def mode_table(analysis, genes_loaded: bool) -> dict:
    """{mode: {applicable, count, blocked, blocked_reason}} for every mode."""
    ctx = analysis.x_context()
    applicable = APPLICABLE.get(ctx.proband_sex or "", ())
    rows = analysis.x_rows()
    out = {}
    for mode in MODES:
        reason = None
        if mode in applicable:
            if mode in NEEDS_REFERENCE and not analysis.has_reference:
                reason = "needs_keep_reference"
            elif mode == "comphet" and not genes_loaded:
                reason = "needs_genes"
        if mode == "comphet":
            count = len(analysis.x_genes()) if mode in applicable else 0
        else:
            count = sum(1 for r in rows if r.mode == mode)
        out[mode] = {
            "applicable": mode in applicable,
            "count": count,
            "blocked": reason is not None,
            "blocked_reason": reason,
        }
    return out


def _candidate_x(out: dict, row, analysis) -> dict:
    ctx = analysis.x_context()
    out["mode"] = row.mode
    out["pattern"] = row.pattern
    out["hemizygous"] = hemizygous_roles(ctx.proband_sex)
    out["father_flag"] = father_flag(row, ctx.proband_sex)
    out["parent_affected"] = parent_affected(row, ctx.affected)
    return out


def _gt(v) -> int:
    return -1 if v is None else int(v)


def _build_xlinked(render: Render) -> dict:
    analysis, limit = render.analysis, render.limit
    ctx = analysis.x_context()
    context = _context_json(analysis)
    if ctx.blocked_reason is not None:
        return {
            "count": 0,
            "blocked": True,
            "blocked_reason": ctx.blocked_reason,
            "truncated": False,
            "needs_genes": False,
            "context": context,
            "modes": {},
            "candidates": [],
            "genes": [],
        }
    modes = mode_table(analysis, render.have["genes"])
    rows = analysis.x_rows()
    shown = rows[:limit]
    base = [
        {
            "chrom": r.chrom,
            "pos": int(r.pos),
            "ref": r.ref,
            "alt": r.alt,
            "proband_gt": int(r.proband_gt),
            "father_gt": _gt(r.father_gt),
            "mother_gt": _gt(r.mother_gt),
            "af": _num(r.af),
        }
        for r in shown
    ]
    outs = decorate_rows(base, shown, render)
    candidates = [
        _candidate_x(o, r, analysis) for o, r in zip(outs, shown, strict=True)
    ]
    genes = analysis.x_genes() if not modes["comphet"]["blocked"] else {}
    runnable = [m for m, v in modes.items() if v["applicable"] and m != "comphet"]
    blocked = all(modes[m]["blocked"] for m in runnable)
    return {
        "count": len(rows),
        "blocked": blocked,
        "blocked_reason": "needs_keep_reference" if blocked else None,
        "truncated": len(rows) > limit,
        "needs_genes": modes["comphet"]["applicable"] and not render.have["genes"],
        "context": context,
        "modes": modes,
        "candidates": candidates,
        "genes": [
            _comphet_entry(sym, genes[sym], analysis.scope)
            for sym in sorted(genes)[:limit]
        ],
    }


def _blocked_text(ctx) -> str:
    if ctx.blocked_reason == "parent_sex_conflict":
        return f"parent sex conflict: {ctx.conflict}"
    if ctx.blocked_reason == "build_conflict":
        return f"build conflict: {ctx.conflict}"
    return _BLOCKED_TEXT.get(ctx.blocked_reason, ctx.blocked_reason or "")


def _summary_xlinked(analysis) -> str | None:
    ctx = analysis.x_context()
    if not ctx.has_x:
        return None
    echo_x_warnings(analysis)
    if ctx.blocked_reason is not None:
        return f"  {'xlinked':10s} {'-':>6}  ({_blocked_text(ctx)})"
    modes = mode_table(analysis, genes_loaded=True)
    shown = [
        f"{m} {v['count']}"
        for m, v in modes.items()
        if v["applicable"]
        and not v["blocked"]
        and m != "comphet"
        and (v["count"] or m in ("recessive", "denovo"))
    ]
    blocked = [
        m
        for m, v in modes.items()
        if v["blocked_reason"] == "needs_keep_reference"
        and m in ("female_het", "denovo")
    ]
    note = ", ".join(shown)
    if blocked:
        note += f"; {', '.join(blocked)} needs --keep-reference"
    n = len(analysis.x_rows())
    return f"  {'xlinked':10s} {n:>6}  ({note})"


def _who(sample: str, sex: str | None, source: str) -> str:
    return f"{sample} {sex or 'unknown'} [{source}]"


def echo_x_warnings(analysis) -> None:
    for w in analysis.x_context().warnings:
        click.echo(f"warning: {w['code']}: {w['message']}", err=True)


def _detail_xlinked(analysis) -> None:
    ctx = analysis.x_context()
    trio = analysis.trio
    echo_x_warnings(analysis)
    if ctx.blocked_reason is not None:
        click.echo(f"\nxlinked: blocked ({_blocked_text(ctx)})")
        return
    s = ctx.sexes
    rows = analysis.x_rows()
    who = ", ".join(
        [
            f"proband={_who(trio.proband, s.proband_sex, s.proband_sex_source)}",
            f"father={_who(trio.father, s.father_sex, s.father_sex_source)}",
            f"mother={_who(trio.mother, s.mother_sex, s.mother_sex_source)}",
        ]
    )
    click.echo(
        f"\nxlinked candidates: {len(rows)}  ({who}; {ctx.build} [{ctx.build_source}])"
    )
    if not analysis.has_reference:
        click.echo(
            "  note: no stored hom-reference calls: X de novo and female_het "
            "need a --keep-reference ingest",
            err=True,
        )
    for r in rows:
        af = "NA" if r.af is None else f"{r.af}"
        fgt = "." if r.father_gt is None else r.father_gt
        mgt = "." if r.mother_gt is None else r.mother_gt
        flag = father_flag(r, ctx.proband_sex)
        click.echo(
            f"  {r.chrom}:{r.pos} {r.ref}>{r.alt}  proband_gt={r.proband_gt} "
            f"father_gt={fgt} mother_gt={mgt}  AF={af}  mode={r.mode} "
            f"pattern={r.pattern}" + (f" father={flag}" if flag else "")
        )
    if ctx.proband_sex == "female" and analysis.has_reference:
        from cli.trio_output import _annotation_status

        genes = analysis.x_genes()
        loaded = _annotation_status()["genes"]
        click.echo(
            f"X compound-het genes: {len(genes)}"
            + (
                ""
                if loaded
                else "  (needs gene annotations: vcfclick annotations load)"
            )
        )
        for sym in sorted(genes):
            click.echo(f"  {sym}")
            for origin in ("paternal", "maternal"):
                for chrom, pos, ref, alt, af in genes[sym][origin]:
                    af_s = "NA" if af is None else f"{af}"
                    click.echo(
                        f"    {origin:8s} {chrom}:{pos} {ref}>{alt}  AF={af_s}"
                        + outside_note(analysis, chrom, pos)
                    )
    excluded = {k: v for k, v in analysis.x_excluded().items() if v}
    if excluded:
        click.echo("  excluded: " + ", ".join(f"{k}={v}" for k, v in excluded.items()))


def _xlinked_scope_sql(trio, gates, scope) -> list[str]:
    """Every statement the X model can run under a scope: both sexes' rows and
    the daughter's compound-het candidates."""
    from cli.xchrom import GRCH38
    from cli.xlinked import comphet_x_sql, xlinked_sql

    return [
        xlinked_sql(trio, gates, "male", GRCH38, scope=scope),
        xlinked_sql(trio, gates, "female", GRCH38, scope=scope),
        comphet_x_sql(trio, gates, GRCH38, scope=scope),
    ]


XLINKED_MODEL = ModelSpec(
    "xlinked",
    lambda trio: True,
    _build_xlinked,
    _summary_xlinked,
    _detail_xlinked,
    scope_builder=_xlinked_scope_sql,
)
