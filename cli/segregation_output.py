"""Text presentation of a segregation document (the JSON is built by
cli/segregation.py; the text is derived from it, so the two cannot differ).

    segregation: family=fam1 father=FA mother=MO affected=S1,S2 unaffected=S3,S4 unknown=S5
      denovo          1  (+2 compatible)
      ...
    chr1:1000000 A>G  FA=0/1 MO=0/0  S1=0/1* S2=0/1* S3=0/0 S4=0/0  consistent

A `*` marks an affected child. A state with no stored genotype prints `-`
(absent_called), `./.` (no_call) or `?` (absent_unknown); a stored row that
failed the gates prints its genotype with `?`; a stored 0/0 that cannot be
trusted (partial-call site, alt reads) prints `0/0!`.
"""

from __future__ import annotations

from collections.abc import Sequence

import click

from cli import segregation_verdict as sv

_GT = {0: "0/0", 1: "0/1", 2: "1/1"}
_BY_STATE = {
    "absent_called": "-",
    "no_call": "./.",
    "absent_unknown": "?",
    "ref_partial_site": "0/0!",
    "ref_alt_reads": "0/0!",
}


def gt_text(state: str, gt: int | None) -> str:
    """One member's genotype as the text listing prints it."""
    if state in _BY_STATE:
        return _BY_STATE[state]
    text = _GT.get(gt if gt is not None else -1, "?")
    return text + "?" if state == "low_quality" else text


def header(doc: dict) -> str:
    fam = doc["family"]
    kids = fam["children"]

    def names(phenotype):
        return ",".join(c["sample"] for c in kids if c["affected"] == phenotype)

    unknown = ",".join(
        c["sample"] for c in kids if c["affected"] not in ("affected", "unaffected")
    )
    line = (
        f"segregation: family={fam['family_id']} father={fam['father']} "
        f"mother={fam['mother']} affected={names('affected')} "
        f"unaffected={names('unaffected')}"
    )
    return line + (f" unknown={unknown}" if unknown else "")


_BLOCKED = {
    "needs_keep_reference": "(needs --keep-reference or --absent called)",
    "needs_genes": "(needs gene annotations: `vcfclick annotations load`)",
}


def _blocked_text(model: str, reason: str) -> str:
    if model == "denovo" and reason == "needs_keep_reference":
        return "(parents must be stored reference: needs --keep-reference)"
    return _BLOCKED.get(reason, f"({reason})")


def summary_line(doc: dict, model: str) -> str:
    m = doc["models"][model]
    extra = ""
    if m["blocked"]:
        extra = "  " + _blocked_text(model, m["blocked_reason"])
    else:
        by = m["by_status"]
        if m["n_compatible"] and not doc["options"]["allow_unresolved"]:
            extra += f"  (+{m['n_compatible']} compatible)"
        if by[sv.REDUCED_PENETRANCE]:
            extra += f"  (incl. {by[sv.REDUCED_PENETRANCE]} reduced penetrance)"
    unit = "  genes" if model == "comphet" else ""
    return f"  {model:10s} {m['count']:>6}{unit}{extra}"


def _calls(calls: Sequence[dict], affected: dict[str, str | None]) -> str:
    return " ".join(
        f"{c['sample']}={gt_text(c['state'], c['gt'])}"
        + ("*" if affected.get(c["sample"]) == "affected" else "")
        for c in calls
    )


def _row_line(doc: dict, c: dict) -> str:
    fam = doc["family"]
    affected = {k["sample"]: k["affected"] for k in fam["children"]}
    parents = (
        f"{fam['father']}={gt_text(c['father_state'], c['father_gt'])} "
        f"{fam['mother']}={gt_text(c['mother_state'], c['mother_gt'])}"
    )
    return (
        f"  {c['chrom']}:{c['pos']} {c['ref']}>{c['alt']}  {parents}  "
        f"{_calls(c['calls'], affected)}  {c['status']}"
    )


def detail(doc: dict, model: str) -> list[str]:
    """The listing of one model."""
    m = doc["models"][model]
    if m["blocked"]:
        return [f"\n{model}: blocked {_blocked_text(model, m['blocked_reason'])}"]
    fam = doc["family"]
    affected = {k["sample"]: k["affected"] for k in fam["children"]}
    if model != "comphet":
        lines = [f"\n{model} candidates: {m['count']}"]
        lines += [_row_line(doc, c) for c in m["candidates"]]
        return lines
    lines = [f"\ncomphet candidate genes: {m['count']}"]
    for g in m["genes"]:
        lines.append(f"  {g['gene']}")
        for origin in ("paternal", "maternal"):
            for v in g[origin]:
                lines.append(
                    f"    {origin:8s} {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}  "
                    + _calls(v["calls"], affected)
                )
        for p in g["pairs"]:
            lines.append(f"    pair {p['paternal']} + {p['maternal']}  {p['status']}")
    return lines


def text_report(doc: dict, category: Sequence[str], header_lines: Sequence[str] = ()):
    """Print one sibship: the header, any preset / scope lines, then a line per
    model (several models) or the listing (one)."""
    click.echo(header(doc))
    for line in header_lines:
        click.echo(line)
    models = list(doc["models"])
    if len(models) == 1:
        for line in detail(doc, models[0]):
            click.echo(line)
    else:
        for model in models:
            click.echo(summary_line(doc, model))
