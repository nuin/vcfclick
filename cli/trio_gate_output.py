"""What a trio document and listing gain when a preset or a strict gate is on.

Nothing here runs for a plain `db trio`: every hook reads `analysis.stats`,
which `db trio` fills only for a run with a preset or a strictness flag, so a
plain run's output is the bytes it always was.

  * the document decorator adds the top-level `preset` and `effective` objects,
    and, once the gnomAD filter has run, how many candidates it found no record
    for (and the warning when that is all of them: a fact about the candidates,
    which a de novo list is expected to have, and not about whether the filter
    binds, which the pre-run probe of the cohort's own variants decides);
  * after the run, for a strict_missing run, how many candidates that pass every
    other gate it dropped only for a missing GQ, DP or allele depth (what the
    same run without strict_missing would list beyond it): a stdout
    WARNING and an entry in `effective.warnings`, so that an empty or short list
    is never quietly the data's missing fields and not the genome;
  * `parent_alt_reads` on each de novo candidate (and a suffix on its line in
    the text listing): the parents' ref / alt read counts at the site, so a
    clinician sees at a glance that a "de novo" has 4 alt reads of 40 in the
    father. One batched lookup per rendered page, not one per candidate.

Importing this module registers the hooks.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from cli.effective_gates import gnomad_record_warning, pct, strict_missing_warning
from cli.trio_output import (
    Render,
    register_candidate_decorator,
    register_document_decorator,
)

_CHUNK = 400  # variants per lookup statement


def _key(row: Sequence) -> tuple[str, int, str, str]:
    return (row[0], int(row[1]), row[2], row[3])


def lookup_parent_reads(analysis, keys: Sequence[tuple]) -> dict:
    """{(chrom, pos, ref, alt): {"father": {"ref", "alt"}, "mother": {...}}}
    for the parents' stored rows at these sites (None counts where the call
    carries no allele depths). A parent with no stored row there is absent
    from the answer's entry; a de novo candidate always has both, because the
    model requires the parents' stored 0/0."""
    from storage import sql_quote_str

    trio = analysis.trio
    names = {trio.father: "father", trio.mother: "mother"}
    listed = ", ".join(sql_quote_str(s) for s in names)
    out: dict = {}
    unique = list(dict.fromkeys(keys))
    for start in range(0, len(unique), _CHUNK):
        tuples = ", ".join(
            f"({sql_quote_str(c)}, {int(p)}, {sql_quote_str(r)}, {sql_quote_str(a)})"
            for c, p, r, a in unique[start : start + _CHUNK]
        )
        sql = (
            "SELECT sample_id, chrom, pos, ref, alt, ad_ref, ad_alt FROM genotypes "
            f"WHERE ingest_id = {sql_quote_str(trio.ingest_id)} "
            f"AND sample_id IN ({listed}) "
            f"AND (chrom, pos, ref, alt) IN ({tuples})"
        )
        result = json.loads(analysis.session.query(sql, "JSONCompact").bytes().decode())
        for sample, chrom, pos, ref, alt, ad_ref, ad_alt in result["data"]:
            entry = out.setdefault((chrom, int(pos), ref, alt), {})
            entry[names[sample]] = {
                "ref": None if ad_ref is None else int(ad_ref),
                "alt": None if ad_alt is None else int(ad_alt),
            }
    return out


# ───────────────────────────── candidate decorator ─────────────────────────────


def _stats(render: Render) -> dict:
    """The run's facts; empty for a document about something that is not a
    TrioAnalysis (another command's rows, a test's stand-in)."""
    return getattr(render.analysis, "stats", {})


def _prefetch(render: Render, rows: Sequence[Sequence]) -> None:
    analysis = render.analysis
    if not _stats(render).get("parent_alt_reads"):
        return
    keys = [_key(r) for r in rows if ("parent_alt", _key(r)) not in render.cache]
    if not keys:
        return
    found = lookup_parent_reads(analysis, keys)
    for key in dict.fromkeys(keys):
        render.cache[("parent_alt", key)] = found.get(key, {})


def _parent_alt_reads(d: dict, row: list, model: str, render: Render) -> dict:
    if _stats(render).get("parent_alt_reads"):
        entry = render.cache.get(("parent_alt", _key(row)), {})
        d["parent_alt_reads"] = {
            side: entry.get(side, {"ref": None, "alt": None})
            for side in ("father", "mother")
        }
    return d


register_candidate_decorator(
    _parent_alt_reads, models=frozenset({"denovo"}), prefetch=_prefetch
)


def show_parent_reads(analysis) -> None:
    """Turn on the parents' read counts for this run's de novo candidates, in
    the JSON (`parent_alt_reads`) and in the text listing. `db trio` does it
    when a preset or --parent-alt-screen is on, and only then."""
    analysis.stats["parent_alt_reads"] = True
    analysis.stats["text_suffix"] = text_suffix(analysis)


def text_suffix(analysis):
    """`stats["text_suffix"]`: for the de novo listing, one string per row
    saying the parents' alt / total read counts (`  parents' alt reads:
    father 4/40, mother 0/38`); nothing for another model."""

    def suffix(model: str, rows: list) -> list[str] | None:
        if model != "denovo" or not rows:
            return None
        found = lookup_parent_reads(analysis, [_key(r) for r in rows])
        out = []
        for row in rows:
            entry = found.get(_key(row), {})
            bits = []
            for side in ("father", "mother"):
                reads = entry.get(side)
                if reads is None or reads["alt"] is None or reads["ref"] is None:
                    bits.append(f"{side} NA")
                else:
                    bits.append(f"{side} {reads['alt']}/{reads['ref'] + reads['alt']}")
            out.append("  parents' alt reads: " + ", ".join(bits))
        return out

    return suffix


# ───────────────────────────── document decorator ─────────────────────────────


def gnomad_reported(analysis) -> bool:
    """Is this run's report about the gnomAD filter? (`db trio` sets it for a run
    that asked for gnomAD and for the whole report; a run that asked only for the
    parent-allele screen is told about that gate alone.)"""
    return bool(analysis.stats.get("gnomad_reported")) and (
        analysis.gnomad_max_af is not None
    )


def gnomad_counts(analysis) -> tuple[int, int]:
    """(candidates the gnomAD filter looked at, those it kept for want of a
    record), over everything this run has rendered so far."""
    stats = analysis.stats
    return len(stats.get("gnomad_seen", ())), len(stats.get("gnomad_absent", ()))


def _late_warning(analysis, effective: dict) -> dict | None:
    """The `gnomad_no_record` warning, once the run knows: only when gnomAD is
    loaded (an unloaded one has its own warning) and was asked for."""
    if analysis.gnomad_max_af is None or not effective["gnomad"]["loaded"]:
        return None
    checked, without = gnomad_counts(analysis)
    sampled = effective["gnomad"].get("cohort_variants_sampled")
    found = effective["gnomad"].get("cohort_variants_with_record")
    return gnomad_record_warning(
        checked=checked,
        without_record=without,
        max_af=analysis.gnomad_max_af,
        info_af_bound=bool(effective["info_af"]["bound"]),
        cohort_rate=found / sampled if sampled else None,
    )


def strict_missing_lost(analysis) -> dict[str, int]:
    """{model: candidates it dropped} for each model `db trio` marked for the check
    (`stats["strict_models"]`) that strict_missing dropped a candidate from: the
    ones that pass every other gate (the gnomAD filter too) and fail only because
    a GQ, DP or allele depth a strict gate needs is missing from their rows. One
    count per model (`TrioAnalysis.dropped_by_strict`), once per analysis."""
    cached = analysis.stats.get("strict_lost")
    if cached is not None:
        return cached
    lost = {}
    for model in analysis.stats.get("strict_models", ()):
        dropped = analysis.dropped_by_strict(model)
        if dropped:
            lost[model] = dropped
    analysis.stats["strict_lost"] = lost
    return lost


def strict_missing_warnings(analysis) -> list[dict]:
    """The `strict_missing_dropped` warnings of a finished run, one per model."""
    return [
        w
        for model, dropped in strict_missing_lost(analysis).items()
        if (w := strict_missing_warning(model, dropped)) is not None
    ]


def strict_missing_footer(analysis) -> list[str]:
    """The `  WARNING:` lines a text run ends with for a strict_missing run."""
    if not analysis.stats.get("strict_models"):
        return []
    return [f"  WARNING: {w['message']}" for w in strict_missing_warnings(analysis)]


def run_footer(analysis) -> list[str]:
    """Everything a text run ends with that only the finished run knows: the
    gnomAD counts and the candidates strict_missing dropped."""
    return gnomad_footer(analysis) + strict_missing_footer(analysis)


def gnomad_footer(analysis) -> list[str]:
    """The lines a text run ends with about the gnomAD filter: a WARNING when
    it found a record for none of the candidates, otherwise how many it could
    not judge. (The `trio:` header is printed before the filter has run.)"""
    effective = analysis.stats.get("effective_block")
    if effective is None or not gnomad_reported(analysis):
        return []
    late = _late_warning(analysis, effective)
    if late is not None:
        return [f"  WARNING: {late['message']}"]
    checked, without = gnomad_counts(analysis)
    if effective["gnomad"]["loaded"] and without:
        sampled = effective["gnomad"].get("cohort_variants_sampled")
        found = effective["gnomad"].get("cohort_variants_with_record")
        reach = (
            f"; it has a record for {pct(found / sampled)} of a sample of the "
            "cohort's variants"
            if sampled
            else ""
        )
        return [
            f"  gnomAD: no record for {without} of {checked} candidates "
            f"(kept as rare){reach}"
        ]
    return []


def _preset_document(doc: dict, render: Render) -> dict:
    stats = _stats(render)
    if "preset_block" in stats:
        doc["preset"] = stats["preset_block"]
    if "effective_block" in stats:
        effective = json.loads(json.dumps(stats["effective_block"]))
        analysis = render.analysis
        # the counts are about a loaded store: over an empty one every candidate is
        # "without a record" and the number says nothing
        if gnomad_reported(analysis) and effective["gnomad"]["loaded"]:
            checked, without = gnomad_counts(analysis)
            effective["gnomad"]["candidates_checked"] = checked
            effective["gnomad"]["candidates_without_record"] = without
            late = _late_warning(analysis, effective)
            if late is not None:
                effective["warnings"].append(late)
                if not effective["info_af"]["bound"]:
                    effective["population_filter"] = False
        if stats.get("strict_models"):
            effective["warnings"].extend(strict_missing_warnings(analysis))
        doc["effective"] = effective
    return doc


register_document_decorator(_preset_document)
