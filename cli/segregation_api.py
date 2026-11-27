"""`db segregate` for callers that are not the command line (the web endpoint
and the MCP tool): the same sibship resolution, gates, AF guard and document,
from plain arguments, on the active database.

It takes no preset and no gene scope (those are command-line features), reads
only, and returns the JSON document `db segregate --format json` prints: one
document when a selector names the sibship, `{"families": [...]}` otherwise.
A refusal that the command line prints as an error raises
`click.ClickException`; the caller turns its message into its own error.
"""

from __future__ import annotations

import click


def segregate_api(
    *,
    family: str | None = None,
    father: str | None = None,
    mother: str | None = None,
    category: str = "all",
    min_gq: int = 20,
    min_dp: int = 10,
    max_af: float | None = None,
    min_ab: float = 0.25,
    max_ab: float = 0.75,
    absent: str = "strict",
    allow_unresolved: bool = False,
    max_unaffected_carriers: int = 0,
    limit: int = 200,
) -> dict:
    from cli import segregation as seg
    from cli.family_checks import trio_af_gate
    from cli.filter_registry import DEFAULT_MAX_AF
    from cli.pedigree_load import load_pedigree_graph
    from cli.trio_analysis import Gates, _has_reference_rows, _sole_ingest_id
    from cli.trio_output import open_render
    from storage import get_readonly_session

    if absent not in seg.ABSENT_CHOICES:
        raise click.ClickException(
            f"absent must be one of {', '.join(seg.ABSENT_CHOICES)}"
        )
    if limit < 1:
        raise click.ClickException("limit must be at least 1")
    models = seg.select_models([t.strip() for t in category.split(",")])
    ingest_id = _sole_ingest_id(None)
    if ingest_id is None:
        raise click.ClickException(
            "this database has zero or multiple ingestions; segregation "
            "analysis needs a single-ingestion cohort"
        )
    sess = get_readonly_session()
    graph = load_pedigree_graph(sess, ingest_id)
    sibships, notes = seg.find_sibships(
        graph, family=family, father=father, mother=mother
    )
    sibships = [seg.with_ingest(s, ingest_id) for s in sibships]
    for sib in sibships:
        seg.check_unaffected_carriers(sib, max_unaffected_carriers)
    threshold = DEFAULT_MAX_AF if max_af is None else max_af
    af_gate, af_message = trio_af_gate(
        sess, ingest_id, threshold, explicit=max_af is not None
    )
    run_notes = [af_message] if af_message else []
    gates = Gates(min_gq, min_dp, threshold, min_ab, max_ab)._replace(af_gate=af_gate)
    opts = seg.SegOptions(
        models=models,
        gates=gates,
        absent=absent,
        allow_unresolved=allow_unresolved,
        max_unaffected_carriers=max_unaffected_carriers,
        limit=limit,
    )
    has_ref = _has_reference_rows(sess)
    caps = seg.probe_caps()
    n_samples = seg.ingest_n_samples(sess, ingest_id)
    docs = []
    with open_render(None, limit) as render:
        for sib in sibships:
            ctx = seg.RunContext(
                sess=sess,
                sib=sib,
                opts=opts,
                has_reference=has_ref,
                caps=caps,
                n_samples=n_samples,
            )
            docs.append(seg.segregation_document(ctx, render, run_notes))
    if notes:
        for doc in docs:
            doc["notes"] = [*doc["notes"], *notes]
    selecting = any((family, father, mother))
    return docs[0] if selecting else {"families": docs}
