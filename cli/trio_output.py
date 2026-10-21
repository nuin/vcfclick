"""Text and JSON presentation of trio inheritance candidates.

Each inheritance model is a `ModelSpec` in a registry: how it appears in a
trio document (`build`), as a count line (`summary`) and as a listing
(`detail`). `_trio_json` and `trio_text` only walk the registry, so a new
model is one `register_model` call and never an edit to the shared
presentation code. A candidate row passes through a chain of decorators
(`register_candidate_decorator`); the first one adds the annotation columns
(genes, gnomAD popmax, ClinVar) and a later feature appends its own keys by
registering another, after the existing ones.

JSON compatibility policy (docs/TRIO.md): additive only. A key that is
already emitted is never removed, renamed or retyped, and the relative
order of the existing keys never changes; new keys are snake_case and are
either appended or placed directly after the key they qualify
(`blocked_reason` follows `blocked`). tests/test_trio_goldens.py pins this
byte-for-byte once the whitelisted new keys are removed.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import click

from cli.trio_analysis import (
    CATEGORIES,
    TrioAnalysis,
    check_blocked_reason,
)


def _annotation_store() -> tuple[dict, object | None]:
    """(which annotation sources hold data, the open store connection).

    Each lookup is skipped when its source is empty. The connection is
    returned so the caller can keep it open: DuckDB shares an open database
    between connections, so every later `get_connection()` is about ten
    times cheaper (0.7 ms against 7 ms with genes loaded) than when the
    last connection has been closed and the file has to be opened again.
    """
    none = {"genes": False, "clinvar": False, "gnomad": False}
    try:
        from annotations.db import get_connection

        conn = get_connection()
    except Exception:  # noqa: BLE001 - no annotation store at all
        return none, None
    have = {}
    for key, table in (
        ("genes", "refseq_genes"),
        ("clinvar", "clinvar_variants"),
        ("gnomad", "gnomad_af"),
    ):
        try:
            have[key] = conn.execute(f"SELECT count(*) FROM {table}").fetchone()[0] > 0
        except Exception:  # noqa: BLE001 - table absent in older stores
            have[key] = False
    return have, conn


def _annotation_status() -> dict:
    """Which annotation sources hold data (each lookup is skipped if not)."""
    have, conn = _annotation_store()
    if conn is not None:
        conn.close()
    return have


def _enrich(chrom: str, pos, ref: str, alt: str, have: dict) -> dict:
    from annotations import clinvar_lookup, gene_at, gnomad_af

    out = {"genes": [], "gnomad_popmax": None, "clinvar": None}
    if have["genes"]:
        out["genes"] = [g.gene_symbol for g in gene_at(chrom, int(pos))]
    if have["gnomad"]:
        g = gnomad_af(chrom, int(pos), ref, alt)
        out["gnomad_popmax"] = None if g is None else g.popmax
    if have["clinvar"]:
        c = clinvar_lookup(chrom, int(pos), ref, alt)
        out["clinvar"] = None if c is None else c.clin_sig
    return out


def _num(v):
    if v is None:
        return None
    f = float(v)
    return int(f) if f.is_integer() else f


@dataclass(frozen=True)
class Render:
    """What building one trio document shares across its models."""

    analysis: TrioAnalysis
    limit: int
    have: dict  # annotation sources that hold data
    store: object | None = None  # open annotation connection, kept alive only
    # Annotation answers prefetched in batches per model (see _prefetch), so
    # a rendered candidate costs a dict lookup, not a store query.
    cache: dict = field(default_factory=dict)


def _prefetch(render: Render, rows: list) -> None:
    """Look up genes, gnomAD and ClinVar for every row about to be rendered
    in one batch per source. A per-candidate store query cost 0.1-1.5 ms
    each (several seconds for a 1,000-candidate model on a slow disk); the
    batch is one join. Rows already cached are skipped."""
    from annotations.db import GeneIndex, clinvar_many, gnomad_af_many

    keys = [(r[0], int(r[1]), r[2], r[3]) for r in rows]
    todo = [k for k in dict.fromkeys(keys) if k not in render.cache]
    if not todo:
        return
    conn = render.store
    genes = gnomad = clinvar = None
    if render.have["genes"]:
        index = render.cache.get("__gene_index__")
        if index is None:
            index = render.cache["__gene_index__"] = GeneIndex(conn)
        genes = index
    if render.have["gnomad"]:
        gnomad = gnomad_af_many(todo, conn=conn)
    if render.have["clinvar"]:
        clinvar = clinvar_many(todo, conn=conn)
    for k in todo:
        g = None if gnomad is None else gnomad.get(k)
        c = None if clinvar is None else clinvar.get(k)
        render.cache[k] = {
            "genes": []
            if genes is None
            else [x.gene_symbol for x in genes.at(k[0], k[1])],
            "gnomad_popmax": None if g is None else g.popmax,
            "clinvar": None if c is None else c.clin_sig,
        }


# A decorator receives the candidate dict built so far, the source row and
# the model name, and returns the dict to emit (usually the same one with
# keys appended).
CandidateDecorator = Callable[[dict, list, str, Render], dict]
_CANDIDATE_DECORATORS: list[tuple[CandidateDecorator, frozenset[str] | None]] = []


def register_candidate_decorator(
    fn: CandidateDecorator, models: frozenset[str] | None = None
) -> CandidateDecorator:
    """Append `fn` to the candidate chain (all models, or only `models`).
    Order of registration is the order keys appear in the output."""
    _CANDIDATE_DECORATORS.append((fn, models))
    return fn


def _candidates(rows: list, model: str, render: Render) -> list[dict]:
    """Rendered candidates for `rows`, annotations fetched in one batch."""
    _prefetch(render, rows)
    return [_candidate(r, model, render) for r in rows]


def _candidate(row: list, model: str, render: Render) -> dict:
    """One 8-column detail row as the JSON candidate object."""
    chrom, pos, ref, alt, pgt, fgt, mgt, af = row
    out: dict = {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "proband_gt": int(pgt),
        "father_gt": int(fgt),
        "mother_gt": int(mgt),
        "af": _num(af),
    }
    for fn, models in _CANDIDATE_DECORATORS:
        if models is None or model in models:
            out = fn(out, row, model, render)
    return out


@register_candidate_decorator
def _annotation_columns(out: dict, row: list, model: str, render: Render) -> dict:
    chrom, pos, ref, alt = row[:4]
    hit = render.cache.get((chrom, int(pos), ref, alt))
    out.update(
        dict(hit) if hit is not None else _enrich(chrom, pos, ref, alt, render.have)
    )
    return out


@dataclass(frozen=True)
class ModelSpec:
    """One inheritance model's place in a trio document and in text output."""

    name: str
    # Does the model exist for this trio? (Later models are conditional.)
    applies: Callable[[object], bool]
    # The JSON object for `models.<name>`.
    build: Callable[[Render], dict]
    # One line for `--category all`, and the listing for `--category <name>`.
    summary: Callable[[TrioAnalysis], str]
    detail: Callable[[TrioAnalysis], None]
    # Needs stored hom-reference calls to prove a parent is 0/0 (the
    # `--keep-reference` note, and the `blocked` flag, key off this).
    needs_reference: bool = False


MODELS: dict[str, ModelSpec] = {}


def register_model(spec: ModelSpec) -> ModelSpec:
    """Add a model after the existing ones: registration order is the order
    of `models` in the document and of the lines in the `all` summary."""
    if spec.name in MODELS:
        raise ValueError(f"model {spec.name!r} is already registered")
    MODELS[spec.name] = spec
    return spec


def _legacy_blocked(needs_reference: bool, analysis: TrioAnalysis) -> bool:
    """The `blocked` flag of the four classic models: no stored hom-ref rows
    anywhere in the database (cli/trio_analysis._has_reference_rows)."""
    return needs_reference and not analysis.has_reference


def _model_reason(blocked: bool, extra: str | None = None) -> str | None:
    """`blocked_reason` for a classic model: the reference requirement when
    the flag says blocked, else whatever other requirement is unmet."""
    return check_blocked_reason("needs_keep_reference" if blocked else extra)


# ───────────────────── de novo / recessive / dominant ─────────────────────


def _row_model(name: str, *, needs_reference: bool) -> ModelSpec:
    def build(render: Render) -> dict:
        analysis, limit = render.analysis, render.limit
        rows = analysis.detail_rows(name)
        blocked = _legacy_blocked(needs_reference, analysis)
        return {
            "count": len(rows),
            "blocked": blocked,
            "blocked_reason": _model_reason(blocked),
            "truncated": len(rows) > limit,
            "candidates": _candidates(rows[:limit], name, render),
        }

    def summary(analysis: TrioAnalysis) -> str:
        n = analysis.count(name)
        note = (
            ""
            if analysis.has_reference or not needs_reference
            else "  (needs --keep-reference)"
        )
        return f"  {name:10s} {n:>6}{note}"

    def detail(analysis: TrioAnalysis) -> None:
        data = analysis.detail_rows(name)
        click.echo(f"\n{name} candidates: {len(data)}")
        for row in data:
            chrom, pos, ref, alt, pgt, fgt, mgt, af = row
            af_s = "NA" if af is None else f"{af}"
            click.echo(
                f"  {chrom}:{pos} {ref}>{alt}  proband_gt={pgt} "
                f"father_gt={fgt} mother_gt={mgt}  AF={af_s}"
            )

    return ModelSpec(
        name,
        lambda trio: True,
        build,
        summary,
        detail,
        needs_reference=needs_reference,
    )


# ───────────────────────────── compound het ─────────────────────────────


def _comphet_entry(sym: str, sides: dict) -> dict:
    """One gene of the comphet model. `paternal` and `maternal` are emitted
    for the sides present; a consumer must treat each as optional (a model
    that cannot tell parent of origin omits them)."""
    entry: dict = {"gene": sym}
    for origin in ("paternal", "maternal"):
        if sides.get(origin) is None:
            continue
        entry[origin] = [
            {"chrom": c, "pos": int(p), "ref": r, "alt": a, "af": _num(af)}
            for c, p, r, a, af in sides[origin]
        ]
    return entry


def _build_comphet(render: Render) -> dict:
    analysis, limit = render.analysis, render.limit
    genes = analysis.comphet_genes()
    entries = [_comphet_entry(sym, genes[sym]) for sym in sorted(genes)[:limit]]
    blocked = _legacy_blocked(True, analysis)
    needs_genes = not render.have["genes"]
    return {
        "count": len(genes),
        "blocked": blocked,
        "blocked_reason": _model_reason(
            blocked, "needs_genes" if needs_genes else None
        ),
        "needs_genes": needs_genes,
        "truncated": len(genes) > limit,
        "genes": entries,
    }


def _summary_comphet(analysis: TrioAnalysis) -> str:
    genes = analysis.comphet_genes()
    note = "" if analysis.has_reference else "  (needs --keep-reference)"
    return f"  {'comphet':10s} {len(genes):>6}{note}  genes"


def _detail_comphet(analysis: TrioAnalysis) -> None:
    genes = analysis.comphet_genes()
    click.echo(f"\ncomphet candidate genes: {len(genes)}")
    for sym in sorted(genes):
        entry = genes[sym]
        click.echo(f"  {sym}")
        for origin in ("paternal", "maternal"):
            for chrom, pos, ref, alt, af in entry[origin]:
                af_s = "NA" if af is None else f"{af}"
                click.echo(f"    {origin:8s} {chrom}:{pos} {ref}>{alt}  AF={af_s}")


register_model(_row_model("denovo", needs_reference=True))
register_model(_row_model("recessive", needs_reference=False))
register_model(_row_model("dominant", needs_reference=True))
register_model(
    ModelSpec(
        "comphet",
        lambda trio: True,
        _build_comphet,
        _summary_comphet,
        _detail_comphet,
        needs_reference=True,
    )
)

# The registry and the shared CATEGORIES constant describe the same models.
if tuple(MODELS) != CATEGORIES:
    raise RuntimeError(f"model registry {tuple(MODELS)} != CATEGORIES {CATEGORIES}")


# ───────────────────────────── the document ─────────────────────────────


def _applicable(analysis: TrioAnalysis, category: str) -> list[ModelSpec]:
    return [
        spec
        for spec in MODELS.values()
        if category in ("all", spec.name) and spec.applies(analysis.trio)
    ]


def _trio_json(analysis: TrioAnalysis, category: str, limit: int) -> dict:
    trio, gates = analysis.trio, analysis.gates
    have, store = _annotation_store()
    render = Render(analysis, limit, have, store)
    try:
        models = {
            spec.name: spec.build(render) for spec in _applicable(analysis, category)
        }
    finally:
        if store is not None:
            store.close()
    return {
        "trio": {
            "ingest_id": trio.ingest_id,
            "proband": trio.proband,
            "father": trio.father,
            "mother": trio.mother,
        },
        "keep_reference": analysis.has_reference,
        "gates": {
            "min_gq": gates.min_gq,
            "min_dp": gates.min_dp,
            "max_af": gates.max_af,
            "min_ab": gates.min_ab,
            "max_ab": gates.max_ab,
            "gnomad_max_af": analysis.gnomad_max_af,
            "af_gate": gates.af_gate,
        },
        "annotations": render.have,
        "models": models,
    }


def trio_text(analysis: TrioAnalysis, category: str) -> None:
    trio = analysis.trio
    proband, father, mother = trio.proband, trio.father, trio.mother
    click.echo(f"trio: proband={proband} father={father} mother={mother}")
    needs_reference = {s.name for s in MODELS.values() if s.needs_reference}
    if not analysis.has_reference and (
        category in needs_reference or category == "all"
    ):
        click.echo(
            "  note: this database has no stored hom-reference calls, so "
            "de-novo/dominant/comphet cannot prove a parent is 0/0 (vs "
            "no-call). Re-ingest with `--keep-reference` for those models.",
            err=True,
        )

    if category == "all":
        for spec in _applicable(analysis, "all"):
            click.echo(spec.summary(analysis))
        return
    MODELS[category].detail(analysis)
