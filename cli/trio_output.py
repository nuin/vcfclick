"""Text and JSON presentation of trio inheritance candidates.

Each inheritance model is a `ModelSpec` in a registry: how it appears in a
trio document (`build`), as a count line (`summary`) and as a listing
(`detail`). `_trio_json` and `trio_text` only walk the registry, so a new
model is one `register_model` call and never an edit to the shared
presentation code.

A candidate row passes through one ordered chain of decorators, and a later
feature appends its own keys by registering another, after the existing
ones. There are two kinds:

  * `register_row_decorator`: a per-variant column that needs only the
    variant's key (chrom, pos, ref, alt). The annotation columns (genes,
    gnomAD popmax, ClinVar) are one. Any command that renders variant rows,
    not only `db trio`, applies exactly these with `decorate_row`.
  * `register_candidate_decorator`: a key that needs the 8-column trio row
    and the model name. Only `db trio` candidates get these.

`register_document_decorator` edits the finished trio document as a whole.

A decorator that needs a lookup should not make one per row (0.1-1.5 ms
each, seconds for a long listing). Both registrars take an optional
`prefetch(render, rows)`, called once per batch of rows before any row of it
is decorated: it does the lookup for the whole batch and leaves the answers
in `render.cache`, and the per-row call only reads them. A row decorator's
prefetch receives the batch's distinct variant keys; a candidate decorator's
receives the trio rows of the one model being rendered (and is skipped for a
model outside its `models`). `db trio` runs these per model; another command
runs them with `decorate_rows` (or `prefetch_rows`, then `decorate_row`).

Both need a `Render`. Another command gets one from `open_render(analysis,
limit)`, a context manager that opens the annotation store and closes it when
the block ends. Do not build a `Render` by hand or call the private
`_annotation_store` from another module: a store left open is a leaked handle.

JSON compatibility policy (docs/TRIO.md): additive only. A key that is
already emitted is never removed, renamed or retyped, and the relative
order of the existing keys never changes; new keys are snake_case and are
either appended or placed directly after the key they qualify
(`blocked_reason` follows `blocked`). tests/test_trio_goldens.py pins this
byte-for-byte once the whitelisted new keys are removed.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
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
    """What building one document shares across its models and rows."""

    # What the document is about. A `TrioAnalysis` for `db trio` (the model
    # builders and candidate decorators rely on it); another command that
    # reuses only the row decorators passes its own object, which nothing in
    # the row chain reads.
    analysis: object
    limit: int
    have: dict  # annotation sources that hold data
    store: object | None = None  # open annotation connection, kept alive only
    # Answers prefetched in batches, so a rendered row costs a dict lookup,
    # not a store query. The annotation columns keep theirs under the bare
    # variant key (see prefetch_annotations); any other decorator keys its
    # entries under its own name, e.g. `render.cache[("parent_alt", key)]`,
    # so that two decorators never share an entry.
    cache: dict = field(default_factory=dict)


@contextmanager
def open_render(analysis: object, limit: int) -> Iterator[Render]:
    """A `Render` over the annotation store, closed when the block ends.

    The way for a command other than `db trio` to get one: hold it open for
    the whole rendering (the store stays open between lookups, see
    `_annotation_store`) and pass it to `decorate_rows`. `analysis` is
    whatever the command's document is about; the row chain never reads it.

        with open_render(my_analysis, limit) as render:
            rows = decorate_rows(outs, variant_rows, render)
    """
    have, store = _annotation_store()
    try:
        yield Render(analysis, limit, have, store)
    finally:
        if store is not None:
            store.close()


def prefetch_annotations(render: Render, rows: Sequence[Sequence]) -> None:
    """Look up genes, gnomAD and ClinVar for every row about to be rendered
    in one batch per source. A per-candidate store query cost 0.1-1.5 ms
    each (several seconds for a 1,000-candidate model on a slow disk); the
    batch is one join. Rows already cached are skipped.

    `rows` need only carry (chrom, pos, ref, alt) as their first four
    columns, so a trio detail row and any other command's variant row both
    work. This is the annotation row decorator's registered prefetch, so the
    batch runs by itself through `prefetch_rows` / `decorate_rows` / `db
    trio`; the decorator reads what is stored in `render.cache`, and without
    a prefetch it falls back to one lookup per row."""
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


_prefetch = prefetch_annotations  # the name before it was public


# One chain, in registration order, which is the order keys appear in the
# output. An entry is either a row decorator (needs only the variant key,
# runs for every kind of row) or a candidate decorator (needs the trio row
# and the model, runs for `db trio` candidates only).
VariantKey = tuple[str, int, str, str]  # (chrom, pos, ref, alt)

# A row decorator receives the dict built so far, the row's variant key and
# the render context, and returns the dict to emit (usually the same one with
# keys appended).
RowDecorator = Callable[[dict, VariantKey, Render], dict]

# A candidate decorator receives the candidate dict built so far, the source
# row and the model name, and returns the dict to emit.
CandidateDecorator = Callable[[dict, list, str, Render], dict]

# A document decorator receives the whole trio document and returns it.
DocumentDecorator = Callable[[dict, Render], dict]

# A prefetch receives the render context and the batch of rows about to be
# decorated, looks up whatever the per-row call will need for all of them at
# once, and stores it in `render.cache`. For a row decorator the rows are the
# batch's distinct variant keys; for a candidate decorator they are the trio
# detail rows of one model. Called once per batch, never for an empty one.
Prefetch = Callable[[Render, Sequence[Sequence]], None]


@dataclass(frozen=True)
class _Decorator:
    fn: RowDecorator | CandidateDecorator
    row_only: bool
    models: frozenset[str] | None = None  # candidate decorators only
    prefetch: Prefetch | None = None


_DECORATORS: list[_Decorator] = []
_DOCUMENT_DECORATORS: list[DocumentDecorator] = []


def register_row_decorator(
    fn: RowDecorator, *, prefetch: Prefetch | None = None
) -> RowDecorator:
    """Append `fn` to the chain as a per-variant column: it sees only the
    variant key (chrom, pos, ref, alt), so it applies to a trio candidate and
    to a row of any other command alike (`decorate_row`). Order of
    registration is the order keys appear in the output.

    `prefetch(render, keys)`, if given, is called once per batch with the
    batch's distinct variant keys, before any row of it is decorated."""
    _DECORATORS.append(_Decorator(fn, row_only=True, prefetch=prefetch))
    return fn


def register_candidate_decorator(
    fn: CandidateDecorator,
    models: frozenset[str] | None = None,
    *,
    prefetch: Prefetch | None = None,
) -> CandidateDecorator:
    """Append `fn` to the chain as a `db trio` candidate key (all models, or
    only `models`). It sees the 8-column trio row, so it is never applied to
    another command's rows. Order of registration is the order keys appear
    in the output.

    `prefetch(render, rows)`, if given, is called once per rendered batch of
    one model's trio rows (only for a model `fn` applies to), before any row
    of it is decorated."""
    _DECORATORS.append(_Decorator(fn, row_only=False, models=models, prefetch=prefetch))
    return fn


def register_document_decorator(fn: DocumentDecorator) -> DocumentDecorator:
    """Append `fn` to the document chain: after every model is built, each
    decorator receives the whole trio document and returns it. It sees the
    nested objects too (`trio`, `gates`, `models`), so it may add a key at
    the top or inside one of them. Run in registration order."""
    _DOCUMENT_DECORATORS.append(fn)
    return fn


def _variant_key(row: Sequence) -> VariantKey:
    return (row[0], int(row[1]), row[2], row[3])


def _decorate(out: dict, row: Sequence, render: Render, model: str | None) -> dict:
    """Run the chain over `out`. `model` is None for a row that is not a trio
    candidate: only the row decorators apply then."""
    key = _variant_key(row)
    for d in _DECORATORS:
        if d.row_only:
            out = d.fn(out, key, render)
        elif model is not None and (d.models is None or model in d.models):
            out = d.fn(out, row, model, render)
    return out


def _prefetch_chain(render: Render, rows: Sequence[Sequence], model: str | None):
    """Run the chain's prefetches once over `rows`, in registration order:
    every row decorator's (given the distinct variant keys) and, when `model`
    is a trio model, those of the candidate decorators that apply to it."""
    if not rows:
        return
    keys = list(dict.fromkeys(_variant_key(r) for r in rows))
    for d in _DECORATORS:
        if d.prefetch is None:
            continue
        if d.row_only:
            d.prefetch(render, keys)
        elif model is not None and (d.models is None or model in d.models):
            d.prefetch(render, rows)


def prefetch_rows(render: Render, rows: Sequence[Sequence]) -> None:
    """Run every row decorator's prefetch once for `rows` (anything that
    starts with (chrom, pos, ref, alt)), so the `decorate_row` calls that
    follow read answers from `render.cache` and make no lookup. Use it when
    the dicts are built row by row; otherwise `decorate_rows` does both."""
    _prefetch_chain(render, rows, None)


def decorate_row(out: dict, row: Sequence, render: Render) -> dict:
    """Append the per-variant columns (the row decorators: genes, gnomAD
    popmax, ClinVar, ...) to `out`, a dict the caller built for `row`.

    `row` need only start with (chrom, pos, ref, alt). Candidate decorators
    are not applied: they are `db trio` keys that read a trio detail row.
    Call `prefetch_rows` on the batch first, or each row costs the lookups
    its decorators make."""
    return _decorate(out, row, render, None)


def decorate_rows(
    outs: Sequence[dict], rows: Sequence[Sequence], render: Render
) -> list[dict]:
    """`decorate_row` for a whole batch: run every row decorator's prefetch
    once over `rows`, then decorate `outs[i]` for `rows[i]`. The entry point
    for a command other than `db trio` that renders many variant rows."""
    if len(outs) != len(rows):
        raise ValueError(f"{len(outs)} dicts for {len(rows)} rows")
    prefetch_rows(render, rows)
    return [decorate_row(o, r, render) for o, r in zip(outs, rows, strict=True)]


def _candidates(rows: list, model: str, render: Render) -> list[dict]:
    """Rendered candidates for `rows`: the prefetches of every decorator that
    applies to `model` run once over the batch first."""
    _prefetch_chain(render, rows, model)
    return [_candidate(r, model, render) for r in rows]


def _candidate_base(row: list) -> dict:
    """The keys of one 8-column trio detail row, before any decorator."""
    chrom, pos, ref, alt, pgt, fgt, mgt, af = row
    return {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "proband_gt": int(pgt),
        "father_gt": int(fgt),
        "mother_gt": int(mgt),
        "af": _num(af),
    }


def _candidate(row: list, model: str, render: Render) -> dict:
    """One 8-column detail row as the JSON candidate object."""
    return _decorate(_candidate_base(row), row, render, model)


def _annotation_columns(out: dict, key: VariantKey, render: Render) -> dict:
    hit = render.cache.get(key)
    out.update(dict(hit) if hit is not None else _enrich(*key, render.have))
    return out


register_row_decorator(_annotation_columns, prefetch=prefetch_annotations)


@dataclass(frozen=True)
class ModelSpec:
    """One inheritance model's place in a trio document and in text output."""

    name: str
    # Does the model exist for this trio? (Later models are conditional.)
    applies: Callable[[object], bool]
    # The JSON object for `models.<name>`.
    build: Callable[[Render], dict]
    # One line for `--category all` (None: no line for this trio), and the
    # listing for `--category <name>`.
    summary: Callable[[TrioAnalysis], str | None]
    detail: Callable[[TrioAnalysis], None]
    # Needs stored hom-reference calls to prove a parent is 0/0 (the
    # `--keep-reference` note, and the `blocked` flag, key off this).
    needs_reference: bool = False
    # Needs the gene table (compound het groups candidates by gene).
    needs_genes: bool = False
    # The SQL the model runs under a scope, for (trio, gates, scope): every
    # statement must carry the scope predicate (cli/scope.py). Registered in
    # cli.scope.SCOPE_AWARE_BUILDERS with the model, and walked by
    # tests/test_scope_models_drift.py.
    scope_builder: Callable | None = None


MODELS: dict[str, ModelSpec] = {}


def register_model(spec: ModelSpec) -> ModelSpec:
    """Add a model after the existing ones: registration order is the order
    of `models` in the document and of the lines in the `all` summary. Its
    scope builder goes into cli.scope.SCOPE_AWARE_BUILDERS; a model without
    one fails the scope drift test."""
    if spec.name in MODELS:
        raise ValueError(f"model {spec.name!r} is already registered")
    MODELS[spec.name] = spec
    if spec.scope_builder is not None:
        from cli.scope import register_scope_builder

        register_scope_builder(spec.name, spec.scope_builder)
    return spec


def _legacy_blocked(needs_reference: bool, analysis: TrioAnalysis) -> bool:
    """The `blocked` flag of the four classic models: no stored hom-ref rows
    anywhere in the database (cli/trio_analysis._has_reference_rows)."""
    return needs_reference and not analysis.has_reference


def model_blockers(
    analysis: TrioAnalysis, spec: ModelSpec, genes_loaded: bool
) -> list[str]:
    """Every `blocked_reason` that stops `spec` from running on this database,
    the reference requirement first (the JSON reports only that one when both
    apply). Empty when the model can run."""
    out = []
    if _legacy_blocked(spec.needs_reference, analysis):
        out.append("needs_keep_reference")
    if spec.needs_genes and not genes_loaded:
        out.append("needs_genes")
    return out


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
            "x_guard_excluded": analysis.x_guard_excluded(name),
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
        # A later feature may append to each listed row (the parents' read
        # counts under --parent-alt-screen): stats["text_suffix"](model, rows)
        # returns one string per row.
        suffix = analysis.stats.get("text_suffix")
        extra = suffix(name, data) if suffix else None
        for i, row in enumerate(data):
            chrom, pos, ref, alt, pgt, fgt, mgt, af = row
            af_s = "NA" if af is None else f"{af}"
            click.echo(
                f"  {chrom}:{pos} {ref}>{alt}  proband_gt={pgt} "
                f"father_gt={fgt} mother_gt={mgt}  AF={af_s}"
                + (extra[i] if extra else "")
            )

    def scope_builder(trio, gates, scope) -> list[tuple[str, str, str]]:
        from cli.trio_analysis import trio_sql

        return [
            (trio_sql(name, trio, gates, count_only=count, scope=scope), "all", "g")
            for count in (False, True)
        ]

    return ModelSpec(
        name,
        lambda trio: True,
        build,
        summary,
        detail,
        needs_reference=needs_reference,
        scope_builder=scope_builder,
    )


# ───────────────────────────── compound het ─────────────────────────────


def _comphet_entry(sym: str, sides: dict, scope=None) -> dict:
    """One gene of the comphet model. `paternal` and `maternal` are emitted
    for the sides present; a consumer must treat each as optional (a model
    that cannot tell parent of origin omits them). With a scope that includes
    genes or regions, each variant carries `in_scope` (a partner outside a
    capture edge is flagged, not dropped)."""
    entry: dict = {"gene": sym}
    flag = scope is not None and scope.flags_variants
    for origin in ("paternal", "maternal"):
        if sides.get(origin) is None:
            continue
        entry[origin] = [
            {"chrom": c, "pos": int(p), "ref": r, "alt": a, "af": _num(af)}
            | ({"in_scope": scope.contains(c, int(p))} if flag else {})
            for c, p, r, a, af in sides[origin]
        ]
    return entry


def outside_note(analysis, chrom: str, pos) -> str:
    """`  (outside the scope)` for a compound-het partner the scope flags."""
    scope = getattr(analysis, "scope", None)
    if scope is None or not scope.flags_variants or scope.contains(chrom, int(pos)):
        return ""
    return "  (outside the scope)"


def _build_comphet(render: Render) -> dict:
    analysis, limit = render.analysis, render.limit
    genes = analysis.comphet_genes()
    entries = [
        _comphet_entry(sym, genes[sym], analysis.scope) for sym in sorted(genes)[:limit]
    ]
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
        "x_guard_excluded": analysis.x_guard_excluded("comphet"),
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
                click.echo(
                    f"    {origin:8s} {chrom}:{pos} {ref}>{alt}  AF={af_s}"
                    + outside_note(analysis, chrom, pos)
                )


def _comphet_scope_sql(trio, gates, scope) -> list[tuple[str, str, str]]:
    from cli.trio_analysis import _comphet_sql

    return [(_comphet_sql(trio, gates, scope=scope), "exclude", "g")]


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
        needs_genes=True,
        scope_builder=_comphet_scope_sql,
    )
)

# The X-linked model: cli/trio_xlinked_output.py registers it here, after the
# four classic models, so it is always the last key of `models`.
from cli.trio_xlinked_output import XLINKED_MODEL  # noqa: E402

register_model(XLINKED_MODEL)

# The registry and the shared CATEGORIES constant describe the same models.
if tuple(MODELS) != CATEGORIES:
    raise RuntimeError(f"model registry {tuple(MODELS)} != CATEGORIES {CATEGORIES}")


# ───────────────────────────── the document ─────────────────────────────


def selected_categories(category: str | Sequence[str]) -> tuple[str, ...]:
    """The models `--category` names, as a tuple of names. `all` anywhere in a
    list means all; a single name (or `all`) behaves exactly as it always has."""
    tokens = (category,) if isinstance(category, str) else tuple(category)
    return ("all",) if "all" in tokens else tokens


def _applicable(
    analysis: TrioAnalysis, category: str | Sequence[str]
) -> list[ModelSpec]:
    wanted = selected_categories(category)
    return [
        spec
        for spec in MODELS.values()
        if ("all" in wanted or spec.name in wanted) and spec.applies(analysis.trio)
    ]


def omitted_models(analysis: TrioAnalysis, category: str | Sequence[str]) -> list[str]:
    """Models a list or `all` names that do not apply to this trio: they are
    left out of the report, and the caller says so. (A single explicit name
    that does not apply is not omitted quietly; that is each model's own rule.)"""
    wanted = selected_categories(category)
    if wanted != ("all",) and len(wanted) == 1:
        return []
    return [
        spec.name
        for spec in MODELS.values()
        if ("all" in wanted or spec.name in wanted) and not spec.applies(analysis.trio)
    ]


def _gates_json(analysis: TrioAnalysis) -> dict:
    gates = analysis.gates
    out = {
        "min_gq": gates.min_gq,
        "min_dp": gates.min_dp,
        "max_af": gates.max_af,
        "min_ab": gates.min_ab,
        "max_ab": gates.max_ab,
        "gnomad_max_af": analysis.gnomad_max_af,
        "af_gate": gates.af_gate,
    }
    # the X model's gates, as the values it applies (cli/xlinked.py)
    from cli.xlinked import male_x_dp, x_max_af

    out["min_male_x_dp"] = male_x_dp(gates)
    out["max_af_x"] = x_max_af(gates)
    out["min_hemi_ab"] = gates.min_hemi_ab
    out["max_parent_alt_reads"] = gates.max_parent_alt_reads
    # additive: a gate is reported only when it is on
    if gates.strict_missing:
        out["strict_missing"] = True
        if gates.strict_exempt:
            out["strict_exempt"] = list(gates.strict_exempt)
    if gates.pass_only:
        out["pass_only"] = True
    if gates.parent_alt_screen:
        out["parent_alt_screen"] = True
        out["max_parent_alt_frac"] = gates.max_parent_alt_frac
    return out


def _trio_json(
    analysis: TrioAnalysis, category: str | Sequence[str], limit: int
) -> dict:
    trio = analysis.trio
    with open_render(analysis, limit) as render:
        models = {
            spec.name: spec.build(render) for spec in _applicable(analysis, category)
        }
        doc = {
            "trio": {
                "ingest_id": trio.ingest_id,
                "proband": trio.proband,
                "father": trio.father,
                "mother": trio.mother,
                "family_id": trio.family_id,
                "kind": trio.kind,
                "siblings": [{"sample": s, "kind": k} for s, k in analysis.siblings],
            },
            "keep_reference": analysis.has_reference,
            "gates": _gates_json(analysis),
            "annotations": render.have,
            "models": models,
            # the X guard on the classic models (cli/xlinked.x_guard_sql)
            "x_guard": analysis.x_guard,
        }
        # Inside the block: the store stays open for a decorator that looks
        # something up.
        for fn in _DOCUMENT_DECORATORS:
            doc = fn(doc, render)
            if not isinstance(doc, dict):
                raise TypeError(
                    f"document decorator {getattr(fn, '__name__', fn)!r} must "
                    f"return the document dict, not {type(doc).__name__}"
                )
    return doc


def _scope_document(doc: dict, render: Render) -> dict:
    """The top-level `scope` object, for a run with genes or regions."""
    scope = getattr(render.analysis, "scope", None)
    if scope is not None:
        from cli.scope import scope_document

        doc["scope"] = scope_document(scope, render.analysis.stats)
    return doc


register_document_decorator(_scope_document)


def trio_text(analysis: TrioAnalysis, category: str | Sequence[str]) -> None:
    trio = analysis.trio
    proband, father, mother = trio.proband, trio.father, trio.mother
    click.echo(f"trio: proband={proband} father={father} mother={mother}")
    # lines a run with a preset adds right under `trio:` (cli/trio_gate_output.py)
    for line in analysis.stats.get("text_header", ()):
        click.echo(line)
    wanted = selected_categories(category)
    needs_reference = {s.name for s in MODELS.values() if s.needs_reference}
    if not analysis.has_reference and (
        "all" in wanted or any(c in needs_reference for c in wanted)
    ):
        click.echo(
            "  note: this database has no stored hom-reference calls, so "
            "de-novo/dominant/comphet cannot prove a parent is 0/0 (vs "
            "no-call). Re-ingest with `--keep-reference` for those models.",
            err=True,
        )

    if len(wanted) == 1 and wanted[0] != "all":
        MODELS[wanted[0]].detail(analysis)  # one model: the listing
    else:
        for spec in _applicable(analysis, wanted):  # `all` or a list: a line each
            line = spec.summary(analysis)
            if line is not None:
                click.echo(line)
    # what only the finished run knows (cli/trio_gate_output.py)
    footer = analysis.stats.get("text_footer")
    for line in footer() if footer else ():
        click.echo(line)
