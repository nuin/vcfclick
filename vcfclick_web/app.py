"""FastAPI app behind `vcfclick web`.

A thin HTTP layer over existing vcfclick internals — it adds no new query
logic. Endpoints reuse storage.get_session (the same path as `db query`
and the MCP run_sql tool), cli.db_trio's SQL builders, and
ingest.combine. Served on localhost only; arbitrary SELECTs are allowed
because the caller already owns the database (same trust level as the
CLI), but write statements are rejected as a guardrail
(storage.sql_guard: an AST check, since a malicious page can POST to a
local server via CSRF / DNS-rebinding).
"""

from __future__ import annotations

import json
import logging
import os
import tempfile
from pathlib import Path

import click
from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from pydantic import BaseModel

from ingest._arrow import (
    GENOTYPES_ARROW_SCHEMA,
    INGESTIONS_ARROW_SCHEMA,
    MISSING_GENOTYPES_ARROW_SCHEMA,
    PEDIGREE_ARROW_SCHEMA,
    POPULATIONS_ARROW_SCHEMA,
    SAMPLES_ARROW_SCHEMA,
    VARIANTS_ARROW_SCHEMA,
    column_names,
)
from ingest.combine import CombineError, combine_vcfs
from storage import get_session, table_columns, table_exists
from cli.trio_analysis import CATEGORIES, DEFAULT_MAX_AF
from storage.sql_guard import is_read_only
from vcfclick_web.page import INDEX_HTML

# Raw exception text is never sent to the browser (it can carry paths and
# internals); the detail goes to the terminal running `vcfclick web`.
_log = logging.getLogger("vcfclick.web")

_TABLES = [
    ("variants", VARIANTS_ARROW_SCHEMA),
    ("genotypes", GENOTYPES_ARROW_SCHEMA),
    ("samples", SAMPLES_ARROW_SCHEMA),
    ("ingestions", INGESTIONS_ARROW_SCHEMA),
    ("pedigree", PEDIGREE_ARROW_SCHEMA),
    ("missing_genotypes", MISSING_GENOTYPES_ARROW_SCHEMA),
    ("populations", POPULATIONS_ARROW_SCHEMA),
]


def _db_name() -> str | None:
    return os.environ.get("VCFCLICK_DB_NAME")


def _run_sql(sql: str) -> dict:
    """Execute SQL via the shared session path; shape it like run_sql."""
    sess = get_session()
    raw = sess.query(sql, "JSONCompact").bytes().decode()
    parsed = json.loads(raw)
    return {
        "sql": sql,
        "columns": [m["name"] for m in parsed.get("meta", [])],
        "rows": parsed.get("data", []),
        "row_count": len(parsed.get("data", [])),
    }


def _scalar_list(sql: str) -> list[str]:
    """Run a one-column query and return the values as strings, or []."""
    try:
        res = _run_sql(sql)
    except Exception:
        return []
    return [str(r[0]) for r in res["rows"] if r and r[0] is not None]


class QueryBody(BaseModel):
    sql: str


class NlBody(BaseModel):
    question: str
    provider: str = "gemini"
    key: str = ""
    model: str = ""


class CombineBody(BaseModel):
    first: str
    second: str
    name1: str = "first"
    name2: str = "second"
    min_callsets: int = 1


app = FastAPI(title="vcfclick web", docs_url=None, redoc_url=None)


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return INDEX_HTML


@app.get("/api/meta")
def meta() -> dict:
    # The live columns, in schema order: a database ingested before a
    # column was added has not been upgraded by the read-only web UI.
    tables = []
    for name, schema in _TABLES:
        if table_exists(name):
            live = set(table_columns(name))
            cols = [c for c in column_names(schema) if c in live]
            tables.append({"name": name, "columns": cols})
    ingest_ids = _scalar_list(
        "SELECT DISTINCT ingest_id FROM ingestions FORMAT JSONCompact"
    )
    samples = _scalar_list(
        "SELECT DISTINCT sample_id FROM samples ORDER BY sample_id LIMIT 1000 FORMAT JSONCompact"
    )
    probands, families = _pedigree_meta(ingest_ids)
    return {
        "db": _db_name(),
        "tables": tables,
        "ingest_ids": ingest_ids,
        "samples": samples,
        "probands": probands,
        "families": families,
    }


_noted: set[str] = set()


def _note_once(key: str, message: str, *args) -> None:
    """Log `message` the first time `key` comes up: /api/meta runs on every
    page load, and the same note each time would bury the log."""
    if key not in _noted:
        _noted.add(key)
        _log.info(message, *args)


def _pedigree_meta(ingest_ids: list[str]) -> tuple[list[str], list[dict]]:
    """(probands, families) for the proband picker.

    `probands` is a flat list of sample ids by the same rule as `db trio`
    and `db families` (cli.pedigree_load builds the graph; see
    storage/pedigree_graph.py): an analysable child not marked unaffected
    who has no children or is affected. `families` is the structure behind
    it: per family its probands and every trio or duo with its siblings.

    Pedigrees are per ingestion. With exactly one, both come from its graph;
    with several (trio analysis does not apply to those databases) the older
    flat list of every child with two declared parents is all there is. A
    database without a pedigree, or one the graph cannot read, gives none.
    """
    from cli.pedigree_load import NO_PEDIGREE

    if not table_exists("pedigree"):
        _note_once(
            "no-pedigree", "pedigree: %s The proband list is empty.", NO_PEDIGREE
        )
        return [], []
    if len(ingest_ids) == 1:
        try:
            from cli.pedigree_load import load_pedigree_graph
            from cli.pedigree_report import web_families

            graph = load_pedigree_graph(get_session(), ingest_ids[0])
            for child in graph.analysable_children():
                if graph.kind_of(child) == "duo":
                    _note_once(
                        f"duo:{child}",
                        "pedigree: %s has one declared parent (a duo); duo "
                        "analysis is not available yet, so it is not offered "
                        "as a proband",
                        child,
                    )
            return graph.probands(), web_families(graph)
        except Exception:
            _log.exception("could not read the pedigree; using the flat list")
    return (
        _scalar_list(
            "SELECT DISTINCT sample_id FROM pedigree "
            "WHERE father_id NOT IN ('0', '') AND mother_id NOT IN ('0', '') "
            "ORDER BY sample_id FORMAT JSONCompact"
        ),
        [],
    )


@app.post("/api/query")
def query(body: QueryBody) -> dict:
    sql = body.sql.strip()
    if not sql:
        return {"error": "empty query"}
    if not is_read_only(sql):
        return {
            "error": "the web UI is read-only — only a single "
            "SELECT/WITH/SHOW/DESCRIBE/EXPLAIN statement is allowed"
        }
    try:
        return _run_sql(sql)
    except Exception:
        _log.exception("query failed")
        return {"error": "query execution failed (details in the vcfclick web log)"}


@app.post("/api/nl")
def nl(body: NlBody) -> dict:
    if not body.question.strip():
        return {"error": "empty question"}
    sql = None
    try:
        from vcfclick_mcp.server import SCHEMA_DESCRIPTION
        from vcfclick_web.llm import LLMError, generate_sql

        try:
            sql = generate_sql(
                body.provider, body.key, body.model, body.question, SCHEMA_DESCRIPTION
            )
        except LLMError as e:
            return {"error": str(e)}
        if not is_read_only(sql):
            return {
                "sql": sql,
                "error": "the model produced a non-read-only statement; not running it",
            }
        result = _run_sql(sql)
        return result
    except Exception:
        _log.exception("natural-language query failed")
        out = {"error": "the generated query failed (details in the vcfclick web log)"}
        if sql:
            out["sql"] = sql
        return out


_KEEP_REFERENCE_NOTE = (
    "No stored hom-reference rows: re-ingest with --keep-reference "
    "for defensible de-novo / dominant results."
)


def _comphet_result(
    sess, trio_ctx, gates, has_reference: bool
) -> tuple[dict, list[str]]:
    """The compound-het model in the same envelope as the row models: one
    table row per candidate variant, grouped under its gene. The model is
    reported per gene (cli.trio_output), so the SQL that ran is the
    candidate query and the rows come from the gene grouping. Returns the
    result and any notes for the user."""
    from cli.trio_analysis import TrioAnalysis, _comphet_sql
    from cli.trio_output import _annotation_status

    analysis = TrioAnalysis(
        session=sess,
        trio=trio_ctx,
        gates=gates,
        has_reference=has_reference,
        gnomad_max_af=None,
    )
    genes = analysis.comphet_genes()
    rows = [
        [sym, origin, chrom, int(pos), ref, alt, af]
        for sym in sorted(genes)
        for origin in ("paternal", "maternal")
        for chrom, pos, ref, alt, af in genes[sym][origin]
    ]
    result = {
        "sql": _comphet_sql(trio_ctx, gates),
        "columns": ["gene", "origin", "chrom", "pos", "ref", "alt", "af"],
        "rows": rows,
        "row_count": len(rows),
        "genes": len(genes),
    }
    notes = []
    if not _annotation_status()["genes"]:
        notes.append(
            "No gene annotations are loaded, so compound hets cannot be "
            "grouped by gene: run `vcfclick annotations load`."
        )
    return result, notes


@app.get("/api/trio")
def trio(
    proband: str,
    category: str = "denovo",
    min_gq: int = 20,
    min_dp: int = 10,
    max_af: float | None = None,
    min_ab: float = 0.25,
    max_ab: float = 0.75,
) -> dict:
    # `max_af` omitted means the default 0.01 *and* that the cohort-derived
    # AF auto-skip may apply (as for `db trio` without --max-af); a value
    # is always honoured.
    if category not in CATEGORIES:
        return {"error": f"unknown category {category!r}"}
    try:
        from cli.db_trio import (
            Gates,
            Trio,
            _has_reference_rows,
            _resolve_parents,
            _sole_ingest_id,
            trio_sql,
        )
        from cli.family_checks import trio_af_gate
        from cli.trio_output import MODELS

        name = _db_name()
        ingest_id = _sole_ingest_id(name)
        if not ingest_id:
            return {
                "error": "this database has zero or multiple ingestions; "
                "trio analysis needs a single-ingestion cohort"
            }
        sess = get_session()
        father, mother = _resolve_parents(sess, ingest_id, proband)
        trio_ctx = Trio(ingest_id, proband, father, mother)
        threshold = DEFAULT_MAX_AF if max_af is None else max_af
        af_gate, af_message = trio_af_gate(
            sess, ingest_id, threshold, explicit=max_af is not None
        )
        gates = Gates(min_gq, min_dp, threshold, min_ab, max_ab)._replace(
            af_gate=af_gate
        )
        needs_reference = MODELS[category].needs_reference
        has_reference = _has_reference_rows(sess) if needs_reference else True
        notes: list[str] = []
        if needs_reference and not has_reference:
            notes.append(_KEEP_REFERENCE_NOTE)
        if category == "comphet":
            result, extra = _comphet_result(sess, trio_ctx, gates, has_reference)
            notes += extra
        else:
            result = _run_sql(trio_sql(category, trio_ctx, gates, count_only=False))
        result["father"] = father
        result["mother"] = mother
        result["af_gate"] = af_gate
        if af_message:
            notes.append(af_message)
        if notes:
            result["note"] = " ".join(notes)
        return result
    except click.ClickException as e:
        # Written for the user (e.g. "no pedigree entry for proband ...").
        return {"error": e.message}
    except Exception:
        _log.exception("trio analysis failed")
        return {"error": "trio analysis failed (details in the vcfclick web log)"}


@app.post("/api/combine")
def combine(body: CombineBody) -> dict:
    try:
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "a.vcf").write_text(body.first)
            (d / "b.vcf").write_text(body.second)
            out = d / "out.vcf"
            combine_vcfs(
                [d / "a.vcf", d / "b.vcf"],
                out,
                names=[body.name1, body.name2],
                min_callsets=body.min_callsets,
            )
            return _parse_combined(out.read_text())
    except CombineError:
        return {"error": "combine failed due to invalid or incompatible VCF input"}
    except Exception:
        return {"error": "internal error while combining VCFs"}


def _parse_combined(text: str) -> dict:
    """Reduce a combined VCF to {samples, fields, records} for the UI."""
    samples: list[str] = []
    fields: list[str] = ["GT"]
    records = []
    for line in text.splitlines():
        if line.startswith("##") or not line.strip():
            continue
        cols = line.split("\t")
        if line.startswith("#CHROM"):
            samples = cols[9:]
            continue
        if len(cols) < 10:
            continue
        info = cols[7]
        set_val = next(
            (kv.split("=", 1)[1] for kv in info.split(";") if kv.startswith("set=")),
            ".",
        )
        fields = cols[8].split(":")
        cells = {s: cols[9 + i] for i, s in enumerate(samples) if 9 + i < len(cols)}
        records.append(
            {
                "chrom": cols[0],
                "pos": int(cols[1]),
                "ref": cols[3],
                "alt": cols[4],
                "set": set_val,
                "cells": cells,
            }
        )
    return {"samples": samples, "fields": fields, "records": records}
