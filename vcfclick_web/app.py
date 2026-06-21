"""FastAPI app behind `vcfclick web`.

A thin HTTP layer over existing vcfclick internals — it adds no new query
logic. Endpoints reuse storage.get_session (the same path as `db query`
and the MCP run_sql tool), cli.db_trio's SQL builders, and
ingest.combine. Served on localhost only; arbitrary SELECTs are allowed
because the caller already owns the database (same trust level as the
CLI), but write statements are rejected as a guardrail.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from pydantic import BaseModel

from ingest._arrow import (
    GENOTYPES_ARROW_SCHEMA,
    INGESTIONS_ARROW_SCHEMA,
    PEDIGREE_ARROW_SCHEMA,
    SAMPLES_ARROW_SCHEMA,
    VARIANTS_ARROW_SCHEMA,
    column_names,
)
from ingest.combine import CombineError, combine_vcfs
from storage import get_session, table_exists
from vcfclick_web.page import INDEX_HTML

_TABLES = [
    ("variants", VARIANTS_ARROW_SCHEMA),
    ("genotypes", GENOTYPES_ARROW_SCHEMA),
    ("samples", SAMPLES_ARROW_SCHEMA),
    ("ingestions", INGESTIONS_ARROW_SCHEMA),
    ("pedigree", PEDIGREE_ARROW_SCHEMA),
]

# Read-only guardrail. The server is localhost-only and the caller owns
# the database, but a malicious page (CSRF / DNS-rebinding) can POST to a
# local server, so writes must be genuinely blocked — not just for queries
# that *start* with a write verb. We allowlist read statements, forbid
# multiple statements, and reject any write construct anywhere (so a
# CTE-prefixed DML or a SELECT ... INTO OUTFILE / COPY ... TO cannot slip
# through). Comments are stripped first so they can neither hide a verb nor
# false-trip the blocklist.
_ALLOWED_START = {"select", "with", "show", "describe", "desc", "explain"}
_COMMENT_RE = re.compile(r"/\*.*?\*/|--[^\n]*", re.DOTALL)
# Verbs forbidden ANYWHERE (catches CTE-prefixed DML like
# `WITH c AS (...) DELETE ...` and file writes). Deliberately excludes
# truncate/replace/merge — those are also SQL *functions*; as statements
# they only appear leading, where the allowlist already rejects them.
_BLOCK_RE = re.compile(
    r"\b(insert|update|delete|drop|alter|create|attach|detach|copy|grant|"
    r"revoke|rename|into\s+outfile|into\s+dumpfile)\b",
    re.IGNORECASE,
)


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


def _is_read_only(sql: str) -> bool:
    """True only for a single read statement with no write construct."""
    bare = _COMMENT_RE.sub(" ", sql).strip()
    body = bare.rstrip(";").strip()
    if not body or ";" in body:  # empty, or more than one statement
        return False
    first = re.match(r"[a-zA-Z]+", body)
    if not first or first.group(0).lower() not in _ALLOWED_START:
        return False
    return _BLOCK_RE.search(body) is None


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
    tables = [
        {"name": name, "columns": column_names(schema)}
        for name, schema in _TABLES
        if table_exists(name)
    ]
    ingest_ids = _scalar_list(
        "SELECT DISTINCT ingest_id FROM ingestions FORMAT JSONCompact"
    )
    samples = _scalar_list(
        "SELECT DISTINCT sample_id FROM samples ORDER BY sample_id LIMIT 1000 FORMAT JSONCompact"
    )
    probands = _scalar_list(
        "SELECT DISTINCT sample_id FROM pedigree "
        "WHERE father_id NOT IN ('0', '') AND mother_id NOT IN ('0', '') "
        "ORDER BY sample_id FORMAT JSONCompact"
    )
    return {
        "db": _db_name(),
        "tables": tables,
        "ingest_ids": ingest_ids,
        "samples": samples,
        "probands": probands,
    }


@app.post("/api/query")
def query(body: QueryBody) -> dict:
    sql = body.sql.strip()
    if not sql:
        return {"error": "empty query"}
    if not _is_read_only(sql):
        return {
            "error": "the web UI is read-only — only a single "
            "SELECT/WITH/SHOW/DESCRIBE/EXPLAIN statement is allowed"
        }
    try:
        return _run_sql(sql)
    except Exception as _e:
        return {"error": "query execution failed"}


@app.post("/api/nl")
def nl(body: NlBody) -> dict:
    if not body.question.strip():
        return {"error": "empty question"}
    try:
        from vcfclick_mcp.server import SCHEMA_DESCRIPTION
        from vcfclick_web.llm import LLMError, generate_sql

        try:
            sql = generate_sql(
                body.provider, body.key, body.model, body.question, SCHEMA_DESCRIPTION
            )
        except LLMError as e:
            return {"error": str(e)}
        if not _is_read_only(sql):
            return {
                "sql": sql,
                "error": "the model produced a non-read-only statement; not running it",
            }
        result = _run_sql(sql)
        return result
    except Exception as e:
        return {"error": str(e)}


@app.get("/api/trio")
def trio(
    proband: str,
    category: str = "denovo",
    min_gq: int = 20,
    min_dp: int = 10,
    max_af: float = 0.01,
    min_ab: float = 0.25,
    max_ab: float = 0.75,
) -> dict:
    if category not in ("denovo", "recessive", "dominant"):
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
        gates = Gates(min_gq, min_dp, max_af, min_ab, max_ab)
        sql = trio_sql(category, trio_ctx, gates, count_only=False)
        result = _run_sql(sql)
        result["father"] = father
        result["mother"] = mother
        if category in ("denovo", "dominant") and not _has_reference_rows(sess):
            result["note"] = (
                "No stored hom-reference rows: re-ingest with --keep-reference "
                "for defensible de-novo / dominant results."
            )
        return result
    except Exception as e:
        return {"error": str(e)}


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
    except CombineError as e:
        return {"error": str(e)}
    except Exception as e:
        return {"error": str(e)}


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
