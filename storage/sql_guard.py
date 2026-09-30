"""AST-based read-only SQL guard shared by `vcfclick web` and MCP `run_sql`.

Both surfaces run caller-supplied SQL (typed by a user or produced by an
LLM) against the cohort database. The web server is localhost-only, but a
malicious page (CSRF / DNS-rebinding) can POST to it, and an LLM can be
prompt-injected, so writes must be genuinely blocked.

Rather than pattern-matching keywords, the statement is parsed with
sqlglot in the active backend's dialect (ClickHouse for chDB, DuckDB for
DuckDB) and the syntax tree is inspected:

  * exactly one statement;
  * the root must be a query (SELECT / WITH ... SELECT / UNION ...),
    DESCRIBE, SHOW, SUMMARIZE, or EXPLAIN of a statement that itself
    passes;
  * no write node anywhere in the tree (INSERT/UPDATE/DELETE/MERGE,
    DDL, COPY, SELECT ... INTO, SET/USE/PRAGMA, ATTACH/DETACH, or any
    unparsed ``Command``).

Anything sqlglot cannot parse is rejected (fail closed). Comments are
handled by the tokenizer, so they can neither hide a verb nor false-trip
the check, and SQL *functions* such as ``replace()`` or ``truncate()``
are not confused with statements.
"""

from __future__ import annotations

import logging
import re

import sqlglot
from sqlglot import exp
from sqlglot.dialects.dialect import Dialect
from sqlglot.errors import SqlglotError
from sqlglot.tokens import TokenType

# Statement roots that only read. exp.Query covers Select and the set
# operations (Union / Intersect / Except), including WITH-prefixed ones.
_READ_ROOTS: tuple[type[exp.Expression], ...] = tuple(
    t
    for t in (
        exp.Query,
        getattr(exp, "Describe", None),
        getattr(exp, "Show", None),
        getattr(exp, "Summarize", None),  # DuckDB: SUMMARIZE <table | query>
    )
    if t is not None
)

# Nodes that write or change session/database state. Looked up by name
# so the guard survives sqlglot versions that lack some of them.
_WRITE_NODES: tuple[type[exp.Expression], ...] = tuple(
    t
    for name in (
        "Insert",
        "Update",
        "Delete",
        "Merge",
        "Create",
        "Drop",
        "Alter",
        "AlterRename",
        "TruncateTable",
        "Copy",
        "Into",
        "LoadData",
        "Set",
        "Use",
        "Pragma",
        "Attach",
        "Detach",
        "Grant",
        "Revoke",
        "Kill",
        "Command",
    )
    if (t := getattr(exp, name, None)) is not None
)

# EXPLAIN variants (ClickHouse and DuckDB) that precede the explained
# statement, e.g. `EXPLAIN SYNTAX SELECT ...`, `EXPLAIN ANALYZE SELECT ...`.
_EXPLAIN_KIND_RE = re.compile(
    r"^\s*(?:AST|SYNTAX|QUERY\s+TREE|PLAN|PIPELINE|ESTIMATE|ANALYZE)\b",
    re.IGNORECASE,
)

# sqlglot logs a warning each time it falls back to an opaque Command
# (e.g. SHOW / EXPLAIN in the ClickHouse dialect); that is expected here.
logging.getLogger("sqlglot").setLevel(logging.ERROR)


def _dialect() -> str:
    """sqlglot dialect for the active storage backend."""
    try:
        from storage.db import backend

        return "duckdb" if backend() == "duckdb" else "clickhouse"
    except Exception:
        return "clickhouse"


def _parse_one(sql: str, dialect: str) -> exp.Expression | None:
    """Parse exactly one statement, or return None (fail closed)."""
    try:
        stmts = [s for s in sqlglot.parse(sql, read=dialect) if s is not None]
    except (SqlglotError, RecursionError):
        # RecursionError: pathologically nested input (e.g. thousands of
        # parentheses) — reject rather than crash the caller.
        return None
    return stmts[0] if len(stmts) == 1 else None


_FILE_WORDS = {"INTO", "OUTFILE", "DUMPFILE"}


def _writes_to_file(text: str, dialect: str) -> bool:
    """True if an unparsed statement tail has an INTO clause (e.g.
    ClickHouse `SHOW TABLES INTO OUTFILE '...'`). Unlexable → True."""
    try:
        tokens = Dialect.get_or_raise(dialect).tokenize(text)
    except (SqlglotError, RecursionError):
        return True
    return any(
        t.token_type == TokenType.INTO
        or (t.token_type != TokenType.STRING and t.text.upper() in _FILE_WORDS)
        for t in tokens
    )


def _check(node: exp.Expression, dialect: str, depth: int = 0) -> bool:
    if isinstance(node, exp.Command):
        verb = str(node.this or "").upper()
        # The tail is a Literal in some dialects and a plain str in others.
        tail = node.expression
        rest = tail.name if isinstance(tail, exp.Expression) else str(tail or "")
        if verb == "SHOW":
            return not _writes_to_file(rest, dialect)
        if verb == "EXPLAIN" and depth == 0:
            inner = _parse_one(_EXPLAIN_KIND_RE.sub("", rest, count=1), dialect)
            return inner is not None and _check(inner, dialect, depth + 1)
        return False
    if not isinstance(node, _READ_ROOTS):
        return False
    return node.find(*_WRITE_NODES) is None


def is_read_only(sql: str, dialect: str | None = None) -> bool:
    """True only for a single statement that cannot write.

    ``dialect`` defaults to the active backend (``clickhouse`` for chDB,
    ``duckdb`` for DuckDB).
    """
    if not sql or not sql.strip():
        return False
    dialect = dialect or _dialect()
    stmt = _parse_one(sql, dialect)
    return stmt is not None and _check(stmt, dialect)
