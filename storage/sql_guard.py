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
    unparsed ``Command``);
  * every FROM / JOIN source is a plain table identifier, at most
    qualified with the default schema (``main`` / ``default``): no table
    function (``file()``, ``url()``, ``read_csv()``, ``numbers()``, ...),
    no string-literal path, no ``system.*`` or other catalog;
  * no scalar function that reads a file or the environment (``file()``,
    ``getenv()``, ``read_text()``, ...);
  * no ``SETTINGS`` clause (a per-query override of the session's limits
    and of ``readonly``), and no ``SHOW`` of ``system`` or the process list.

Anything sqlglot cannot parse is rejected (fail closed). Comments are
handled by the tokenizer, so they can neither hide a verb nor false-trip
the check, and SQL *functions* such as ``replace()`` or ``truncate()``
are not confused with statements.

The guard is the first of two layers. `storage.get_readonly_session`
also locks the session itself (DuckDB ``enable_external_access = false``
plus ``lock_configuration``, chDB ``readonly = 1``), so a statement the
parser misreads still cannot reach a file or change a setting.
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

# Table-source nodes sqlglot builds for DuckDB file readers and generators
# when it recognises them (an unknown function becomes `Anonymous`), plus
# `Schema`, a table with a column list: in a read it only arises when the
# parser takes a table function's arguments for columns
# (`DESCRIBE TABLE file('/x', 'CSV')`).
_EXTERNAL_NODES: tuple[type[exp.Expression], ...] = tuple(
    t
    for name in ("ReadCSV", "ReadParquet", "ReadJSON", "GenerateSeries", "Schema")
    if (t := getattr(exp, name, None)) is not None
)

# Scalar functions that read a file, a URL or the environment. Table
# functions are refused wholesale by `_table_ok`; this list is for the same
# names in a SELECT list or WHERE clause (ClickHouse `file(path)` returns
# the file's content as a string).
_EXTERNAL_FUNCTIONS = frozenset(
    {
        "file",
        "url",
        "getenv",
        "glob",
        "read_text",
        "read_blob",
        "read_csv",
        "read_csv_auto",
        "read_parquet",
        "read_json",
        "read_json_auto",
        "read_ndjson",
        "read_ndjson_auto",
        "sniff_csv",
    }
)

# A FROM source must be a bare table name: the cohort tables, or a CTE.
_TABLE_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

# The only schema qualifier accepted: each backend's default schema.
_DEFAULT_SCHEMAS = frozenset({"main", "default"})

# Words in an unparsed `SHOW ...` tail that reach server state.
_SYSTEM_WORDS = frozenset({"SYSTEM", "INFORMATION_SCHEMA", "PROCESSLIST"})

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


def _show_tail_ok(text: str, dialect: str) -> bool:
    """False if an unparsed `SHOW` tail has an INTO clause (ClickHouse
    `SHOW TABLES INTO OUTFILE '...'` writes a file) or names server state
    (`SHOW TABLES FROM system`, `SHOW PROCESSLIST`). Unlexable → False."""
    try:
        tokens = Dialect.get_or_raise(dialect).tokenize(text)
    except (SqlglotError, RecursionError):
        return False
    for t in tokens:
        if t.token_type == TokenType.INTO:
            return False
        if t.token_type == TokenType.STRING:
            continue
        word = t.text.upper()
        if word in _FILE_WORDS or word in _SYSTEM_WORDS:
            return False
    return True


def _table_ok(table: exp.Table) -> bool:
    """A FROM / JOIN source is a bare table name (a cohort table or a CTE),
    at most qualified with the default schema. Everything else (a table
    function, a quoted path, another catalog, `system.*`) is refused."""
    if not isinstance(table.this, exp.Identifier):
        return False
    name = table.this.name
    # A dialect that lacks a table function parses `file('/x')` as a table
    # named `file` with alias columns; no cohort table carries such a name.
    if not _TABLE_NAME_RE.match(name) or name.lower() in _EXTERNAL_FUNCTIONS:
        return False
    if table.args.get("catalog"):
        return False
    db = table.args.get("db")
    return db is None or db.name.lower() in _DEFAULT_SCHEMAS


def _show_ok(show: exp.Show) -> bool:
    """A parsed SHOW may list the cohort's own tables and columns, never
    another schema's (`SHOW TABLES FROM system`)."""
    return not (show.args.get("from_") or show.args.get("db"))


def _tree_ok(node: exp.Expression) -> bool:
    """No write, no external source, no settings override below `node`."""
    if node.find(*_WRITE_NODES, *_EXTERNAL_NODES) is not None:
        return False
    if any(not _table_ok(t) for t in node.find_all(exp.Table)):
        return False
    if any(not _show_ok(s) for s in node.find_all(exp.Show)):
        return False
    if any(f.name.lower() in _EXTERNAL_FUNCTIONS for f in node.find_all(exp.Anonymous)):
        return False
    return not any(s.args.get("settings") for s in node.find_all(exp.Select))


def _check(node: exp.Expression, dialect: str, depth: int = 0) -> bool:
    if isinstance(node, exp.Command):
        verb = str(node.this or "").upper()
        # The tail is a Literal in some dialects and a plain str in others.
        tail = node.expression
        rest = tail.name if isinstance(tail, exp.Expression) else str(tail or "")
        if verb == "SHOW":
            return _show_tail_ok(rest, dialect)
        if verb == "EXPLAIN" and depth == 0:
            inner = _parse_one(_EXPLAIN_KIND_RE.sub("", rest, count=1), dialect)
            return inner is not None and _check(inner, dialect, depth + 1)
        return False
    if not isinstance(node, _READ_ROOTS):
        return False
    return _tree_ok(node)


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
