"""DuckDB backend.

Mirrors the surface of the chDB session that the rest of the codebase
depends on so call sites stay unchanged. The two engines differ in:

  * SQL dialect (ClickHouse vs DuckDB) — handled by backend-specific
    helpers in storage.db (`rollback_where_sql`, `parquet_file_expr`,
    `count_expr`) and by per-backend schema directories.

  * Result delivery — chDB returns a result handle with `.bytes()`
    plus a chDB-specific FORMAT clause. DuckDB returns Python rows via
    `cursor.fetchall()`. The `DuckDBSession` class below renders
    fetchall rows into the formats vcfclick actually uses
    (JSONCompact, TSV, CSV, default), all decode-to-str at call
    sites that look like `sess.query(...).bytes().decode()`.

  * Concurrency model — chDB has the EmbeddedServer / async-load race
    handled by `_open_session_with_retry`. DuckDB uses POSIX file
    locking on its own file handle; multi-process access to the same
    DuckDB file is single-writer enforced by DuckDB itself. The
    storage layer's `ingest_id_lock` (fcntl) still serialises
    same-(DB, ingest_id) ingests at the application level.
"""

from __future__ import annotations

import csv
import io
import json
import logging
import re
from pathlib import Path

import duckdb

log = logging.getLogger(__name__)


# Map `FORMAT <name>` suffix on a chDB-style query to a renderer.
# Only the formats vcfclick actually uses are supported; anything
# else falls back to a default Pretty rendering for `db query`.
_SUPPORTED_FORMATS = {"JSONCompact", "TSV", "CSV", "Vertical", "PrettyCompact"}

# Strip an arbitrary trailing `FORMAT <name>` clause off a SQL statement.
# chDB accepts the clause as a query terminator. DuckDB does not — the
# format is requested out-of-band via the renderer.
_FORMAT_CLAUSE_RE = re.compile(r"\s+FORMAT\s+(\w+)\s*;?\s*$", re.IGNORECASE)


class _DuckDBQueryResult:
    """Mimics the result of `chdb.session.Session.query()` enough for
    vcfclick call sites. Specifically:

        sess.query(sql).bytes().decode()
        sess.query(sql, "CSV").bytes().decode()

    both work. `.bytes()` returns the rendered output for the format
    that was requested (or the default if none).
    """

    def __init__(self, payload: bytes) -> None:
        self._payload = payload

    def bytes(self) -> bytes:
        return self._payload

    def __bytes__(self) -> bytes:
        return self._payload

    def __str__(self) -> str:
        return self._payload.decode()


def _render_jsoncompact(columns: list[str], rows: list[tuple]) -> bytes:
    """chDB's JSONCompact format: `{"meta": [...], "data": [[...], ...], "rows": N}`."""
    payload = {
        "meta": [{"name": c} for c in columns],
        "data": [list(r) for r in rows],
        "rows": len(rows),
    }
    return json.dumps(payload, default=str).encode()


def _fmt_cell(v: object) -> str:
    """Render one cell of a TSV/CSV row in a way that matches chDB output:

    * NULL → literal `\\N`.
    * float → `.6g` to cap at Float32-ish precision (DuckDB returns
      REAL values widened to Python double, whose default `str()` form
      shows trailing ULP noise like `0.9200000166893005`; chDB hides
      that at the renderer level).
    * everything else → `str(v)`.
    """
    if v is None:
        return "\\N"
    if isinstance(v, float):
        # `.6g` drops trailing zeros, prints integers without a decimal
        # point, and caps at six significant digits — visually identical
        # to chDB's Float32 rendering for the value ranges vcfclick uses.
        return f"{v:.6g}"
    return str(v)


def _render_tsv(rows: list[tuple]) -> bytes:
    """Match chDB's TabSeparated: raw tab-joined fields, NULL → `\\N`."""
    out = ["\t".join(_fmt_cell(v) for v in r) for r in rows]
    return ("\n".join(out) + ("\n" if out else "")).encode()


def _render_csv(rows: list[tuple]) -> bytes:
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    for r in rows:
        writer.writerow([_fmt_cell(v) for v in r])
    return buf.getvalue().encode()


def _render_vertical(columns: list[str], rows: list[tuple]) -> bytes:
    """One row per block; one column per line. Matches chDB Vertical."""
    width = max((len(c) for c in columns), default=0)
    out = []
    for i, r in enumerate(rows):
        out.append(f"Row {i + 1}:")
        out.append("──────")
        for c, v in zip(columns, r, strict=False):
            out.append(f"{c:<{width}}: {'' if v is None else v}")
        out.append("")
    return ("\n".join(out)).encode()


def _render_pretty(columns: list[str], rows: list[tuple]) -> bytes:
    """Lightweight Pretty/PrettyCompact replacement. chDB's box-drawing
    output isn't replicated 1:1; users get a readable table with column
    headers."""
    widths = [len(c) for c in columns]
    for r in rows:
        for i, v in enumerate(r):
            widths[i] = max(widths[i], len("" if v is None else str(v)))
    line = " │ ".join(c.ljust(widths[i]) for i, c in enumerate(columns))
    sep = "─┼─".join("─" * w for w in widths)
    body = []
    for r in rows:
        body.append(
            " │ ".join(
                ("" if v is None else str(v)).ljust(widths[i]) for i, v in enumerate(r)
            )
        )
    return ("\n".join([line, sep, *body])).encode()


class DuckDBSession:
    """Backend-specific session wrapping a duckdb.DuckDBPyConnection.

    The chDB-shaped `query(sql, format=None)` API is preserved so the
    bulk of vcfclick code (which does `sess.query(...).bytes().decode()`)
    works without modification.
    """

    def __init__(self, path: Path) -> None:
        self._path = Path(path)
        self._path.mkdir(parents=True, exist_ok=True)
        # Single .duckdb file per database directory. Keeps DB content
        # co-located with locks/ and any other vcfclick working files,
        # matching the chDB directory model.
        self._db_file = self._path / "store.duckdb"
        self._conn = duckdb.connect(str(self._db_file))

    def query(self, sql: str, format: str | None = None) -> _DuckDBQueryResult:
        """Run a SQL statement. Format defaults to a Pretty render."""
        # Strip any `FORMAT <name>` suffix the caller appended chDB-style.
        # An inline FORMAT clause is explicit and wins over the `format`
        # kwarg, which is often a default ("PrettyCompact") rather than
        # a deliberate caller choice.
        m = _FORMAT_CLAUSE_RE.search(sql)
        if m:
            sql = _FORMAT_CLAUSE_RE.sub("", sql).rstrip().rstrip(";")
            format = m.group(1)

        cur = self._conn.execute(sql)

        # DDL/DML — no rows to render. Return empty payload.
        if cur.description is None:
            return _DuckDBQueryResult(b"")

        columns = [d[0] for d in cur.description]
        rows = cur.fetchall()

        fmt = (format or "PrettyCompact").strip().lower()
        # chDB has overlapping format aliases; map them to our four
        # canonical renderers so callers can use either name.
        if fmt in ("jsoncompact", "json"):
            return _DuckDBQueryResult(_render_jsoncompact(columns, rows))
        if fmt in ("tsv", "tabseparated", "tabseparatedwithnames"):
            return _DuckDBQueryResult(_render_tsv(rows))
        if fmt in ("csv", "csvwithnames"):
            return _DuckDBQueryResult(_render_csv(rows))
        if fmt == "vertical":
            return _DuckDBQueryResult(_render_vertical(columns, rows))
        # Default + PrettyCompact + anything unknown — Pretty render.
        return _DuckDBQueryResult(_render_pretty(columns, rows))

    def close(self) -> None:
        try:
            self._conn.close()
        except Exception as exc:
            log.warning("DuckDB connection close failed for %s: %s", self._db_file, exc)
