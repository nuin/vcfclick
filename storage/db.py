"""Embedded chDB session — the single storage backend for vcfclick.

Why chDB and not the ClickHouse server:
  - Single binary deployment. `pip install vcfclick && vcfclick demo`
    instead of "install Docker, pull ClickHouse, configure ports."
  - Same SQL surface, same MergeTree engines, same projections.
  - One less process, one less port, one less Gatekeeper dialog.

Phase 1 limitation made explicit:
  chDB's persistent session is a directory on disk, opened exclusively
  by one process at a time. The ingester and the MCP server cannot run
  concurrently against the same DB path. Run one at a time, or use the
  Parquet export path for "looking at the data while ingesting" needs.

  Phase 2 option: switch to chdb-server (chDB's TCP server mode) for
  concurrent reader/writer access. Same SQL, modest operational cost.

Override the DB location via the VCFCLICK_DB env var (default: ./.chdb).
"""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from chdb import session as _session

DB_PATH = Path(os.environ.get("VCFCLICK_DB", Path.cwd() / ".chdb"))


_cached_session: _session.Session | None = None


def get_session() -> _session.Session:
    """Open (and cache) the chDB session at DB_PATH."""
    global _cached_session
    if _cached_session is None:
        DB_PATH.mkdir(parents=True, exist_ok=True)
        _cached_session = _session.Session(str(DB_PATH))
    return _cached_session


_TABLE_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def insert_via_parquet(table: str, schema: pa.Schema, rows: list[dict]) -> None:
    """Bulk-insert rows into a chDB table via a Parquet staging file.

    Uses the same Parquet-staging path the variants/genotypes loaders use,
    just exposed for small writes (samples, ingestions catalog) so we never
    string-interpolate user-controlled data into SQL. Sample IDs from a
    hostile VCF header could otherwise inject — Parquet's typed columnar
    serialisation removes that surface entirely.

    `table` is validated against a strict identifier regex; it must be
    caller-controlled (never user input from a VCF or HTTP request).
    """
    if not rows:
        return
    if not _TABLE_NAME_RE.match(table):
        raise ValueError(f"Unsafe table name: {table!r}")

    sess = get_session()
    with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
        tmp_path = Path(f.name)
    try:
        arrays = [
            pa.array([r.get(field.name) for r in rows], type=field.type)
            for field in schema
        ]
        pq.write_table(pa.Table.from_arrays(arrays, schema=schema), tmp_path)
        # tmp_path is from tempfile.NamedTemporaryFile — controlled, safe
        # to interpolate into SQL (no quotes, no shell metacharacters).
        sess.query(f"INSERT INTO {table} SELECT * FROM file('{tmp_path}', 'Parquet')")
    finally:
        tmp_path.unlink(missing_ok=True)


def apply_schema(schema_dir: Path) -> None:
    """Apply every .sql file in schema_dir in name order.

    Strips `--` line comments before splitting on `;` so multi-statement
    files (with comments) apply cleanly. Idempotent: CREATE TABLE
    statements should be IF NOT EXISTS-friendly OR caller should run
    this against a fresh DB.
    """
    sess = get_session()
    for f in sorted(Path(schema_dir).glob("*.sql")):
        text = f.read_text()
        # Strip line comments before splitting on `;` — comment fragments
        # otherwise contaminate statement boundaries.
        clean = "\n".join(re.sub(r"--.*$", "", ln) for ln in text.splitlines())
        for stmt in clean.split(";"):
            stmt = stmt.strip()
            if stmt:
                sess.query(stmt)
