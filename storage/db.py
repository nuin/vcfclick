"""Embedded chDB storage with named-database multi-tenancy.

Storage model
─────────────
vcfclick keeps each cohort / study / VCF in its own small database. Every
named database is a self-contained chDB session directory under
`$VCFCLICK_HOME/dbs/<name>/`. Defaults:

    $VCFCLICK_HOME  →  ~/.vcfclick

This makes vcfclick easy to share, dump, sync, and discard at the unit a
researcher actually thinks about (a cohort), instead of forcing everything
into one mega-store.

Selecting the active database
─────────────────────────────
`get_session(name)` takes an explicit DB name. If `name` is None, the
session falls back to (in order):

  1. the VCFCLICK_DB_NAME env var (used by the CLI, ingest workers, and
     the MCP server so subprocess children inherit the same DB),
  2. the VCFCLICK_DB env var (legacy: explicit path to a chDB directory),
  3. `./.chdb` (legacy: the original single-DB layout, kept for
     backward compatibility).

The legacy paths exist because earlier versions of vcfclick had no notion
of multiple named databases. New code should pass `name` explicitly or
set VCFCLICK_DB_NAME.
"""

from __future__ import annotations

import os
import re
import shutil
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from chdb import session as _session


VCFCLICK_HOME = Path(os.environ.get("VCFCLICK_HOME", Path.home() / ".vcfclick"))
DB_ROOT = VCFCLICK_HOME / "dbs"

LEGACY_DEFAULT_PATH = Path.cwd() / ".chdb"

_TABLE_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_DB_NAME_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9_-]{0,62}$")

_sessions: dict[str, _session.Session] = {}


def _resolve_name(name: str | None) -> str | None:
    """Return the name to use, after checking the env-var fallback."""
    if name:
        return name
    return os.environ.get("VCFCLICK_DB_NAME") or None


def db_path(name: str | None = None) -> Path:
    """On-disk directory for a database (named or legacy)."""
    resolved = _resolve_name(name)
    if resolved is None:
        # Legacy fallbacks for pre-multi-DB layouts.
        if env := os.environ.get("VCFCLICK_DB"):
            return Path(env)
        return LEGACY_DEFAULT_PATH
    if not _DB_NAME_RE.match(resolved):
        raise ValueError(
            f"Unsafe DB name {resolved!r}; allowed: letter then letters, "
            "digits, underscores, or hyphens (max 63 chars)."
        )
    return DB_ROOT / resolved


def get_session(name: str | None = None) -> _session.Session:
    """Open (and cache) the chDB session for a named DB."""
    resolved = _resolve_name(name)
    cache_key = resolved or "_legacy"
    if cache_key not in _sessions:
        path = db_path(name)
        path.mkdir(parents=True, exist_ok=True)
        _sessions[cache_key] = _session.Session(str(path))
    return _sessions[cache_key]


def list_dbs() -> list[str]:
    """All named DBs under VCFCLICK_HOME/dbs/, sorted."""
    if not DB_ROOT.exists():
        return []
    return sorted(p.name for p in DB_ROOT.iterdir() if p.is_dir())


def drop_db(name: str) -> None:
    """Remove a named DB and any cached session pointing at it."""
    path = db_path(name)
    if not path.exists():
        raise FileNotFoundError(f"db {name!r} does not exist at {path}")
    if name in _sessions:
        # chDB sessions can be GC'd implicitly; explicit close would be
        # preferable but the API is stable enough that pop() suffices here.
        _sessions.pop(name)
    shutil.rmtree(path)


def db_disk_size(name: str | None = None) -> int:
    """Total bytes on disk for a named DB."""
    path = db_path(name)
    if not path.exists():
        return 0
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def insert_via_parquet(
    table: str, schema: pa.Schema, rows: list[dict]
) -> None:
    """Bulk-insert rows into the active DB via a Parquet staging file.

    Same Parquet-staging path the variants/genotypes loaders use, just
    exposed for small writes (samples, ingestions catalog) so we never
    string-interpolate user-controlled data into SQL.
    """
    if not rows:
        return
    if not _TABLE_NAME_RE.match(table):
        raise ValueError(f"Unsafe table name: {table!r}")

    sess = get_session()
    with tempfile.NamedTemporaryFile(
        suffix=".parquet", delete=False
    ) as f:
        tmp_path = Path(f.name)
    try:
        arrays = [
            pa.array([r.get(field.name) for r in rows], type=field.type)
            for field in schema
        ]
        pq.write_table(
            pa.Table.from_arrays(arrays, schema=schema), tmp_path
        )
        sess.query(
            f"INSERT INTO {table} SELECT * FROM file('{tmp_path}', 'Parquet')"
        )
    finally:
        tmp_path.unlink(missing_ok=True)


def apply_schema(schema_dir: Path | None = None) -> None:
    """Apply every .sql file in schema_dir in name order.

    Strips `--` line comments before splitting on `;` so multi-statement
    files (with comments) apply cleanly. Idempotent: CREATE TABLE
    statements should be IF NOT EXISTS-friendly OR caller should run
    this against a fresh DB.
    """
    if schema_dir is None:
        # Default: the repo's schema/ directory next to this module.
        schema_dir = Path(__file__).parent.parent / "schema"

    sess = get_session()
    for f in sorted(Path(schema_dir).glob("*.sql")):
        text = f.read_text()
        clean = "\n".join(re.sub(r"--.*$", "", ln) for ln in text.splitlines())
        for stmt in clean.split(";"):
            stmt = stmt.strip()
            if stmt:
                sess.query(stmt)
