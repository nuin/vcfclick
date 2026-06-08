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

import contextlib
import importlib
import logging
import os
import re
import shutil
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

log = logging.getLogger(__name__)


# `VCFCLICK_BACKEND` picks the variant store:
#   - "chdb"    → chDB (ClickHouse engine, requires the [chdb] extra)
#   - "duckdb"  → DuckDB (no extra needed; always available)
# If unset, auto-detect: prefer chdb (legacy default; existing user
# DBs are chdb-format) when the import succeeds; fall back to duckdb.
# A future major release may flip the default once duckdb is the
# common case via bioconda. Until then no user gets surprised by a
# silent backend change.
_VALID_BACKENDS = ("chdb", "duckdb")


def backend() -> str:
    """Return the active storage backend name."""
    raw = os.environ.get("VCFCLICK_BACKEND")
    if raw:
        b = raw.strip().lower()
        if b not in _VALID_BACKENDS:
            raise ValueError(
                f"VCFCLICK_BACKEND={raw!r} not recognised; "
                f"expected one of {_VALID_BACKENDS}"
            )
        return b
    # Auto-detect: chdb wins if importable (legacy default), else duckdb.
    try:
        import chdb  # noqa: F401

        return "chdb"
    except ImportError:
        return "duckdb"


# Paths are computed at every call rather than cached at import time.
# Tests (and any in-process flow that changes VCFCLICK_HOME after import)
# depend on this — module-level caching was a footgun: a test that
# monkeypatched the env after `import storage` would still see the
# old path.
def _vcfclick_home() -> Path:
    return Path(os.environ.get("VCFCLICK_HOME", str(Path.home() / ".vcfclick")))


def _db_root() -> Path:
    return _vcfclick_home() / "dbs"


# Backward-compatible attribute access — code that imports
# `storage.VCFCLICK_HOME` / `storage.DB_ROOT` still works, but the
# value is recomputed on every attribute lookup via __getattr__ below.
LEGACY_DEFAULT_PATH = Path.cwd() / ".chdb"


def __getattr__(name: str):
    if name == "VCFCLICK_HOME":
        return _vcfclick_home()
    if name == "DB_ROOT":
        return _db_root()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


_TABLE_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_DB_NAME_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9_-]{0,62}$")

_sessions: dict[str, object] = {}


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
    return _db_root() / resolved


def _open_session(path: Path):
    """Open a backend-specific session at `path`."""
    if backend() == "duckdb":
        from storage._duckdb import DuckDBSession

        return DuckDBSession(path)
    from storage._chdb import open_session as _open_chdb

    return _open_chdb(path)


def get_session(name: str | None = None):
    """Open (and cache) the storage session for a named DB."""
    resolved = _resolve_name(name)
    # Cache key includes the backend so the same Python process
    # juggling both engines (tests, in-process backend switch)
    # doesn't hand back the wrong handle.
    cache_key = f"{backend()}::{resolved or '_legacy'}"
    if cache_key not in _sessions:
        path = db_path(name)
        path.mkdir(parents=True, exist_ok=True)
        _sessions[cache_key] = _open_session(path)
    return _sessions[cache_key]


def list_dbs() -> list[str]:
    """All named DBs under VCFCLICK_HOME/dbs/, sorted."""
    root = _db_root()
    if not root.exists():
        return []
    return sorted(p.name for p in root.iterdir() if p.is_dir())


def drop_db(name: str) -> None:
    """Remove a named DB and any cached session pointing at it."""
    path = db_path(name)
    if not path.exists():
        raise FileNotFoundError(f"db {name!r} does not exist at {path}")
    cache_key = f"{backend()}::{_resolve_name(name) or '_legacy'}"
    session = _sessions.pop(cache_key, None)
    if hasattr(session, "close"):
        session.close()
    shutil.rmtree(path)


def db_disk_size(name: str | None = None) -> int:
    """Total bytes on disk for a named DB."""
    path = db_path(name)
    if not path.exists():
        return 0
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


# Ingest IDs come from user CLI input. Validating them up front lets us
# safely interpolate the value into the rollback DELETE statements below
# without quote-escaping — chDB's parameterised binding does not extend
# to ALTER TABLE ... DELETE expressions.
_INGEST_ID_RE = re.compile(r"^[A-Za-z0-9_.\-]+$")


def validate_ingest_id(ingest_id: str) -> None:
    """Reject ingest IDs that aren't safe to interpolate into SQL.

    Allowed: ASCII letters, digits, underscore, dot, hyphen. This covers
    UUIDs, batch labels like `batch_a`, dated tags like `2026.05.31`,
    and the slugs users actually type. Anything else (quotes, spaces,
    semicolons) is rejected.
    """
    if not _INGEST_ID_RE.fullmatch(ingest_id):
        raise ValueError(
            f"invalid ingest_id {ingest_id!r}: only ASCII letters, "
            f"digits, underscore, dot, hyphen allowed"
        )


# Tables that carry per-ingestion rows — these are the ones rollback
# has to scrub when an ingest fails mid-stream. Keep in sync with the
# schema/*.sql files.
_INGESTION_SCOPED_TABLES = ("variants", "genotypes", "samples", "ingestions")


def sql_quote_str(s: str) -> str:
    """ClickHouse string literal — backslashes AND quotes both escaped.

    ClickHouse / chDB recognises two escape forms inside string
    literals: `''` is a single quote, and `\\'` is also a single
    quote. Plain quote-doubling (the portable-SQL form) is therefore
    NOT enough on its own: a payload like `\\'; DROP TABLE …` would,
    after quote-doubling, become `\\''; DROP TABLE …` — the engine
    reads the first quote as backslash-escaped, the second as the
    string terminator, and arbitrary SQL after it. Codex round 9
    caught exactly this.

    Order of operations matters: escape backslashes FIRST (each `\\`
    → `\\\\`), then quotes (each `'` → `''`). Doing it the other way
    would double the just-added backslashes again.
    """
    escaped = s.replace("\\", "\\\\").replace("'", "''")
    return "'" + escaped + "'"


def _fcntl_module():
    """Load POSIX fcntl only where the ingest lock is used."""
    return importlib.import_module("fcntl")


@contextlib.contextmanager
def ingest_id_lock(ingest_id: str):
    """Serialize concurrent ingests under the same `(DB, ingest_id)`.

    Without this, two `vcfclick db ingest` invocations sharing an
    ingest_id (across subprocesses, or two threads, or a workflow
    runner firing them in parallel) would race on the staging
    directory, the rollback_ingest() call, and the bulk-import
    glob — one run could delete or overwrite the other's Parquets
    and import a mixed set.

    Uses an exclusive `fcntl.flock` on a per-id lockfile under
    `<DB>/locks/`. Blocking — the second caller waits until the
    first releases. Auto-released when the context exits (process
    death, exception, or normal completion all close the FD).

    Unix-only: fcntl is a Linux + macOS API. Windows isn't a
    supported target for vcfclick today.
    """
    validate_ingest_id(ingest_id)
    lock_dir = db_path() / "locks"
    lock_dir.mkdir(parents=True, exist_ok=True)
    lock_path = lock_dir / f"{ingest_id}.lock"

    fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o600)
    fcntl = _fcntl_module()
    try:
        # LOCK_EX blocks until acquired. If the holding process dies
        # mid-ingest, the kernel releases the lock automatically and
        # the next caller proceeds.
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)


def rollback_ingest(ingest_id: str) -> None:
    """Delete every row carrying `ingest_id` from the active DB.

    Used by both ingest paths to leave the DB clean when ingestion
    fails partway — without this, a half-loaded ingest_id stays
    queryable and the user has to `db rm` the whole cohort to recover.

    Backend-specific:
      * chDB uses `ALTER TABLE … DELETE WHERE` with
        `mutations_sync = 2` so the deletion is observable before
        the next INSERT runs.
      * DuckDB uses standard `DELETE FROM … WHERE`; deletions are
        synchronous and observable in the same transaction.
    """
    validate_ingest_id(ingest_id)
    sess = get_session()
    for table in _INGESTION_SCOPED_TABLES:
        sess.query(delete_where_sql(table, f"ingest_id = '{ingest_id}'"))


def delete_where_sql(table: str, where: str) -> str:
    """Return the backend-specific synchronous DELETE statement."""
    if backend() == "duckdb":
        return f"DELETE FROM {table} WHERE {where}"
    return f"ALTER TABLE {table} DELETE WHERE {where} SETTINGS mutations_sync = 2"


def parquet_file_expr(path: str) -> str:
    """Return the backend-specific SQL fragment for reading a Parquet file.

      * chDB:    `file('/abs/path.parquet', 'Parquet')`
      * DuckDB:  `read_parquet('/abs/path.parquet')`

    Pass the value already wrapped in a quoted string literal — the
    caller is responsible for using `sql_quote_str` on the path.
    """
    if backend() == "duckdb":
        return f"read_parquet({sql_quote_str(path)})"
    return f"file({sql_quote_str(path)}, 'Parquet')"


def count_expr() -> str:
    """ClickHouse permits `count()`; DuckDB requires `count(*)`."""
    return "count(*)" if backend() == "duckdb" else "count()"


def table_exists(name: str) -> bool:
    """Return True if `name` is a table in the active database.

    chDB exposes a `system.tables` table; DuckDB exposes the SQL-standard
    `information_schema.tables`. Wrap the dialect difference here.
    """
    if not _TABLE_NAME_RE.match(name):
        raise ValueError(f"Unsafe table name: {name!r}")
    sess = get_session()
    if backend() == "duckdb":
        sql = (
            "SELECT count(*) FROM information_schema.tables "
            f"WHERE table_name = '{name}'"
        )
    else:
        sql = (
            "SELECT count() FROM system.tables "
            f"WHERE database = currentDatabase() AND name = '{name}'"
        )
    result = sess.query(sql, "CSV").bytes().decode().strip()
    # Strip a possible CSV header on the DuckDB path (chDB CSV is headerless).
    last = result.splitlines()[-1].strip() if result else "0"
    return int(last) > 0


def insert_via_parquet(table: str, schema: pa.Schema, rows: list[dict]) -> None:
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
    with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
        tmp_path = Path(f.name)
    try:
        arrays = [
            pa.array([r.get(field.name) for r in rows], type=field.type)
            for field in schema
        ]
        pq.write_table(pa.Table.from_arrays(arrays, schema=schema), tmp_path)
        # Explicit column list derived from the Arrow schema so the
        # INSERT is immune to either engine shifting Parquet handling
        # from name-based to positional mapping. Same rationale as the
        # variants/genotypes loaders in ingest.vcf_load / ingest.parallel.
        # Double-quoted identifiers work on both chDB and DuckDB.
        cols = ", ".join(f'"{f.name}"' for f in schema)
        sess.query(
            f"INSERT INTO {table} ({cols}) "
            f"SELECT {cols} FROM {parquet_file_expr(str(tmp_path))}"
        )
    finally:
        tmp_path.unlink(missing_ok=True)


def schema_dir_for_backend() -> Path:
    """Return the SQL-file directory matching the active backend.

      * chDB    → repo/schema/*.sql       (the original chDB-flavoured DDL)
      * DuckDB  → repo/schema/duckdb/*.sql

    The two directories carry the same set of CREATE TABLE statements
    and the same column names and column ORDER; only the type names
    and engine clauses differ.
    """
    base = Path(__file__).parent.parent / "schema"
    if backend() == "duckdb":
        return base / "duckdb"
    return base


def apply_schema(schema_dir: Path | None = None) -> None:
    """Apply every .sql file in `schema_dir` in name order.

    Strips `--` line comments before splitting on `;` so multi-statement
    files apply cleanly. Idempotent against a fresh DB; CREATE TABLE
    statements are not IF NOT EXISTS, so callers are expected to invoke
    this once at db_create.
    """
    if schema_dir is None:
        schema_dir = schema_dir_for_backend()

    sess = get_session()
    for f in sorted(Path(schema_dir).glob("*.sql")):
        text = f.read_text()
        clean = "\n".join(re.sub(r"--.*$", "", ln) for ln in text.splitlines())
        for stmt in clean.split(";"):
            stmt = stmt.strip()
            if stmt:
                sess.query(stmt)
