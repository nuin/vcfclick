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
import fcntl
import logging
import os
import re
import shutil
import tempfile
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from chdb import session as _session

log = logging.getLogger(__name__)


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
    return _db_root() / resolved


# Substrings that flag the known chDB EmbeddedServer async-load race
# (a `recursive_mutex` invalid-state failure during session init). The
# race is intermittent — retrying with a small backoff has been enough
# in observation to bring per-open success to effectively 100%. If
# the exception message contains NONE of these substrings, we surface
# it immediately so genuine errors aren't masked by the retry loop.
_CHDB_RACE_MARKERS = (
    "BAD_ARGUMENTS",
    "recursive_mutex lock failed",
    "ASYNC_LOAD_WAIT_FAILED",
)

# Backoff schedule. Total worst-case retry cost is ~1.3 s on cold open;
# well under the chDB session-init baseline so retried opens still feel
# instant in normal CLI use.
_CHDB_RETRY_DELAYS_S = (0.1, 0.3, 0.9)


def _open_session_with_retry(path: str) -> _session.Session:
    """Open a chDB session, retrying through the known EmbeddedServer
    async-load race. Non-race exceptions propagate immediately."""
    last_exc: Exception | None = None
    for attempt in range(len(_CHDB_RETRY_DELAYS_S) + 1):
        try:
            return _session.Session(path)
        except Exception as exc:
            msg = str(exc)
            if not any(marker in msg for marker in _CHDB_RACE_MARKERS):
                raise
            last_exc = exc
            if attempt < len(_CHDB_RETRY_DELAYS_S):
                delay = _CHDB_RETRY_DELAYS_S[attempt]
                log.warning(
                    "[storage] chDB session open hit the async-load race "
                    "(attempt %d/%d), retrying in %.1fs",
                    attempt + 1,
                    len(_CHDB_RETRY_DELAYS_S) + 1,
                    delay,
                )
                time.sleep(delay)
    assert last_exc is not None  # narrows for type-checkers
    raise last_exc


def get_session(name: str | None = None) -> _session.Session:
    """Open (and cache) the chDB session for a named DB."""
    resolved = _resolve_name(name)
    cache_key = resolved or "_legacy"
    if cache_key not in _sessions:
        path = db_path(name)
        path.mkdir(parents=True, exist_ok=True)
        _sessions[cache_key] = _open_session_with_retry(str(path))
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

    Runs synchronously (`mutations_sync = 2`) so callers can rely on
    the cleanup having completed before they re-raise the original
    error.
    """
    validate_ingest_id(ingest_id)
    sess = get_session()
    for table in _INGESTION_SCOPED_TABLES:
        sess.query(
            f"ALTER TABLE {table} DELETE WHERE ingest_id = '{ingest_id}' "
            f"SETTINGS mutations_sync = 2"
        )


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
        # INSERT is immune to chDB shifting Parquet handling from
        # name-based to positional mapping. Same rationale as the
        # variants/genotypes loaders in ingest.vcf_load / ingest.parallel.
        cols = ", ".join(f"`{f.name}`" for f in schema)
        sess.query(
            f"INSERT INTO {table} ({cols}) "
            f"SELECT {cols} FROM file({sql_quote_str(str(tmp_path))}, 'Parquet')"
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
