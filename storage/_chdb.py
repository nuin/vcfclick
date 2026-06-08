"""chDB backend.

Wraps `chdb.session.Session` with the EmbeddedServer async-load
race retry that previously lived in `storage.db._open_session_with_retry`.
Lazy-imports `chdb` so installations without the `[chdb]` extra can
still run on the DuckDB backend without import errors.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

log = logging.getLogger(__name__)


# Substrings that flag the known chDB EmbeddedServer async-load race
# (a `recursive_mutex` invalid-state failure during session init).
_CHDB_RACE_MARKERS = (
    "BAD_ARGUMENTS",
    "recursive_mutex lock failed",
    "ASYNC_LOAD_WAIT_FAILED",
)

# Backoff schedule. Total worst-case retry cost is ~1.3 s on cold open;
# well under the chDB session-init baseline.
_CHDB_RETRY_DELAYS_S = (0.1, 0.3, 0.9)


def open_session(path: Path):
    """Open a chDB session, retrying through the known EmbeddedServer
    async-load race. Lazy-imports chdb so DuckDB-only installs don't
    error at module-import time.
    """
    try:
        from chdb import session as _session
    except ImportError as e:
        raise RuntimeError(
            "VCFCLICK_BACKEND=chdb requested but chdb is not installed. "
            "Run `pip install vcfclick[chdb]` to add the chDB backend, "
            "or set VCFCLICK_BACKEND=duckdb."
        ) from e

    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    last_exc: Exception | None = None
    for attempt in range(len(_CHDB_RETRY_DELAYS_S) + 1):
        try:
            return _session.Session(str(path))
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
    assert last_exc is not None
    raise last_exc
