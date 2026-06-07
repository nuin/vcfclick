"""Tests for the two storage hardening pieces codex round 8 caught:

1. `sql_quote_str` — SQL-standard quote-doubling for any string
   interpolated into a chDB query. Defence against single quotes in
   paths (the public `staging_dir` parameter to `ingest_parallel`
   accepts arbitrary callers).

2. `ingest_id_lock` — fcntl-based per-(DB, ingest_id) file lock that
   serialises concurrent ingests sharing an ingest_id. Without it,
   two `vcfclick db ingest` calls running in parallel would race on
   the staging dir and the rollback/import path.
"""

from __future__ import annotations

import multiprocessing
import os
import shutil
import subprocess
import time
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")


# ─────────────────────── sql_quote_str ───────────────────────


def test_sql_quote_str_plain_string():
    from storage import sql_quote_str

    assert sql_quote_str("plain") == "'plain'"


def test_sql_quote_str_doubles_embedded_quote():
    """SQL standard: a literal single quote inside a string is two
    consecutive single quotes."""
    from storage import sql_quote_str

    assert sql_quote_str("o'brien") == "'o''brien'"


def test_sql_quote_str_path_with_quote_doesnt_break_out():
    """The path `/tmp/x'.parquet` interpolated into
    `file('{path}', 'Parquet')` would otherwise close the file
    literal at position 9 and let the rest of the path appear as
    raw SQL — classic injection vector."""
    from storage import sql_quote_str

    path = "/tmp/x'; DROP TABLE variants; --"
    quoted = sql_quote_str(path)
    # The quoted form is one balanced single-quoted literal.
    assert quoted.startswith("'") and quoted.endswith("'")
    # All embedded single quotes are doubled.
    body = quoted[1:-1]
    assert "'" not in body.replace("''", "")


def test_sql_quote_str_classic_injection_payload():
    """Make sure the famous injection pattern would be safely
    contained as a string literal."""
    from storage import sql_quote_str

    quoted = sql_quote_str("' OR 1=1 --")
    # When pasted into `WHERE x = {quoted}`, the result is
    # `WHERE x = ''' OR 1=1 --'` — three single quotes at the
    # start (empty string + doubled), then literal content, then
    # closing quote. The OR is inside the string literal.
    assert quoted == "''' OR 1=1 --'"


# ─────────────────────── ingest_id_lock ───────────────────────


def test_ingest_id_lock_releases_on_context_exit(vcfclick_home, monkeypatch):
    """A single use of the lock acquires and releases cleanly. The
    second use in the same process must not block forever."""
    monkeypatch.chdir(REPO)  # ensure storage's paths resolve correctly
    # Initialise a DB so db_path() resolves to a real directory.
    subprocess.run(
        [VCFCLICK_BIN, "db", "create", "smoke"],
        env={**os.environ, "VCFCLICK_HOME": str(vcfclick_home)},
        check=True,
        capture_output=True,
    )
    os.environ["VCFCLICK_HOME"] = str(vcfclick_home)
    os.environ["VCFCLICK_DB_NAME"] = "smoke"

    from storage import ingest_id_lock

    with ingest_id_lock("batch_a"):
        pass  # acquire + release

    # Second acquisition must succeed quickly (lock was released).
    started = time.time()
    with ingest_id_lock("batch_a"):
        pass
    assert time.time() - started < 1.0, "lock did not release"


def _hold_lock(
    home_str: str,
    ingest_id: str,
    hold_seconds: float,
    acquired_evt,
) -> None:
    """Helper: child process acquires the lock, signals the parent
    via `acquired_evt`, then holds the lock for `hold_seconds`."""
    os.environ["VCFCLICK_HOME"] = home_str
    os.environ["VCFCLICK_DB_NAME"] = "smoke"
    from storage import ingest_id_lock

    with ingest_id_lock(ingest_id):
        acquired_evt.set()  # tells the parent we're in the critical section
        time.sleep(hold_seconds)


def test_ingest_id_lock_blocks_concurrent_holder(vcfclick_home):
    """A second process trying to acquire the same lock must block
    until the first releases. The child signals the parent via an
    `Event` when it actually holds the lock; the parent then tries
    to acquire and times its wait."""
    subprocess.run(
        [VCFCLICK_BIN, "db", "create", "smoke"],
        env={**os.environ, "VCFCLICK_HOME": str(vcfclick_home)},
        check=True,
        capture_output=True,
    )

    # Use spawn so the child runs a fresh Python — no inherited
    # in-process locks, mirrors the real-world subprocess case.
    ctx = multiprocessing.get_context("spawn")
    acquired = ctx.Event()
    hold = 0.6
    child = ctx.Process(
        target=_hold_lock, args=(str(vcfclick_home), "batch_a", hold, acquired)
    )
    child.start()

    # Wait until the child confirms it holds the lock. 5s ceiling
    # catches the case where the child crashed before acquiring.
    assert acquired.wait(timeout=5.0), "child never acquired the lock"

    os.environ["VCFCLICK_HOME"] = str(vcfclick_home)
    os.environ["VCFCLICK_DB_NAME"] = "smoke"
    from storage import ingest_id_lock

    started = time.time()
    with ingest_id_lock("batch_a"):
        held_for = time.time() - started

    child.join(timeout=5.0)
    # Parent must have waited at least until the child slept-and-
    # released. Give scheduling slack but ensure it blocked at all.
    assert held_for >= 0.3, (
        f"parent acquired lock too fast: waited only {held_for:.2f}s "
        f"while child was supposed to hold for {hold}s — the lock is not "
        f"actually blocking concurrent acquisition"
    )


def test_different_ingest_ids_do_not_block_each_other(vcfclick_home):
    """The lock is per-ingest_id. Two ingests under different
    ingest_ids must run concurrently — assert acquiring `batch_a`
    while another process holds `batch_b` doesn't block."""
    subprocess.run(
        [VCFCLICK_BIN, "db", "create", "smoke"],
        env={**os.environ, "VCFCLICK_HOME": str(vcfclick_home)},
        check=True,
        capture_output=True,
    )

    ctx = multiprocessing.get_context("spawn")
    acquired = ctx.Event()
    child = ctx.Process(
        target=_hold_lock, args=(str(vcfclick_home), "batch_b", 1.0, acquired)
    )
    child.start()
    assert acquired.wait(timeout=5.0), "child never acquired the lock"

    os.environ["VCFCLICK_HOME"] = str(vcfclick_home)
    os.environ["VCFCLICK_DB_NAME"] = "smoke"
    from storage import ingest_id_lock

    started = time.time()
    with ingest_id_lock("batch_a"):
        waited = time.time() - started
    child.join(timeout=5.0)

    # batch_a and batch_b are independent locks — acquisition should
    # be near-instant despite child holding batch_b.
    assert waited < 0.3, f"different-id lock blocked unexpectedly ({waited:.2f}s)"
