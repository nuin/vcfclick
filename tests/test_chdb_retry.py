"""Unit tests for the chDB session-init retry wrapper.

The retry path is impossible to exercise reliably against real chDB —
the race is statistical. These tests mock `chdb.session.Session` to
return controlled failure sequences and assert the wrapper's behaviour:

  - Race-marker errors trigger retry (up to 3 retries by default).
  - Non-race errors propagate on the first attempt with no retry.
  - When retries are exhausted, the last race error is re-raised.
  - The backoff schedule is what the module documents.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest


def test_non_race_error_propagates_without_retry():
    """A garden-variety exception (e.g., disk full) must not get
    retried — that would mask real failures behind the same wrapper."""
    import chdb.session as cdb_session

    import storage._chdb as sdb

    call_count = {"n": 0}

    def fake_session(path):
        call_count["n"] += 1
        raise RuntimeError("Disk full, no retry please")

    with patch.object(cdb_session, "Session", side_effect=fake_session):
        with pytest.raises(RuntimeError, match="Disk full"):
            sdb.open_session("/tmp/whatever")

    assert call_count["n"] == 1, "non-race error must NOT be retried"


def test_race_error_retried_until_success():
    """A race-marker error twice in a row, then success: wrapper
    returns the success-attempt session."""
    import chdb.session as cdb_session

    import storage._chdb as sdb

    calls = {"n": 0}

    def fake_session(path):
        calls["n"] += 1
        if calls["n"] <= 2:
            raise RuntimeError(
                "Code: 36. DB::Exception: Error initializing EmbeddedServer: "
                "Code: 695. … recursive_mutex lock failed: Invalid argument. "
                "(ASYNC_LOAD_WAIT_FAILED) (BAD_ARGUMENTS)"
            )
        return object()  # opaque success-sentinel

    # Zero out the backoff delays so the test doesn't sleep 1.3s.
    with patch.object(sdb, "_CHDB_RETRY_DELAYS_S", (0, 0, 0)):
        with patch.object(cdb_session, "Session", side_effect=fake_session):
            result = sdb.open_session("/tmp/test")

    assert calls["n"] == 3, "should have retried twice before success"
    assert result is not None


def test_race_error_exhausts_retries_then_raises():
    """If every attempt loses the race, the last exception surfaces
    with the original message so the caller can see what happened."""
    import chdb.session as cdb_session

    import storage._chdb as sdb

    def fake_session(path):
        raise RuntimeError(
            "Code: 36. DB::Exception: Error initializing EmbeddedServer: "
            "recursive_mutex lock failed: Invalid argument. (BAD_ARGUMENTS)"
        )

    with patch.object(sdb, "_CHDB_RETRY_DELAYS_S", (0, 0, 0)):
        with patch.object(cdb_session, "Session", side_effect=fake_session):
            with pytest.raises(RuntimeError, match="recursive_mutex"):
                sdb.open_session("/tmp/test")


def test_race_markers_cover_observed_chdb_error_patterns():
    """Sanity-pin the marker list against the exact phrasings we've
    seen on CI so a chDB version bump silently dropping one marker
    doesn't disable the retry."""
    import storage._chdb as sdb

    expected_markers = {
        "BAD_ARGUMENTS",
        "recursive_mutex lock failed",
        "ASYNC_LOAD_WAIT_FAILED",
    }
    assert set(sdb._CHDB_RACE_MARKERS) == expected_markers


def test_backoff_is_capped_at_reasonable_total():
    """The retry loop shouldn't be able to burn an open-ended amount
    of time — keep the cumulative budget tight enough that a CLI
    invocation hitting persistent failure still surfaces it quickly."""
    import storage._chdb as sdb

    total = sum(sdb._CHDB_RETRY_DELAYS_S)
    assert total < 2.0, f"retry budget {total}s is too long for interactive use"
    # And at least non-trivial — instant retries don't give the
    # async-load scheduler time to settle.
    assert sdb._CHDB_RETRY_DELAYS_S[0] >= 0.05
