"""Storage session cache lifecycle tests."""

from __future__ import annotations

import logging


def test_drop_db_removes_backend_qualified_session_cache(vcfclick_home, monkeypatch):
    """drop_db must remove the same backend-qualified key get_session writes."""
    monkeypatch.setenv("VCFCLICK_BACKEND", "duckdb")

    import storage.db as sdb

    sdb._sessions.clear()
    sdb.get_session("cachetest")
    assert "duckdb::cachetest" in sdb._sessions

    sdb.drop_db("cachetest")

    assert "duckdb::cachetest" not in sdb._sessions


def test_duckdb_close_logs_close_failures(tmp_path, caplog):
    """DuckDB close failures are intentionally non-fatal but must be visible."""
    from storage._duckdb import DuckDBSession

    class BrokenConnection:
        def close(self) -> None:
            raise RuntimeError("close broke")

    session = DuckDBSession(tmp_path / "db")
    session._conn = BrokenConnection()

    with caplog.at_level(logging.WARNING):
        session.close()

    assert "DuckDB connection close failed" in caplog.text
