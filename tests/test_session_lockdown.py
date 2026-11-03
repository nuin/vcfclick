"""`storage.get_readonly_session`: the engine-level lock behind the SQL
guard for the web UI and MCP `run_sql`.

DuckDB: `enable_external_access = false` and `lock_configuration = true`
on the shared connection, so no statement on that session can open a
file, fetch a URL or change a setting again. chDB: `readonly = 1` on the
session, which refuses writes, table functions and SETTINGS clauses.
The DuckDB check runs in-process; the chDB check runs in a subprocess
(one embedded server per process)."""

from __future__ import annotations

import pytest

from tests.conftest import run_cli


def test_duckdb_readonly_session_blocks_external_access(
    vcfclick_home, tmp_path, monkeypatch
):
    monkeypatch.setenv("VCFCLICK_BACKEND", "duckdb")
    import storage.db as sdb
    from storage import get_readonly_session, get_session

    secret = tmp_path / "secret.txt"
    secret.write_text("secretline\n")
    sess = get_readonly_session("lockdb")
    try:
        # The cohort itself stays readable and writable through the storage
        # layer (ingest, schema upgrades): the lock is about leaving it.
        sess.query("CREATE TABLE t (x INTEGER)")
        sess.query("INSERT INTO t VALUES (1)")
        assert sess.query("SELECT x FROM t", "CSV").bytes().decode().strip() == "1"
        with pytest.raises(Exception, match="disabled"):
            sess.query(f"SELECT content FROM read_text('{secret}')", "CSV")
        with pytest.raises(Exception, match="locked"):
            sess.query("SET enable_external_access = true")
        with pytest.raises(Exception, match="locked"):
            sess.query("SET lock_configuration = false")
        # Idempotent, and the same object every other caller gets.
        assert get_readonly_session("lockdb") is sess
        assert get_session("lockdb") is sess
    finally:
        sdb._sessions.pop("duckdb::lockdb").close()


@pytest.mark.parametrize("backend", ["chdb"])
def test_chdb_readonly_session_refuses_files_writes_and_settings(
    vcfclick_home, tmp_path, run_python, backend
):
    run_cli(vcfclick_home, backend, "db", "create", "lockdb")
    secret = tmp_path / "secret.txt"
    secret.write_text("secretline\n")
    code = f"""
from storage import get_readonly_session
sess = get_readonly_session("lockdb")
print("plain", sess.query("SELECT 1", "CSV").bytes().decode().strip())
for label, sql in [
    ("file", "SELECT * FROM file({str(secret)!r}, 'LineAsString')"),
    ("url", "SELECT * FROM url('http://127.0.0.1:9/x', CSV)"),
    ("insert", "INSERT INTO ingestions (ingest_id) VALUES ('x')"),
    ("drop", "DROP TABLE ingestions"),
    ("set", "SET readonly = 0"),
    ("settings", "SELECT 1 SETTINGS max_threads = 1"),
]:
    try:
        sess.query(sql, "CSV")
        print(label, "RAN")
    except Exception as exc:
        print(label, "READONLY" if "READONLY" in str(exc) else "OTHER " + str(exc)[:80])
print("same", get_readonly_session("lockdb") is sess)
"""
    out = run_python(vcfclick_home, code, VCFCLICK_BACKEND=backend)
    lines = dict(line.split(" ", 1) for line in out.strip().splitlines())
    assert lines["plain"] == "1"
    for label in ("file", "url", "insert", "drop", "set", "settings"):
        assert lines[label] == "READONLY", (label, lines[label])
    assert lines["same"] == "True"


def test_web_and_mcp_run_sql_use_the_locked_session(monkeypatch):
    """Both surfaces fetch their session through get_readonly_session."""
    import vcfclick_web.app as web
    from vcfclick_mcp import server

    calls: list[str] = []

    class FakeResult:
        def bytes(self):
            return b'{"meta": [{"name": "n"}], "data": [[1]]}'

    class FakeSession:
        def query(self, sql, fmt=None):
            return FakeResult()

    def fake_readonly_session(name=None):
        calls.append("readonly")
        return FakeSession()

    def plain_session(name=None):
        raise AssertionError("run_sql must use get_readonly_session")

    monkeypatch.setattr(web, "get_readonly_session", fake_readonly_session)
    monkeypatch.setattr(web, "get_session", plain_session)
    monkeypatch.setattr(server, "get_readonly_session", fake_readonly_session)
    monkeypatch.setattr(server, "get_session", plain_session)

    assert web._run_sql("SELECT 1")["rows"] == [[1]]
    assert server.run_sql("SELECT 1")["rows"] == [[1]]
    assert calls == ["readonly", "readonly"]
