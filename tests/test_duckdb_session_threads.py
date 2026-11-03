"""DuckDBSession under concurrent callers.

`vcfclick web` serves its sync endpoints from a threadpool, and every
endpoint reads the one cached session, so two requests run on the same
DuckDB connection at once. A DuckDBPyConnection is not thread-safe: with
a shared cursor, one thread's fetchall() returns another thread's rows,
or nothing. Each query must run on a cursor of its own."""

from __future__ import annotations

import json
import threading

from storage._duckdb import DuckDBSession

THREADS = 4
QUERIES = 200


def test_concurrent_queries_each_get_their_own_rows(tmp_path):
    sess = DuckDBSession(tmp_path / "db")
    sess.query("CREATE TABLE t AS SELECT range AS x FROM range(1000)")
    wrong: list[tuple[int, int, object]] = []
    start = threading.Barrier(THREADS)

    def worker(tid: int) -> None:
        start.wait()
        for i in range(QUERIES):
            tag = tid * QUERIES + i
            raw = sess.query(
                f"SELECT {tag} AS tag, (SELECT count(*) FROM t) AS n", "JSONCompact"
            ).bytes()
            try:
                rows = json.loads(raw)["data"]
            except Exception:
                rows = raw
            if rows != [[tag, 1000]]:
                wrong.append((tid, i, rows))

    threads = [threading.Thread(target=worker, args=(t,)) for t in range(THREADS)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    sess.close()
    assert not wrong, f"{len(wrong)} of {THREADS * QUERIES} results wrong: {wrong[:5]}"


def test_concurrent_arrow_queries_each_get_their_own_rows(tmp_path):
    sess = DuckDBSession(tmp_path / "db")
    wrong: list[tuple[int, int, object]] = []
    start = threading.Barrier(THREADS)

    def worker(tid: int) -> None:
        start.wait()
        for i in range(QUERIES):
            tag = tid * QUERIES + i
            try:
                got = sess.query_arrow(f"SELECT {tag} AS tag").column("tag").to_pylist()
            except Exception as exc:
                got = repr(exc)
            if got != [tag]:
                wrong.append((tid, i, got))

    threads = [threading.Thread(target=worker, args=(t,)) for t in range(THREADS)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    sess.close()
    assert not wrong, f"{len(wrong)} of {THREADS * QUERIES} results wrong: {wrong[:5]}"


def test_get_session_opens_one_session_under_concurrent_first_calls(
    vcfclick_home, monkeypatch
):
    monkeypatch.setenv("VCFCLICK_BACKEND", "duckdb")
    import storage.db as sdb

    opened: list[object] = []
    real_open = sdb._open_session

    def slow_open(path):
        s = real_open(path)
        opened.append(s)
        return s

    monkeypatch.setattr(sdb, "_open_session", slow_open)
    got: list[object] = []
    start = threading.Barrier(THREADS)

    def worker() -> None:
        start.wait()
        got.append(sdb.get_session("racedb"))

    threads = [threading.Thread(target=worker) for _ in range(THREADS)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    try:
        assert len(opened) == 1
        assert all(s is opened[0] for s in got)
    finally:
        for s in opened:
            s.close()
        sdb._sessions.pop("duckdb::racedb", None)
