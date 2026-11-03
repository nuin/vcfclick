"""Reference-data loaders (GENCODE genes, ClinVar, etc.) → DuckDB."""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager

import duckdb


@contextmanager
def store_transaction() -> Iterator[duckdb.DuckDBPyConnection]:
    """The annotation store open for writing inside one transaction: committed
    when the block ends, rolled back when it raises, closed either way. A
    loader parses its input before it enters the block, so a bad file never
    reaches the store and the store is held only for the swap."""
    from annotations.db import get_connection

    conn = get_connection()
    try:
        conn.execute("BEGIN TRANSACTION")
        try:
            yield conn
        except BaseException:
            conn.execute("ROLLBACK")
            raise
        conn.execute("COMMIT")
    finally:
        conn.close()
