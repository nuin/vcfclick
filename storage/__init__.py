"""Public storage API.

Functions are re-exported eagerly because they don't change between
calls. The two PATH attributes (`DB_ROOT`, `VCFCLICK_HOME`) are
forwarded via `__getattr__` instead — they depend on
`os.environ["VCFCLICK_HOME"]` at lookup time, and importing them
eagerly here would snapshot the env at import time and silently
ignore any subsequent change. Codex round 9 caught this exact gap.
"""

from storage.db import (
    apply_schema,
    backend,
    count_expr,
    db_disk_size,
    db_path,
    delete_where_sql,
    drop_db,
    get_session,
    ingest_id_lock,
    insert_via_parquet,
    list_dbs,
    parquet_file_expr,
    rollback_ingest,
    schema_dir_for_backend,
    sql_quote_str,
    table_exists,
    validate_ingest_id,
)


def __getattr__(name: str):
    """Lazy forward for env-dependent module attributes."""
    if name in ("DB_ROOT", "VCFCLICK_HOME"):
        from storage import db as _db

        return getattr(_db, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "get_session",
    "apply_schema",
    "insert_via_parquet",
    "rollback_ingest",
    "sql_quote_str",
    "validate_ingest_id",
    "ingest_id_lock",
    "db_path",
    "list_dbs",
    "drop_db",
    "db_disk_size",
    "backend",
    "schema_dir_for_backend",
    "delete_where_sql",
    "parquet_file_expr",
    "count_expr",
    "table_exists",
    "DB_ROOT",
    "VCFCLICK_HOME",
]
