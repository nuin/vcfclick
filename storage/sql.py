"""SQL fragments for the supported variant-store dialects."""

import re

_TABLE_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def sql_quote_str(s: str) -> str:
    """ClickHouse string literal — backslashes AND quotes both escaped.

    ClickHouse / chDB recognises two escape forms inside string
    literals: `''` is a single quote, and `\\'` is also a single
    quote. Plain quote-doubling (the portable-SQL form) is therefore
    NOT enough on its own: a payload like `\\'; DROP TABLE …` would,
    after quote-doubling, become `\\''; DROP TABLE …` — the engine
    reads the first quote as backslash-escaped, the second as the
    string terminator, and arbitrary SQL after it.

    Order of operations matters: escape backslashes FIRST (each `\\`
    → `\\\\`), then quotes (each `'` → `''`). Doing it the other way
    would double the just-added backslashes again.
    """
    escaped = s.replace("\\", "\\\\").replace("'", "''")
    return "'" + escaped + "'"


def delete_where_sql(table: str, where: str) -> str:
    """Return the backend-specific synchronous DELETE statement."""

    from storage.db import backend

    if backend() == "duckdb":
        return f"DELETE FROM {table} WHERE {where}"
    return f"ALTER TABLE {table} DELETE WHERE {where} SETTINGS mutations_sync = 2"


def parquet_file_expr(path: str) -> str:
    """Return the backend-specific SQL fragment for reading a Parquet file.

      * chDB:    `file('/abs/path.parquet', 'Parquet')`
      * DuckDB:  `read_parquet('/abs/path.parquet')`

    The raw path is escaped and quoted with `sql_quote_str` here.
    """

    from storage.db import backend

    if backend() == "duckdb":
        return f"read_parquet({sql_quote_str(path)})"
    return f"file({sql_quote_str(path)}, 'Parquet')"


def count_expr() -> str:
    """ClickHouse permits `count()`; DuckDB requires `count(*)`."""

    from storage.db import backend

    return "count(*)" if backend() == "duckdb" else "count()"


def typed_columns_sql(table: str) -> str:
    """SQL listing a table's columns as (name, type, is_flag).

    `is_flag` marks a column whose "populated" test is `!= 0` rather than
    `IS NOT NULL`: in ClickHouse those are the non-`Nullable` flag columns;
    in DuckDB the same columns are nullable but declared `DEFAULT 0`, so
    `IS NOT NULL` would count every row.
    """

    from storage.db import backend

    if not _TABLE_NAME_RE.match(table):
        raise ValueError(f"Unsafe table name: {table!r}")
    if backend() == "duckdb":
        return (
            "SELECT column_name, data_type, "
            "CASE WHEN column_default = '0' THEN 1 ELSE 0 END "
            "FROM information_schema.columns "
            f"WHERE table_name = '{table}' ORDER BY ordinal_position"
        )
    # '%Nullable%', not 'Nullable%': ClickHouse wraps types, so a nullable
    # column can read `LowCardinality(Nullable(String))`. A prefix match would
    # misclassify those as flags and emit `!= 0` against a String.
    return (
        "SELECT name, type, CASE WHEN type LIKE '%Nullable%' THEN 0 ELSE 1 END "
        "FROM system.columns "
        f"WHERE table = '{table}' AND database = currentDatabase() "
        "ORDER BY position"
    )


def populated_expr(col: str, is_flag: bool) -> str:
    """Per-column "is populated" aggregate. ClickHouse has `countIf`; DuckDB
    uses the SQL-standard `count(*) FILTER (WHERE ...)`."""

    from storage.db import backend

    test = f'"{col}" != 0' if is_flag else f'"{col}" IS NOT NULL'
    if backend() == "duckdb":
        return f'count(*) FILTER (WHERE {test}) AS "{col}"'
    test = test.replace('"', "`")
    return f"countIf({test}) AS `{col}`"


def map_keys_from(table: str, map_col: str) -> str:
    """A FROM-clause exposing one row per Map key as `k`.

    ClickHouse flattens a Map with `ARRAY JOIN mapKeys(m)`; DuckDB uses
    `unnest(map_keys(m))` in a subquery.
    """

    from storage.db import backend

    if not _TABLE_NAME_RE.match(table) or not _TABLE_NAME_RE.match(map_col):
        raise ValueError(f"Unsafe identifier: {table!r}/{map_col!r}")
    if backend() == "duckdb":
        return f"(SELECT unnest(map_keys({map_col})) AS k FROM {table})"
    return f"{table} ARRAY JOIN mapKeys({map_col}) AS k"
