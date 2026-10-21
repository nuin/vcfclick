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


# ───────────────────────── shared query fragments ─────────────────────────
#
# Small builders every genotype query repeats. They exist so a lesson learned
# the hard way (a LEFT JOIN that reads an absent parent as 0/0 on chDB, an
# unsigned add that overflows on DuckDB) is encoded once, with a test on both
# backends, instead of being rediscovered by each feature.

_IDENT_PATH_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(\.[A-Za-z_][A-Za-z0-9_]*)?$")


def _ident(name: str) -> str:
    if not _IDENT_PATH_RE.match(name):
        raise ValueError(f"Unsafe identifier: {name!r}")
    return name


def str_list(values) -> str:
    """`'a', 'b'` for an `IN (...)` list; every value through sql_quote_str."""
    return ", ".join(sql_quote_str(v) for v in values)


def autosome_filter(col: str = "chrom") -> str:
    """`chrom IN (<every autosome spelling>)`. The names come from
    storage.contigs, the one definition of the contig alias sets."""
    from storage.contigs import autosome_names

    return f"{_ident(col)} IN ({str_list(autosome_names())})"


def sample_join(
    alias: str, sample_id: str, *, base: str = "g", table: str = "genotypes"
) -> str:
    """`INNER JOIN genotypes <alias> ...` bringing in SAMPLE_ID's call at the
    same (ingest_id, chrom, pos, ref, alt) as `base`.

    Always an INNER join, on purpose. A LEFT JOIN would hand back, for a
    sample with no stored row at that site, a default instead of NULL on
    chDB (an unmatched row reads gt = 0), so absence would look like a
    proven 0/0. With an INNER join a site where the sample has no row simply
    does not match. Use `genotype_pivot_sql` when absence must be visible.
    """
    a, b = _ident(alias), _ident(base)
    key = " AND ".join(
        f"{a}.{c} = {b}.{c}" for c in ("ingest_id", "chrom", "pos", "ref", "alt")
    )
    return (
        f"INNER JOIN {_ident(table)} {a} ON {key} AND "
        f"{a}.sample_id = {sql_quote_str(sample_id)}"
    )


def to_int64(expr: str) -> str:
    """NULL-safe widening to a 64-bit integer. chDB's CAST(x AS BIGINT)
    throws on a NULL in a Nullable column; toInt64 propagates the NULL."""
    from storage.db import backend

    if backend() == "duckdb":
        return f"CAST({expr} AS BIGINT)"
    return f"toInt64({expr})"


def to_double(expr: str) -> str:
    from storage.db import backend

    if backend() == "duckdb":
        return f"CAST({expr} AS DOUBLE)"
    return f"toFloat64({expr})"


def allele_fraction_sql(alt_col: str, ref_col: str) -> str:
    """alt / (ref + alt) as a double, NULL when either count is NULL.

    The read counts are 16-bit unsigned columns. DuckDB adds two USMALLINTs
    in USMALLINT and raises on overflow (65535 + 1), so both are widened to
    BIGINT before they are added. A zero denominator is the caller's to
    guard (`ref + alt > 0`): both engines return nan or inf, not NULL.
    """
    alt, ref = to_int64(_ident(alt_col)), to_int64(_ident(ref_col))
    return f"({to_double(alt)} / {to_double(f'{alt} + {ref}')})"


def genotype_pivot_sql(
    ingest_id: str, samples: dict[str, str], fields: tuple[str, ...] = ("gt",)
) -> str:
    """A subquery with one row per site and, per ROLE in `samples`
    ({"p": "CHILD", "f": "DAD"}), one column per field named `<role>_<field>`
    (`p_gt`, `f_gt`, ...), plus ingest_id, chrom, pos, ref, alt.

    A sample with no stored row at the site reads NULL, never 0: the value
    is `max(CASE WHEN sample_id = ... THEN field END)`, and an aggregate
    over no matching row is NULL on both backends. That is the null-safe
    way to see an absent parent; the equivalent LEFT JOIN reads gt = 0 on
    chDB. Only sites where at least one listed sample has a row appear.
    Wrap it as `FROM (<this>) x` and alias the columns you select.
    """
    if not samples:
        raise ValueError("genotype_pivot_sql needs at least one sample")
    cols = []
    for role, sample in samples.items():
        for field in fields:
            cols.append(
                f"max(CASE WHEN sample_id = {sql_quote_str(sample)} "
                f"THEN {_ident(field)} END) AS {_ident(role)}_{_ident(field)}"
            )
    return (
        "SELECT ingest_id, chrom, pos, ref, alt, "
        + ", ".join(cols)
        + f" FROM genotypes WHERE ingest_id = {sql_quote_str(ingest_id)} "
        f"AND sample_id IN ({str_list(samples.values())}) "
        "GROUP BY ingest_id, chrom, pos, ref, alt"
    )
