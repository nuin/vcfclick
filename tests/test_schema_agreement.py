"""Drift guard: every pyarrow schema and SQL DDL must agree column-by-
column.

The Arrow schemas in `ingest._arrow` drive Parquet column order, the
INSERT column lists, and the row-builder column lists in
`ingest.vcf_load`. The SQL DDL in `schema/*.sql` drives chDB's table
shape. If those two ever drift, chDB will silently misalign columns
on any future release where Parquet imports become positional (rather
than name-matched as they are today). The unit tests still pass in
that scenario — they only check specific values, not column-by-column
agreement.

This module is the explicit drift guard. Any column added on the SQL
side must also be added in the same position on the Arrow side, and
vice versa.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

SCHEMA_DIR = Path(__file__).resolve().parent.parent / "schema"


def _columns_from_sql(sql_path: Path) -> list[str]:
    """Parse a CREATE TABLE statement and return the column names in
    declaration order. Skips the `ingested_at` housekeeping column,
    which chDB fills from the DEFAULT and the ingester does not write."""
    text = sql_path.read_text()

    # CREATE TABLE <name> ( <body> ) ENGINE = ...
    m = re.search(
        r"CREATE\s+TABLE\s+\w+\s*\((.*?)\)\s*ENGINE",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    if not m:
        raise AssertionError(f"no CREATE TABLE found in {sql_path}")

    columns: list[str] = []
    for raw_line in m.group(1).splitlines():
        line = raw_line.strip()
        # Skip blank lines, full-line comments, and trailing-comma artefacts.
        if not line or line.startswith("--"):
            continue
        # Strip inline comments before tokenising.
        if "--" in line:
            line = line.split("--", 1)[0].strip()
        if not line:
            continue
        # `name <type>...,` — pull the first whitespace-delimited token.
        first = line.split()[0].rstrip(",")
        # Filter anything that isn't a plain identifier (catches stray
        # keywords if a future schema rewrite changes the body format).
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", first):
            continue
        # `ingested_at` is filled by DEFAULT; ingester does not write it.
        if first == "ingested_at":
            continue
        columns.append(first)
    return columns


@pytest.mark.parametrize(
    "table,sql_file,arrow_attr",
    [
        ("variants", "01_variants.sql", "VARIANTS_ARROW_SCHEMA"),
        ("genotypes", "02_genotypes.sql", "GENOTYPES_ARROW_SCHEMA"),
    ],
)
def test_arrow_and_sql_columns_agree_in_order(table, sql_file, arrow_attr):
    """For each ingest-target table, the Arrow schema must list the
    same columns in the same order as the SQL DDL declares them.

    Failure mode this catches: someone adds a column on the SQL side
    and forgets to update the Arrow schema (or vice versa). Today's
    chDB does name-based Parquet import so the row data still lands
    in the right column, but the Parquet column ORDER would then
    diverge from the DDL ORDER and a future chDB change could
    silently misalign every typed INFO/FORMAT column.
    """
    import ingest._arrow as arrow_mod

    sql_columns = _columns_from_sql(SCHEMA_DIR / sql_file)
    arrow_columns = [f.name for f in getattr(arrow_mod, arrow_attr)]

    # The Arrow schema is the one we control and re-order; show the
    # exact diff if they disagree.
    if sql_columns != arrow_columns:
        first_diff = next(
            (
                i
                for i, (a, b) in enumerate(
                    zip(sql_columns, arrow_columns, strict=False)
                )
                if a != b
            ),
            min(len(sql_columns), len(arrow_columns)),
        )
        raise AssertionError(
            f"{table}: Arrow schema diverges from SQL DDL at position "
            f"{first_diff}.\n  SQL:   {sql_columns[first_diff : first_diff + 3]} ...\n"
            f"  Arrow: {arrow_columns[first_diff : first_diff + 3]} ...\n"
            f"  full SQL list:   {sql_columns}\n"
            f"  full Arrow list: {arrow_columns}"
        )


def test_samples_arrow_matches_sql():
    """samples.sql column order vs SAMPLES_ARROW_SCHEMA. Same contract."""
    from ingest._arrow import SAMPLES_ARROW_SCHEMA

    sql_columns = _columns_from_sql(SCHEMA_DIR / "03_samples.sql")
    arrow_columns = [f.name for f in SAMPLES_ARROW_SCHEMA]
    assert sql_columns == arrow_columns, (
        f"samples Arrow diverges from SQL DDL.\n"
        f"  SQL:   {sql_columns}\n  Arrow: {arrow_columns}"
    )
