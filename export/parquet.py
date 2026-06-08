"""Export chDB tables to Parquet for interop with the wider data-lake
ecosystem (DuckDB, Snowflake, BigQuery, Spark, Iceberg-on-Parquet).

vcfclick is a complete embedded system in itself, but the Parquet
export means users are never locked in — every byte of their data
lands as a portable open-format file with a single SQL statement.

Usage:
    python -m export.parquet variants /path/out.parquet
    python -m export.parquet variants /path/out.parquet \\
        --where "chrom='chr17' AND pos BETWEEN 43044295 AND 43125483"

Multiple tables in one call:
    python -m export.parquet --all /path/output_dir/
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

from storage import get_session

log = logging.getLogger(__name__)


TABLES = ["variants", "genotypes", "samples", "ingestions"]


def export_table(table: str, out_path: Path, where: str | None = None) -> None:
    """Materialise a table (or filtered slice) as a Parquet file.

    Note on safety:
      - `table` is checked against an allowlist (above).
      - `out_path` is NEVER interpolated into the SQL. chDB writes to
        a tempfile under our control, then we atomic-move to out_path.
      - `where` IS interpolated. It is treated as a trusted SQL
        fragment passed in by the local operator (via `vcfclick db dump
        --where ...`). Do NOT expose this entry point to untrusted
        callers (HTTP, etc.) without first switching to a parameterised
        query builder.
    """
    from storage import backend

    if table not in TABLES:
        raise ValueError(f"Unknown table {table}; expected one of {TABLES}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sess = get_session()
    with tempfile.NamedTemporaryFile(
        suffix=".parquet", delete=False, dir=out_path.parent
    ) as f:
        staging = Path(f.name)
    try:
        select_sql = f"SELECT * FROM {table}"
        if where:
            select_sql += f" WHERE {where}"
        # staging is from tempfile (path is controlled, no quotes).
        if backend() == "duckdb":
            # COPY ... TO uses single-quoted file path; TO overwrites by default.
            sql = f"COPY ({select_sql}) TO '{staging}' (FORMAT 'parquet')"
        else:
            # chDB: TRUNCATE lets it overwrite the empty file tempfile
            # created for us — otherwise FILE_ALREADY_EXISTS.
            sql = f"{select_sql} INTO OUTFILE '{staging}' TRUNCATE FORMAT Parquet"
        sess.query(sql)
        staging.replace(out_path)
    except Exception:
        staging.unlink(missing_ok=True)
        raise

    size = out_path.stat().st_size
    log.info("[export] %s → %s (%s bytes)", table, out_path, f"{size:,}")


def export_all(out_dir: Path) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for t in TABLES:
        export_table(t, out_dir / f"{t}.parquet")


# Library module — invoke via `vcfclick db dump <name>`.
# The public CLI lives in cli/main.py.
