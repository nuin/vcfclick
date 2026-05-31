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

import argparse
import tempfile
from pathlib import Path

from storage import get_session


TABLES = ["variants", "genotypes", "samples", "ingestions"]


def export_table(table: str, out_path: Path, where: str | None = None) -> None:
    """Materialise a table (or filtered slice) as a Parquet file.

    Note on safety:
      - `table` is checked against an allowlist (above).
      - `out_path` is NEVER interpolated into the SQL. chDB writes to
        a tempfile under our control, then we atomic-move to out_path.
      - `where` IS interpolated. It's a CLI argument provided by the
        operator and treated as trusted SQL fragment input. Do NOT
        expose this entry point to untrusted callers (HTTP, etc.)
        without first switching to a parameterised query builder.
    """
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
        sql = f"SELECT * FROM {table}"
        if where:
            sql += f" WHERE {where}"
        # staging is from tempfile (path is controlled, no quotes).
        # TRUNCATE lets chDB overwrite the empty file tempfile created
        # for us — otherwise FILE_ALREADY_EXISTS.
        sql += f" INTO OUTFILE '{staging}' TRUNCATE FORMAT Parquet"
        sess.query(sql)
        staging.replace(out_path)
    except Exception:
        staging.unlink(missing_ok=True)
        raise

    size = out_path.stat().st_size
    print(f"[export] {table} → {out_path} ({size:,} bytes)")


def export_all(out_dir: Path) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for t in TABLES:
        export_table(t, out_dir / f"{t}.parquet")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("table", nargs="?", help="Table name (variants/genotypes/samples/ingestions)")
    ap.add_argument("out_path", nargs="?", help="Output Parquet path")
    ap.add_argument("--where", default=None, help="Optional WHERE clause")
    ap.add_argument(
        "--all", action="store_true",
        help="Export all tables to <out_path>/<table>.parquet",
    )
    args = ap.parse_args()

    if args.all:
        if not args.table:
            ap.error("--all requires a directory argument")
        export_all(Path(args.table))
    else:
        if not args.table or not args.out_path:
            ap.error("table and out_path required (or use --all <dir>)")
        export_table(args.table, Path(args.out_path), args.where)


if __name__ == "__main__":
    main()
