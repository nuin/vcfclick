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
from pathlib import Path

from storage import get_session


TABLES = ["variants", "genotypes", "samples", "ingestions"]


def export_table(table: str, out_path: Path, where: str | None = None) -> None:
    """Materialise a table (or filtered slice) as a Parquet file."""
    if table not in TABLES:
        raise ValueError(f"Unknown table {table}; expected one of {TABLES}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sess = get_session()
    sql = f"SELECT * FROM {table}"
    if where:
        sql += f" WHERE {where}"
    sql += f" INTO OUTFILE '{out_path}' FORMAT Parquet"

    sess.query(sql)
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
