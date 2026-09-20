"""Schema-population stats for `vcfclick db`."""

from __future__ import annotations

import click

from cli.main import _set_db, db

_STATS_SKIP = {
    "variants": {
        "ingest_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "vcf_id",
        "qual",
        "filter",
        "info_extra",
        "ingested_at",
    },
    "genotypes": {
        "ingest_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "sample_id",
        "gt",
        "phased",
        "format_extra",
        "ingested_at",
    },
}


def _list_typed_columns(sess, table: str) -> list[tuple[str, bool]]:
    """Typed, non-structural columns for `table` as (name, is_flag)."""
    from storage import typed_columns_sql

    out = sess.query(typed_columns_sql(table), "TSV").bytes().decode().strip()
    skip = _STATS_SKIP.get(table, set())
    cols: list[tuple[str, bool]] = []
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name, _type, is_flag = parts[0], parts[1], parts[2]
        if name not in skip:
            cols.append((name, is_flag.strip() in ("1", "true", "True")))
    return cols


def _query_population(sess, table: str, columns: list[tuple[str, bool]]) -> dict:
    """Return `{column_name: populated_row_count}` for typed columns."""
    from storage import populated_expr

    if not columns:
        return {}
    exprs = ", ".join(populated_expr(n, flag) for n, flag in columns)
    out = sess.query(f"SELECT {exprs} FROM {table}", "TSV").bytes().decode().strip()
    values = [int(v) for v in out.split("\t")]
    return dict(zip([n for n, _ in columns], values))


def _query_map_keys(
    sess, table: str, map_col: str, top: int
) -> tuple[list[tuple[str, int]], int]:
    """Return top Map keys and total distinct key count."""
    from storage import count_expr, map_keys_from

    src = map_keys_from(table, map_col)
    n_distinct = int(
        sess.query(f"SELECT count(DISTINCT k) FROM {src}", "TSV")
        .bytes()
        .decode()
        .strip()
        or "0"
    )
    if n_distinct == 0:
        return [], 0
    out = (
        sess.query(
            f"SELECT k, {count_expr()} AS n FROM {src} "
            f"GROUP BY k ORDER BY n DESC, k LIMIT {int(top)}",
            "TSV",
        )
        .bytes()
        .decode()
        .strip()
    )
    rows: list[tuple[str, int]] = []
    for line in out.splitlines():
        k, n = line.split("\t", 1)
        rows.append((k, int(n)))
    return rows, n_distinct


def _pct(num: int, denom: int) -> str:
    if denom == 0:
        return "  0.0%"
    return f"{100.0 * num / denom:>5.1f}%"


def _query_rows(sess, sql: str) -> list[list[str]]:
    out = sess.query(sql, "TSV").bytes().decode().strip()
    return [line.split("\t") for line in out.splitlines() if line]


def _stats_payload(sess, top: int) -> dict:
    from storage import backend, count_expr

    cnt = count_expr()
    counts = {
        t: int(
            sess.query(f"SELECT {cnt} FROM {t}", "TSV").bytes().decode().strip() or "0"
        )
        for t in ("variants", "genotypes", "samples", "ingestions")
    }
    variants_cols = _list_typed_columns(sess, "variants")
    genotypes_cols = _list_typed_columns(sess, "genotypes")
    # ClickHouse counts distinct tuples directly; DuckDB needs the pair
    # concatenated into one expression.
    distinct_pair = (
        "count(DISTINCT (ingest_id, sample_id))"
        if backend() != "duckdb"
        else "count(DISTINCT (ingest_id || '' || sample_id))"
    )
    return {
        "counts": counts,
        "cohorts": _query_rows(
            sess,
            f"SELECT cohort, {distinct_pair} AS n "
            "FROM samples GROUP BY cohort ORDER BY n DESC, cohort",
        ),
        "contigs": _query_rows(
            sess,
            f"SELECT chrom, {cnt} AS n FROM variants "
            "GROUP BY chrom ORDER BY n DESC, chrom",
        ),
        "variants_pop": _query_population(sess, "variants", variants_cols),
        "info_extra": _query_map_keys(sess, "variants", "info_extra", top),
        "genotypes_pop": _query_population(sess, "genotypes", genotypes_cols),
        "format_extra": _query_map_keys(sess, "genotypes", "format_extra", top),
    }


def _render_population(title: str, total: int, population: dict[str, int]) -> None:
    click.echo()
    click.echo(f"{title} (of {total:,} rows):")
    for col, n in sorted(population.items(), key=lambda kv: (-kv[1], kv[0])):
        click.echo(f"  {col:<32} {n:>10,}  ({_pct(n, total)})")


def _render_map_keys(
    title: str,
    rows: list[tuple[str, int]],
    n_distinct: int,
    total: int,
    top: int,
) -> None:
    click.echo()
    header = (
        f"{title} - overflow keys (top {min(top, n_distinct)} of {n_distinct})"
        if n_distinct
        else f"{title} - overflow keys: (none)"
    )
    click.echo(header)
    for key, n in rows:
        click.echo(f"  {key:<32} {n:>10,}  ({_pct(n, total)})")


@db.command(name="stats")
@click.argument("name")
@click.option(
    "--top",
    "top",
    type=int,
    default=20,
    show_default=True,
    help="Show at most TOP overflow-Map keys per table.",
)
def db_stats(name: str, top: int) -> None:
    """Schema-population stats for an ingested cohort."""
    from storage import db_disk_size, db_path, get_session

    path = db_path(name)
    if not path.exists():
        raise click.ClickException(f"db {name!r} not found")

    size_mb = db_disk_size(name) / 1_000_000
    _set_db(name)
    stats = _stats_payload(get_session(name), top)
    counts = stats["counts"]

    click.echo(f"db:        {name}")
    click.echo(f"path:      {path}")
    click.echo(f"size:      {size_mb:.1f} MB")
    click.echo()
    click.echo("counts:")
    for t in ("variants", "genotypes", "samples", "ingestions"):
        click.echo(f"  {t:<12} {counts[t]:>10,}")

    if stats["cohorts"]:
        click.echo()
        click.echo("cohorts:")
        for cohort, n in stats["cohorts"]:
            click.echo(f"  {cohort:<20} {int(n):>10,} samples")

    if stats["contigs"]:
        click.echo()
        click.echo("contigs:")
        for chrom, n in stats["contigs"]:
            click.echo(f"  {chrom:<20} {int(n):>10,} variants")

    v_total = counts["variants"]
    g_total = counts["genotypes"]
    _render_population(
        "variants - typed INFO column population", v_total, stats["variants_pop"]
    )
    info_rows, info_n = stats["info_extra"]
    _render_map_keys("variants.info_extra", info_rows, info_n, v_total, top)
    _render_population(
        "genotypes - typed column population", g_total, stats["genotypes_pop"]
    )
    fmt_rows, fmt_n = stats["format_extra"]
    _render_map_keys("genotypes.format_extra", fmt_rows, fmt_n, g_total, top)
