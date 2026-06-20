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


def _list_typed_columns(sess, table: str) -> list[tuple[str, str]]:
    """Return typed, non-structural columns for `table`."""
    out = (
        sess.query(
            f"SELECT name, type FROM system.columns "
            f"WHERE table = '{table}' AND database = currentDatabase() "
            f"ORDER BY position",
            "TSV",
        )
        .bytes()
        .decode()
        .strip()
    )
    skip = _STATS_SKIP.get(table, set())
    cols: list[tuple[str, str]] = []
    for line in out.splitlines():
        name, type_str = line.split("\t", 1)
        if name not in skip:
            cols.append((name, type_str))
    return cols


def _population_expr(col_name: str, col_type: str) -> str:
    """Build the per-column aggregation fragment."""
    if "Nullable" in col_type:
        return f"countIf(`{col_name}` IS NOT NULL) AS `{col_name}`"
    return f"countIf(`{col_name}` != 0) AS `{col_name}`"


def _query_population(
    sess, table: str, columns: list[tuple[str, str]]
) -> dict[str, int]:
    """Return `{column_name: populated_row_count}` for typed columns."""
    if not columns:
        return {}
    exprs = ", ".join(_population_expr(n, t) for n, t in columns)
    out = sess.query(f"SELECT {exprs} FROM {table}", "TSV").bytes().decode().strip()
    values = [int(v) for v in out.split("\t")]
    return dict(zip([n for n, _ in columns], values))


def _query_map_keys(
    sess, table: str, map_col: str, top: int
) -> tuple[list[tuple[str, int]], int]:
    """Return top Map keys and total distinct key count."""
    n_distinct = int(
        sess.query(
            f"SELECT count(DISTINCT k) FROM {table} ARRAY JOIN mapKeys({map_col}) AS k",
            "TSV",
        )
        .bytes()
        .decode()
        .strip()
        or "0"
    )
    if n_distinct == 0:
        return [], 0
    out = (
        sess.query(
            f"SELECT k, count() AS n FROM {table} "
            f"ARRAY JOIN mapKeys({map_col}) AS k "
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
    counts = {
        t: int(
            sess.query(f"SELECT count() FROM {t}", "TSV").bytes().decode().strip()
            or "0"
        )
        for t in ("variants", "genotypes", "samples", "ingestions")
    }
    variants_cols = _list_typed_columns(sess, "variants")
    genotypes_cols = _list_typed_columns(sess, "genotypes")
    return {
        "counts": counts,
        "cohorts": _query_rows(
            sess,
            "SELECT cohort, count(DISTINCT (ingest_id, sample_id)) AS n "
            "FROM samples GROUP BY cohort ORDER BY n DESC, cohort",
        ),
        "contigs": _query_rows(
            sess,
            "SELECT chrom, count() AS n FROM variants "
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
    from storage import backend, db_disk_size, db_path, get_session

    if backend() == "duckdb":
        # The current implementation depends on chDB-specific SQL
        # (system.columns, countIf, ARRAY JOIN mapKeys). Porting to
        # DuckDB requires SQL-dialect helpers (information_schema,
        # FILTER clause, UNNEST(map_keys)); that lands as a follow-up.
        raise click.ClickException(
            "db stats is not yet implemented on the DuckDB backend. "
            "Set VCFCLICK_BACKEND=chdb (and `pip install vcfclick[chdb]`) "
            "to run stats today."
        )

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
