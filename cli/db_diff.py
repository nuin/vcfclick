"""Cohort-diff `vcfclick db` command."""

from __future__ import annotations

import json

import click

from cli.main import _set_db, db
from storage.sql import sql_quote_str


@db.command(name="diff")
@click.argument("name")
@click.option(
    "--cohort-a", required=True, help="First cohort name (as stored on samples.cohort)."
)
@click.option("--cohort-b", required=True, help="Second cohort name.")
@click.option(
    "--top",
    "top",
    type=int,
    default=50,
    show_default=True,
    help="Limit to the top N variants by absolute AF difference. 0 = no limit.",
)
@click.option(
    "--format",
    "fmt",
    default="PrettyCompact",
    show_default=True,
    help="chDB output format (PrettyCompact, JSON, CSV, TSV, ...).",
)
def db_diff(name: str, cohort_a: str, cohort_b: str, top: int, fmt: str) -> None:
    """Per-variant allele-frequency comparison across two cohorts."""
    from storage import db_path, get_session

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} not found")

    _set_db(name)
    sess = get_session(name)

    a = sql_quote_str(cohort_a)
    b = sql_quote_str(cohort_b)

    # JSONCompact, not TSV: TSV escapes a backslash in a name (`c\\ne`), so
    # the stored name would never equal the one typed.
    rows = json.loads(
        sess.query(
            f"SELECT cohort FROM samples WHERE cohort IN ({a}, {b}) GROUP BY cohort",
            "JSONCompact",
        )
        .bytes()
        .decode()
    )["data"]
    found = {r[0] for r in rows}
    missing = sorted({cohort_a, cohort_b} - found)
    if missing:
        raise click.ClickException(
            f"unknown cohort(s): {', '.join(missing)}. "
            f"Available: {', '.join(sorted(found)) or '(none)'}"
        )

    limit_clause = f"LIMIT {int(top)}" if top > 0 else ""
    sql = f"""
WITH cohort_sizes AS (
    SELECT cohort, count(DISTINCT (ingest_id, sample_id)) * 2 AS an
    FROM samples
    WHERE cohort IN ({a}, {b})
    GROUP BY cohort
),
calls AS (
    SELECT
        g.chrom, g.pos, g.ref, g.alt,
        COALESCE(sum(g.gt) FILTER (WHERE s.cohort = {a}), 0) AS ac_a,
        COALESCE(sum(g.gt) FILTER (WHERE s.cohort = {b}), 0) AS ac_b
    FROM genotypes g
    INNER JOIN samples s
        ON s.ingest_id = g.ingest_id AND s.sample_id = g.sample_id
    WHERE s.cohort IN ({a}, {b})
    GROUP BY g.chrom, g.pos, g.ref, g.alt
)
SELECT
    chrom, pos, ref, alt,
    ac_a, an_a, round(ac_a / an_a, 4) AS af_a,
    ac_b, an_b, round(ac_b / an_b, 4) AS af_b,
    round(ac_a / an_a - ac_b / an_b, 4) AS af_diff
FROM calls
CROSS JOIN (SELECT an AS an_a FROM cohort_sizes WHERE cohort = {a}) ca
CROSS JOIN (SELECT an AS an_b FROM cohort_sizes WHERE cohort = {b}) cb
WHERE an_a > 0 AND an_b > 0
ORDER BY abs(af_diff) DESC, chrom, pos
{limit_clause}
""".strip()

    out = sess.query(sql, fmt).bytes().decode()
    click.echo(out)
