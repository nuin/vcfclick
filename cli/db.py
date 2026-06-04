"""`vcfclick db ...` subcommands.

Registered against the ``db`` group defined in :mod:`cli.main` at import
time — the module is imported for its side effects from cli.main.
"""

from __future__ import annotations

import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path

import click

from cli.main import _set_db, db


# Parquet files a vcfclick bundle is expected to contain.
BUNDLE_TABLES = ("variants", "genotypes", "samples", "ingestions")

# Reasonable default for `vcfclick db ingest`.
DEFAULT_WORKERS = 4


def _quote_str(s: str) -> str:
    """SQL-standard single-quoted string literal.

    chDB's `session.query()` doesn't take parameter bindings, so we
    interpolate strings into SQL by hand. Doubling embedded single
    quotes is the portable escape and matches what ClickHouse parses.
    """
    return "'" + s.replace("'", "''") + "'"


# ─── create ─────────────────────────────────────────────────────────


@db.command(name="create")
@click.argument("name")
def db_create(name: str) -> None:
    """Create a new empty database with the vcfclick schema applied."""
    from storage import apply_schema, db_path, get_session

    path = db_path(name)
    if path.exists() and any(path.iterdir()):
        raise click.ClickException(
            f"db {name!r} already exists at {path}. Use `vcfclick db rm {name}` "
            "first if you want to start over."
        )

    _set_db(name)
    get_session(name)  # creates the directory + chDB session
    apply_schema()  # applies schema/*.sql against the new session
    click.echo(f"created  {name}  →  {path}")


# ─── list ───────────────────────────────────────────────────────────


@db.command(name="list")
def db_list() -> None:
    """List all named databases."""
    from storage import db_disk_size, list_dbs

    names = list_dbs()
    if not names:
        click.echo("(no databases yet — `vcfclick db create <name>` to start)")
        return

    click.echo(f"{'NAME':32s}  {'SIZE':>10s}")
    for n in names:
        size = db_disk_size(n)
        click.echo(f"{n:32s}  {size / 1_000_000:>7.1f} MB")


# ─── info ───────────────────────────────────────────────────────────


@db.command(name="info")
@click.argument("name")
def db_info(name: str) -> None:
    """Show metadata about a database (row counts, ingestions, size)."""
    from storage import db_disk_size, db_path, get_session

    path = db_path(name)
    if not path.exists():
        raise click.ClickException(f"db {name!r} not found at {path}")

    _set_db(name)
    sess = get_session(name)

    def scalar(sql: str) -> str:
        try:
            return sess.query(sql, "CSV").bytes().decode().strip()
        except Exception:
            return "(schema not yet applied)"

    click.echo(f"db:        {name}")
    click.echo(f"path:      {path}")
    click.echo(f"size:      {db_disk_size(name) / 1_000_000:.1f} MB")
    click.echo(f"variants:  {scalar('SELECT count() FROM variants')}")
    click.echo(f"genotypes: {scalar('SELECT count() FROM genotypes')}")
    click.echo(f"samples:   {scalar('SELECT count() FROM samples')}")
    click.echo(f"ingestions:{scalar('SELECT count() FROM ingestions')}")


# ─── ingest ─────────────────────────────────────────────────────────


@db.command(name="ingest")
@click.argument("name")
@click.argument("vcf_path", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--cohort",
    default="default",
    show_default=True,
    help="Cohort label this VCF belongs to.",
)
@click.option(
    "--ingest-id",
    default=None,
    help="Stable upload identifier (UUID4 if omitted). Reuse to replace prior data.",
)
@click.option(
    "--workers",
    default=DEFAULT_WORKERS,
    type=int,
    show_default=True,
    help="Parallel worker processes for ingestion.",
)
@click.option(
    "--serial",
    is_flag=True,
    help="Use the single-process serial ingester instead of parallel.",
)
def db_ingest(
    name: str,
    vcf_path: str,
    cohort: str,
    ingest_id: str | None,
    workers: int,
    serial: bool,
) -> None:
    """Ingest a (normalised) VCF into a named database."""
    from storage import db_path

    if not db_path(name).exists():
        raise click.ClickException(
            f"db {name!r} does not exist. Run `vcfclick db create {name}` first."
        )

    _set_db(name)

    if serial:
        from ingest.vcf_load import ingest as ingest_serial

        ingest_serial(vcf_path, cohort, ingest_id)
    else:
        from ingest.parallel import ingest_parallel

        ingest_parallel(
            vcf_path,
            cohort,
            ingest_id=ingest_id,
            workers=workers,
        )


# ─── query ──────────────────────────────────────────────────────────


@db.command(name="query")
@click.argument("name")
@click.argument("sql")
@click.option(
    "--format",
    "fmt",
    default="PrettyCompact",
    show_default=True,
    help="chDB output format (PrettyCompact, JSON, CSV, TSV, Vertical, ...).",
)
def db_query(name: str, sql: str, fmt: str) -> None:
    """Run a SQL query against a named database and print the result."""
    from storage import db_path, get_session

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} not found")

    _set_db(name)
    sess = get_session(name)
    out = sess.query(sql, fmt).bytes().decode()
    click.echo(out)


# ─── dump ───────────────────────────────────────────────────────────


@db.command(name="dump")
@click.argument("name")
@click.option(
    "--out",
    default=None,
    type=click.Path(file_okay=False),
    help="Output directory (default: ./<name>-dump/).",
)
def db_dump(name: str, out: str | None) -> None:
    """Export all tables from a named database to Parquet files."""
    from export.parquet import export_all
    from storage import db_path, get_session

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} not found")

    _set_db(name)
    get_session(name)  # ensure session live
    out_dir = Path(out) if out else Path(f"{name}-dump")
    export_all(out_dir)
    click.echo(f"\ndumped to {out_dir.resolve()}")


# ─── rm ─────────────────────────────────────────────────────────────


@db.command(name="rm")
@click.argument("name")
@click.confirmation_option(prompt="really delete this database? (cannot be undone)")
def db_rm(name: str) -> None:
    """Permanently delete a named database and all its data."""
    from storage import db_path

    path = db_path(name)
    if not path.exists():
        raise click.ClickException(f"db {name!r} not found")
    shutil.rmtree(path)
    click.echo(f"removed  {name}  ←  {path}")


# ─── path ───────────────────────────────────────────────────────────


@db.command(name="path")
@click.argument("name")
def db_path_cmd(name: str) -> None:
    """Print the on-disk path of a named database (no checks)."""
    from storage import db_path

    click.echo(db_path(name))


# ─── push ───────────────────────────────────────────────────────────


@db.command(name="push")
@click.argument("name")
@click.argument("out_path", type=click.Path(dir_okay=False))
def db_push(name: str, out_path: str) -> None:
    """Dump a database and bundle it as a portable tar.gz file.

    The bundle is a flat archive of variants.parquet, genotypes.parquet,
    samples.parquet, ingestions.parquet — exactly what `db pull` consumes.
    Upload the resulting file anywhere (S3, HTTPS, scp, USB) and the
    receiver runs `vcfclick db pull <name> <url|path>` to restore it.
    """
    from export.parquet import export_all
    from storage import db_path, get_session

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} not found")

    out = Path(out_path).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    _set_db(name)
    get_session(name)

    with tempfile.TemporaryDirectory(prefix="vcfclick_push_") as tmp:
        tmp_dir = Path(tmp)
        export_all(tmp_dir)
        with tarfile.open(out, "w:gz") as tar:
            for f in sorted(tmp_dir.iterdir()):
                tar.add(f, arcname=f.name)

    size_mb = out.stat().st_size / 1_000_000
    click.echo(f"\npushed   {name}  →  {out}  ({size_mb:.1f} MB)")


# ─── pull ───────────────────────────────────────────────────────────


def _acquire_tarball(source: str, tmp_dir: Path) -> Path:
    """Resolve `source` (URL or local path) to a tarball path on disk."""
    if source.startswith(("http://", "https://")):
        tarball = tmp_dir / "bundle.tar.gz"
        click.echo(f"downloading  {source}")
        urllib.request.urlretrieve(source, tarball)
        return tarball
    tarball = Path(source).expanduser().resolve()
    if not tarball.exists():
        raise click.ClickException(f"{tarball} not found")
    return tarball


def _import_table(sess, table: str, parquets: list[Path], extract_dir: Path) -> None:
    """Find <table>.parquet in `parquets`, validate safety, INSERT it."""
    match = next((p for p in parquets if p.stem == table), None)
    if match is None:
        click.echo(f"  (no {table}.parquet in bundle, skipping)")
        return

    extract_root = extract_dir.resolve()
    resolved = match.resolve()
    if not resolved.is_relative_to(extract_root):
        raise click.ClickException(
            f"refusing to import {match}: resolves outside the extraction directory"
        )

    safe_path = extract_dir / f"{table}.parquet"
    if match != safe_path:
        match.rename(safe_path)

    # Triple-belt: rename + extract_root check above already make
    # `safe_path` SQL-safe (tempfile dirs cannot contain quotes on
    # any supported OS). The SQL-standard quote-doubling here is
    # defensive against a future regression that might point
    # tempfile at a user-controlled location.
    sql_path = str(safe_path).replace("'", "''")
    sess.query(f"INSERT INTO {table} SELECT * FROM file('{sql_path}', 'Parquet')")
    count = sess.query(f"SELECT count() FROM {table}", "CSV").bytes().decode().strip()
    click.echo(f"  imported  {table:11s}  {count:>12s} rows")


@db.command(name="pull")
@click.argument("name")
@click.argument("source")
def db_pull(name: str, source: str) -> None:
    """Restore a database from a tar.gz bundle (HTTPS URL or local file).

    Creates a new named DB called <name>, applies the schema, and imports
    each Parquet file in the bundle. Refuses to overwrite an existing
    DB — run `vcfclick db rm <name>` first if you want to replace.

    Source can be:
      · an https://… or http://… URL (downloaded to a temp file)
      · a local .tar.gz / .tgz path
    """
    from storage import apply_schema, db_path, get_session

    target = db_path(name)
    if target.exists() and any(target.iterdir()):
        raise click.ClickException(
            f"db {name!r} already exists at {target}. "
            f"Run `vcfclick db rm {name}` first."
        )

    with tempfile.TemporaryDirectory(prefix="vcfclick_pull_") as tmp:
        tmp_dir = Path(tmp)
        tarball = _acquire_tarball(source, tmp_dir)

        # Extract with the safest standard filter (Python 3.12+).
        extract_dir = tmp_dir / "extract"
        extract_dir.mkdir()
        with tarfile.open(tarball, "r:gz") as tar:
            tar.extractall(extract_dir, filter="data")

        parquets = list(extract_dir.rglob("*.parquet"))
        if not parquets:
            raise click.ClickException("no parquet files in bundle")

        _set_db(name)
        sess = get_session(name)
        apply_schema()

        # SAFETY: the path interpolated into `INSERT ... FROM file('{p}')`
        # MUST be fully controlled by this code; a single quote in the
        # filename would break out of the SQL string literal. _import_table
        # rejects anything resolving outside extract_dir (defence-in-depth
        # even though extractall(filter='data') already prevents path
        # traversal) and renames to a canonical `<table>.parquet` before
        # any SQL interpolation.
        for table in BUNDLE_TABLES:
            _import_table(sess, table, parquets, extract_dir)

    click.echo(f"\npulled   {name}  ←  {source}")


# ─── diff ───────────────────────────────────────────────────────────


@db.command(name="diff")
@click.argument("name")
@click.option("--cohort-a", required=True, help="First cohort name (as stored on samples.cohort).")
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
    """Per-variant allele-frequency comparison across two cohorts.

    Computes AC / AN / AF for cohorts A and B (and their difference)
    for every variant with at least one non-reference call in either
    cohort. AN is `2 × distinct samples in the cohort`; samples
    absent from the sparse genotypes table at a given variant are
    assumed 0/0 (the sparse-storage convention). Results are sorted
    by absolute AF difference, descending.

    Multi-ingestion semantics: a cohort spans every (ingest_id, sample_id)
    where samples.cohort matches. Samples ingested twice under
    different ingest_ids count as two observations — by design,
    matching the schema's cross-ingestion non-merging rule.
    """
    from storage import db_path, get_session

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} not found")

    _set_db(name)
    sess = get_session(name)

    a = _quote_str(cohort_a)
    b = _quote_str(cohort_b)

    # Pre-check: both cohorts must have at least one sample. Without
    # this the main query happily returns an empty result and the user
    # can't tell whether it's "no variants differ" or "you typoed a
    # cohort name".
    rows = (
        sess.query(
            f"SELECT cohort FROM samples WHERE cohort IN ({a}, {b}) GROUP BY cohort",
            "TSV",
        )
        .bytes()
        .decode()
        .strip()
    )
    found = set(rows.splitlines()) if rows else set()
    missing = sorted({cohort_a, cohort_b} - found)
    if missing:
        raise click.ClickException(
            f"unknown cohort(s): {', '.join(missing)}. "
            f"Available: {', '.join(sorted(found)) or '(none)'}"
        )

    limit_clause = f"LIMIT {int(top)}" if top > 0 else ""

    # `sumIf(gt, cohort = X)` counts alt alleles in each cohort.
    # `cohort_sizes` provides the denominator (2 × distinct samples).
    # CROSS JOIN over a single-row table is the chDB-friendly way to
    # bring a scalar into a per-row SELECT — works around the lack of
    # uncorrelated scalar subqueries in older ClickHouse / chDB.
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
        sumIf(g.gt, s.cohort = {a}) AS ac_a,
        sumIf(g.gt, s.cohort = {b}) AS ac_b
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
