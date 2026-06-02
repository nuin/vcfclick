"""`vcfclick` — small VCF databases, one per cohort.

Usage:
    vcfclick db create <name>
    vcfclick db list
    vcfclick db ingest <name> <vcf_path>
    vcfclick db query <name> "<sql>"
    vcfclick db info <name>
    vcfclick db dump <name> [--out <dir>]
    vcfclick db rm <name>

Each named DB is a self-contained chDB session at
$VCFCLICK_HOME/dbs/<name>/ (default ~/.vcfclick/dbs/<name>/). Setting the
VCFCLICK_DB_NAME env var inside a command propagates the choice to
ingest worker subprocesses and to the MCP server.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import click


# Reasonable default for `vcfclick db ingest`.
DEFAULT_WORKERS = 4


def _set_db(name: str) -> None:
    """Make the named DB the active target for storage + subprocess workers."""
    os.environ["VCFCLICK_DB_NAME"] = name


@click.group()
@click.version_option(package_name="vcfclick")
def cli() -> None:
    """vcfclick — small VCF databases."""


@cli.group()
def db() -> None:
    """Manage named VCF databases."""


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
    get_session(name)        # creates the directory + chDB session
    apply_schema()           # applies schema/*.sql against the new session
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
@click.option("--cohort", default="default", show_default=True,
              help="Cohort label this VCF belongs to.")
@click.option("--ingest-id", default=None,
              help="Stable upload identifier (UUID4 if omitted). Reuse to replace prior data.")
@click.option("--workers", default=DEFAULT_WORKERS, type=int, show_default=True,
              help="Parallel worker processes for ingestion.")
@click.option("--serial", is_flag=True,
              help="Use the single-process serial ingester instead of parallel.")
def db_ingest(
    name: str, vcf_path: str, cohort: str, ingest_id: str | None,
    workers: int, serial: bool,
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
            vcf_path, cohort, ingest_id=ingest_id, workers=workers,
        )


# ─── query ──────────────────────────────────────────────────────────

@db.command(name="query")
@click.argument("name")
@click.argument("sql")
@click.option("--format", "fmt", default="PrettyCompact", show_default=True,
              help="chDB output format (PrettyCompact, JSON, CSV, TSV, Vertical, ...).")
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
@click.option("--out", default=None, type=click.Path(file_okay=False),
              help="Output directory (default: ./<name>-dump/).")
def db_dump(name: str, out: str | None) -> None:
    """Export all tables from a named database to Parquet files."""
    from storage import db_path, get_session
    from export.parquet import export_all

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} not found")

    _set_db(name)
    get_session(name)                                # ensure session live
    out_dir = Path(out) if out else Path(f"{name}-dump")
    export_all(out_dir)
    click.echo(f"\ndumped to {out_dir.resolve()}")


# ─── rm ─────────────────────────────────────────────────────────────

@db.command(name="rm")
@click.argument("name")
@click.confirmation_option(
    prompt="really delete this database? (cannot be undone)"
)
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


if __name__ == "__main__":
    cli()
