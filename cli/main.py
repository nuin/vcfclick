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
import tarfile
import tempfile
import urllib.request
from pathlib import Path

import click


# Parquet files a vcfclick bundle is expected to contain.
BUNDLE_TABLES = ("variants", "genotypes", "samples", "ingestions")


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


# ─── annotations group ─────────────────────────────────────────────

@cli.group()
def annotations() -> None:
    """Manage the embedded annotation reference store (DuckDB).

    Annotations are reference data shared across all named databases:
    gene coordinates today (GENCODE), transcript / exon / CDS depth and
    ClinVar in Phase 2. They live at `annotations/annotations.duckdb`
    inside the installed package directory.
    """


@annotations.command(name="load")
@click.option(
    "--gff", type=click.Path(exists=True, dir_okay=False), default=None,
    help="Local GENCODE GFF3 (.gff3 or .gff3.gz). Downloads v45 from EBI if omitted.",
)
@click.option(
    "--keep-existing", is_flag=True,
    help="Don't truncate refseq_genes before loading (default: replace).",
)
def annotations_load(gff: str | None, keep_existing: bool) -> None:
    """Populate the gene-coordinates table from a GENCODE GFF3.

    Required once after `pip install vcfclick` for the MCP server's
    `position_for_gene` tool to resolve symbols → coordinates. Re-run
    when GENCODE releases a new annotation version (yearly-ish).
    """
    from annotations.loaders.gencode_genes import load
    from pathlib import Path as _Path

    n = load(_Path(gff) if gff else None, replace=not keep_existing)
    click.echo(f"\nloaded   {n:,} genes into the annotation store")


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
    from storage import db_path, get_session
    from export.parquet import export_all

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

        # Acquire the tarball.
        if source.startswith(("http://", "https://")):
            tarball = tmp_dir / "bundle.tar.gz"
            click.echo(f"downloading  {source}")
            urllib.request.urlretrieve(source, tarball)
        else:
            tarball = Path(source).expanduser().resolve()
            if not tarball.exists():
                raise click.ClickException(f"{tarball} not found")

        # Extract with the safest standard filter (Python 3.12+).
        extract_dir = tmp_dir / "extract"
        extract_dir.mkdir()
        with tarfile.open(tarball, "r:gz") as tar:
            tar.extractall(extract_dir, filter="data")

        parquets = list(extract_dir.rglob("*.parquet"))
        if not parquets:
            raise click.ClickException("no parquet files in bundle")

        # Create the target DB and apply the schema.
        _set_db(name)
        sess = get_session(name)
        apply_schema()

        # Import each table by basename.
        #
        # SAFETY: the path interpolated into `INSERT ... FROM file('{p}')`
        # MUST be fully controlled by this code; a single quote in the
        # filename would break out of the SQL string literal. Two layers
        # of defence:
        #   1. resolve() + is_relative_to() rejects anything pointing
        #      outside extract_dir (defence-in-depth even though
        #      extractall(filter='data') already prevents path
        #      traversal).
        #   2. rename to a canonical `<table>.parquet` under extract_dir
        #      before referencing in SQL — extract_dir is a tempfile
        #      (no quotes/control chars), <table> is from BUNDLE_TABLES
        #      (hardcoded literals), so the resulting path is guaranteed
        #      to be SQL-safe.
        extract_root = extract_dir.resolve()
        for table in BUNDLE_TABLES:
            match = next((p for p in parquets if p.stem == table), None)
            if match is None:
                click.echo(f"  (no {table}.parquet in bundle, skipping)")
                continue

            resolved = match.resolve()
            if not resolved.is_relative_to(extract_root):
                raise click.ClickException(
                    f"refusing to import {match}: resolves outside the "
                    f"extraction directory"
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
            sess.query(
                f"INSERT INTO {table} "
                f"SELECT * FROM file('{sql_path}', 'Parquet')"
            )
            count = (
                sess.query(f"SELECT count() FROM {table}", "CSV")
                .bytes()
                .decode()
                .strip()
            )
            click.echo(f"  imported  {table:11s}  {count:>12s} rows")

    click.echo(f"\npulled   {name}  ←  {source}")


if __name__ == "__main__":
    cli()
