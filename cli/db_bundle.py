"""Bundle-oriented `vcfclick db` commands."""

from __future__ import annotations

import tarfile
import tempfile
import urllib.request
from pathlib import Path

import click

from cli.main import _set_db, db

BUNDLE_TABLES = ("variants", "genotypes", "samples", "ingestions")


@db.command(name="push")
@click.argument("name")
@click.argument("out_path", type=click.Path(dir_okay=False))
def db_push(name: str, out_path: str) -> None:
    """Dump a database and bundle it as a portable tar.gz file."""
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
    click.echo(f"\npushed   {name}  ->  {out}  ({size_mb:.1f} MB)")


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

    from storage import count_expr, parquet_file_expr

    sess.query(f"INSERT INTO {table} SELECT * FROM {parquet_file_expr(str(safe_path))}")
    raw = (
        sess.query(f"SELECT {count_expr()} FROM {table}", "CSV")
        .bytes()
        .decode()
        .strip()
    )
    # CSV from DuckDB carries a header line; chDB CSV is headerless.
    last = [ln for ln in raw.splitlines() if ln.strip()]
    count = last[-1] if last else "0"
    click.echo(f"  imported  {table:11s}  {count:>12s} rows")


@db.command(name="pull")
@click.argument("name")
@click.argument("source")
def db_pull(name: str, source: str) -> None:
    """Restore a database from a tar.gz bundle (HTTPS URL or local file)."""
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

        for table in BUNDLE_TABLES:
            _import_table(sess, table, parquets, extract_dir)

    click.echo(f"\npulled   {name}  <-  {source}")
