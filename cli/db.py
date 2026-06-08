"""`vcfclick db ...` subcommands.

Registered against the ``db`` group defined in :mod:`cli.main` at import
time — the module is imported for its side effects from cli.main.
"""

from __future__ import annotations

import csv
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from dataclasses import dataclass
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


# ─── ingest-parquet ─────────────────────────────────────────────────


@db.command(name="ingest-parquet")
@click.argument("name")
@click.argument("dump_dir", type=click.Path(exists=True, file_okay=False))
@click.option(
    "--cohort",
    default="default",
    show_default=True,
    help="Cohort label to assign to the imported data.",
)
@click.option(
    "--ingest-id",
    default=None,
    help="Stable upload identifier (UUID4 if omitted). Reuse to replace prior data.",
)
def db_ingest_parquet(
    name: str,
    dump_dir: str,
    cohort: str,
    ingest_id: str | None,
) -> None:
    """Ingest a Parquet dump (produced by `db dump`) into NAME.

    DUMP_DIR is a directory containing variants.parquet (required) and
    optionally genotypes.parquet / samples.parquet. The ingest_id and
    cohort columns of the source files are overridden with --ingest-id
    and --cohort. Reuse an --ingest-id to replace prior data under it.
    """
    from storage import db_path

    if not db_path(name).exists():
        raise click.ClickException(
            f"db {name!r} does not exist. Run `vcfclick db create {name}` first."
        )

    _set_db(name)

    from ingest.parquet_load import ingest_from_parquet

    ingest_from_parquet(dump_dir, cohort, ingest_id)


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


# ─── ingest-batch ───────────────────────────────────────────────────


@dataclass(frozen=True)
class _BatchEntry:
    ingest_id: str
    vcf_path: Path
    cohort: str


# VCF extensions stripped when deriving an ingest_id from a filename.
# Order matters: ".vcf.gz" before ".vcf" so the bgzipped form wins.
_VCF_EXTS = (".vcf.gz", ".vcf.bgz", ".vcf")


def _derive_ingest_id(vcf_path: Path) -> str:
    """`/path/HG00096.vcf.gz` → `HG00096`."""
    name = vcf_path.name
    for ext in _VCF_EXTS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return vcf_path.stem


def _scan_dir_for_vcfs(directory: Path, cohort: str) -> list[_BatchEntry]:
    """List `*.vcf.gz` and `*.vcf.bgz` under `directory`.

    Non-recursive: per-sample VCFs from one batch typically live in a
    single directory; descending into subdirectories would surprise
    users who have unrelated VCFs nearby.
    """
    vcfs = sorted(directory.glob("*.vcf.gz")) + sorted(directory.glob("*.vcf.bgz"))
    if not vcfs:
        raise click.ClickException(f"no *.vcf.gz files found under {directory}")
    return [
        _BatchEntry(
            ingest_id=_derive_ingest_id(v),
            vcf_path=v.resolve(),
            cohort=cohort,
        )
        for v in vcfs
    ]


def _parse_manifest(path: Path, default_cohort: str) -> list[_BatchEntry]:
    """Parse a TSV manifest.

    Required column: `vcf_path` (resolved relative to the manifest file
    if not absolute, matching nf-core / Snakemake convention).

    Optional columns:
      * `sample_id` or `ingest_id` — defaults to the filename-derived
        ingest_id from the vcf_path basename.
      * `cohort` — defaults to the `--cohort` flag value passed to the
        CLI, so a manifest can mix cohorts on a per-row basis or
        delegate cohort assignment to the CLI flag entirely.
    """
    entries: list[_BatchEntry] = []
    seen_ids: set[str] = set()

    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None or "vcf_path" not in reader.fieldnames:
            raise click.ClickException(
                f"manifest must have a `vcf_path` column (found: {reader.fieldnames})"
            )

        for line_no, row in enumerate(reader, start=2):  # 1 = header
            raw = (row.get("vcf_path") or "").strip()
            if not raw:
                continue  # blank rows allowed
            vcf = Path(raw)
            if not vcf.is_absolute():
                vcf = (path.parent / vcf).resolve()
            if not vcf.exists():
                raise click.ClickException(
                    f"manifest line {line_no}: VCF not found at {vcf}"
                )

            ingest_id = (
                row.get("ingest_id") or row.get("sample_id") or _derive_ingest_id(vcf)
            ).strip()
            if ingest_id in seen_ids:
                raise click.ClickException(
                    f"manifest line {line_no}: duplicate ingest_id "
                    f"{ingest_id!r} (must be unique within a manifest)"
                )
            seen_ids.add(ingest_id)

            cohort = (row.get("cohort") or default_cohort).strip()
            if not cohort:
                raise click.ClickException(
                    f"manifest line {line_no}: empty cohort and no "
                    f"--cohort default; supply one or the other"
                )

            entries.append(_BatchEntry(ingest_id, vcf, cohort))

    if not entries:
        raise click.ClickException(f"manifest has no entries: {path}")
    return entries


@db.command(name="ingest-batch")
@click.argument("name")
@click.option(
    "--from-dir",
    "from_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Ingest every *.vcf.gz under DIR. ingest_id is the filename stem.",
)
@click.option(
    "--manifest",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="TSV with required `vcf_path` column; optional `sample_id`/"
    "`ingest_id` and `cohort` columns override the defaults.",
)
@click.option(
    "--cohort",
    default=None,
    help="Default cohort label. Required with --from-dir; used as the "
    "fallback for manifest rows that don't carry their own `cohort`.",
)
def db_ingest_batch(
    name: str,
    from_dir: Path | None,
    manifest: Path | None,
    cohort: str | None,
) -> None:
    """Ingest many per-sample VCFs into NAME as one cohort.

    Per-sample VCFs are what clinical pipelines (DRAGEN, GATK's
    `-ERC GVCF` joint-calling workflow, etc.) emit by default. Each
    input file becomes its own `ingest_id`; samples table rows
    accumulate under a shared cohort label.

    Files are ingested SEQUENTIALLY using the existing single-file
    serial loader. The atomic guarantee from `db ingest` applies
    per-file: a failure on one VCF rolls back only that file's
    writes and the batch continues. Exit code is non-zero if any
    file failed, so a shell pipeline of
        vcfclick db ingest-batch ... && echo all-clean
    still works.
    """
    from storage import db_path

    # XOR validation — pick exactly one source.
    if from_dir and manifest:
        raise click.ClickException(
            "--from-dir and --manifest are mutually exclusive; pick one"
        )
    if not from_dir and not manifest:
        raise click.ClickException(
            "provide --from-dir or --manifest to point at the VCFs to ingest"
        )
    if from_dir and not cohort:
        raise click.ClickException(
            "--cohort is required with --from-dir (no per-row source for it)"
        )

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} not found")
    _set_db(name)

    if from_dir:
        entries = _scan_dir_for_vcfs(from_dir.resolve(), cohort or "")
    else:
        assert manifest is not None  # narrows for type-checkers
        entries = _parse_manifest(manifest, default_cohort=cohort or "")

    click.echo(f"[batch] {len(entries)} VCFs to ingest into {name}")

    # Sequential single-file ingest via library call so the chDB
    # session is shared across iterations. Spawning a fresh subprocess
    # per file would amortise chDB startup × N which dominates for
    # small per-sample VCFs. The existing rollback inside _ingest
    # ensures any failure on one file is contained.
    from ingest.vcf_load import ingest as _ingest

    successes: list[_BatchEntry] = []
    failures: list[tuple[_BatchEntry, str]] = []

    for i, entry in enumerate(entries, start=1):
        click.echo(
            f"[batch] [{i:>4}/{len(entries)}] {entry.ingest_id} (cohort={entry.cohort})"
        )
        try:
            _ingest(
                str(entry.vcf_path),
                cohort=entry.cohort,
                ingest_id=entry.ingest_id,
            )
            successes.append(entry)
        except Exception as e:
            # _ingest already ran its rollback; this branch just
            # records the file as failed and lets the loop continue.
            short_err = str(e).splitlines()[0][:200]
            click.echo(f"[batch]   skip — {short_err}", err=True)
            failures.append((entry, short_err))

    click.echo()
    click.echo("Batch ingest summary:")
    click.echo(f"  total:    {len(entries)}")
    click.echo(f"  ingested: {len(successes)}")
    click.echo(f"  failed:   {len(failures)}")
    if failures:
        click.echo()
        click.echo("Failed:")
        for entry, err in failures:
            click.echo(f"  {entry.ingest_id} ({entry.vcf_path.name}): {err}")
        sys.exit(1)


# ─── stats ──────────────────────────────────────────────────────────


# Structural / Map columns that aren't "typed data" — keys, sort columns,
# the overflow Maps themselves, ingest timestamps. We exclude these from
# population stats because "what % of rows have chrom?" is uninteresting.
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
    """Return [(name, type), ...] for typed (non-structural) columns
    on `table`, in schema-declaration order."""
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
        if name in skip:
            continue
        cols.append((name, type_str))
    return cols


def _population_expr(col_name: str, col_type: str) -> str:
    """Build the per-column aggregation fragment. Nullable columns
    count `IS NOT NULL`; non-nullable flag-style columns (UInt8 default
    0) count non-zero — both align with the human intuition for "is
    this field actually used in this data?"."""
    if "Nullable" in col_type:
        return f"countIf(`{col_name}` IS NOT NULL) AS `{col_name}`"
    # Non-nullable: defaults to 0 in our schema. Count the non-default
    # values, which for flag columns matches "% rows where the flag is set".
    return f"countIf(`{col_name}` != 0) AS `{col_name}`"


def _query_population(
    sess, table: str, columns: list[tuple[str, str]]
) -> dict[str, int]:
    """One aggregation query for every typed column on `table`. Returns
    `{column_name: populated_row_count}`."""
    if not columns:
        return {}
    exprs = ", ".join(_population_expr(n, t) for n, t in columns)
    out = sess.query(f"SELECT {exprs} FROM {table}", "TSV").bytes().decode().strip()
    values = [int(v) for v in out.split("\t")]
    return dict(zip([n for n, _ in columns], values))


def _query_map_keys(
    sess, table: str, map_col: str, top: int
) -> tuple[list[tuple[str, int]], int]:
    """Top-N most frequent keys in `map_col`. Returns ([(key, n), ...],
    total_distinct_keys)."""
    # Distinct key count first — drives the "top X of Y" header.
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
    """Schema-population stats for an ingested cohort.

    Reports row counts, cohort and contig breakdowns, the populated
    fraction of every typed column, and the most frequent
    overflow-Map keys. Use this when you want to know "what fields
    are actually in this cohort's data?" — the discover command
    answers the same question for a VCF before ingest; `db stats`
    answers it for the stored data after.
    """
    from storage import db_disk_size, db_path, get_session

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} not found")

    path = db_path(name)
    size_mb = db_disk_size(name) / 1_000_000

    _set_db(name)
    sess = get_session(name)

    counts = {
        t: int(
            sess.query(f"SELECT count() FROM {t}", "TSV").bytes().decode().strip()
            or "0"
        )
        for t in ("variants", "genotypes", "samples", "ingestions")
    }

    cohorts_out = (
        sess.query(
            "SELECT cohort, count(DISTINCT (ingest_id, sample_id)) AS n "
            "FROM samples GROUP BY cohort ORDER BY n DESC, cohort",
            "TSV",
        )
        .bytes()
        .decode()
        .strip()
    )
    cohorts = [line.split("\t") for line in cohorts_out.splitlines() if line]

    contigs_out = (
        sess.query(
            "SELECT chrom, count() AS n FROM variants "
            "GROUP BY chrom ORDER BY n DESC, chrom",
            "TSV",
        )
        .bytes()
        .decode()
        .strip()
    )
    contigs = [line.split("\t") for line in contigs_out.splitlines() if line]

    variants_cols = _list_typed_columns(sess, "variants")
    variants_pop = _query_population(sess, "variants", variants_cols)
    info_extra_top, info_extra_n = _query_map_keys(sess, "variants", "info_extra", top)

    genotypes_cols = _list_typed_columns(sess, "genotypes")
    genotypes_pop = _query_population(sess, "genotypes", genotypes_cols)
    fmt_extra_top, fmt_extra_n = _query_map_keys(sess, "genotypes", "format_extra", top)

    # ─── render ──────────────────────────────────────────────────────

    click.echo(f"db:        {name}")
    click.echo(f"path:      {path}")
    click.echo(f"size:      {size_mb:.1f} MB")
    click.echo()
    click.echo("counts:")
    for t in ("variants", "genotypes", "samples", "ingestions"):
        click.echo(f"  {t:<12} {counts[t]:>10,}")

    if cohorts:
        click.echo()
        click.echo("cohorts:")
        for cohort, n in cohorts:
            click.echo(f"  {cohort:<20} {int(n):>10,} samples")

    if contigs:
        click.echo()
        click.echo("contigs:")
        for chrom, n in contigs:
            click.echo(f"  {chrom:<20} {int(n):>10,} variants")

    v_total = counts["variants"]
    click.echo()
    click.echo(f"variants — typed INFO column population (of {v_total:,} rows):")
    # Sort by population descending so populated columns surface first.
    sorted_v = sorted(variants_pop.items(), key=lambda kv: (-kv[1], kv[0]))
    for col, n in sorted_v:
        click.echo(f"  {col:<32} {n:>10,}  ({_pct(n, v_total)})")

    click.echo()
    header = (
        f"variants.info_extra — overflow keys "
        f"(top {min(top, info_extra_n)} of {info_extra_n})"
        if info_extra_n
        else "variants.info_extra — overflow keys: (none)"
    )
    click.echo(header)
    for key, n in info_extra_top:
        click.echo(f"  {key:<32} {n:>10,}  ({_pct(n, v_total)})")

    g_total = counts["genotypes"]
    click.echo()
    click.echo(f"genotypes — typed column population (of {g_total:,} rows):")
    sorted_g = sorted(genotypes_pop.items(), key=lambda kv: (-kv[1], kv[0]))
    for col, n in sorted_g:
        click.echo(f"  {col:<32} {n:>10,}  ({_pct(n, g_total)})")

    click.echo()
    header = (
        f"genotypes.format_extra — overflow keys "
        f"(top {min(top, fmt_extra_n)} of {fmt_extra_n})"
        if fmt_extra_n
        else "genotypes.format_extra — overflow keys: (none)"
    )
    click.echo(header)
    for key, n in fmt_extra_top:
        click.echo(f"  {key:<32} {n:>10,}  ({_pct(n, g_total)})")
