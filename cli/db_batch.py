"""Batch-ingest `vcfclick db` command."""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import click

from cli.main import _set_db, db


@dataclass(frozen=True)
class _BatchEntry:
    ingest_id: str
    vcf_path: Path
    cohort: str


_VCF_EXTS = (".vcf.gz", ".vcf.bgz", ".vcf")


def _derive_ingest_id(vcf_path: Path) -> str:
    """Derive an ingest_id from a VCF filename."""
    name = vcf_path.name
    for ext in _VCF_EXTS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return vcf_path.stem


def _scan_dir_for_vcfs(directory: Path, cohort: str) -> list[_BatchEntry]:
    """List batch VCFs under `directory` without recursing."""
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
    """Parse a TSV manifest with required `vcf_path` column."""
    entries: list[_BatchEntry] = []
    seen_ids: set[str] = set()

    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None or "vcf_path" not in reader.fieldnames:
            raise click.ClickException(
                f"manifest must have a `vcf_path` column (found: {reader.fieldnames})"
            )

        for line_no, row in enumerate(reader, start=2):
            raw = (row.get("vcf_path") or "").strip()
            if not raw:
                continue
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
    """Ingest many per-sample VCFs into NAME as one cohort."""
    from storage import db_path

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
        assert manifest is not None
        entries = _parse_manifest(manifest, default_cohort=cohort or "")

    click.echo(f"[batch] {len(entries)} VCFs to ingest into {name}")

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
            short_err = str(e).splitlines()[0][:200]
            click.echo(f"[batch]   skip - {short_err}", err=True)
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
