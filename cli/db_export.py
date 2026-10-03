"""`vcfclick db export` — a database, or a slice of it, back to VCF."""

from __future__ import annotations

from dataclasses import dataclass

import click

from cli.main import _set_db, db
from cli.options import command_options


@dataclass(frozen=True)
class ExportOptions:
    name: str
    out: str
    ingest_id: str | None
    samples: str | None
    regions: tuple[str, ...]
    regions_bed: str | None
    gene: str | None
    flank: int
    where: str | None
    pass_only: bool
    sites_only: bool
    absent_as: str


@db.command(name="export")
@click.argument("name")
@click.option(
    "-o",
    "--out",
    required=True,
    help="Output path: .vcf.gz (BGZF, tabix-indexable), .vcf, or - for stdout.",
)
@click.option(
    "--ingest-id",
    default=None,
    help="Ingestion to export (required when the database has several).",
)
@click.option("--samples", default=None, help="Comma-separated sample IDs to keep.")
@click.option(
    "--region",
    "regions",
    multiple=True,
    help="chr1, chr1:1000 or chr1:1000-2000 (repeatable).",
)
@click.option(
    "--regions-bed",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="BED file of regions.",
)
@click.option(
    "--gene",
    default=None,
    help="Only variants in this gene (needs `vcfclick annotations load`).",
)
@click.option(
    "--flank",
    type=int,
    default=0,
    show_default=True,
    help="Extend --gene by this many bp each side.",
)
@click.option(
    "--where",
    default=None,
    help='Extra SQL condition on the variants table, e.g. "info_AF < 0.01".',
)
@click.option(
    "--pass-only", is_flag=True, help="Only variants whose FILTER is PASS or unset."
)
@click.option(
    "--sites-only",
    is_flag=True,
    help="No genotypes: write the first eight columns only.",
)
@click.option(
    "--absent-as",
    type=click.Choice(["ref", "nocall"]),
    default="ref",
    show_default=True,
    help="How to write samples with no stored genotype (ignored for --keep-reference ingestions, "
    "where absent always means no-call).",
)
@command_options(ExportOptions)
def db_export(options: ExportOptions) -> None:
    """Export an ingestion to VCF, optionally sliced by region, gene, samples or SQL."""
    from importlib.metadata import PackageNotFoundError, version

    from export.vcf import ExportError, export_vcf, parse_region, read_bed
    from storage import db_path, get_session

    if not db_path(options.name).exists():
        raise click.ClickException(f"db {options.name!r} does not exist.")
    try:
        region_list = [parse_region(r) for r in options.regions]
        if options.regions_bed:
            region_list += read_bed(options.regions_bed)
        if options.gene:
            from annotations.db import position_for_gene

            g = position_for_gene(options.gene)
            if g is None:
                raise ExportError(
                    f"gene {options.gene.upper()!r} not found; load gene coordinates with `vcfclick annotations load`"
                )
            region_list.append(
                (
                    g.chrom,
                    max(1, int(g.start_pos) - options.flank),
                    int(g.end_pos) + options.flank,
                )
            )
        try:
            ver = version("vcfclick")
        except PackageNotFoundError:
            ver = ""
        _set_db(options.name)
        summary = export_vcf(
            get_session(options.name),
            options.name,
            options.out,
            ingest_id=options.ingest_id,
            samples=[s.strip() for s in options.samples.split(",") if s.strip()]
            if options.samples
            else None,
            regions=region_list or None,
            where=options.where,
            pass_only=options.pass_only,
            sites_only=options.sites_only,
            absent_as=options.absent_as,
            version=ver,
        )
    except ExportError as e:
        raise click.ClickException(str(e)) from None
    dest = "stdout" if options.out == "-" else options.out
    click.echo(
        f"exported {summary['variants']:,} variants x {summary['samples']:,} samples "
        f"from {summary['ingest_id']} to {dest}",
        err=True,
    )
