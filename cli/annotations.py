"""`vcfclick annotations ...` subcommands.

Registered against the ``annotations`` group defined in :mod:`cli.main`
at import time.
"""

from __future__ import annotations

from pathlib import Path

import click

from cli.main import annotations


@annotations.command(name="load")
@click.option(
    "--gff",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Local GENCODE GFF3 (.gff3 or .gff3.gz). Downloads v45 from EBI if omitted.",
)
@click.option(
    "--keep-existing",
    is_flag=True,
    help="Don't truncate refseq_genes before loading (default: replace).",
)
def annotations_load(gff: str | None, keep_existing: bool) -> None:
    """Populate the gene-coordinates table from a GENCODE GFF3.

    Required once after `pip install vcfclick` for the MCP server's
    `position_for_gene` tool to resolve symbols → coordinates. Re-run
    when GENCODE releases a new annotation version (yearly-ish).
    """
    from annotations.loaders.gencode_genes import load

    n = load(Path(gff) if gff else None, replace=not keep_existing)
    click.echo(f"\nloaded   {n:,} genes into the annotation store")


@annotations.command(name="load-clinvar")
@click.option(
    "--vcf",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Local ClinVar VCF (.vcf.gz). Downloads the current NCBI weekly release if omitted.",
)
@click.option(
    "--keep-existing",
    is_flag=True,
    help="Don't truncate clinvar_variants before loading (default: replace).",
)
def annotations_load_clinvar(vcf: str | None, keep_existing: bool) -> None:
    """Populate the ClinVar significance table from the NCBI ClinVar VCF.

    Required for the MCP server's `clinvar_lookup` tool to return real
    significance calls. The NCBI VCF refreshes weekly; re-run monthly
    (or before any clinically-adjacent demo) to stay current. Bare
    numeric contigs are normalised to `chr`-prefixed during load so
    lookups against sample data (which uses chr-style) compose.
    """
    from annotations.loaders.clinvar import load

    n = load(Path(vcf) if vcf else None, replace=not keep_existing)
    click.echo(f"\nloaded   {n:,} ClinVar variants into the annotation store")


@annotations.command(name="load-gnomad")
@click.argument("vcf", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--replace",
    is_flag=True,
    help="Truncate gnomad_af before loading (default: append, so several "
    "per-chromosome slices can be loaded incrementally).",
)
def annotations_load_gnomad(vcf: str, replace: bool) -> None:
    """Load gnomAD allele frequencies from a gnomAD sites VCF.

    gnomAD is too large to bundle, so pass a VCF you supply — a region
    slice, an af-only file, or a per-chromosome sites VCF. A small region
    can be pulled with tabix-over-HTTPS from the public gnomAD bucket; see
    docs/MCP.md. Powers the `gnomad_lookup` MCP tool and the
    `db trio --gnomad-max-af` rarity filter.
    """
    from annotations.loaders.gnomad import load

    n = load(Path(vcf), replace=replace)
    click.echo(f"\nloaded   {n:,} gnomAD allele frequencies into the annotation store")
