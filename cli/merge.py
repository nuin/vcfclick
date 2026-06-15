"""`vcfclick merge` — combine per-sample VCFs into one joint VCF.

    vcfclick merge proband.vcf.gz father.vcf.gz mother.vcf.gz -o trio.vcf.gz

The joint output is what trio analysis needs: all samples under one
ingest_id so a single query compares their genotypes at each site.

    vcfclick merge p.vcf.gz f.vcf.gz m.vcf.gz -o trio.vcf.gz
    vcfclick db create fam1
    vcfclick db ingest fam1 trio.vcf.gz --cohort trio --ingest-id fam1
    vcfclick db ped fam1 fam1.ped
    vcfclick db trio fam1 --proband PROBAND
"""

from __future__ import annotations

import click

from cli.main import cli


@cli.command(name="merge")
@click.argument("vcfs", nargs=-1, required=True, type=click.Path(dir_okay=False))
@click.option(
    "-o",
    "--out",
    "output",
    required=True,
    type=click.Path(dir_okay=False),
    help="Output joint VCF path (.vcf.gz).",
)
@click.option(
    "--multiallelic",
    default="none",
    show_default=True,
    help="bcftools merge -m mode. 'none' keeps records decomposed so the "
    "output stays ingest-ready (vcfclick rejects multi-allelic sites).",
)
def merge(vcfs: tuple[str, ...], output: str, multiallelic: str) -> None:
    """Merge two or more per-sample VCFs into a joint multi-sample VCF.

    Wraps `bcftools merge`. Inputs must be bgzip-compressed with disjoint
    sample names; missing indexes are created automatically. A sample
    absent at a site another sample calls becomes ./. (missing), not 0/0
    — a variant-only VCF does not assert hom-reference where it is silent.
    """
    from ingest.merge import MergeError, merge_vcfs

    if len(vcfs) < 2:
        raise click.ClickException("merge needs at least two input VCFs.")
    try:
        out = merge_vcfs(list(vcfs), output, multiallelic=multiallelic)
    except MergeError as e:
        raise click.ClickException(str(e)) from e
    click.echo(f"merged → {out}")
