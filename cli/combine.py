"""`vcfclick combine` — merge multiple VCF call sets with provenance.

The GATK3 CombineVariants functionality GATK4 dropped. Unlike `merge`
(disjoint samples → one joint VCF via bcftools), `combine` unions call
sets that may *share* samples — two callers over the same cohort, or a
pre/post-filter pair — resolving overlaps by priority and recording
where each variant came from:

    vcfclick combine gatk.vcf.gz deepvariant.vcf.gz -o consensus.vcf.gz
    vcfclick combine a.vcf.gz b.vcf.gz c.vcf.gz -o all.vcf.gz --min-callsets 2

Input order is priority: a sample called in more than one input takes
its genotype from the first input that has a non-missing call. Each
output record carries set= naming the inputs it appears in
(Intersection when all). --min-callsets keeps only sites seen in at
least N inputs (consensus calling).
"""

from __future__ import annotations

import click

from cli.main import cli


@cli.command(name="combine")
@click.argument("vcfs", nargs=-1, required=True, type=click.Path(dir_okay=False))
@click.option(
    "-o",
    "--out",
    "output",
    required=True,
    type=click.Path(dir_okay=False),
    help="Output combined VCF path (.vcf or .vcf.gz).",
)
@click.option(
    "--name",
    "names",
    multiple=True,
    help="Set name for an input (repeat once per input, in order). "
    "Default: derived from each filename.",
)
@click.option(
    "--min-callsets",
    type=int,
    default=1,
    show_default=True,
    help="Keep only sites present in at least this many inputs "
    "(consensus filter).",
)
def combine(
    vcfs: tuple[str, ...],
    output: str,
    names: tuple[str, ...],
    min_callsets: int,
) -> None:
    """Combine two or more VCF call sets into one, with set= provenance.

    Inputs are unioned by (chrom, pos, ref, alt) and must be decomposed
    (one ALT per record). A sample shared across inputs is resolved by
    PRIORITY — input order is highest-first. Output carries GT plus a
    set= INFO field; per-sample FORMAT fields (GQ/DP/AD) are not yet
    propagated.
    """
    from ingest.combine import CombineError, combine_vcfs

    if len(vcfs) < 2:
        raise click.ClickException("combine needs at least two input VCFs.")
    try:
        out = combine_vcfs(
            list(vcfs),
            output,
            names=list(names) or None,
            min_callsets=min_callsets,
        )
    except CombineError as e:
        raise click.ClickException(str(e)) from e
    click.echo(f"combined → {out}")
