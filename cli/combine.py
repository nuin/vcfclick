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

from dataclasses import dataclass

import click

from cli.main import cli
from cli.options import command_options


@dataclass(frozen=True)
class CombineOptions:
    vcfs: tuple[str, ...]
    output: str
    names: tuple[str, ...]
    min_callsets: int
    pass_only: bool
    count_by: str
    reference: str | None
    atomize: bool
    carry_info: bool
    genotype_merge: str
    merge_info_with_max_ac: bool
    gatk3: bool
    filtered_are_uncalled: bool
    priority: str | None


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
    help="Keep only sites present in at least this many inputs (consensus filter).",
)
@click.option(
    "--pass-only",
    is_flag=True,
    default=False,
    help="Count only PASS calls toward --min-callsets; a filtered input is still "
    "named filterIn<name> in set= (GATK convention). FILTER '.' counts as PASS.",
)
@click.option(
    "--count-by",
    type=click.Choice(["allele", "site"]),
    default="allele",
    show_default=True,
    help="Count inputs by exact allele (default) or by position. 'site' keeps "
    "every allele at a position counted by enough inputs (GATK --minimumN).",
)
@click.option(
    "--reference",
    type=click.Path(dir_okay=False),
    help="Reference FASTA. Split multi-allelics and left-align/trim internally "
    "(the bcftools norm -m - -f equivalent) instead of refusing. Needs pyfaidx "
    "(vcfclick[benchmark]).",
)
@click.option(
    "--atomize",
    is_flag=True,
    default=False,
    help="Split complex alleles into primitives (SNPs + one indel) after "
    "left-alignment, so a caller packing substitutions into one record "
    "converges with one that emits them separately. Requires --reference.",
)
@click.option(
    "--carry-info",
    is_flag=True,
    default=False,
    help="Carry QUAL, FILTER, and INFO (plus their header lines) from the "
    "highest-priority input that called each allele.",
)
@click.option(
    "--genotype-merge",
    type=click.Choice(["prioritize", "unsorted"]),
    default="prioritize",
    show_default=True,
    help="Genotype of a sample present in several inputs. 'prioritize' takes the "
    "first input with a non-missing call (a ./. does not win). 'unsorted' is GATK3 "
    "UNSORTED: the first input that has the sample wins, even with ./. "
    "(identical to GATK3 PRIORITIZE in input order).",
)
@click.option(
    "--merge-info-with-max-ac",
    is_flag=True,
    default=False,
    help="Take INFO (and with the normalising engine QUAL and FILTER) from the "
    "input with the highest INFO AC for the allele; the first input wins a tie "
    "(GATK3 --mergeInfoWithMaxAC). Implies --carry-info without --gatk3.",
)
@click.option(
    "--gatk3",
    is_flag=True,
    default=False,
    help="Reproduce GATK 3.8-1 CombineVariants record for record: no "
    "normalisation (records merge by position and exact alleles), records "
    "counted and split by type as GATK does, GATK's QUAL/FILTER/ID/INFO/"
    "AC/AF/AN/set rules and header merge. --min-callsets is then GATK "
    "--minimumN (records at the position, minus filtered ones). --reference only "
    "supplies contig order and ##contig lines.",
)
@click.option(
    "--filtered-are-uncalled",
    is_flag=True,
    default=False,
    help="With --gatk3: treat filtered records as uncalled (GATK "
    "-filteredAreUncalled): they are left out of the merge. Filtered records "
    "never count toward --min-callsets, with or without this flag.",
)
@click.option(
    "--priority",
    help="With --gatk3 and --genotype-merge prioritize: comma-separated input "
    "names, highest priority first (GATK -priority). Default: input order.",
)
@command_options(CombineOptions)
def combine(options: CombineOptions) -> None:
    """Combine two or more VCF call sets into one, with set= provenance.

    Inputs are unionized by (chrom, pos, ref, alt); a sample shared across
    inputs is resolved by PRIORITY (input order, highest first). Output
    carries GT + GQ/DP/AD from the priority source, plus a set= INFO field.

    Opt-in GATK CombineVariants parity (defaults keep current behaviour):
    --pass-only, --count-by site, --reference (split + normalize internally,
    so raw multi-allelic/padded inputs combine with no bcftools step), and
    --carry-info (carry QUAL/FILTER/INFO from the priority input).

    --gatk3 reproduces GATK 3.8-1 CombineVariants exactly (see docs/COMBINE.md
    for the equivalent of `--genotypemergeoption UNSORTED
    --mergeInfoWithMaxAC --minimumN 2`).
    """
    from ingest.combine import CombineError, combine_vcfs

    if len(options.vcfs) < 2:
        raise click.ClickException("combine needs at least two input VCFs.")
    try:
        out = combine_vcfs(
            list(options.vcfs),
            options.output,
            names=list(options.names) or None,
            min_callsets=options.min_callsets,
            pass_only=options.pass_only,
            count_by=options.count_by,
            reference=options.reference,
            atomize=options.atomize,
            carry_info=options.carry_info,
            genotype_merge=options.genotype_merge,
            merge_info_with_max_ac=options.merge_info_with_max_ac,
            gatk3=options.gatk3,
            filtered_are_uncalled=options.filtered_are_uncalled,
            priority=options.priority.split(",") if options.priority else None,
        )
    except CombineError as e:
        raise click.ClickException(str(e)) from e
    click.echo(f"combined → {out}")
