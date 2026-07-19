"""`vcfclick benchmark` — normalized genotype concordance of a query VCF against
a truth VCF over a reference FASTA.

Thin click verb: validate flags, delegate to `benchmark.pipeline.run_benchmark`,
surface BenchmarkError/UnsupportedFeatureError as clean ClickExceptions.

    vcfclick benchmark --truth truth.vcf.gz --query query.vcf.gz \\
        --ref hg38.fa --regions conf.bed -o out/
"""

from __future__ import annotations

import click

from cli.main import cli

_ALL_FORMATS = ("csv", "json", "html")


@cli.command(name="benchmark")
@click.option(
    "--truth",
    required=True,
    type=click.Path(dir_okay=False),
    help="Truth (gold-standard) VCF.",
)
@click.option(
    "--query",
    required=True,
    type=click.Path(dir_okay=False),
    help="Query (call set under test) VCF.",
)
@click.option(
    "--ref",
    required=True,
    type=click.Path(dir_okay=False),
    help="Reference FASTA (indexed .fai).",
)
@click.option(
    "--regions",
    type=click.Path(dir_okay=False),
    help="Confident-region BED. Omit to treat every call as confident.",
)
@click.option(
    "-o",
    "--output",
    "output",
    required=True,
    type=click.Path(file_okay=False),
    help="Output directory for reports.",
)
@click.option(
    "--engine",
    type=click.Choice(["normalized", "haplotype"]),
    default="normalized",
    show_default=True,
    help="Reconciliation engine (haplotype is a P2 feature).",
)
@click.option(
    "--report-formats",
    default="all",
    show_default=True,
    help="Comma-separated subset of csv,json,html (or 'all').",
)
@click.option(
    "--on-ref-mismatch",
    type=click.Choice(["error", "skip"]),
    default="error",
    show_default=True,
    help="Behaviour when a REF allele disagrees with the reference.",
)
def benchmark(
    truth: str,
    query: str,
    ref: str,
    regions: str | None,
    output: str,
    engine: str,
    report_formats: str,
    on_ref_mismatch: str,
) -> None:
    """Benchmark a query VCF against a truth VCF (normalized genotype concordance)."""
    from benchmark.pipeline import run_benchmark
    from benchmark.reconcile import UnsupportedFeatureError
    from benchmark.reference import BenchmarkError

    if report_formats.strip().lower() == "all":
        formats = list(_ALL_FORMATS)
    else:
        formats = [f.strip() for f in report_formats.split(",") if f.strip()]
        unknown = [f for f in formats if f not in _ALL_FORMATS]
        if unknown:
            raise click.ClickException(
                f"unknown report format(s): {', '.join(unknown)}"
            )

    try:
        res = run_benchmark(
            truth,
            query,
            ref,
            output,
            regions=regions,
            engine=engine,
            report_formats=formats,
            on_ref_mismatch=on_ref_mismatch,
        )
    except (BenchmarkError, UnsupportedFeatureError) as e:
        raise click.ClickException(str(e)) from e

    for row in res["summary"]:
        click.echo(
            f"{row['Type']:5} {row['Filter']:4}  "
            f"recall={row['recall']:.4f} precision={row['precision']:.4f} "
            f"f1={row['f1']:.4f}"
        )
    click.echo(f"reports → {output}")
