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

_ALL_FORMATS = ("csv", "json", "parquet", "html")


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
    type=click.Choice(["normalized", "haplotype", "exact"]),
    default="normalized",
    show_default=True,
    help="Reconciliation engine. 'haplotype' adds local-haplotype rescue; "
    "'exact' is a diagnostic that skips reference left-alignment.",
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
@click.option(
    "--conf-containment",
    type=click.Choice(["start", "full"]),
    default="start",
    show_default=True,
    help="Confident-region membership: variant start, or whole ref span.",
)
@click.option(
    "--decompose-mnp/--no-decompose-mnp",
    default=False,
    show_default=True,
    help="Atomize MNPs into per-position SNPs (loses phase; off by default).",
)
@click.option(
    "--strict",
    is_flag=True,
    default=False,
    help="Promote warnings (e.g. missing --regions) to hard errors.",
)
@click.option(
    "--pass-only/--all",
    "pass_only",
    default=None,
    help="Show only the PASS (or only the ALL) filter view in the headline.",
)
@click.option(
    "--stratify",
    default=None,
    help="Comma-separated concordance stratification axes (gnomad,clinvar,gene) "
    "joined against the annotation store; writes stratified_<axis>.csv.",
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
    conf_containment: str,
    decompose_mnp: bool,
    strict: bool,
    pass_only: bool | None,
    stratify: str | None,
) -> None:
    """Benchmark a query VCF against a truth VCF (normalized genotype concordance)."""
    from benchmark.pipeline import run_benchmark
    from benchmark.reconcile import UnsupportedFeatureError
    from benchmark.reference import BenchmarkError

    strat = [a.strip() for a in stratify.split(",") if a.strip()] if stratify else None

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
            conf_containment=conf_containment,
            decompose_mnp=decompose_mnp,
            strict=strict,
            stratify=strat,
        )
    except (BenchmarkError, UnsupportedFeatureError) as e:
        raise click.ClickException(str(e)) from e

    rows = res["summary"]
    if pass_only is True:
        rows = [r for r in rows if r["Filter"] == "PASS"]
    elif pass_only is False:
        rows = [r for r in rows if r["Filter"] == "ALL"]
    for row in rows:
        click.echo(
            f"{row['Type']:5} {row['Filter']:4}  "
            f"recall={row['recall']:.4f} precision={row['precision']:.4f} "
            f"f1={row['f1']:.4f}"
        )
    click.echo(f"reports → {output}")


@cli.command(name="benchmark-cohort")
@click.option(
    "--truth", required=True, type=click.Path(dir_okay=False), help="Truth VCF."
)
@click.option(
    "--ref", required=True, type=click.Path(dir_okay=False), help="Reference FASTA."
)
@click.option(
    "--regions", type=click.Path(dir_okay=False), help="Confident-region BED."
)
@click.option(
    "--caller",
    "callers",
    multiple=True,
    required=True,
    metavar="NAME=QUERY.vcf",
    help="A named query call set; repeat for each caller.",
)
@click.option(
    "-o", "--output", "output", required=True, type=click.Path(file_okay=False)
)
@click.option(
    "--engine",
    type=click.Choice(["normalized", "haplotype", "exact"]),
    default="haplotype",
    show_default=True,
)
@click.option(
    "--on-ref-mismatch", type=click.Choice(["error", "skip"]), default="error"
)
@click.option(
    "--history",
    type=click.Path(dir_okay=False),
    help="Append per-caller metrics to this history CSV.",
)
@click.option(
    "--label",
    default="run",
    show_default=True,
    help="Label for the history row (e.g. a pipeline version).",
)
def benchmark_cohort(
    truth: str,
    ref: str,
    regions: str | None,
    callers: tuple[str, ...],
    output: str,
    engine: str,
    on_ref_mismatch: str,
    history: str | None,
    label: str,
) -> None:
    """Benchmark several callers against one truth; report per-caller concordance
    and each caller's relative blind spots (variants others catch)."""
    import csv
    import os

    from benchmark.cohort import (
        append_run,
        benchmark_callers,
        per_caller_metrics,
        variants_missed_by,
    )
    from benchmark.reference import BenchmarkError

    caller_map: dict[str, str] = {}
    for spec in callers:
        if "=" not in spec:
            raise click.ClickException(f"--caller must be NAME=path, got {spec!r}")
        name, path = spec.split("=", 1)
        caller_map[name.strip()] = path.strip()

    os.makedirs(output, exist_ok=True)
    try:
        frame = benchmark_callers(
            truth,
            ref,
            caller_map,
            regions=regions,
            engine=engine,
            on_ref_mismatch=on_ref_mismatch,
        )
    except BenchmarkError as e:
        raise click.ClickException(str(e)) from e

    metrics = per_caller_metrics(frame)
    with open(os.path.join(output, "per_caller_metrics.csv"), "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=list(metrics[0].keys())
            if metrics
            else [
                "caller",
                "vtype",
                "truth_tp",
                "truth_fn",
                "query_tp",
                "query_fp",
                "recall",
                "precision",
            ],
        )
        w.writeheader()
        w.writerows(metrics)

    for row in metrics:
        click.echo(
            f"{row['caller']:12} {row['vtype']:5}  "
            f"recall={row['recall']:.4f} precision={row['precision']:.4f}"
        )
    for name in caller_map:
        missed = variants_missed_by(frame, name, min_others=1)
        if missed:
            click.echo(f"{name}: misses {len(missed)} variant(s) other callers catch")

    if history:
        append_run(history, label, metrics)
        click.echo(f"history → {history} (label={label})")
    click.echo(f"reports → {output}")
