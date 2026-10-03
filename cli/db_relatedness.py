"""`vcfclick db relatedness` — pairwise kinship between samples (KING-robust).

For every pair of samples in an ingestion, the KING-robust kinship
coefficient (Manichaikul et al. 2010, Bioinformatics 26:2867):

    phi_ij = (N_het,het - 2 * N_ibs0) / (N_het_i + N_het_j)

counted over biallelic SNVs where both samples are called. N_ibs0 is the
number of sites where one sample is hom-ref and the other hom-alt; it is
~0 for parent/child and clearly non-zero for full siblings, which is how
the two first-degree relationships are told apart.

Expected values: duplicate / MZ twin 0.5, first degree 0.25, second 0.125,
third 0.0625, unrelated ~0. Pairs are binned with KING's thresholds.

Genotypes come from the sparse `genotypes` table, so how an absent row is
read matters:
- ingested with --keep-reference: 0/0 calls are stored, so absent means a
  no-call and the site is skipped for that pair;
- otherwise only non-reference calls are stored and absent is read as 0/0
  (a no-call cannot be told apart). Missing calls then look like hom-ref,
  which slightly lowers kinship. Joint-called VCFs have few, so the
  estimate is still reliable for relationship checks.

The computation is dense matrix products (numpy), one ingestion at a time,
since only samples genotyped together share a site list.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import click

from cli.main import _set_db, db
from cli.options import command_options
from cli.relatedness_analysis import (
    _THIRD,
    _rows,
)
from cli.relatedness_analysis import classify as classify
from cli.relatedness_analysis import king_robust as king_robust
from cli.relatedness_analysis import relatedness as relatedness


@dataclass(frozen=True)
class RelatednessOptions:
    name: str
    ingest_id: str | None
    min_kinship: float
    include_all: bool
    max_sites: int
    force: bool
    out_format: str


@db.command(name="relatedness")
@click.argument("name")
@click.option(
    "--ingest-id",
    default=None,
    help="Only this ingestion (default: each ingestion separately).",
)
@click.option(
    "--min-kinship",
    type=float,
    default=_THIRD,
    show_default=True,
    help="Report pairs at or above this kinship (third degree by default).",
)
@click.option(
    "--all",
    "include_all",
    is_flag=True,
    help="Report every pair, including unrelated ones.",
)
@click.option(
    "--max-sites",
    type=int,
    default=200000,
    show_default=True,
    help="Thin to at most this many SNVs (evenly spread) to bound memory.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Classify pairs even when the data covers too small a region.",
)
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["table", "json"]),
    default="table",
    show_default=True,
)
@command_options(RelatednessOptions)
def db_relatedness(options: RelatednessOptions) -> None:
    """Pairwise kinship (KING-robust): duplicates, relatives, pedigree errors."""
    from storage import db_path, get_session

    if not db_path(options.name).exists():
        raise click.ClickException(f"db {options.name!r} does not exist.")
    _set_db(options.name)
    sess = get_session(options.name)

    ingests = [
        r[0]
        for r in _rows(
            sess, "SELECT ingest_id FROM ingestions ORDER BY ingested_at, ingest_id"
        )
    ]
    if options.ingest_id is not None:
        if options.ingest_id not in ingests:
            raise click.ClickException(
                f"no ingestion {options.ingest_id!r} in {options.name!r}; have: {', '.join(ingests) or 'none'}"
            )
        ingests = [options.ingest_id]
    results = [
        relatedness(
            sess,
            i,
            max_sites=options.max_sites,
            include_all=options.include_all,
            min_kinship=options.min_kinship,
            force=options.force,
        )
        for i in ingests
    ]

    if options.out_format == "json":
        click.echo(json.dumps(results, indent=2))
        return
    for res in results:
        _print_table(res, options.include_all)


def _print_table(res: dict, include_all: bool) -> None:
    click.echo(
        f"ingestion {res['ingest_id']}: {res['n_samples']:,} samples, "
        f"{res['sites_used']:,} SNVs"
        + (
            f" (thinned from {res['n_sites']:,})"
            if res["sites_used"] < res["n_sites"]
            else ""
        )
    )
    if res["warning"]:
        click.echo(f"  WARNING: {res['warning']}")
    if res["mode"] == "sparse":
        click.echo(
            "  absent genotypes read as 0/0 (this ingestion stored only non-reference calls)"
        )
    if not res["pairs"]:
        click.echo("  no related pairs" if not include_all else "  no pairs")
    else:
        click.echo(
            f"  {'sample_a':<16}{'sample_b':<16}{'kinship':>9}{'ibs0':>8}{'sites':>9}  relationship"
        )
        for p in res["pairs"]:
            k = "n/a" if p["kinship"] is None else f"{p['kinship']:.3f}"
            b = "n/a" if p["ibs0"] is None else f"{p['ibs0']:.3f}"
            click.echo(
                f"  {p['sample_a']:<16}{p['sample_b']:<16}{k:>9}{b:>8}{p['n_sites']:>9,}  {p['relationship']}"
            )
    if res["pedigree_checks"]:
        click.echo("  pedigree:")
        for c in res["pedigree_checks"]:
            k = "n/a" if c["kinship"] is None else f"{c['kinship']:.3f}"
            mark = "" if c["verdict"] == "ok" else "  <-- " + c["verdict"]
            click.echo(
                f"    {c['parent']} ({c['role']}) of {c['child']}: kinship {k}, "
                f"observed {c['observed'] or 'n/a'}{mark}"
            )
    click.echo("")
