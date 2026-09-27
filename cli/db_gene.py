"""`vcfclick db gene` — every variant in a gene, with carrier counts."""

from __future__ import annotations

import json

import click

from cli.main import _set_db, db


@db.command(name="gene")
@click.argument("name")
@click.argument("symbol")
@click.option(
    "--flank",
    type=int,
    default=0,
    show_default=True,
    help="Extend the gene by this many bp on each side.",
)
@click.option(
    "--limit",
    type=int,
    default=500,
    show_default=True,
    help="At most this many variants.",
)
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["table", "json", "tsv"]),
    default="table",
    show_default=True,
)
def db_gene(name: str, symbol: str, flank: int, limit: int, out_format: str) -> None:
    """Variants in a gene (by HGNC symbol), with carrier and hom-alt counts."""
    from storage import db_path, get_session
    from storage.gene_query import GeneNotFound, gene_variants

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} does not exist.")
    _set_db(name)
    try:
        res = gene_variants(get_session(name), symbol, flank=flank, limit=limit)
    except GeneNotFound:
        raise click.ClickException(
            f"gene {symbol.upper()!r} not found in the annotation store. "
            "Load GENCODE gene coordinates first with `vcfclick annotations load`."
        ) from None

    if out_format == "json":
        click.echo(json.dumps(res, indent=2))
        return
    cols = (
        [c for c in res["columns"] if c != "ingest_id"]
        if len({r["ingest_id"] for r in res["variants"]}) <= 1
        else res["columns"]
    )
    if out_format == "tsv":
        click.echo("\t".join(cols))
        for r in res["variants"]:
            click.echo("\t".join("" if r[c] is None else str(r[c]) for c in cols))
        return

    g, reg = res["gene"], res["region"]
    click.echo(
        f"{g['gene_symbol']} {reg['chrom']}:{reg['start']}-{reg['end']}"
        + (f" (gene {g['start_pos']}-{g['end_pos']}, +/-{flank} bp)" if flank else "")
        + f": {res['row_count']:,} variant{'s' if res['row_count'] != 1 else ''}"
        + (f", first {limit:,} shown" if res["truncated"] else "")
    )
    if not res["variants"]:
        return
    click.echo(
        f"{'chrom':<7}{'pos':>12}  {'ref':<8}{'alt':<8}{'qual':>8}  {'filter':<8}{'af':>9}{'carriers':>10}{'hom_alt':>9}"
    )
    for r in res["variants"]:
        af = "" if r["info_AF"] is None else f"{r['info_AF']:.4g}"
        qual = "" if r["qual"] is None else f"{r['qual']:.1f}"
        click.echo(
            f"{r['chrom']:<7}{r['pos']:>12,}  {r['ref'][:8]:<8}{r['alt'][:8]:<8}{qual:>8}  "
            f"{(r['filter'] or '.'):<8}{af:>9}{r['carriers']:>10,}{r['hom_alt']:>9,}"
        )
