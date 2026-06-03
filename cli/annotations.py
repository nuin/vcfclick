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
