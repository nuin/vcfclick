"""`vcfclick db panel` — load a sample -> population panel.

    vcfclick db panel kg integrated_call_samples_v3.20130502.ALL.panel
    vcfclick db panel kg labels.csv --sample-col IID --pop-col group

The panel feeds `vcfclick db popgen`, which groups samples by
`population` (or `super_population`). See ingest/panel.py for the
accepted formats.
"""

from __future__ import annotations

from dataclasses import dataclass

import click

from cli.main import _set_db, db
from cli.options import command_options

_SHOW = 10  # sample ids listed per "not matched" report line


@dataclass(frozen=True)
class PanelOptions:
    name: str
    panel_path: str
    ingest_id: str | None
    sample_col: str | None
    pop_col: str | None
    super_pop_col: str | None
    sex_col: str | None


def _preview(items: list[str]) -> str:
    head = ", ".join(items[:_SHOW])
    more = len(items) - _SHOW
    return head + (f", ... (+{more} more)" if more > 0 else "")


@db.command(name="panel")
@click.argument("name")
@click.argument("panel_path", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--ingest-id",
    default=None,
    help="Apply the panel to this ingestion only (default: every "
    "ingestion that contains each sample).",
)
@click.option("--sample-col", default=None, help="Sample-id column name.")
@click.option("--pop-col", default=None, help="Population column name.")
@click.option(
    "--super-pop-col", default=None, help="Super-population column name (optional)."
)
@click.option("--sex-col", default=None, help="Sex/gender column name (optional).")
@command_options(PanelOptions)
def db_panel(options: PanelOptions) -> None:
    """Load a sample -> population panel into NAME.

    Reads the 1000 Genomes panel format (sample, pop, super_pop, gender;
    tab-separated with a header) as is, or any TSV/CSV with a header —
    columns are matched by name, or named with the --*-col options.
    Re-loading replaces the labels of the samples it lists.
    """
    from ingest.panel import PanelColumns, PanelError, load_panel, parse_panel
    from storage import db_path

    if not db_path(options.name).exists():
        raise click.ClickException(f"db {options.name!r} does not exist.")
    _set_db(options.name)

    try:
        rows = parse_panel(
            options.panel_path,
            PanelColumns(
                sample=options.sample_col,
                population=options.pop_col,
                super_population=options.super_pop_col,
                sex=options.sex_col,
            ),
        )
        report = load_panel(rows, options.ingest_id)
    except (PanelError, ValueError) as e:
        raise click.ClickException(str(e)) from e

    if report.loaded == 0:
        raise click.ClickException(
            f"none of the {len(rows)} panel samples is in the database "
            f"(first: {_preview(report.panel_not_in_db)})"
        )
    for ingest_id, n in sorted(report.by_ingest.items()):
        click.echo(f"loaded panel: {n} samples under ingest_id={ingest_id}")
    if report.panel_not_in_db:
        click.echo(
            f"panel samples not in the database ({len(report.panel_not_in_db)}): "
            f"{_preview(report.panel_not_in_db)}"
        )
    if report.db_not_in_panel:
        missing = [f"{s} ({ing})" for ing, s in report.db_not_in_panel]
        click.echo(
            f"database samples missing from the panel ({len(missing)}): "
            f"{_preview(missing)}"
        )
