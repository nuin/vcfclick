"""Register the optional `vcfclick tui` command."""

from __future__ import annotations

import click

from cli.main import cli


@cli.command(name="tui")
@click.option("--db", "db_name", default=None, help="Database to open initially.")
def tui_cmd(db_name: str | None) -> None:
    """Launch the optional terminal UI."""
    try:
        from tui.app import VcfclickTuiApp
    except ModuleNotFoundError as exc:
        if exc.name != "textual":
            raise
        raise click.ClickException(
            'Install the TUI extra first: pip install "vcfclick[tui]"'
        ) from exc

    VcfclickTuiApp(initial_db=db_name).run()
