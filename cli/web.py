"""`vcfclick web` — optional local browser UI for a cohort database.

    vcfclick web epilepsy_2026

Starts a localhost server (the `[web]` extra) and opens a browser to a
SQL explorer, a natural-language→SQL box, and trio / combine panels over
the named database. Everything runs on your machine; nothing is hosted.
"""

from __future__ import annotations

import threading
import webbrowser

import click

from cli.main import _set_db, cli


@cli.command(name="web")
@click.argument("name")
@click.option("--port", default=8765, type=int, show_default=True, help="Port to bind.")
@click.option(
    "--host",
    default="127.0.0.1",
    show_default=True,
    help="Interface to bind. Defaults to localhost; this is a single-user "
    "local tool with no authentication — only change this if you understand "
    "the exposure.",
)
@click.option("--no-browser", is_flag=True, help="Do not open a browser automatically.")
def web_cmd(name: str, port: int, host: str, no_browser: bool) -> None:
    """Launch the optional web UI for database NAME."""
    try:
        import uvicorn

        from vcfclick_web.app import app
    except ModuleNotFoundError as exc:
        if exc.name not in ("fastapi", "uvicorn", "starlette", "pydantic"):
            raise
        raise click.ClickException(
            'Install the web extra first: pip install "vcfclick[web]" '
            '(or: uv tool install "vcfclick[web]").'
        ) from exc

    from storage import db_path

    if not db_path(name).exists():
        raise click.ClickException(
            f"database {name!r} not found. Create it first: vcfclick db create {name}"
        )

    _set_db(name)
    url = f"http://{host}:{port}"
    click.echo(f"vcfclick web → {url}  (database: {name}, Ctrl-C to stop)")
    if not no_browser:
        # Open after a short delay so the server is accepting connections.
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    uvicorn.run(app, host=host, port=port, log_level="warning")
