from __future__ import annotations

from textual.app import ComposeResult
from textual.containers import Vertical
from textual.widgets import Button, Input, Static

from tui import services


class OperationsPane(Vertical):
    """Database operations for the TUI."""

    def compose(self) -> ComposeResult:
        yield Static("Operations", id="mode-title")
        yield Static(self._database_text(), id="operations-body")
        yield Input(placeholder="database name", id="db-name-input")
        yield Button("Use Database", id="use-db")
        yield Button("Show Stats", id="show-stats")
        yield Static(
            "Ingest launcher: choose a DB first, then use existing CLI for v1.",
            id="ingest-note",
        )

    def _database_text(self) -> str:
        names = services.list_database_names()
        if not names:
            return "No databases found."
        return "Databases:\n" + "\n".join(f"- {name}" for name in names)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        body = self.query_one("#operations-body", Static)
        name = self.query_one("#db-name-input", Input).value.strip()

        if event.button.id == "use-db":
            try:
                validated_name = services.validate_database(name)
                summary = services.database_summary(validated_name)
            except services.TuiServiceError as exc:
                body.update(str(exc))
                return
            self.app.active_db = validated_name
            body.update(
                f"Active DB: {summary.name}\n"
                f"Path: {summary.path}\n"
                f"Variants: {summary.variants}\n"
                f"Genotypes: {summary.genotypes}\n"
                f"Samples: {summary.samples}\n"
                f"Ingestions: {summary.ingestions}"
            )

        if event.button.id == "show-stats":
            active = getattr(self.app, "active_db", None)
            if active is None:
                body.update("Select a database first.")
                return
            try:
                stats = services.stats_summary(active)
            except services.TuiServiceError as exc:
                body.update(str(exc))
                return
            body.update(services.render_stats(stats))
