from __future__ import annotations

from textual.app import ComposeResult
from textual.containers import Vertical
from textual.widgets import Button, DataTable, Input, Static

from tui import services


class LocusPane(Vertical):
    """Gene/region summary home view."""

    def compose(self) -> ComposeResult:
        yield Static("Locus", id="mode-title")
        yield Input(placeholder="BRCA1 or chr17:43044295-43125483", id="locus-input")
        yield Static("", id="locus-status")
        yield Button("Open SQL", id="open-locus-sql", disabled=True)
        yield DataTable(id="locus-preview")

    def _reset_handoff(self) -> None:
        self.query_one("#open-locus-sql", Button).disabled = True
        self.query_one("#locus-preview", DataTable).clear(columns=True)
        if hasattr(self, "_last_sql"):
            del self._last_sql

    def on_input_submitted(self, event: Input.Submitted) -> None:
        self._reset_handoff()
        active_db = getattr(self.app, "active_db", None)
        if active_db is None:
            self.query_one("#locus-status", Static).update("Select a database first.")
            return

        try:
            parsed = services.parse_locus_input(event.value)
            resolved = services.resolve_locus(parsed)
            summary = services.build_locus_summary(active_db, resolved)
        except services.TuiServiceError as exc:
            self.query_one("#locus-status", Static).update(str(exc))
            return

        self.query_one("#locus-status", Static).update(
            f"{summary.locus.label}: {summary.locus.chrom}:"
            f"{summary.locus.start_pos}-{summary.locus.end_pos}"
        )
        table = self.query_one("#locus-preview", DataTable)
        table.clear(columns=True)
        for column in summary.preview.columns:
            table.add_column(column)
        for row in summary.preview.rows:
            table.add_row(*[str(cell) for cell in row])
        self._last_sql = summary.preview.sql
        self.query_one("#open-locus-sql", Button).disabled = False

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "open-locus-sql" and hasattr(self, "_last_sql"):
            self.app.open_sql(self._last_sql)
