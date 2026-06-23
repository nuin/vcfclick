from __future__ import annotations

from textual.app import ComposeResult
from textual.containers import Vertical
from textual.widgets import Button, DataTable, Static, TextArea

from tui import services


class SqlPane(Vertical):
    """SQL editor and result table."""

    BINDINGS = [("ctrl+r", "run", "Run query")]

    def compose(self) -> ComposeResult:
        yield Static("SQL", id="mode-title")
        yield TextArea("", id="sql-editor")
        yield Button("Run (Ctrl+R)", id="sql-run")
        yield Static("", id="sql-status")
        yield DataTable(id="sql-results")

    def set_sql(self, sql: str) -> None:
        self.query_one("#sql-editor", TextArea).text = sql

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "sql-run":
            self.action_run()

    def action_run(self) -> None:
        status = self.query_one("#sql-status", Static)
        active_db = getattr(self.app, "active_db", None)
        if active_db is None:
            status.update("Select a database first (Operations → Use Database).")
            return

        sql = self.query_one("#sql-editor", TextArea).text
        try:
            result = services.execute_sql(active_db, sql)
        except services.TuiServiceError as exc:
            status.update(str(exc))
            return

        status.update(f"{result.row_count} row(s)")
        table = self.query_one("#sql-results", DataTable)
        table.clear(columns=True)
        for column in result.columns:
            table.add_column(column)
        for row in result.rows:
            table.add_row(*[str(cell) for cell in row])
