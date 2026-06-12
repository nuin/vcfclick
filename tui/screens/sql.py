from __future__ import annotations

from textual.app import ComposeResult
from textual.containers import Vertical
from textual.widgets import DataTable, Static, TextArea


class SqlPane(Vertical):
    """SQL editor and result table."""

    def compose(self) -> ComposeResult:
        yield Static("SQL", id="mode-title")
        yield TextArea("", id="sql-editor")
        yield DataTable(id="sql-results")

    def set_sql(self, sql: str) -> None:
        self.query_one("#sql-editor", TextArea).text = sql
