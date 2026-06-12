from __future__ import annotations

from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical
from textual.widgets import Button, Footer, Header, Static

from tui.screens.locus import LocusPane
from tui.screens.operations import OperationsPane
from tui.screens.sql import SqlPane


class VcfclickTuiApp(App):
    """Genomics-first terminal UI for vcfclick."""

    CSS = """
    #nav {
        width: 18;
        border-right: solid $surface;
    }
    #main {
        width: 1fr;
    }
    Button {
        width: 100%;
    }
    """

    BINDINGS = [
        ("1", "show_locus", "Locus"),
        ("2", "show_operations", "Operations"),
        ("3", "show_sql", "SQL"),
    ]

    def __init__(self, initial_db: str | None = None) -> None:
        super().__init__()
        self.active_db = initial_db

    def compose(self) -> ComposeResult:
        yield Header()
        with Horizontal():
            with Vertical(id="nav"):
                yield Static("vcfclick")
                yield Button("Locus", id="nav-locus")
                yield Button("Operations", id="nav-operations")
                yield Button("SQL", id="nav-sql")
            with Vertical(id="main"):
                yield LocusPane()
        yield Footer()

    def action_show_locus(self) -> None:
        main = self.query_one("#main", Vertical)
        main.remove_children()
        main.mount(LocusPane())

    def action_show_operations(self) -> None:
        main = self.query_one("#main", Vertical)
        main.remove_children()
        main.mount(OperationsPane())

    def action_show_sql(self) -> None:
        main = self.query_one("#main", Vertical)
        main.remove_children()
        main.mount(SqlPane())

    def open_sql(self, sql: str) -> None:
        self.action_show_sql()
        self.call_after_refresh(self.query_one(SqlPane).set_sql, sql)
