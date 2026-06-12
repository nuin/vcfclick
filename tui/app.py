"""Minimal Textual application shell for `vcfclick tui`."""

from __future__ import annotations

from textual.app import App, ComposeResult
from textual.widgets import Static


class VcfclickTuiApp(App):
    """Small app shell for the optional terminal UI."""

    def __init__(self, initial_db: str | None = None) -> None:
        super().__init__()
        self.initial_db = initial_db
        self.active_db = initial_db

    def compose(self) -> ComposeResult:
        yield Static("vcfclick")
