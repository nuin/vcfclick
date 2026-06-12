from __future__ import annotations

import asyncio

import pytest

pytest.importorskip("textual")

from tui.app import VcfclickTuiApp
from textual.widgets import Static


def test_app_starts_with_no_active_db():
    async def run() -> None:
        app = VcfclickTuiApp()
        async with app.run_test():
            assert app.active_db is None
            assert app.query_one("#mode-title", Static).content == "Locus"

    asyncio.run(run())


def test_open_sql_sets_editor_text():
    async def run() -> None:
        app = VcfclickTuiApp()
        async with app.run_test() as pilot:
            app.open_sql("SELECT 1")
            await pilot.pause()
            editor = app.query_one("#sql-editor")
            assert editor.text == "SELECT 1"

    asyncio.run(run())


def test_locus_submit_renders_missing_db_error():
    async def run() -> None:
        app = VcfclickTuiApp()
        async with app.run_test() as pilot:
            await pilot.click("#locus-input")
            await pilot.press("B", "R", "C", "A", "1")
            await pilot.press("enter")
            await pilot.pause()
            assert "Select a database" in app.query_one("#locus-status", Static).content

    asyncio.run(run())
