from __future__ import annotations

import asyncio

import pytest

pytest.importorskip("textual")

from tui import services
from tui.app import VcfclickTuiApp
from textual.widgets import Button, Input, Static


def test_app_starts_with_no_active_db():
    async def run() -> None:
        app = VcfclickTuiApp()
        async with app.run_test():
            assert app.active_db is None
            assert app.query_one("#mode-title", Static).content == "Locus"

    asyncio.run(run())


def test_locus_error_clears_prior_sql_handoff(monkeypatch):
    def resolve_locus(parsed: services.ParsedLocus) -> services.ResolvedLocus:
        assert parsed.gene_symbol == "BRCA1"
        return services.ResolvedLocus(
            label="BRCA1",
            chrom="chr17",
            start_pos=43044295,
            end_pos=43125483,
            gene_symbol="BRCA1",
            source="gene",
        )

    def build_locus_summary(
        name: str, locus: services.ResolvedLocus
    ) -> services.LocusSummary:
        assert name == "smoke"
        return services.LocusSummary(
            locus=locus,
            counts=services.QueryResult(sql="", columns=[], rows=[], row_count=0),
            cohorts=services.QueryResult(sql="", columns=[], rows=[], row_count=0),
            quality=services.QueryResult(sql="", columns=[], rows=[], row_count=0),
            preview=services.QueryResult(
                sql="SELECT stale",
                columns=["chrom"],
                rows=[["chr17"]],
                row_count=1,
            ),
        )

    monkeypatch.setattr(services, "resolve_locus", resolve_locus)
    monkeypatch.setattr(services, "build_locus_summary", build_locus_summary)

    async def run() -> None:
        app = VcfclickTuiApp(initial_db="smoke")
        async with app.run_test() as pilot:
            await pilot.click("#locus-input")
            locus_input = app.query_one("#locus-input", Input)
            locus_input.value = "BRCA1"
            await pilot.press("enter")
            await pilot.pause()

            open_sql = app.query_one("#open-locus-sql", Button)
            assert open_sql.disabled is False

            locus_input.value = "chr1"
            await pilot.press("enter")
            await pilot.pause()

            assert open_sql.disabled is True
            assert not hasattr(app.query_one("#locus-input").parent, "_last_sql")

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


def test_operations_screen_lists_databases(vcfclick_home):
    (vcfclick_home / "dbs" / "smoke").mkdir(parents=True)

    async def run() -> None:
        app = VcfclickTuiApp()
        async with app.run_test() as pilot:
            app.action_show_operations()
            await pilot.pause()
            assert "smoke" in app.query_one("#operations-body", Static).content

    asyncio.run(run())


def test_operations_use_db_keeps_active_db_when_summary_fails(monkeypatch):
    monkeypatch.setattr(services, "list_database_names", lambda: ["broken"])
    monkeypatch.setattr(services, "validate_database", lambda name: name)

    def fail_database_summary(name: str) -> services.DatabaseSummary:
        raise services.TuiServiceError("summary failed")

    monkeypatch.setattr(services, "database_summary", fail_database_summary)

    async def run() -> None:
        app = VcfclickTuiApp(initial_db="stable")
        async with app.run_test() as pilot:
            app.action_show_operations()
            await pilot.pause()
            app.query_one("#db-name-input", Input).value = "broken"
            await pilot.click("#use-db")
            await pilot.pause()

            assert app.active_db == "stable"
            assert "summary failed" in app.query_one("#operations-body", Static).content

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
