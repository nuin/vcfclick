# Genomics-First TUI Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an optional `vcfclick tui` terminal UI whose default workflow is gene/region summary, with operations and SQL available as adjacent modes.

**Architecture:** Add a small `tui/` package with pure service functions plus Textual screens. Keep Textual optional and lazy-loaded from the `vcfclick tui` command so the existing CLI remains lightweight. Service functions return structured data and SQL strings, making the core locus/query behavior testable without a terminal UI.

**Tech Stack:** Python 3.11+, Click, Textual optional extra, pytest, existing `storage` and `annotations` APIs.

---

## File Structure

- Create `tui/__init__.py`: package marker and version-independent exports.
- Create `tui/services.py`: dataclasses, user-facing service errors, DB metadata, locus parsing/resolution, SQL execution, summary query generation.
- Create `tui/app.py`: Textual app shell and mode switching.
- Create `tui/screens/__init__.py`: screen package marker.
- Create `tui/screens/locus.py`: Locus mode widgets and generated SQL handoff.
- Create `tui/screens/operations.py`: Operations mode widgets for DB list/info/path/stats/ingest entry point.
- Create `tui/screens/sql.py`: SQL editor/result table mode.
- Create `cli/tui.py`: Click command registration for `vcfclick tui`; imports Textual lazily.
- Modify `cli/main.py`: side-effect import `cli.tui`.
- Modify `pyproject.toml`: add optional `tui` dependency extra and include `tui` in wheel/sdist packages.
- Create `tests/test_tui_services.py`: service-layer tests that do not require Textual.
- Create `tests/test_tui_cli.py`: CLI lazy-import and missing-dependency behavior.
- Create `tests/test_tui_app.py`: Textual smoke tests skipped when Textual is unavailable.

---

### Task 1: Locus Parser And Service Types

**Files:**
- Create: `tui/__init__.py`
- Create: `tui/services.py`
- Test: `tests/test_tui_services.py`

- [ ] **Step 1: Write failing parser and error tests**

Add this initial test file:

```python
from __future__ import annotations

import pytest

from tui.services import LocusInputError, ParsedLocus, parse_locus_input


def test_parse_range_with_commas():
    locus = parse_locus_input("chr17:43,044,295-43,125,483")
    assert locus == ParsedLocus(
        kind="range",
        label="chr17:43044295-43125483",
        chrom="chr17",
        start_pos=43044295,
        end_pos=43125483,
        gene_symbol=None,
    )


def test_parse_gene_symbol():
    locus = parse_locus_input("BRCA1")
    assert locus == ParsedLocus(
        kind="gene",
        label="BRCA1",
        chrom=None,
        start_pos=None,
        end_pos=None,
        gene_symbol="BRCA1",
    )


@pytest.mark.parametrize(
    "raw",
    ["", "chr17", "chr17:431-430", "chr17:start-end", "BRCA1 BRCA2"],
)
def test_parse_invalid_locus(raw):
    with pytest.raises(LocusInputError):
        parse_locus_input(raw)
```

- [ ] **Step 2: Run parser tests to verify failure**

Run:

```bash
uv run pytest tests/test_tui_services.py -q
```

Expected: FAIL during import because `tui.services` does not exist.

- [ ] **Step 3: Add service dataclasses and parser**

Create `tui/__init__.py`:

```python
"""Terminal UI support for vcfclick."""
```

Create `tui/services.py`:

```python
"""Pure service layer for the vcfclick Textual UI.

This module intentionally imports no Textual symbols. The TUI and tests both
use these functions so parsing, SQL generation, and backend behavior stay
testable without a terminal event loop.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import re
from typing import Any, Literal


class TuiServiceError(Exception):
    """Recoverable user-facing service error."""

    code = "service_error"


class LocusInputError(TuiServiceError):
    """Raised when a gene/range input cannot be interpreted."""

    code = "invalid_locus"


class DatabaseError(TuiServiceError):
    """Raised for missing or invalid active database state."""

    code = "database_error"


class AnnotationUnavailableError(TuiServiceError):
    """Raised when annotation lookup cannot answer a request."""

    code = "annotation_unavailable"


class UnsupportedFeatureError(TuiServiceError):
    """Raised when a backend does not support a TUI operation yet."""

    code = "unsupported_feature"


@dataclass(frozen=True)
class ParsedLocus:
    kind: Literal["gene", "range"]
    label: str
    chrom: str | None
    start_pos: int | None
    end_pos: int | None
    gene_symbol: str | None


@dataclass(frozen=True)
class ResolvedLocus:
    label: str
    chrom: str
    start_pos: int
    end_pos: int
    gene_symbol: str | None = None
    source: Literal["gene", "range"] = "range"


@dataclass(frozen=True)
class QueryResult:
    sql: str
    columns: list[str]
    rows: list[list[Any]]
    row_count: int


@dataclass(frozen=True)
class DatabaseSummary:
    name: str
    path: str
    size_bytes: int
    variants: int | None
    genotypes: int | None
    samples: int | None
    ingestions: int | None


@dataclass(frozen=True)
class LocusSummary:
    locus: ResolvedLocus
    counts: QueryResult
    cohorts: QueryResult
    quality: QueryResult
    preview: QueryResult


_RANGE_RE = re.compile(
    r"^(?P<chrom>[A-Za-z0-9_.-]+):(?P<start>[0-9,]+)-(?P<end>[0-9,]+)$"
)
_GENE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


def parse_locus_input(raw: str) -> ParsedLocus:
    """Parse a user-entered gene symbol or `chrom:start-end` range."""
    text = raw.strip()
    if not text:
        raise LocusInputError("Enter a gene symbol or chrom:start-end range.")

    range_match = _RANGE_RE.fullmatch(text)
    if range_match:
        chrom = range_match.group("chrom")
        start_pos = int(range_match.group("start").replace(",", ""))
        end_pos = int(range_match.group("end").replace(",", ""))
        if start_pos < 1 or end_pos < start_pos:
            raise LocusInputError("Range must use positive coordinates with start <= end.")
        label = f"{chrom}:{start_pos}-{end_pos}"
        return ParsedLocus(
            kind="range",
            label=label,
            chrom=chrom,
            start_pos=start_pos,
            end_pos=end_pos,
            gene_symbol=None,
        )

    if _GENE_RE.fullmatch(text):
        symbol = text.upper()
        return ParsedLocus(
            kind="gene",
            label=symbol,
            chrom=None,
            start_pos=None,
            end_pos=None,
            gene_symbol=symbol,
        )

    raise LocusInputError("Use a gene symbol or range like chr17:43044295-43125483.")
```

- [ ] **Step 4: Run parser tests to verify pass**

Run:

```bash
uv run pytest tests/test_tui_services.py -q
```

Expected: PASS for the parser tests.

- [ ] **Step 5: Commit parser service foundation**

```bash
git add tui/__init__.py tui/services.py tests/test_tui_services.py
git commit -m "Add TUI locus parser services"
```

---

### Task 2: Database Metadata And SQL Execution Services

**Files:**
- Modify: `tui/services.py`
- Modify: `tests/test_tui_services.py`

- [ ] **Step 1: Add failing DB metadata and SQL execution tests**

Append to `tests/test_tui_services.py`:

```python
from storage import apply_schema, get_session
from tui.services import (
    DatabaseError,
    database_summary,
    execute_sql,
    list_database_names,
    validate_database,
)


def test_list_database_names_uses_vcfclick_home(vcfclick_home):
    (vcfclick_home / "dbs" / "alpha").mkdir(parents=True)
    (vcfclick_home / "dbs" / "beta").mkdir(parents=True)
    assert list_database_names() == ["alpha", "beta"]


def test_validate_database_rejects_missing(vcfclick_home):
    with pytest.raises(DatabaseError):
        validate_database("missing")


def test_execute_sql_returns_structured_rows(vcfclick_home):
    name = "smoke"
    validate_home = vcfclick_home / "dbs" / name
    validate_home.mkdir(parents=True)
    get_session(name)
    apply_schema()

    result = execute_sql(name, "SELECT count() AS n FROM variants")
    assert result.sql == "SELECT count() AS n FROM variants"
    assert result.columns == ["n"]
    assert result.rows == [[0]]
    assert result.row_count == 1


def test_database_summary_counts_empty_schema(vcfclick_home):
    name = "smoke"
    (vcfclick_home / "dbs" / name).mkdir(parents=True)
    get_session(name)
    apply_schema()

    summary = database_summary(name)
    assert summary.name == name
    assert summary.variants == 0
    assert summary.genotypes == 0
    assert summary.samples == 0
    assert summary.ingestions == 0
```

- [ ] **Step 2: Run DB service tests to verify failure**

Run:

```bash
VCFCLICK_BACKEND=duckdb uv run pytest tests/test_tui_services.py -q
```

Expected: FAIL because the imported service functions are not defined.

- [ ] **Step 3: Implement DB metadata and SQL execution**

Append these functions to `tui/services.py`:

```python
def list_database_names() -> list[str]:
    """Return local database names sorted by the storage layer."""
    from storage import list_dbs

    return list_dbs()


def validate_database(name: str) -> str:
    """Return `name` when it exists, otherwise raise a recoverable error."""
    from storage import db_path

    path = db_path(name)
    if not path.exists():
        raise DatabaseError(f"Database {name!r} does not exist.")
    return name


def _query_json(name: str, sql: str) -> QueryResult:
    """Execute SQL and normalize chDB/DuckDB JSONCompact output."""
    from storage import get_session

    validate_database(name)
    sess = get_session(name)
    raw = sess.query(sql, "JSONCompact").bytes().decode()
    parsed = json.loads(raw)
    columns = [m["name"] for m in parsed.get("meta", [])]
    rows = parsed.get("data", [])
    return QueryResult(
        sql=sql,
        columns=columns,
        rows=rows,
        row_count=len(rows),
    )


def execute_sql(name: str, sql: str) -> QueryResult:
    """Run user-entered SQL against a named database."""
    cleaned = sql.strip()
    if not cleaned:
        raise TuiServiceError("Enter a SQL query.")
    return _query_json(name, cleaned)


def _scalar_count(name: str, table: str) -> int | None:
    from storage import count_expr

    result = _query_json(name, f"SELECT {count_expr()} AS n FROM {table}")
    if not result.rows:
        return None
    return int(result.rows[0][0])


def database_summary(name: str) -> DatabaseSummary:
    """Return basic DB metadata for the Operations mode."""
    from storage import db_disk_size, db_path

    validate_database(name)
    path = db_path(name)
    counts: dict[str, int | None] = {}
    for table in ("variants", "genotypes", "samples", "ingestions"):
        try:
            counts[table] = _scalar_count(name, table)
        except Exception:
            counts[table] = None

    return DatabaseSummary(
        name=name,
        path=str(path),
        size_bytes=db_disk_size(name),
        variants=counts["variants"],
        genotypes=counts["genotypes"],
        samples=counts["samples"],
        ingestions=counts["ingestions"],
    )
```

- [ ] **Step 4: Run DB service tests to verify pass**

Run:

```bash
VCFCLICK_BACKEND=duckdb uv run pytest tests/test_tui_services.py -q
```

Expected: PASS.

- [ ] **Step 5: Commit DB service layer**

```bash
git add tui/services.py tests/test_tui_services.py
git commit -m "Add TUI database services"
```

---

### Task 3: Gene Resolution And Locus Summary Services

**Files:**
- Modify: `tui/services.py`
- Modify: `tests/test_tui_services.py`

- [ ] **Step 1: Add failing locus resolution and summary tests**

Append to `tests/test_tui_services.py`:

```python
from annotations import GeneRange
from tui.services import (
    AnnotationUnavailableError,
    ResolvedLocus,
    build_locus_summary,
    resolve_locus,
)


def test_resolve_range_without_annotations():
    parsed = parse_locus_input("chr1:100-200")
    assert resolve_locus(parsed) == ResolvedLocus(
        label="chr1:100-200",
        chrom="chr1",
        start_pos=100,
        end_pos=200,
        gene_symbol=None,
        source="range",
    )


def test_resolve_gene_uses_annotation_lookup(monkeypatch):
    def fake_position_for_gene(symbol: str):
        assert symbol == "BRCA1"
        return GeneRange("BRCA1", "chr17", 43044295, 43125483, "-")

    monkeypatch.setattr("annotations.position_for_gene", fake_position_for_gene)

    resolved = resolve_locus(parse_locus_input("BRCA1"))
    assert resolved == ResolvedLocus(
        label="BRCA1",
        chrom="chr17",
        start_pos=43044295,
        end_pos=43125483,
        gene_symbol="BRCA1",
        source="gene",
    )


def test_resolve_gene_not_found(monkeypatch):
    monkeypatch.setattr("annotations.position_for_gene", lambda symbol: None)
    with pytest.raises(AnnotationUnavailableError):
        resolve_locus(parse_locus_input("NOPE1"))


def test_build_locus_summary_returns_sql(vcfclick_home):
    name = "smoke"
    (vcfclick_home / "dbs" / name).mkdir(parents=True)
    get_session(name)
    apply_schema()

    summary = build_locus_summary(name, ResolvedLocus("chr1:1-1000", "chr1", 1, 1000))

    assert summary.locus.chrom == "chr1"
    assert "FROM variants" in summary.counts.sql
    assert "FROM samples" in summary.cohorts.sql
    assert "gq" in summary.quality.sql
    assert "LIMIT 50" in summary.preview.sql
```

- [ ] **Step 2: Run locus summary tests to verify failure**

Run:

```bash
VCFCLICK_BACKEND=duckdb uv run pytest tests/test_tui_services.py -q
```

Expected: FAIL because `resolve_locus` and `build_locus_summary` are missing.

- [ ] **Step 3: Implement gene resolution and summary queries**

Append to `tui/services.py`:

```python
def resolve_locus(parsed: ParsedLocus) -> ResolvedLocus:
    """Resolve parsed user input into concrete coordinates."""
    if parsed.kind == "range":
        assert parsed.chrom is not None
        assert parsed.start_pos is not None
        assert parsed.end_pos is not None
        return ResolvedLocus(
            label=parsed.label,
            chrom=parsed.chrom,
            start_pos=parsed.start_pos,
            end_pos=parsed.end_pos,
            gene_symbol=None,
            source="range",
        )

    assert parsed.gene_symbol is not None
    try:
        import annotations

        gene = annotations.position_for_gene(parsed.gene_symbol)
    except Exception as exc:
        raise AnnotationUnavailableError(
            f"Gene annotations are unavailable: {exc}"
        ) from exc

    if gene is None:
        raise AnnotationUnavailableError(f"Gene {parsed.gene_symbol!r} was not found.")

    return ResolvedLocus(
        label=parsed.gene_symbol,
        chrom=gene.chrom,
        start_pos=int(gene.start_pos),
        end_pos=int(gene.end_pos),
        gene_symbol=gene.gene_symbol,
        source="gene",
    )


def _locus_where(locus: ResolvedLocus, alias: str | None = None) -> str:
    from storage import sql_quote_str

    prefix = f"{alias}." if alias else ""
    return (
        f"{prefix}chrom = {sql_quote_str(locus.chrom)} "
        f"AND {prefix}pos BETWEEN {int(locus.start_pos)} AND {int(locus.end_pos)}"
    )


def _quality_sql(locus: ResolvedLocus) -> str:
    where = _locus_where(locus)
    return (
        "SELECT "
        "count() AS genotype_rows, "
        "count(gq) AS rows_with_gq, "
        "count(dp) AS rows_with_dp "
        f"FROM genotypes WHERE {where}"
    )


def build_locus_summary(name: str, locus: ResolvedLocus) -> LocusSummary:
    """Run the v1 summary query set for a resolved locus."""
    where_v = _locus_where(locus, "v")
    where_g = _locus_where(locus, "g")

    counts_sql = (
        "SELECT "
        "count(DISTINCT (v.ingest_id, v.chrom, v.pos, v.ref, v.alt)) AS variants, "
        "count(DISTINCT (g.ingest_id, g.sample_id)) AS carrier_samples, "
        "count(g.sample_id) AS non_ref_genotype_rows "
        "FROM variants v "
        "LEFT JOIN genotypes g "
        "ON g.ingest_id = v.ingest_id "
        "AND g.chrom = v.chrom AND g.pos = v.pos "
        "AND g.ref = v.ref AND g.alt = v.alt "
        f"WHERE {where_v}"
    )
    cohorts_sql = (
        "SELECT cohort, count(DISTINCT (ingest_id, sample_id)) AS samples "
        "FROM samples GROUP BY cohort ORDER BY samples DESC, cohort"
    )
    preview_sql = (
        "SELECT chrom, pos, ref, alt, vcf_id, qual, filter, info_AF, info_AC "
        "FROM variants "
        f"WHERE {_locus_where(locus)} "
        "ORDER BY chrom, pos, ref, alt LIMIT 50"
    )

    return LocusSummary(
        locus=locus,
        counts=_query_json(name, counts_sql),
        cohorts=_query_json(name, cohorts_sql),
        quality=_query_json(name, _quality_sql(locus)),
        preview=_query_json(name, preview_sql),
    )
```

- [ ] **Step 4: Run locus summary tests to verify pass**

Run:

```bash
VCFCLICK_BACKEND=duckdb uv run pytest tests/test_tui_services.py -q
```

Expected: PASS.

- [ ] **Step 5: Commit locus summary services**

```bash
git add tui/services.py tests/test_tui_services.py
git commit -m "Add TUI locus summary services"
```

---

### Task 4: Optional Dependency And `vcfclick tui` Command

**Files:**
- Modify: `pyproject.toml`
- Create: `cli/tui.py`
- Modify: `cli/main.py`
- Create: `tests/test_tui_cli.py`

- [ ] **Step 1: Add failing CLI registration tests**

Create `tests/test_tui_cli.py`:

```python
from __future__ import annotations

from click.testing import CliRunner

from cli.main import cli


def test_tui_command_is_registered():
    runner = CliRunner()
    result = runner.invoke(cli, ["tui", "--help"])
    assert result.exit_code == 0
    assert "Launch the optional terminal UI" in result.output
    assert "--db" in result.output


def test_tui_command_reports_missing_textual(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "tui.app":
            raise ModuleNotFoundError("No module named 'textual'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    runner = CliRunner()
    result = runner.invoke(cli, ["tui"])

    assert result.exit_code != 0
    assert "Install the TUI extra" in result.output
```

- [ ] **Step 2: Run CLI tests to verify failure**

Run:

```bash
uv run pytest tests/test_tui_cli.py -q
```

Expected: FAIL because `tui` is not a registered command.

- [ ] **Step 3: Add packaging metadata**

Modify `pyproject.toml`:

```toml
[project.optional-dependencies]
tui = [
    "textual>=0.89",
]
```

Also add `"tui"` to `[tool.hatch.build.targets.wheel].packages` and add `"tui"` to `[tool.hatch.build.targets.sdist].include`.

- [ ] **Step 4: Add lazy Click command**

Create `cli/tui.py`:

```python
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
        raise click.ClickException(
            'Install the TUI extra first: pip install "vcfclick[tui]"'
        ) from exc

    VcfclickTuiApp(initial_db=db_name).run()
```

Modify the module registration tuple in `cli/main.py`:

```python
for _module in (
    "cli.db",
    "cli.annotations",
    "cli.discover",
    "cli.tui",
):
    importlib.import_module(_module)
```

- [ ] **Step 5: Run CLI tests to verify pass**

Run:

```bash
uv run pytest tests/test_tui_cli.py -q
```

Expected: PASS.

- [ ] **Step 6: Commit optional command wiring**

```bash
git add pyproject.toml cli/tui.py cli/main.py tests/test_tui_cli.py
git commit -m "Add optional TUI command"
```

---

### Task 5: Textual App Shell And SQL Mode

**Files:**
- Create: `tui/app.py`
- Create: `tui/screens/__init__.py`
- Create: `tui/screens/sql.py`
- Create: `tests/test_tui_app.py`

- [ ] **Step 1: Add failing Textual smoke tests**

Create `tests/test_tui_app.py`:

```python
from __future__ import annotations

import pytest

pytest.importorskip("textual")

from tui.app import VcfclickTuiApp


async def test_app_starts_with_no_active_db():
    app = VcfclickTuiApp()
    async with app.run_test() as pilot:
        assert app.active_db is None
        assert app.query_one("#mode-title").renderable == "Locus"


async def test_open_sql_sets_editor_text():
    app = VcfclickTuiApp()
    async with app.run_test() as pilot:
        app.open_sql("SELECT 1")
        await pilot.pause()
        editor = app.query_one("#sql-editor")
        assert editor.text == "SELECT 1"
```

- [ ] **Step 2: Run app smoke tests to verify failure**

Run:

```bash
uv run pytest tests/test_tui_app.py -q
```

Expected when Textual is installed: FAIL because `tui.app` is missing. Expected when Textual is not installed: SKIP.

- [ ] **Step 3: Add app shell and SQL screen**

Create `tui/screens/__init__.py`:

```python
"""Textual screens for the vcfclick TUI."""
```

Create `tui/screens/sql.py`:

```python
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
```

Create `tui/app.py`:

```python
from __future__ import annotations

from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical
from textual.widgets import Button, Footer, Header, Static

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
                yield Static("Locus", id="mode-title")
                yield Static("Enter a gene or range to begin.", id="empty-state")
        yield Footer()

    def action_show_locus(self) -> None:
        main = self.query_one("#main", Vertical)
        main.remove_children()
        main.mount(Static("Locus", id="mode-title"))
        main.mount(Static("Enter a gene or range to begin.", id="empty-state"))

    def action_show_operations(self) -> None:
        main = self.query_one("#main", Vertical)
        main.remove_children()
        main.mount(Static("Operations", id="mode-title"))
        main.mount(Static("Database operations will appear here.", id="empty-state"))

    def action_show_sql(self) -> None:
        main = self.query_one("#main", Vertical)
        main.remove_children()
        main.mount(SqlPane())

    def open_sql(self, sql: str) -> None:
        self.action_show_sql()
        self.call_after_refresh(self.query_one(SqlPane).set_sql, sql)
```

- [ ] **Step 4: Run app smoke tests**

Run:

```bash
uv run pytest tests/test_tui_app.py -q
```

Expected when Textual is installed: PASS. Expected when Textual is not installed: SKIP.

- [ ] **Step 5: Commit app shell**

```bash
git add tui/app.py tui/screens/__init__.py tui/screens/sql.py tests/test_tui_app.py
git commit -m "Add TUI app shell"
```

---

### Task 6: Locus Screen Integration

**Files:**
- Create: `tui/screens/locus.py`
- Modify: `tui/app.py`
- Modify: `tests/test_tui_app.py`

- [ ] **Step 1: Add failing locus screen smoke test**

Append to `tests/test_tui_app.py`:

```python
async def test_locus_submit_renders_missing_db_error():
    app = VcfclickTuiApp()
    async with app.run_test() as pilot:
        await pilot.click("#locus-input")
        await pilot.press("B", "R", "C", "A", "1")
        await pilot.press("enter")
        await pilot.pause()
        assert "Select a database" in app.query_one("#locus-status").renderable
```

- [ ] **Step 2: Run locus screen test to verify failure**

Run:

```bash
uv run pytest tests/test_tui_app.py -q
```

Expected when Textual is installed: FAIL because `#locus-input` is missing. Expected when Textual is not installed: SKIP.

- [ ] **Step 3: Add locus screen**

Create `tui/screens/locus.py`:

```python
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

    def on_input_submitted(self, event: Input.Submitted) -> None:
        app = self.app
        if getattr(app, "active_db", None) is None:
            self.query_one("#locus-status", Static).update("Select a database first.")
            return

        try:
            parsed = services.parse_locus_input(event.value)
            resolved = services.resolve_locus(parsed)
            summary = services.build_locus_summary(app.active_db, resolved)
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
```

Modify `tui/app.py` to import and mount `LocusPane`:

```python
from tui.screens.locus import LocusPane
```

Replace both places that mount the Locus empty state with:

```python
yield LocusPane()
```

in `compose`, and:

```python
main.mount(LocusPane())
```

in `action_show_locus`.

- [ ] **Step 4: Run locus screen tests**

Run:

```bash
uv run pytest tests/test_tui_app.py -q
```

Expected when Textual is installed: PASS. Expected when Textual is not installed: SKIP.

- [ ] **Step 5: Commit locus screen**

```bash
git add tui/app.py tui/screens/locus.py tests/test_tui_app.py
git commit -m "Add TUI locus screen"
```

---

### Task 7: Operations Screen With DB List, Info, Stats Guard, And Ingest Entry Point

**Files:**
- Create: `tui/screens/operations.py`
- Modify: `tui/app.py`
- Modify: `tui/services.py`
- Modify: `tests/test_tui_services.py`
- Modify: `tests/test_tui_app.py`

- [ ] **Step 1: Add failing stats guard service test**

Append to `tests/test_tui_services.py`:

```python
from tui.services import stats_summary


def test_stats_summary_reports_duckdb_unsupported(monkeypatch, vcfclick_home):
    name = "smoke"
    (vcfclick_home / "dbs" / name).mkdir(parents=True)
    monkeypatch.setenv("VCFCLICK_BACKEND", "duckdb")

    with pytest.raises(UnsupportedFeatureError):
        stats_summary(name)
```

- [ ] **Step 2: Add failing operations screen smoke test**

Append to `tests/test_tui_app.py`:

```python
async def test_operations_screen_lists_databases(vcfclick_home):
    (vcfclick_home / "dbs" / "smoke").mkdir(parents=True)
    app = VcfclickTuiApp()
    async with app.run_test() as pilot:
        app.action_show_operations()
        await pilot.pause()
        assert "smoke" in app.query_one("#operations-body").renderable
```

- [ ] **Step 3: Run tests to verify failure**

Run:

```bash
VCFCLICK_BACKEND=duckdb uv run pytest tests/test_tui_services.py tests/test_tui_app.py -q
```

Expected: service test FAIL because `stats_summary` is missing; app test FAIL/SKIP depending on Textual availability.

- [ ] **Step 4: Implement stats guard and operations screen**

Append to `tui/services.py`:

```python
def stats_summary(name: str, top: int = 20) -> dict[str, Any]:
    """Return stats payload where supported by the active backend."""
    from storage import backend, get_session

    validate_database(name)
    if backend() == "duckdb":
        raise UnsupportedFeatureError("Stats are not implemented on DuckDB yet.")

    from cli.db_stats import _stats_payload

    return _stats_payload(get_session(name), top)
```

Create `tui/screens/operations.py`:

```python
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
        yield Static("Ingest launcher: choose a DB first, then use existing CLI for v1.", id="ingest-note")

    def _database_text(self) -> str:
        names = services.list_database_names()
        if not names:
            return "No databases found."
        return "Databases:\\n" + "\\n".join(f"- {name}" for name in names)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        body = self.query_one("#operations-body", Static)
        name = self.query_one("#db-name-input", Input).value.strip()

        if event.button.id == "use-db":
            try:
                self.app.active_db = services.validate_database(name)
                summary = services.database_summary(name)
            except services.TuiServiceError as exc:
                body.update(str(exc))
                return
            body.update(
                f"Active DB: {summary.name}\\n"
                f"Path: {summary.path}\\n"
                f"Variants: {summary.variants}\\n"
                f"Genotypes: {summary.genotypes}\\n"
                f"Samples: {summary.samples}\\n"
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
            body.update(str(stats))
```

Modify `tui/app.py` to import and mount `OperationsPane`:

```python
from tui.screens.operations import OperationsPane
```

Replace the body of `action_show_operations` with:

```python
def action_show_operations(self) -> None:
    main = self.query_one("#main", Vertical)
    main.remove_children()
    main.mount(OperationsPane())
```

- [ ] **Step 5: Run operations tests**

Run:

```bash
VCFCLICK_BACKEND=duckdb uv run pytest tests/test_tui_services.py tests/test_tui_app.py -q
```

Expected: PASS for service tests; app tests PASS when Textual is installed or SKIP when it is not.

- [ ] **Step 6: Commit operations mode**

```bash
git add tui/services.py tui/screens/operations.py tui/app.py tests/test_tui_services.py tests/test_tui_app.py
git commit -m "Add TUI operations screen"
```

---

### Task 8: Final Verification And Documentation Touches

**Files:**
- Modify: `README.md`
- Modify: `docs/superpowers/specs/2026-06-11-vcfclick-tui-design.md` only if implementation intentionally diverged.

- [ ] **Step 1: Add README TUI section**

Add this section after the existing installation/storage backend material in `README.md`:

```markdown
### Optional terminal UI

vcfclick also has an optional Textual terminal UI for local exploration:

```bash
pip install "vcfclick[tui]"
vcfclick tui --db my-cohort
```

The TUI starts with a genomics-first Locus view: enter a gene symbol such as
`BRCA1` or a range such as `chr17:43044295-43125483`, inspect the summary, then
open the generated SQL when you want to edit or rerun it.
```
```

- [ ] **Step 2: Run focused test suite**

Run:

```bash
VCFCLICK_BACKEND=duckdb uv run pytest tests/test_tui_services.py tests/test_tui_cli.py tests/test_tui_app.py -q
```

Expected: PASS for service/CLI tests; Textual tests PASS if Textual is installed or SKIP if it is not.

- [ ] **Step 3: Run existing CLI smoke tests**

Run:

```bash
VCFCLICK_BACKEND=duckdb uv run pytest tests/test_cli.py -q
```

Expected: PASS. This checks that lazy TUI registration did not break existing `vcfclick db ...` behavior.

- [ ] **Step 4: Run package check**

Run:

```bash
uv build
```

Expected: source distribution and wheel build successfully, with `tui` included in both.

- [ ] **Step 5: Inspect git history and status**

Run:

```bash
git log --oneline -5
git status --short
```

Expected: recent commits have project-focused messages only, and the only remaining untracked files are unrelated pre-existing files such as `assets/`.

- [ ] **Step 6: Commit README and final adjustments**

```bash
git add README.md docs/superpowers/specs/2026-06-11-vcfclick-tui-design.md
git commit -m "Document optional TUI"
```

If the spec did not change, run:

```bash
git add README.md
git commit -m "Document optional TUI"
```

---

## Self-Review Notes

- Spec coverage: the plan covers optional packaging, `vcfclick tui [--db NAME]`, active DB state, gene/range parsing, summary-first locus flow, variant preview, generated SQL handoff, operations DB list/info/path/stats guard, ingest entry-point note, recoverable service errors, and service/Textual test split.
- Deferred by design: full ingest progress UI, natural-language query generation, Harlequin plugin integration, multi-query notebooks, and edit workflows remain out of scope.
- Type consistency: `ParsedLocus`, `ResolvedLocus`, `QueryResult`, `DatabaseSummary`, and `LocusSummary` are introduced in Task 1 and reused consistently in later tasks.
