# Terminal UI

The optional TUI is a Textual interface for exploring local vcfclick
databases without switching between many CLI commands.

Install it with the `tui` extra:

```bash
uv tool install "vcfclick[tui]"
```

From a source checkout:

```bash
uv sync --extra tui --group dev
uv run vcfclick tui --help
```

If Textual is missing, `vcfclick tui` prints an install hint instead of
failing with a Python traceback.

## Start The TUI

Open a specific database:

```bash
vcfclick tui --db demo
```

Or start without an active database and choose one from Operations:

```bash
vcfclick tui
```

Use the same environment variables as the CLI:

```bash
VCFCLICK_HOME=/data/vcfclick VCFCLICK_BACKEND=chdb vcfclick tui --db study
```

## Panes

The current TUI has three panes.

| Key | Pane | Use |
|---|---|---|
| `1` | Locus | Start with a gene symbol or genomic range. |
| `2` | Operations | List databases, select the active database, inspect counts, and run stats where supported. |
| `3` | SQL | Inspect or edit SQL handed off from another pane. |

### Locus

Enter either a gene symbol:

```text
BRCA1
```

or a coordinate range:

```text
chr17:43044295-43125483
```

The pane resolves the locus, runs a summary query against the active
database, and shows a small result preview. Use **Open SQL** to hand the
generated query to the SQL pane.

Gene symbols require the annotation store to be loaded:

```bash
vcfclick annotations load
```

Coordinate ranges do not require annotations.

### Operations

Operations lists databases under `VCFCLICK_HOME`, lets you choose the
active database, and shows basic row counts.

`Show Stats` uses `vcfclick db stats` behavior. It is currently
available for chDB and reports a clear unsupported-feature message on
DuckDB.

Ingest is still CLI-first in the current TUI release:

```bash
vcfclick db create study
vcfclick db ingest study cohort.vcf.gz --cohort cases
```

### SQL

The SQL pane is the handoff point for generated SQL. It currently shows
the query text and result table surface; use the CLI for production
query execution, exports, and scripted workflows:

```bash
vcfclick db query demo "SELECT count() FROM variants"
```

## Troubleshooting

If the TUI says no databases were found, check:

```bash
vcfclick db list
echo "$VCFCLICK_HOME"
```

If gene lookup fails, load annotations:

```bash
vcfclick annotations load
```

If stats fail on DuckDB, switch to chDB for that database or use
`vcfclick db info` for basic counts.
