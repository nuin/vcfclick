# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

vcfclick is a Python CLI that turns VCF cohorts into local SQL databases (one named DB per cohort) on an embedded chDB (ClickHouse) or DuckDB backend, with trio/QC/popgen/combine/benchmark analyses, a Textual TUI, a FastAPI web UI and an MCP server. Research preview, not clinical.

## Commands

```bash
uv sync --extra tui --extra web --extra benchmark --group dev   # dev setup (CI uses web+benchmark+dev)
uv run vcfclick --help

uv run pytest tests/                       # full suite (several minutes; chDB tests run in subprocesses)
uv run pytest tests/test_cli.py            # one file
uv run pytest tests/test_trio.py::test_denovo_excludes_no_call_parent_site   # one test
uv run pytest tests/ -k discover           # by name
VCFCLICK_BACKEND=duckdb uv run pytest tests/   # DuckDB backend (CI runs both)
uv run pytest tests/ -m "not timing"     # skip the wall-clock gates (CI runs them in one `timing` job)
uv run pytest tests/ -m "not timing" -n auto   # in parallel (pytest-xdist), as CI does; keep timing runs serial

uvx ruff@0.15.16 format --check .          # CI pins this ruff version
uvx ruff@0.15.16 check .
uvx ruff@0.15.16 format .
```

Tests need `bgzip`, `tabix` and `bcftools` on `$PATH` (`brew install htslib bcftools`); some tests skip without them.

Run the CLI via the entry-point script (`uv run vcfclick` / `.venv/bin/vcfclick`), never `python -m cli.main`: the module gets loaded twice (`__main__` and `cli.main`), subcommands attach to the other instance's groups, and every command reports "No such command".

## Architecture

Top-level packages (each is a wheel package listed in `pyproject.toml`; add new ones to both the wheel and sdist lists):

- `cli/` — Click entry point `cli.main:cli`. Commands register by **side-effect import**: `cli/main.py` imports `cli.db`, `cli.annotations`, …; `cli/db.py` in turn imports the `cli/db_*.py` modules, each adding `@db.command`s. A new subcommand module must be added to one of those import lists. Commands are kept in small focused modules; heavy analysis logic lives in `*_analysis.py` / `popgen/`, not the command file.
- Filters and presets: `cli/filter_registry.py` is the only place a `db trio` filter is declared (file key, Click parameter, flag, type, range, default, owner, commands); the Click options come from it (`gate_options`) and `tests/test_filter_registry_drift.py` fails when an option is added without being registered or listed in `NEVER_PRESET`. A step that adds a filter or a model registers it in its own change; a built-in preset (`cli/builtin_presets.py`) only changes with a revision bump, a history entry and a CHANGELOG line, and a new inheritance model fails `tests/test_builtin_presets.py` until each built-in runs or excludes it. A value that came from a preset is never "explicit" (`is_explicit`), and a run with a preset prints which gates actually bind (`cli/effective_gates.py`; a GQ / DP / AD gate binds on the CANDIDATE ROWS of each selected model (`collect_candidate_stats`: one grouped count over the join `trio_sql` runs, from the shared `MODEL_PREDICATES` plus the SQL rarity and FILTER gates (`site_gate_sql`: a site they remove is no candidate; the Python gnomAD filter cannot narrow the pool, `strict_missing_dropped` backstops it), before any quality gate; of those rows the genotype classes the gate reads, `gate_reads`, the weakest model deciding), not on the share over all of a sample's rows or of one genotype class, which fall back only for a model with no candidate rows; gnomAD binds only when a probe of the cohort's own variants finds records, not when the store merely has a row, and `gnomad.bound` is null for a store nobody asked to consult; `pass_only` binds only where some variant has a FILTER label; a requested gate that no selected model reads is `bound: false` with a `gate_not_applicable` warning (never a refusal); the INFO/AF gate below 100% of the variants is a `gate_partial`; `thin_coverage` is the median depth of the proband's non-reference calls; a selected model that cannot run is a `model_blocked` warning and a refusal for `require_bound`). `strict_missing` is held against a gate only when that gate is asked for (`Gates.strict_for`), and a `--no-require-bound-gates` run exempts the gates that cannot bind (`Gates.strict_exempt`) instead of failing every call, and after the run the candidates it dropped that pass every other gate, only for a missing field, are counted and warned (`strict_missing_lost`: the tolerant count less the strict one); a name is always checked with `name_problem` before it touches the file system (`find_preset`, `user_path`), and a stored preset is called what its file is called (`read_preset_file` rewrites a declared name that differs, so every writer writes the file name). A report limited to some gates (`build_report(only=...)`, a typed `--parent-alt-screen` or `--pass-only`) says `null`, not true or false, for the `bound` flags and `population_filter` it did not examine, and speaks for a model that cannot run only when the typed gate is about it (`_GATE_MODELS`: `--parent-alt-screen` and de novo).
- `storage/` — backend abstraction. `storage.backend()` picks `chdb` or `duckdb` from `VCFCLICK_BACKEND` (auto-detects chDB if importable). All backend-specific SQL goes through helpers in `storage/sql.py` (`count_expr`, `delete_where_sql`, `parquet_file_expr`, `sql_quote_str`, …) — never hard-code ClickHouse- or DuckDB-only syntax in callers. `storage/_chdb.py` / `_duckdb.py` hold the session implementations. Paths (`VCFCLICK_HOME`, `DB_ROOT`) are recomputed on every access via module `__getattr__` so env changes after import take effect; keep it that way.
- `schema/*.sql` (ClickHouse DDL) and `schema/duckdb/*.sql` (DuckDB DDL) must stay in sync with each other **and** with the Arrow schemas in `ingest/_arrow.py`, column-by-column and in order — `tests/test_schema_agreement.py` enforces this. Schema evolution of existing DBs goes through `storage.db.upgrade_schema()` (`_ADDED_VARIANT_COLUMNS`, `_ADDED_TABLE_FILES`): additive only, never rewrite or drop data.
- `ingest/` — VCF → cyvcf2 rows → Parquet batches → `INSERT … SELECT FROM <parquet_file_expr>`. Serial (`vcf_load.py`) and parallel (`parallel*.py`: workers only write Parquet, have no storage dependency; main process does the import) share that path. `ingest/routing.py` is the single source of truth for which INFO/FORMAT fields get typed columns vs. the `info_extra` / `format_extra` Map overflow. Ingests are atomic per `ingest_id` (lock + `rollback_ingest`). Input is expected to be multi-allelic-split (`bcftools norm -m -`).
- `annotations/` — shared DuckDB reference store (GENCODE genes/transcripts, ClinVar, gnomAD) at `annotations/annotations.duckdb` inside the package dir, independent of the variant backend.
- `popgen/`, `benchmark/`, `export/` — analysis and output layers built on SQL over the cohort tables.
- `tui/`, `vcfclick_web/`, `vcfclick_mcp/` — front-ends. Web and MCP both run caller/LLM-supplied SQL only after `storage.sql_guard.is_read_only` (sqlglot AST check, fails closed); keep any new SQL-executing surface behind it.

### Data-model rules that matter when writing SQL

- `genotypes` is **sparse**: hom-ref (`0/0`) calls are not stored; full no-calls go to `missing_genotypes`. Per-site called counts are `variants.n_called` / `an_called` / `ac_called` (NULL on DBs ingested before those columns existed). See `docs/SCHEMA.md#common-query-patterns` before writing AF or hom-ref logic.
- Keys include `ingest_id`: `variants` is one row per `(ingest_id, chrom, pos, ref, alt)`, `samples` per `(ingest_id, sample_id)`.

## Testing constraints

- Every test gets an isolated `VCFCLICK_HOME` (`vcfclick_home` fixture) and, via a session fixture, an empty annotation store (`VCFCLICK_ANNOTATIONS_DB`), never the developer's `annotations/annotations.duckdb`; tests that seed annotations use `isolated_annotation_db`. Tests must stay independent: CI runs them in parallel.
- **Never open a chDB session inside the pytest process.** chDB allows one embedded server per process bound to the first path; an autouse fixture in `tests/conftest.py` fails any test that leaks one. Exercise storage through the CLI (`run_cli(home, backend, *args)` in conftest) or the `run_python` fixture (fresh subprocess).
- chDB native crashes/hangs are a known CI flake; pytest is configured to rerun only those signatures (`pyproject.toml` addopts) with a 180s per-test timeout. A real assertion failure still fails.
- Fixtures live in `tests/fixtures/` (`tiny.vcf.gz` is the 5-variant/3-sample baseline; `bgzip_vcf()` builds indexed VCFs from plain-text fixtures).

## Dependency pins (don't loosen casually)

- macOS chDB is pinned to `chdb==4.2.0` + `chdb-core==26.5.0` because other builds fail to load on macOS 26+/Darwin 27 (misaligned LINKEDIT); see `docs/BACKENDS.md`.
- `duckdb<2`: DuckDB 2.0 files can't be read by 1.x, which would break shared `db push`/`db pull` bundles.

## Conventions

- Commit messages: noun-first, present tense, brief subject; body wrapped at ~72 chars.
- Releases: update `CHANGELOG.md` (`[Unreleased]` → version), bump `version` in `pyproject.toml` and the matching line in `uv.lock`, tag `v0.x.y` from `main`; `release.yml` publishes to PyPI and creates the GitHub release. If a newer local uv bumps the lockfile format revision, keep only the version line. Don't touch the bioconda recipe — BiocondaBot updates it after each PyPI release.
- This is a public repo. Never commit private working docs (`PLAN.md`, strategy/roadmap/notes drafts, `workspace-*.md`, `paper/`, `web/`); several are gitignored for that reason. Check what's staged before pushing.
- The desktop GUIs (SwiftUI, Avalonia) live in a separate repo, `vcfclick-desktop`, and drive this CLI as a subprocess — keep these command interfaces and their JSON stable: `db query --format JSON`, `db create`, `db ingest` (incl. `--serial`, `--keep-reference`), `db ped --proband/--father/--mother`, `db families --format json`, `db trio --format json`, `db relatedness --format json`, `db qc --format json`, `db gene --format json`, `db export`.
- User docs live in `docs/` (one page per feature, `docs/CLI.md` covers every command) — update them alongside CLI changes.
- Out of scope: structural variants, non-stdio MCP transports.
