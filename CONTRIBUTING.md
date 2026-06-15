# Contributing to vcfclick

Thanks for considering a contribution. The bar is "the change makes
the project better and doesn't break the parts that work" — small
patches are welcome.

## Setup

```bash
git clone https://github.com/nuin/vcfclick.git
cd vcfclick
uv sync --group dev
```

This installs the runtime + the dev extras (`pytest`). vcfclick's
native deps (`cyvcf2`, `chdb`, `duckdb`, `pyarrow`) ship as prebuilt
wheels for macOS arm64 and Linux x86_64; other platforms build from
source — `cyvcf2` needs `htslib` headers on `$PATH`.

You also need `bgzip` and `tabix` on `$PATH` — used by the smoke
tests to build fixture VCFs and by `bench/run.sh`. Get them via
`brew install htslib` (macOS) or `apt install tabix` (Debian/Ubuntu).

## Running tests

```bash
uv run pytest tests/             # full suite, ~20s
uv run pytest tests/test_cli.py  # one file
uv run pytest tests/ -k discover # filtered by name
```

The suite is hermetic — every test uses an isolated `VCFCLICK_HOME`
(via the `vcfclick_home` fixture in `tests/conftest.py`) so it can't
disturb your real `~/.vcfclick/dbs/` data. The DuckDB annotation
store is similarly isolated via `isolated_annotation_db` for tests
that touch the gene / ClinVar tables.

## Formatting

The project uses ruff for formatting and linting. CI pins the ruff
version, so use the same commands locally when checking a PR:

```bash
uvx ruff@0.15.16 format --check .
uvx ruff@0.15.16 check .
```

To apply formatting:

```bash
uvx ruff@0.15.16 format .
```

Run the formatter locally before opening a PR if the check fails.

## Pull requests

1. Fork, branch from `main`, make the change.
2. Add tests for new behavior. Bug fixes should include a test that
   would have failed before the fix.
3. `uv run pytest tests/` must be green.
4. Open the PR. CI will run the suite across Ubuntu and macOS,
   Python 3.11–3.13.
5. Commit messages follow the existing style: noun-first, present
   tense, brief. Wrap the body at ~72 chars. Example:

   ```
   Fix FORMAT field routing: PL was dropped, DP/AD silently NULL

   Two bugs surfaced while writing routing tests. ...
   ```

## What we're looking for

- Bug fixes (with tests).
- New typed-column routings — promoting fields from the
  `info_extra` / `format_extra` overflow Maps into typed columns
  via `ingest/routing.py` + the schema SQL.
- New annotation loaders (ClinVar variant of variant, dbSNP, gnomAD
  population frequencies) under `annotations/loaders/`.
- Better SCHEMA_DESCRIPTION prompt for `vcfclick_mcp/server.py` —
  ideally accompanied by transcripts showing the change steers an
  LLM client away from a specific failure mode.
- Performance improvements to ingest, ideally with `bench/run.sh`
  numbers before/after.

## What we're not looking for right now

- Structural variant (SV) support — out of audience scope (clinical
  rather than research-cohort focused).
- New MCP transports (SSE, HTTP) — stdio covers Claude Desktop and
  the test harness; no current user has asked for the others.

## Release process

vcfclick releases are tagged from `main`:

1. Update `CHANGELOG.md` — move items from `[Unreleased]` to a new
   versioned section.
2. Bump `version` in `pyproject.toml`.
3. Commit + push.
4. Tag: `git tag -a v0.1.x -m "v0.1.x: …" && git push origin v0.1.x`.
5. The `release.yml` workflow publishes to PyPI on tag push and
   creates a GitHub release with notes derived from `CHANGELOG.md`.
