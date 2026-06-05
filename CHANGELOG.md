# Changelog

All notable changes to vcfclick are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning
follows [SemVer](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- chDB session-init retry wrapper in `storage.db._open_session_with_retry`.
  Catches the known intermittent EmbeddedServer async-load race
  (`recursive_mutex lock failed: Invalid argument` → `BAD_ARGUMENTS` /
  `ASYNC_LOAD_WAIT_FAILED`) and retries up to 3× with 0.1 / 0.3 / 0.9 s
  backoff. Non-race exceptions propagate immediately so real failures
  aren't masked. Logs a warning when a retry fires. 5 unit tests in
  `tests/test_chdb_retry.py` lock the contract.
- DRAGEN-specific INFO routing — `FractionInformativeReads`, `HAPCOMP`,
  `HAPDOM`, `DragenSnvHardQUAL`, `DragenIndelHardQUAL` now promote to
  typed columns on `variants` instead of landing in `info_extra`.
  DRAGEN cohorts can query QC metrics with the same typed-column
  ergonomics as GATK cohorts. Per-sample DRAGEN FORMAT fields
  (`ML_PROB` etc.) stay in `format_extra` until a Float FORMAT
  routing path is added.
- `vcfclick db ingest-batch <name>` — multi-file ingest for per-sample
  VCFs (DRAGEN, GATK `-ERC GVCF`, clinical-pipeline output). Takes
  either `--from-dir <path>` (scans `*.vcf.gz`, derives ingest_id
  from filename) or `--manifest <tsv>` (nf-core/Snakemake-style;
  required `vcf_path` column, optional `sample_id`/`ingest_id` and
  `cohort` overrides). Files are ingested sequentially via library
  call (shared chDB session). Per-file atomic — a failure rolls back
  only that file and the batch continues; the summary lists what
  failed and exit code is non-zero if anything did.
- `vcfclick db diff <name> --cohort-a A --cohort-b B [--top N]` —
  per-variant allele-frequency comparison across two cohorts in the
  same database. Computes AC, AN, AF, and AF difference for every
  variant called in either cohort, sorted by absolute AF diff
  descending. Honest about sparse-table semantics (absent samples
  count as 0/0) and cross-ingestion accumulation (same sample under
  two ingest_ids counts as two observations). Cohort names go
  through SQL-standard quote-doubling — embedded quotes don't break
  out of the literal.
- CHANGELOG.md, CONTRIBUTING.md, and a release-on-tag GitHub Actions
  workflow that runs tests, builds, publishes to PyPI, and creates
  a GitHub release in one tag push. Requires a `PYPI_API_TOKEN`
  repository secret.
- Atomic ingest — failed ingests (multi-allelic record, batch flush
  error, Ctrl-C) now roll back every row written under the failed
  `ingest_id` so the database is left in its pre-call state, instead
  of forcing a `vcfclick db rm` and full re-ingest.
- `validate_ingest_id()` on the storage layer — rejects ingest IDs
  that aren't safe to interpolate into the rollback `ALTER DELETE`
  statement (chDB has no parameter binding for ALTER paths).
- GitHub Actions CI — runs the test suite on every push and PR
  across Ubuntu and macOS, Python 3.11–3.13 (Linux × 3.13 excluded
  pending a chDB upstream fix).
- `tests/test_mcp_server.py` — 11 integration tests for the MCP layer:
  tool-registry drift, JSON-Schema generation, SCHEMA_DESCRIPTION
  invariants, DuckDB-backed tool semantics, and a stdio-transport
  smoke test that spawns the actual server.
- `tests/test_ingest_atomic.py` — 4 tests covering rollback behavior
  and `validate_ingest_id` rejection.
- README CI + PyPI badges.

### Changed
- **Schema added five new typed columns on `variants`** for the DRAGEN
  fields above. Existing databases created on 0.1.1 don't have these
  columns; ingest against them will fail until the database is
  re-created (`vcfclick db rm <name>` then `vcfclick db create <name>`
  and re-ingest). Fresh databases on this version onward pick up the
  new schema automatically.
- `INFO_SCALAR`/`_PAIR`/`_FLAG`, `FORMAT_SCALAR`/`_PAIR`/`_TRIPLE`,
  and `classify_header()` moved from `ingest/vcf_load.py` into a
  dedicated `ingest/routing.py` module. `vcf_load`, `parallel`,
  `cli/discover` all import from the new location. Single source of
  truth for schema routing.
- README now leads with `uv tool install vcfclick` / `pipx install
  vcfclick` for CLI use, with `pip install vcfclick` documented as
  the library-embedding form.

## [0.1.1] — 2026-06-02

### Added
- `vcfclick discover <vcf>` — previews which INFO/FORMAT fields land
  in typed columns vs the `info_extra` / `format_extra` overflow Maps,
  with suggested promotion DDL. No ingestion needed.
- `vcfclick annotations load-clinvar` — pulls the NCBI ClinVar weekly
  VCF and populates the embedded DuckDB. Replaces the stubbed
  `clinvar_lookup` MCP tool with a real one. Normalises bare numeric
  contigs to `chr`-prefixed during load.
- `bench/run.sh` — reproducible ingest benchmark against the 1000G
  chr17:40-50M slice. Caches the slice locally, runs all three
  configurations (serial, parallel-4, parallel-8) from cold.
- `examples/brca1-cohort.md` — worked MCP session with verbatim SQL
  and chDB outputs against the demo bundle. Five canonical questions,
  the tool calls each one makes, and the actual result.
- `tests/test_cli.py`, `tests/test_ingest_routing.py`,
  `tests/test_discover.py`, `tests/test_clinvar_loader.py` —
  29 tests covering CLI lifecycle, schema routing, the discover
  command, and the ClinVar loader.

### Fixed
- FORMAT field routing — `FORMAT/PL` and `FORMAT/GL` were silently
  dropped during ingest (the schema reserved the columns but
  `build_genotype_rows` never populated them). Re-ingest after
  upgrading to recover the lost data.
- cyvcf2 `gt_depths` / `gt_quals` / `gt_ref_depths` / `gt_alt_depths`
  shortcuts return `-1` when FORMAT column ordering varies between
  records in the same VCF — even when the field is clearly present.
  Switched to `variant.format(field)` which reads the right column
  unconditionally.

## [0.1.0] — 2026-05-31

Initial PyPI release.

### Added
- Embedded chDB (ClickHouse engine) for sample data; embedded DuckDB
  for reference annotations.
- Multi-database layout at `$VCFCLICK_HOME/dbs/<name>/` — one chDB
  session per cohort, each independently dumpable and shareable.
- `vcfclick db` CLI with `create`, `list`, `info`, `ingest`, `query`,
  `dump`, `push`, `pull`, `rm`, `path` subcommands.
- Parallel ingestion (`--workers N`) using a tabix `.tbi`-driven
  variant-count-aware splitter. ~2,500 v/s at parallel-4, 3,400 v/s
  at parallel-8 on the 1000G chr17 benchmark.
- Push / pull bundle workflow — tar.gz of Parquet exports, HTTPS URLs
  accepted for `pull`. Powers the GitHub-Releases demo distribution.
- MCP server (`vcfclick_mcp.server`) exposing `get_schema`, `run_sql`,
  `position_for_gene`, `gene_at`, `clinvar_lookup` tools to any
  MCP-capable client (Claude Desktop, etc.).
- GENCODE v45 gene-coordinates loader (`vcfclick annotations load`).
- Apache 2.0 license.

[Unreleased]: https://github.com/nuin/vcfclick/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/nuin/vcfclick/releases/tag/v0.1.1
[0.1.0]: https://github.com/nuin/vcfclick/releases/tag/v0.1.0
