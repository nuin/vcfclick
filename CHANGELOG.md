# Changelog

All notable changes to vcfclick are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning
follows [SemVer](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.3] — 2026-06-07

### Fixed
- `storage.sql_quote_str` now escapes BACKSLASHES as well as single
  quotes. ClickHouse/chDB recognises both `''` and `\\'` as quote
  escapes inside string literals; a path containing `\\'` would,
  under plain quote-doubling, become `\\''` — the engine reads the
  first quote as backslash-escaped, the second as the string
  terminator, and arbitrary SQL after it. Closes a remaining
  injection bypass codex round 9 caught after the round-8 fix.
  Round-trip test through chDB itself (FORMAT JSONCompact to
  bypass TSV escaping) locks the contract for backslash + quote
  payloads including the classic `\\'; DROP TABLE …` shape.
- `storage.DB_ROOT` / `storage.VCFCLICK_HOME` now lazy at the
  package level too. The round-8 fix made them lazy in
  `storage.db` via module `__getattr__`, but `storage/__init__.py`
  had already imported them eagerly — so reads through the public
  `storage.DB_ROOT` API kept seeing the import-time value. Replaced
  the eager imports with package-level `__getattr__` forwarding.
  New `test_storage_db_root_reflects_env_change_at_call_time`
  test locks the lazy contract end-to-end.

### Added
- `storage.sql_quote_str(s)` — SQL-standard single-quoted literal with
  embedded quotes doubled. Used at every site where we interpolate a
  string into a chDB query (parquet file paths in
  `INSERT … FROM file('…')`, etc.). Defence against a single quote
  in the `staging_dir` parameter to `ingest_parallel` closing the
  `file()` literal and letting raw SQL through. Plain unit tests in
  `tests/test_ingest_concurrency.py` lock the contract including the
  classic `' OR 1=1 --` injection payload.
- `storage.ingest_id_lock(ingest_id)` — `fcntl.flock`-based exclusive
  per-`(DB, ingest_id)` file lock. Both ingest paths
  (`ingest.vcf_load.ingest`, `ingest.parallel.ingest_parallel`)
  acquire the lock at the start of the call and release on exit.
  Two concurrent `vcfclick db ingest` invocations sharing an
  `ingest_id` (workflow runner firing parallel jobs, retry script
  racing the previous run) used to race on the staging dir,
  rollback, and bulk-import glob, producing a mixed corrupt state.
  Blocking acquisition means the second call waits until the first
  finishes. Unix-only (Windows isn't a supported target).
- `tests/test_ingest_concurrency.py` — 7 tests covering both
  hardening pieces, including a multiprocessing test that asserts
  the lock blocks concurrent holders across processes and a
  same-DB-but-different-`ingest_id` test that asserts the locks
  don't over-serialize.

### Changed
- `storage.VCFCLICK_HOME` and `storage.DB_ROOT` are now computed at
  attribute-access time (module-level `__getattr__`) rather than
  cached at import. Tests and in-process flows that change
  `VCFCLICK_HOME` after `import storage` now see the new path
  immediately. No effect on existing callers.

### Fixed
- Allele-frequency briefing in `SCHEMA_DESCRIPTION` and `docs/SCHEMA.md`
  now teaches the correct denominator pattern: compute cohort size from
  `samples` ALONE as a CTE/subquery, then bring it into the genotypes
  aggregation via CROSS JOIN. The previous wording — added when
  cohort_sizes_mv was dropped — implied counting `DISTINCT samples`
  *through* the join to `genotypes`, which silently shrinks the
  denominator to the non-reference set (the genotypes table is sparse)
  and produces inflated AF. `vcfclick db diff` already uses the
  correct pattern; the briefing now matches it explicitly. Test in
  `tests/test_mcp_server.py::test_briefing_includes_non_obvious_invariants`
  locks the contract — the briefing must mention `cohort_size` and
  `CROSS JOIN`. Caught by codex in the sixth review pass.

### Removed
- **`cohort_sizes_mv` materialized view** is no longer created on
  fresh databases. `SummingMergeTree` never decrements on DELETE, so
  rolling back a samples insert (failed ingest, `--ingest-id`
  replacement) left the view's count permanently inflated and any
  AF query that used it as the denominator produced wrong numbers
  after retries. The MCP `SCHEMA_DESCRIPTION` briefing and
  `docs/SCHEMA.md` now teach the direct
  `count(DISTINCT (s.ingest_id, s.sample_id))` pattern instead;
  `vcfclick db diff` already used it. At realistic cohort scale the
  `count(DISTINCT)` runs in microseconds, so this is correctness
  with no measurable cost. Existing 0.1.2 databases still have the
  old view (`apply_schema` uses `CREATE … IF NOT EXISTS` semantics
  and won't drop it) — harmless but stale; recreate the DB if you
  want a clean 0.1.3 schema. Spotted in a codex CLI review.

### Fixed
- Parallel ingest staging dir is now wiped before workers start, so a
  retry under the same `ingest_id` after a failed parallel parse can't
  pick up stale Parquet files from the previous run via the bulk-import
  glob. Closes a corruption path codex flagged after the stage-then-
  commit restructure.
- Serial ingest staging now happens on the same volume as the
  destination DB (`$VCFCLICK_HOME/dbs/`) rather than the system `/tmp`.
  With stage-then-commit the whole VCF's Parquet is held on disk
  before Phase 2; on systems where `/tmp` is a small ramdisk or a
  separate partition smaller than the DB volume, a fresh large ingest
  could fail before commit purely from staging-disk exhaustion. The
  DB volume necessarily has room for the destination data, so it has
  room for the staging equivalent. Also flagged by codex.

### Changed
- Both ingest paths now use a **two-phase stage-then-commit** flow.
  Phase 1 reads the VCF, validates per-record invariants
  (multi-allelic check fires here), and writes Parquet batches to a
  tempdir — **no chDB writes happen during parsing**. Phase 2 runs only
  if Phase 1 completes successfully: rollback prior rows under the
  ingest_id, insert samples, bulk-import the staged Parquets in one
  go, then write the ingestions catalog. A `commit_started` flag
  controls the except arm: failures during Phase 1 (multi-allelic,
  malformed body, corrupt header, Ctrl-C mid-parse) skip the rollback
  so prior data under that ingest_id is preserved.

### Fixed
- Mid-stream re-ingest failures no longer wipe prior data. Closes the
  third codex CLI finding: a re-ingest under an existing ingest_id
  whose new VCF opens fine but fails inside the variant loop (e.g.,
  multi-allelic record encountered) used to delete the prior good
  rows and leave no replacement. The stage-then-commit restructure
  above means the multi-allelic check now raises in Phase 1, before
  rollback ever runs. New test
  `test_reingest_with_multi_allelic_preserves_prior_data` locks the
  contract.
- Atomic `--ingest-id` replacement against early failures. The previous
  fix called `rollback_ingest()` before opening the new VCF, which
  silently wiped prior data on every failed re-ingest (corrupt header,
  unreadable file, classification failure). The rollback now runs
  inside the `try` block — only after the new VCF opens and classifies
  successfully. A bad header raises before any deletion happens, so a
  re-ingest under an existing id whose new VCF fails to read leaves
  the prior rows intact. New test
  `test_reingest_with_corrupt_vcf_preserves_prior_data` locks the
  contract. Docstring updated to be honest about the residual gap:
  failures during the variant loop still wipe prior data; for full
  atomicity against mid-stream failures, use a fresh ingest_id.
  Spotted by codex CLI in a second review pass.
- Arrow ↔ SQL column-order drift in `ingest/_arrow.py`. The Arrow
  schemas now match `schema/01_variants.sql` and
  `schema/02_genotypes.sql` byte-for-byte (previously `info_AD_ref/_alt`
  and the flag block were in the wrong position on the variants side,
  and `mq/ft/ps/pq` were misplaced on the genotypes side). All four
  `INSERT … SELECT … FROM file('*.parquet')` sites in `ingest/vcf_load.py`,
  `ingest/parallel.py`, and `storage/db.py` now use explicit
  `INSERT INTO t (col1, col2, ...) SELECT col1, col2, ... FROM file(...)`
  column lists derived from the Arrow schema, so the import is immune
  to any future chDB change shifting Parquet handling from name-based
  to positional. A new drift-guard test in
  `tests/test_schema_agreement.py` parses each schema/*.sql file and
  asserts column-by-column agreement with the matching Arrow schema —
  fails the build on any future drift. Spotted as P1 in the codex
  CLI review.
- `--ingest-id` reuse now truly replaces prior data, matching what the
  CLI help text always claimed. Previously the ingest path relied on
  `ReplacingMergeTree` dedup on the sorting key, so re-ingesting under
  an existing `ingest_id` was an upsert: rows from the prior ingest
  that weren't in the new VCF stayed queryable. Both ingest paths
  (`ingest.vcf_load.ingest`, `ingest.parallel.ingest_parallel`) now
  call `rollback_ingest(ingest_id)` before writing, so the documented
  "replace" semantics are literal. No-op on fresh IDs. Test in
  `tests/test_ingest_atomic.py::test_reingest_same_id_truly_replaces_prior_data`
  locks the contract. Spotted in an external code review by codex CLI.

## [0.1.2] — 2026-06-05

### Added
- `docs/SCHEMA.md` — flat reference for every column on every table
  (`variants`, `genotypes`, `samples`, `ingestions`), plus the three
  conventions SQL writers have to internalise (sparse genotypes,
  cross-ingestion non-merging, the GQ/DP NULL silent-failure trap)
  and the four-step promotion recipe for moving an overflow field
  to a typed column. README links to it.
- `vcfclick db stats <name>` — schema-population stats for an ingested
  cohort. Reports row counts (variants / genotypes / samples /
  ingestions), per-cohort and per-contig breakdowns, the populated
  fraction of every typed column on `variants` and `genotypes`
  (sorted by population descending so the actually-used fields
  surface first), and the most frequent `info_extra` / `format_extra`
  overflow Map keys. `--top N` caps the overflow listings. Closes
  the "I don't know what's in this DB" gap that the `discover`
  command only closed for VCFs pre-ingest.
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

[Unreleased]: https://github.com/nuin/vcfclick/compare/v0.1.3...HEAD
[0.1.3]: https://github.com/nuin/vcfclick/releases/tag/v0.1.3
[0.1.2]: https://github.com/nuin/vcfclick/releases/tag/v0.1.2
[0.1.1]: https://github.com/nuin/vcfclick/releases/tag/v0.1.1
[0.1.0]: https://github.com/nuin/vcfclick/releases/tag/v0.1.0
