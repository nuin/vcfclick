# Changelog

All notable changes to vcfclick are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning
follows [SemVer](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Runnable demo notebook** (`examples/vcfclick-demo.ipynb`, Colab-ready).
  Self-contained: it installs vcfclick, downloads its own example VCFs from
  the repo, and walks the headline flows — ingest → SQL, GIAB-validated trio
  de novo, compound-het + gnomAD rarity filtering, sample QC, and `combine` —
  with real captured outputs. Doubles as a release smoke-test.

## [0.7.0] — 2026-06-17

### Added
- **gnomAD population-allele-frequency annotation.** A `gnomad_af` table
  in the annotation store, loaded from a gnomAD sites VCF with
  `vcfclick annotations load-gnomad <vcf>` (gnomAD is too large to bundle,
  so you supply a region slice — a small one pulls in seconds via
  tabix-over-HTTPS from the public gnomAD bucket). Exposed as the
  `gnomad_lookup` MCP tool and a `db trio --gnomad-max-af` filter that
  drops candidates whose gnomAD popmax AF exceeds the threshold —
  defensible rare-variant filtering on true population frequency rather
  than the cohort's own (often missing) `info_AF`. Variants absent from
  the loaded slice are treated as rare, never AF 0.
- **`vcfclick db qc` — per-sample QC metrics.** One pass over the cohort
  (both backends): heterozygous / homozygous-alt counts and ratio,
  transition/transversion ratio (SNVs), and a chromosome-X heterozygosity
  sex check. When a pedigree is loaded the inferred sex is compared to the
  declared sex and a mismatch is flagged — the classic sample-swap signal.
  Honest about the sparse store: call-rate / missingness is not reported
  because no-calls aren't stored. `--format json` for pipelines. See
  [docs/QC.md](docs/QC.md).
- **Compound-heterozygous trio model** — `vcfclick db trio --category
  comphet` reports genes where the proband carries two rare hets in
  *trans* (one inherited from each parent), the recessive mechanism
  per-variant filters miss. Candidate variants are rare proband hets
  with a clear parent-of-origin (needs `--keep-reference`); they are
  grouped per gene via the annotation store (`vcfclick annotations
  load`), since genes can't be SQL-joined to the cohort. Reported per
  gene, not per variant.
- **`VCFCLICK_ANNOTATIONS_DB`** — environment override for the annotation
  store path, so a custom/shared store (or a test store reachable by
  subprocesses) can be used instead of the bundled default.

### Tests
- **GIAB benchmark validation for trio analysis** (`docs/VALIDATION.md`).
  The inheritance models are tested against real genotypes from the
  Genome in a Bottle Ashkenazi trio (HG002/HG003/HG004, NIST v4.2.1
  GRCh38), not only hand-built fixtures — recessive, dominant, and
  compound-het (CFTR) reproduce on published data. De novo is validated
  against an **independent** ground truth: GIAB's high-confidence BED. In
  a 5 Mb chr20 region 50 HG002-only variants look de novo to the naive
  rule, but only 4 have both parents confidently hom-ref per the BED;
  vcfclick's `--keep-reference` de novo recovers exactly the 4 confident
  sites; the checked-in test exercises a 7-site subset (the 4 plus 3 of
  the no-call sites) and `bcftools +mendelian2` independently flags the
  same 4 (asserted in CI). A documented one-time `slivar` run on the same
  fixtures flags the same 4 de-novo and 2 recessive sites, with the same
  no-call exclusions. Source URLs are in the fixture headers.

## [0.6.0] — 2026-06-16

### Added
- **`vcfclick web` — optional local browser UI** (the `[web]` extra,
  FastAPI + uvicorn). Starts a localhost server and serves a single-page
  app over a cohort database with four panels: a SQL explorer (read-only
  guard), a natural-language→SQL box (same MCP schema briefing, BYO
  Gemini/Anthropic key kept in the browser), a trio panel (reusing the
  `db trio` SQL builders), and a combine panel (reusing `combine_vcfs`).
  It is a thin HTTP layer over existing internals — the SQL path is the
  same `get_session().query()` as `db query` and MCP `run_sql` — binds to
  `127.0.0.1` only, and ships the UI embedded in the wheel (no JS build).
  See [docs/WEB.md](docs/WEB.md).

## [0.5.1] — 2026-06-16

### Added
- **GATK3 equivalence test for `combine`.** A frozen golden file of real
  GATK 3.8-1 `CombineVariants` output
  (`tests/fixtures/gatk3/combine_priority.gatk3.vcf`) plus a test
  asserting `vcfclick combine` reproduces its site union, `set=`
  provenance, and PRIORITIZE genotype resolution. `docs/COMBINE.md` gains
  a side-by-side GATK3-vs-vcfclick comparison.

### Documentation
- README feature list and docs now advertise `combine` and trio/family
  analysis; the User Guide gains "Combine Call Sets" and "Trio / Family
  Analysis" sections and the FAQ a `merge` vs `combine` entry. Corrected
  the stale "does not auto-load pedigree" limit (now `db ped`).

## [0.5.0] — 2026-06-16

### Added
- **`vcfclick combine` — multi-callset merge with provenance.** A native
  reimplementation of GATK3's `CombineVariants` (dropped in GATK4 and not
  replaced by `bcftools`). Unlike `merge` (disjoint samples → one joint
  VCF), `combine` unions call sets that may *share* samples — two callers
  over the same cohort, or a pre/post-filter pair:
  - sites are unioned by `(chrom, pos, ref, alt)`;
  - each output record carries a `set=` INFO field naming the inputs it
    appears in (`Intersection` when present in all);
  - a sample shared across inputs is resolved by PRIORITY (input order):
    its genotype comes from the first input with a non-missing call;
  - `--min-callsets N` keeps only sites seen in at least N inputs
    (consensus calling — the "variants present in all / a fraction of
    the call sets" feature GATK4 specifically lost).
  - the `GQ`/`DP`/`AD` FORMAT fields are carried through from the same
    priority-source record that supplied each genotype, so depth and
    allele balance describe the call that was kept and a combined VCF
    feeds straight into the `db trio` quality gates. Output FORMAT lists
    only the passthrough fields some input actually carries.
- **Trio / family analysis.** A four-part pipeline for Mendelian
  inheritance-model candidate filtering:
  - `vcfclick merge a.vcf.gz b.vcf.gz … -o joint.vcf.gz` — combine
    per-sample VCFs into one joint multi-sample VCF (wraps
    `bcftools merge` with disjoint-sample validation, auto-indexing,
    and `-m none` to stay ingest-ready). Trios are usually delivered
    as separate per-sample VCFs; this produces the single joint VCF
    that trio analysis needs.
  - `pedigree` table + `vcfclick db ped <name> <file.ped>` — load
    standard PED/FAM relationships (family/father/mother/sex/affected),
    validated against the cohort's samples, round-tripped through
    `db dump` and `db push/pull`.
  - `vcfclick db ingest --keep-reference` — additionally store
    confident hom-reference (gt=0) calls at variant sites (no-calls
    still dropped), so de-novo analysis can PROVE a parent is 0/0
    rather than guessing from absence. Opt-in; normal cohorts stay
    sparse.
  - `vcfclick db trio <name> --proband <id>` — report candidate
    variants under de-novo, recessive, and dominant models with
    slivar-style genotype quality gates (GQ, depth, het allele
    balance) and population-AF rarity. A no-call parent is correctly
    excluded from de novo (the INNER JOIN to a provable parent 0/0
    finds no row), so it reports defensible candidates, not the naive
    "neither parent carries it" false positives.
  - The MCP schema briefing teaches the LLM the pedigree table, the
    inheritance models, and the sparse de-novo caveat, so trio
    questions can be asked in natural language with the SQL shown.

  Honest scope: candidate FILTERING, not variant calling. Quality
  gating is only as strong as the FORMAT fields present; PL-based
  Bayesian de-novo refinement is future work.

## [0.4.1] — 2026-06-12

### Fixed
- `vcfclick db pull` now restores older bundle archives whose Parquet
  files lack newer nullable/defaulted columns. Bundle restore imports
  matching columns by name instead of relying on positional `SELECT *`,
  which fixes the published 1000 Genomes BRCA1 demo bundle on the
  current DuckDB backend.

## [0.4.0] — 2026-06-12

### Added
- Optional Textual terminal UI via `vcfclick tui [--db NAME]`,
  installable with the `tui` extra. The first screen is a
  genomics-first Locus view for gene symbols and `chrom:start-end`
  ranges, with generated SQL available for inspection/editing.
- TUI service layer for locus parsing/resolution, database metadata,
  SQL execution, locus summaries, and backend-aware stats handling.
- Operations and SQL TUI panes. Operations lists local databases,
  switches the active database after metadata validation, shows
  basic counts/path information, and reports unsupported stats on
  DuckDB as a recoverable UI error.

## [0.3.2] — 2026-06-10

### Added
- Zenodo integration. `.zenodo.json` ships rich DOI metadata
  (title, description, keywords, license, author ORCID
  0000-0003-4915-664X, related-identifier links to GitHub and
  PyPI) so the DOI mints clean instead of with Zenodo's
  GitHub-fallback defaults. New `docs/CITATION.md` documents the
  citation format and the one-time OAuth setup steps. The first
  Zenodo webhook fires on this release tag.

### Fixed
- Parallel ingest works on the DuckDB backend. The docstring already
  claimed it was backend-agnostic but two latent splitter bugs
  prevented it from ever actually running through the workers:
  - `ingest._tabix.variant_density()` set the FINAL 16Kb
    linear-index bucket's byte_cost to 0, then filtered out zero-cost
    buckets, dropping the trailing position bucket on every contig.
    On 1000G phase 3 chr21 that silently lost the 134 variants past
    position 48,100,000. Now uses placeholder cost 1 so the final
    bucket survives.
  - `ingest.parallel.ingest_parallel()` checked `if regions is None`
    to fall back from the tabix splitter to the cyvcf2 pre-pass
    splitter. The tabix splitter actually returns `[]` (empty list,
    not None) when the linear index is too sparse to balance — typical
    for VCFs under a few hundred variants. Empty-list silently meant
    "0 workers, 0 variants ingested." Switched to `if not regions`.
- New `tests/test_ingest_parallel.py` locks both fixes in: the
  end-to-end test exercises parallel ingest on the 5-variant tiny
  fixture (catches the empty-list bug — without the fix, 0 rows
  land); a direct unit test on `variant_density()` catches the
  trailing-bucket bug.

### Changed
- `ingest/parallel.py` docstring + in-body comments reframed away
  from chDB-specific language now that the path is exercised on both
  backends. Architecture diagram updated to show
  `parquet_file_expr()` instead of the literal `file('p', 'Parquet')`
  fragment.

### Performance
- Parallel ingest on 1000 Genomes phase 3 chr21 (1.11M variants,
  2,504 samples) under DuckDB: **163 s with 8 workers vs 944 s
  serial — 5.8× speedup.** Per-region throughput is ~1,000 v/s per
  worker, dominated by Python row construction inside the workers.

## [0.3.1] — 2026-06-08

### Fixed
- Fixed the Bioconda recipe smoke-test database name so it satisfies
  vcfclick's database-name validator.
- Fixed `storage.drop_db()` to evict backend-qualified session-cache
  entries and close cached sessions before removing the database
  directory.
- Logged DuckDB connection close failures instead of silently
  swallowing them.

## [0.3.0] — 2026-06-08

### Added
- **DuckDB storage backend.** vcfclick now runs on either chDB
  (ClickHouse engine, the original) or DuckDB. Selection is via
  `VCFCLICK_BACKEND=chdb|duckdb`; the default auto-detects (chDB
  wins when importable, otherwise DuckDB). The DuckDB backend is
  what unblocks distribution through bioconda, where chDB is not
  available as a conda package.
  - New: `storage/_duckdb.py` (session wrapper that mimics
    `chdb.session.Session.query(sql, format)` so most call sites
    are backend-agnostic — JSONCompact / TabSeparated / CSV /
    Vertical / Pretty renderers reproduce chDB's output shapes
    including the `\N` NULL marker and Float32 precision).
  - New: `storage/_chdb.py` (extracted chDB session open with the
    EmbeddedServer async-load retry, lazy-imports `chdb` so
    DuckDB-only installs don't error at import time).
  - New: `schema/duckdb/{01_variants,02_genotypes,03_samples}.sql`
    — DuckDB-flavoured DDL with identical column names and order
    to the chDB schemas. `LowCardinality(X)` becomes `VARCHAR`,
    `Nullable(X)` collapses to `X`, `ReplacingMergeTree` engine
    config drops entirely (the application-level
    `rollback_ingest()` already enforces the same idempotent-
    replace semantics across both backends).
  - New: `storage.backend()`, `storage.parquet_file_expr()`,
    `storage.delete_where_sql()`, `storage.count_expr()`,
    `storage.table_exists()`, `storage.schema_dir_for_backend()`
    — dialect helpers used by ingest, export, CLI, and MCP code
    so SQL emission picks the right form per backend without
    spreading `if backend() == "duckdb"` everywhere.
- **Bioconda recipe** at `packaging/bioconda/meta.yaml`. Lists
  cyvcf2, pyarrow, duckdb, mcp, click as run requirements;
  omits chdb (unavailable on conda-forge as of 2026-06; see the
  closed upstream issue chdb-io/chdb#189 for why). Auto-detect
  picks DuckDB when chDB isn't installed, so `conda install
  -c bioconda vcfclick` followed by `vcfclick db create` works
  with no extra environment-variable setup.
- CI matrix gains a `VCFCLICK_BACKEND=duckdb` cell so every push
  exercises both backends.

### Changed
- `ingest/_arrow.py` column-list helper switched from backtick
  identifier quoting to double quotes. DuckDB rejects backticks;
  chDB accepts both. Same INSERT contracts otherwise.
- `cli/db_diff.py` cohort allele-frequency SQL switched from
  chDB's `sumIf()` to standard-SQL `sum() FILTER (WHERE …)` with
  `COALESCE(..., 0)` so cohorts with no matching rows produce 0
  rather than NULL on either engine.
- `export/parquet.py` (`db dump`) emits `COPY … TO 'path' (FORMAT
  'parquet')` on DuckDB and the existing `INTO OUTFILE 'path'
  TRUNCATE FORMAT Parquet` on chDB.
- `cli/db_stats.py` raises a clear ClickException on the DuckDB
  backend (the chDB-specific SQL — `system.columns`, `countIf`,
  `ARRAY JOIN mapKeys`, type-string `Nullable` inspection — has
  not been ported yet). The 10 `test_stats.py` tests are skipped
  on the DuckDB cell of the CI matrix until the port lands.

### Internal
- CLI commands split out of the monolithic `cli/db.py` into
  focused modules (`cli/db_basic.py`, `cli/db_batch.py`,
  `cli/db_bundle.py`, `cli/db_diff.py`, `cli/db_stats.py`).
  `cli/db.py` is now an importlib shim that loads them by side
  effect, matching the pattern the project's quality gate
  expects for per-file size. No user-visible behaviour change.
- Ingest module similarly split: `ingest/parallel_split.py`
  (variant-count-aware region splitter) and `ingest/vcf_rows.py`
  (row builders) factored out of `ingest/parallel.py` and
  `ingest/vcf_load.py`.

## [0.2.0] — 2026-06-08

### Added
- `vcfclick db ingest-parquet <db> <dump_dir>` — the symmetric
  inverse of `db dump`. Reads variants/genotypes/samples Parquet
  files matching the locked Arrow schemas in `ingest/_arrow.py`
  and lands them in chDB under a new (cohort, ingest_id) label.
  The full round-trip is `db dump → ingest-parquet` — the same
  three Parquet files that fall out of a dump are the wire format
  that comes back in.
  - The source `ingest_id` and `cohort` columns are rewritten to
    the caller's `--ingest-id` and `--cohort` via the SELECT list
    of an `INSERT … SELECT FROM file()`, so moving a dump between
    cohorts/labels is a single command.
  - Samples are imported as-is when `samples.parquet` is present,
    derived via `SELECT DISTINCT sample_id` against
    `genotypes.parquet` when not, or skipped entirely for a
    variants-only cohort-summary ingest.
  - Phase-1 schema validation rejects files with the wrong column
    set BEFORE the per-ingest_id file lock is acquired, so a bad
    re-ingest under a stable id leaves the prior data intact —
    same atomicity contract as the VCF path.
  - The `ingested_at` `DEFAULT now()` ReplacingMergeTree version
    column is tolerated on input but not carried through; chDB
    re-defaults it on the new INSERT.
  - 9 new tests in `tests/test_ingest_parquet.py`: round-trip
    counts preservation, ingest_id/cohort override, replay-replace
    semantics, ingestions-catalog provenance, variants-only valid,
    samples-derived-from-genotypes, missing-variants rejection,
    schema-mismatch rejection before chDB touch, ingest_id
    validation.

### Changed
- `docs/SCHEMA.md` adds a "Parquet as a public interchange format"
  section documenting the dump-ingest round-trip, the column-set
  contract, server-default column handling, and how external tools
  (DuckDB, polars, Spark) can produce conforming Parquet.

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

[Unreleased]: https://github.com/nuin/vcfclick/compare/v0.7.0...HEAD
[0.7.0]: https://github.com/nuin/vcfclick/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/nuin/vcfclick/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/nuin/vcfclick/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/nuin/vcfclick/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/nuin/vcfclick/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/nuin/vcfclick/releases/tag/v0.4.0
[0.3.2]: https://github.com/nuin/vcfclick/releases/tag/v0.3.2
[0.3.1]: https://github.com/nuin/vcfclick/releases/tag/v0.3.1
[0.3.0]: https://github.com/nuin/vcfclick/releases/tag/v0.3.0
[0.2.0]: https://github.com/nuin/vcfclick/releases/tag/v0.2.0
[0.1.3]: https://github.com/nuin/vcfclick/releases/tag/v0.1.3
[0.1.2]: https://github.com/nuin/vcfclick/releases/tag/v0.1.2
[0.1.1]: https://github.com/nuin/vcfclick/releases/tag/v0.1.1
[0.1.0]: https://github.com/nuin/vcfclick/releases/tag/v0.1.0
