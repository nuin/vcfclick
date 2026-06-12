# vcfclick

[![test](https://github.com/nuin/vcfclick/actions/workflows/test.yml/badge.svg)](https://github.com/nuin/vcfclick/actions/workflows/test.yml)
[![PyPI](https://img.shields.io/pypi/v/vcfclick.svg)](https://pypi.org/project/vcfclick/)
[![PyPI downloads](https://img.shields.io/pypi/dm/vcfclick.svg)](https://pypistats.org/packages/vcfclick)
[![GitHub stars](https://img.shields.io/github/stars/nuin/vcfclick.svg?style=flat)](https://github.com/nuin/vcfclick/stargazers)

A modern VCF database for research labs and bioinformatics teams.
Embedded chDB (ClickHouse engine, serverless) for sample data, embedded
DuckDB for reference annotations, and a NL query layer
that turns plain English into SQL for demonstration.

Single binary. `uv run vcfclick`.  no port, no server, no
Gatekeeper dialog. The headline demo runs from a clean `git clone`.

Status: research preview. Architecture validated against real 1000 Genomes data.
Next gnomAD.

## Why

Two complaints heard repeatedly in research bioinformatics:

1. *"My cohort grew and `bcftools | pandas` stopped scaling."* When
   you have 500+ samples, ad-hoc cohort correlation queries become
   painfully slow. The standard answer is "go install Hail," which is
   correct and expensive.

2. *"I can write the SQL, but I shouldn't have to type the boilerplate
   every time — and when it's written for me, I want to see it."*
   Bioinformaticians don't want SQL hidden. They want it generated and
   visible, because trust comes from being able to read what ran.

vcfclick closes both:

- **chDB** (ClickHouse embedded as a library) handles cohort scale.
  We've measured ~963 variants/sec single-process ingest, 6% sparse
  compression vs dense, in-process Native query speed.
- The **MCP server** lets any LLM client translate plain English into
  the SQL underneath. The generated SQL is shown alongside the result —
  it's *part* of the answer, not a debug trace.

## Architecture

```
┌────────────────────────────────────┐
│  Tiny web UI (separate repo)       │   English in → SQL + result out
└────────────────┬───────────────────┘
                 │
┌────────────────▼───────────────────┐
│  MCP server (Python)               │   Composes the two embedded stores
│  Tools: get_schema, run_sql,       │
│    position_for_gene, gene_at,     │
│    clinvar_lookup                  │
└────┬─────────────────────────┬─────┘
     │                         │
┌────▼──────────────┐  ┌───────▼────────────┐
│  chDB             │  │  DuckDB            │
│  (embedded)       │  │  (embedded)        │
│  sample data      │  │  reference data    │
│  - variants       │  │  - genes (RefSeq)  │
│  - genotypes      │  │  - clinvar_*       │
│  - samples        │  │                    │
│  - ingestions     │  │                    │
└───────────────────┘  └────────────────────┘
```

Two embedded stores, distinct purposes:

- **chDB** holds sample data: wide pre-declared schema for VCF 4.3
  reserved + common GATK INFO/FORMAT fields, with
  `Map(String, String)` overflow for anything else. **Same SQL surface,
  same MergeTree engines, same projections as full ClickHouse — serverless
  .** Persistent on disk under `.chdb/`.
- **DuckDB** holds reference data: RefSeq genes, ClinVar. Embedded,
  swappable, monthly refresh. Never touches sample data.

The MCP server composes across them at query time. Annotation lookups
happen first (DuckDB), then their results parameterise the sample
query (chDB). The chain of reasoning is visible in the UI.

## Installation

```bash
uv tool install vcfclick     # recommended for CLI use
# or
pipx install vcfclick
```

Both install an isolated environment behind the scenes and expose the
`vcfclick` command on `$PATH`. Upgrade with `uv tool upgrade vcfclick`
(or `pipx upgrade vcfclick`).

If you're embedding vcfclick as a library inside your own Python
project (e.g., importing `vcfclick_mcp.server`), use the project
form instead:

```bash
pip install vcfclick         # inside your own venv
uv pip install vcfclick      # on a uv environment
```

vcfclick is pure-Python; its native dependencies (`cyvcf2`, `chdb`,
`duckdb`) ship as prebuilt wheels for macOS arm64 and Linux x86_64.
Other platforms build from source — `cyvcf2` needs `htslib` headers
on `$PATH`.

Listing: <https://pypi.org/project/vcfclick/>.

### Storage backend

vcfclick runs on either of two embedded SQL engines: **chDB**
(ClickHouse) or **DuckDB**. The choice is per-process via
`VCFCLICK_BACKEND=chdb|duckdb`; if unset, auto-detect picks chDB
when its import succeeds and otherwise falls back to DuckDB.

| Install path | What you get |
|---|---|
| `pip install vcfclick` / `uv tool install vcfclick` | Both backends installed; default to chDB |
| `conda install -c bioconda vcfclick` | DuckDB only (chDB is not packageable for conda; see `packaging/bioconda/README.md`); auto-detect picks DuckDB |
| `pip install vcfclick && pip uninstall chdb` | DuckDB only on pip; auto-detect picks DuckDB |

The two backends use distinct on-disk formats (chDB writes a
ClickHouse data directory; DuckDB writes a single `store.duckdb`
file), so a database created under one backend is NOT readable by
the other. Move data between them with `db dump` → `db ingest-parquet`
under the destination backend.

### Optional terminal UI

vcfclick also has an optional Textual terminal UI for local exploration:

```bash
uv tool install "vcfclick[tui]"   # or: pipx install "vcfclick[tui]"
vcfclick tui --db my-cohort
```

The TUI starts with a genomics-first Locus view: enter a gene symbol such as
`BRCA1` or a range such as `chr17:43044295-43125483`, inspect the summary, then
open the generated SQL when you want to inspect or edit it.

## 30-second demo

A pre-built 1000 Genomes Phase 3 BRCA1 cohort (3,014 variants × 3,202
samples) included as an asset. Three commands from a clean machine:

```bash
uv tool install vcfclick

vcfclick db pull demo \
    https://github.com/nuin/vcfclick/releases/download/v0.1.0/1000g-brca1-demo.tar.gz

vcfclick db query demo \
    "SELECT count(DISTINCT (ingest_id, sample_id)) FROM genotypes
     WHERE chrom='chr17' AND pos BETWEEN 43044295 AND 43170245"
```

## Using vcfclick on your own data

Each cohort / study / VCF lives in its own small database under
`~/.vcfclick/dbs/<name>/`. The `vcfclick` CLI manages them.

```bash
# Normalise the VCF (one-time per file)
bcftools norm -m - input.vcf.gz | bgzip > normalised.vcf.gz

# Preview which INFO/FORMAT fields will be stored in typed columns vs the
# overflow Maps — and what DDL would promote an overflow field to typed
vcfclick discover normalised.vcf.gz

# Create a database for this cohort
vcfclick db create my-cohort

# Ingest the VCF into it
vcfclick db ingest my-cohort normalised.vcf.gz \
    --cohort demo --ingest-id batch_a

# Or ingest many per-sample VCFs (DRAGEN, GATK -ERC GVCF, etc.) as one
# cohort — each file becomes its own ingest_id, atomic per-file:
vcfclick db ingest-batch my-cohort \
    --from-dir per_sample_vcfs/ --cohort study1
# ...or with an nf-core/Snakemake-style manifest (TSV with vcf_path
# column; optional sample_id and cohort columns):
vcfclick db ingest-batch my-cohort \
    --manifest samples.tsv --cohort fallback

# Inspect what's in it
vcfclick db info my-cohort

# Run SQL directly
vcfclick db query my-cohort "SELECT count() FROM variants"

# Export the whole database as Parquet (interop with DuckDB,
# Snowflake, BigQuery, Spark, Iceberg)
vcfclick db dump my-cohort --out my-cohort-export/

# Show what fields are actually populated in this cohort — the analog
# of `vcfclick discover` (which works on a VCF) but operating on the
# stored data after ingest.
vcfclick db stats my-cohort

# Compare allele frequencies across two cohorts in the same DB
# (case-vs-control, drug-vs-placebo, population-vs-population)
vcfclick db diff my-cohort --cohort-a case --cohort-b control --top 50

# Bundle a database as a single tar.gz for sharing
vcfclick db push my-cohort /path/to/my-cohort.tar.gz

# Restore from a bundle — local file or HTTPS URL
vcfclick db pull other-cohort https://example.com/other-cohort.tar.gz

# List, remove
vcfclick db list
vcfclick db rm my-cohort
```

Each database is a self-contained chDB session — the on-disk format is
byte-identical to a full ClickHouse server. Multiple databases sit side
by side; each is cheap to create, dump, share, or delete.

The ingester prints a classification of the VCF's INFO/FORMAT fields
on startup — what is stored in typed columns vs. the overflow Maps. That
log line is the "adapts to any VCF" claim made literally visible.

**Per-ingestion identity inside a database.** Every row carries
`ingest_id`. Rows are NOT merged across uploads — the same
`(chrom, pos, ref, alt)` observed in two different VCFs is two rows,
because annotations and QC origin can differ. Re-running with the same
`--ingest-id` is idempotent (silently replaces prior rows via
`ReplacingMergeTree`). Using a new `--ingest-id` appends.

**Parallel ingestion** is the default; pass `--serial` to force the
single-process loader. The parallel splitter does a single-pass count
of variants per 100Kb position bucket via the tabix `.tbi` index (~1 ms)
and greedy-splits each contig into ranges of approximately equal
variant count — so dense subregions (gene panels, exomes) don't leave
N–1 workers idle.

### Pointing the MCP server at a specific database

In your Claude Desktop / MCP-client config, set `VCFCLICK_DB_NAME` to
the database you want the LLM to talk to:

```jsonc
"vcfclick": {
  "command": "/path/to/vcfclick/.venv/bin/python",
  "args": ["-m", "vcfclick_mcp.server"],
  "cwd": "/path/to/vcfclick",
  "env": {
    "PYTHONPATH": "/path/to/vcfclick",
    "VCFCLICK_DB_NAME": "my-cohort"
  }
}
```

Register multiple `vcfclick-<dbname>` entries if you want the LLM to be
able to switch between cohorts in a single Claude Desktop session.

Worked example with real SQL and real outputs: see
[`examples/brca1-cohort.md`](examples/brca1-cohort.md) — five
canonical questions against the demo bundle, the MCP tools the LLM
calls for each, the SQL it generates, and verbatim chDB results.

### Annotation reference store

The MCP server's annotation tools (`position_for_gene`, `gene_at`,
`clinvar_lookup`) read from the embedded DuckDB. Two one-time loads
after installing:

```bash
# Gene coordinates (GENCODE v45 — ~60 MB, ~61,000 genes).
# Required for position_for_gene / gene_at.
vcfclick annotations load

# Pathogenic / benign variant calls (NCBI ClinVar weekly release —
# ~80 MB compressed, ~3M variants). Required for clinvar_lookup.
vcfclick annotations load-clinvar
```

GENCODE updates yearly; ClinVar updates weekly. Re-run either command
to refresh. Both default to downloading the canonical source; pass
`--gff` or `--vcf` to load from a local file instead.

## Schema reference

Every column on `variants`, `genotypes`, `samples`, `ingestions`,
plus the three conventions that matter at query time (sparse
genotypes, cross-ingestion non-merging, the GQ/DP NULL trap):
[`docs/SCHEMA.md`](docs/SCHEMA.md).

## Layout

- `schema/` — ClickHouse DDL (chDB applies it unchanged).
- `storage/db.py` — chDB session singleton; `apply_schema()` helper.
- `ingest/vcf_load.py` — serial cyvcf2-based ingester.
- `ingest/parallel.py` — multi-process variant; Parquet staging.
- `ingest/_arrow.py` — pyarrow schemas matching the ClickHouse tables.
- `export/parquet.py` — table → Parquet export CLI.
- `annotations/db.py` — DuckDB annotation API (gene, ClinVar).
- `annotations/transcripts.py` — transcript/exon/CDS API stubs (Phase 2).
- `vcfclick_mcp/server.py` — MCP server (chDB + DuckDB tool surface).
  Renamed from `mcp/` so the directory does not shadow the upstream
  `mcp` Python SDK.
- `data/` — VCF inputs (gitignored).

## Validated against real data

| Workload | Vars | Samples | Calls stored | Throughput |
|---|---|---|---|---|
| BRCA1 region (1000G 30x) | 1,863 | 3,202 | 369,776 | small-VCF baseline |
| 10 Mb chr17 (1000G 30x) — serial | 235,768 | 3,202 | **44,986,737** | 952 v/s |
| 10 Mb chr17 (1000G 30x) — parallel 4 workers | 235,768 | 3,202 | **44,986,737** | 1,983 v/s (2.1×) |
| 10 Mb chr17 (1000G 30x) — parallel 8 workers | 235,768 | 3,202 | **44,986,737** | 2,466 v/s (2.6×) |

Parallel speedup comes from the variant-count-aware splitter — each
worker gets approximately equal work regardless of where the data
actually lives along the chromosome. Sparse-table compression
empirically 6.2% of dense theoretical max.

## Design comparison with TileDB-VCF

vcfclick and TileDB-VCF have different design centres, not different
points on the same axis. The categorical differences below are
intrinsic to what each tool is built for:

| | vcfclick | TileDB-VCF |
|---|---|---|
| Intended input | joint VCF | per-sample VCFs |
| Joint VCF support | native | not currently supported (`Combined VCFs are currently not supported` runtime error) |
| Pre-processing for joint-VCF input | none | `bcftools +split` per sample + `tabix` × N |
| Pre-processing disk overhead | none | per-sample VCFs duplicate headers; on the 235k-variant 1000G slice, 114 MB joint expands to ~97 GB across 3,202 per-sample files after `bcftools +split` |
| Storage model | chDB MergeTree (ClickHouse engine) | TileDB 2D sparse array |
| Query surface | SQL via chDB | `tiledbvcf-cli export` to VCF stream |
| Cross-cohort comparison | `samples.cohort` JOIN in SQL | per-array; application-level |
| Primary audience | joint-VCF cohort analysis | per-sample clinical pipelines |

Neither shape is universally correct. Joint VCFs are the output of
population-scale variant calling (1000G, gnomAD); per-sample VCFs are
the output of single-patient clinical pipelines. The pre-processing
row above is a *consequence* of the input-shape difference, not
TileDB-VCF being slow.

vcfclick's own ingest configuration sweep is in
[`bench/BENCHMARK.md`](bench/BENCHMARK.md).

## License

Apache License 2.0. Full text in [`LICENSE`](LICENSE); rationale in
[`LICENSING.md`](LICENSING.md).

## TODO

- Phase 2: transcript / exon / CDS hierarchy + corresponding MCP tools.
- LLM-prompt stress-testing of `SCHEMA_DESCRIPTION` against a real
  client. The MCP transport + tool wiring is covered by
  `tests/test_mcp_server.py`; whether the prompt actually steers a
  model away from common mistakes (NULL GQ/DP traps, sparse-table
  joins) needs real LLM runs to confirm.
