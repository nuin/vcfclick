# vcfclick

A modern VCF database for research labs and bioinformatics teams.
Embedded chDB (ClickHouse engine, no server) for sample data, embedded
DuckDB for reference annotations, and a natural-language query layer
that turns plain English into SQL you can read.

Single binary. `uv run vcfclick`. No Docker, no port, no server, no
Gatekeeper dialog. The wedge demo runs from a clean `git clone`.

Designed as a portfolio + distribution piece for
[Bioinformat](https://bioinformat.org).

Status: Phase 0 — architecture validated against real 1000 Genomes data.

## Why

Two complaints heard repeatedly in research bioinformatics:

1. *"My cohort grew and `bcftools | pandas` stopped scaling."* When
   you have 500+ samples, ad-hoc cohort correlation queries become
   painfully slow. The standard answer is "go install Hail," which is
   correct and operationally expensive.

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
  same MergeTree engines, same projections as full ClickHouse — no
  server.** Persistent on disk under `.chdb/`.
- **DuckDB** holds reference data: RefSeq genes, ClinVar. Embedded,
  swappable, monthly refresh. Never touches sample data.

The MCP server composes across them at query time. Annotation lookups
happen first (DuckDB), then their results parameterise the sample
query (chDB). The chain of reasoning is visible in the UI.

## Ingestion

```bash
# 1. Normalise multi-allelic sites (one-time per VCF)
bcftools norm -m - input.vcf.gz | bgzip > normalised.vcf.gz

# 2. Load
uv run python -m ingest.vcf_load normalised.vcf.gz \
    --cohort demo \
    --ingest-id batch_2026Q2_a
```

The ingester prints a classification of the VCF's INFO/FORMAT fields
on startup — what landed in typed columns vs. the overflow Maps. That
log line is the "adapts to any VCF" claim made literally visible.

**Per-ingestion identity.** Every row carries `ingest_id`. Rows are
NOT merged across uploads — the same `(chrom, pos, ref, alt)` observed
in two different VCFs is two rows, because annotations and QC origin
can differ. Re-running with the same `--ingest-id` is idempotent
(silently replaces prior rows via `ReplacingMergeTree`). Using a new
`--ingest-id` appends.

### Parallel ingest

```bash
uv run python -m ingest.parallel normalised.vcf.gz \
    --cohort demo --ingest-id batch --workers 4
```

Workers parse VCF regions in parallel, write Parquet to staging,
main process bulk-imports via
`INSERT INTO ... SELECT * FROM file('staging/*.parquet')`.
Use `--keep-staging` to retain the worker Parquet files as exports.

The splitter does a single-pass count of variants per 100Kb position
bucket and greedy-splits each contig into ranges of approximately
equal variant count — so dense subregions (gene panels, exomes)
don't leave N-1 workers idle. Workers flush Parquet in batches
during parsing, keeping per-worker memory bounded.

## Export

Parquet is the open-format export — usable from DuckDB, Snowflake,
BigQuery, Spark, or as the storage layer for an Iceberg table:

```bash
# Single table
uv run python -m export.parquet variants /path/out.parquet

# With a filter
uv run python -m export.parquet variants /path/brca1.parquet \
    --where "chrom='chr17' AND pos BETWEEN 43044295 AND 43125483"

# Everything
uv run python -m export.parquet --all /path/output_dir/
```

The parallel ingester's `--keep-staging` flag gives you Parquet
exports as a side effect of ingestion.

## Layout

- `schema/` — ClickHouse DDL (chDB applies it unchanged).
- `storage/db.py` — chDB session singleton; `apply_schema()` helper.
- `ingest/vcf_load.py` — serial cyvcf2-based ingester.
- `ingest/parallel.py` — multi-process variant; Parquet staging.
- `ingest/_arrow.py` — pyarrow schemas matching the ClickHouse tables.
- `export/parquet.py` — table → Parquet export CLI.
- `annotations/db.py` — DuckDB annotation API (gene, ClinVar).
- `annotations/transcripts.py` — transcript/exon/CDS API stubs (Phase 2).
- `mcp/server.py` — MCP server (chDB + DuckDB tool surface).
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

## TileDB-VCF comparison

End-to-end on the same 235k-variant / 3,202-sample workload, native
arm64 (vcfclick) vs Rosetta-emulated linux/amd64 (TileDB-VCF Docker):

| | vcfclick | TileDB-VCF |
|---|---|---|
| Source VCF format | joint VCF ingested directly | per-sample VCFs only ("Combined VCFs are currently not supported") |
| Pre-processing | none | bcftools +split + tabix × 3,202 ≈ 8+ min |
| Source VCF disk | 114 MB | 15.1 GB (132× inflation) |
| Ingest, best stable config | **69 s** (parallel-8) | **~79 min** projected (single-thread, multi-thread failed) |
| End-to-end | **~1 min** | **~87 min** |

Full methodology, caveats (including the Rosetta penalty), and
reproduction commands: [`bench/BENCHMARK.md`](bench/BENCHMARK.md).

## Licensing

Single open-source license, all features included. See
[LICENSING.md](LICENSING.md) for the choice (AGPL vs BSL, TBD) and the
Bioinformat business model around the OSS — hosted SaaS, support
contracts, consulting. No feature gating.

## Open work

- VCF schema auto-discovery utility (`vcf-discover`).
- ClinVar VCF loader under `annotations/loaders/` (the GENCODE gene
  loader is in; ClinVar significance lookup is still stubbed).
- Phase 2: transcript / exon / CDS hierarchy + corresponding MCP tools.
- End-to-end MCP integration test with a real LLM client — the
  `SCHEMA_DESCRIPTION` prompt is theoretical until it's stress-tested.
