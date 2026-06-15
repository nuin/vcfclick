# User Guide

This page describes the main CLI workflows for local vcfclick databases.

## Database Names And Home Directory

Each database is a named directory under:

```text
~/.vcfclick/dbs/<name>/
```

Set `VCFCLICK_HOME` to isolate projects or avoid your real home
directory during testing:

```bash
VCFCLICK_HOME=/tmp/vcfclick-demo vcfclick db list
```

Database names are managed with:

```bash
vcfclick db create my-cohort
vcfclick db list
vcfclick db info my-cohort
vcfclick db path my-cohort
vcfclick db rm my-cohort
```

## Prepare VCF Input

vcfclick stores one ALT allele per row. Decompose multi-allelic records
before ingest:

```bash
bcftools norm -m - input.vcf.gz | bgzip > normalised.vcf.gz
tabix -p vcf normalised.vcf.gz
```

The `.tbi` index is required for parallel ingestion and useful for
range-aware splitting.

Preview field routing:

```bash
vcfclick discover normalised.vcf.gz
```

`discover` reports which INFO/FORMAT fields land in typed columns and
which go into `info_extra` or `format_extra` maps.

## Ingest One VCF

```bash
vcfclick db create my-cohort

vcfclick db ingest my-cohort normalised.vcf.gz \
  --cohort case \
  --ingest-id batch_a
```

Important options:

| Option | Meaning |
|---|---|
| `--cohort` | logical cohort label stored on `samples.cohort` |
| `--ingest-id` | stable upload identifier; reuse it to replace prior data |
| `--workers` | parallel worker count, default `4` |
| `--serial` | force the single-process ingester |

Rows are not merged across ingestions. The same variant in two
different `ingest_id`s remains two variant rows because QC and INFO
origin can differ.

## Ingest Many Per-Sample VCFs

Directory mode:

```bash
vcfclick db ingest-batch my-cohort \
  --from-dir per_sample_vcfs/ \
  --cohort study1
```

Manifest mode:

```bash
vcfclick db ingest-batch my-cohort \
  --manifest samples.tsv \
  --cohort fallback
```

Manifest columns:

| Column | Required | Meaning |
|---|---:|---|
| `vcf_path` | yes | path to the VCF |
| `sample_id` | no | sample override |
| `ingest_id` | no | upload identifier override |
| `cohort` | no | cohort override |

## Query A Database

Run SQL directly:

```bash
vcfclick db query my-cohort \
  "SELECT count(*) AS variants FROM variants"
```

Use explicit limits when returning rows:

```bash
vcfclick db query my-cohort \
  "SELECT chrom, pos, ref, alt, qual
   FROM variants
   WHERE chrom = 'chr17'
     AND pos BETWEEN 43000000 AND 43200000
   ORDER BY pos
   LIMIT 20"
```

For allele frequency and hom-ref patterns, read
[Schema Reference: Common Query Patterns](SCHEMA.md#common-query-patterns).

## Inspect Stored Data

Basic counts and ingestions:

```bash
vcfclick db info my-cohort
```

Field population and overflow keys:

```bash
vcfclick db stats my-cohort
```

`stats` is the post-ingest complement to `discover`: it tells you what
actually appeared in stored data. It is currently implemented for chDB;
on DuckDB, use `vcfclick db info` for basic counts.

## Compare Cohorts

If two cohort labels exist in the same database:

```bash
vcfclick db diff my-cohort \
  --cohort-a case \
  --cohort-b control \
  --top 50
```

This reports per-variant allele-frequency differences using cohort
sizes from `samples`, not from sparse `genotypes`.

## Export And Import Parquet

Export all tables:

```bash
vcfclick db dump my-cohort --out my-cohort-export/
```

Import into another database:

```bash
vcfclick db create imported
vcfclick db ingest-parquet imported my-cohort-export/ \
  --cohort imported \
  --ingest-id import_a
```

Parquet export is the bridge to DuckDB, Spark, BigQuery, Snowflake,
Iceberg, and backend migration.

## Share A Database Bundle

Create a portable `.tar.gz` bundle:

```bash
vcfclick db push my-cohort my-cohort.tar.gz
```

Restore from a local file:

```bash
vcfclick db pull restored ./my-cohort.tar.gz
```

Restore from HTTPS:

```bash
vcfclick db pull restored \
  https://example.org/path/to/my-cohort.tar.gz
```

Bundles contain Parquet exports, not backend-specific database files,
so they are suitable for sharing and backend migration.

## Work With The TUI

```bash
vcfclick tui --db my-cohort
```

The TUI is optional and requires the `vcfclick[tui]` extra. See
[Terminal UI](TUI.md).

## Use Natural-Language SQL Through MCP

MCP is optional. It lets an MCP client ask vcfclick for schema context,
gene coordinates, and SQL execution while keeping generated SQL visible.

See [MCP and annotations](MCP.md).
