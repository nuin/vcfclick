# Getting Started

This guide gets you from a clean machine to a queryable demo database.

For the browser-only experience, use the live demo instead:
<https://nuin.github.io/vcfclick-demo/>.

## 1. Install The CLI

Recommended:

```bash
uv tool install vcfclick
vcfclick --help
```

Alternative:

```bash
pipx install vcfclick
vcfclick --help
```

For the optional terminal UI:

```bash
uv tool install "vcfclick[tui]"
vcfclick tui --help
```

From a source checkout:

```bash
git clone https://github.com/nuin/vcfclick.git
cd vcfclick
uv sync --extra tui --group dev
uv run vcfclick --help
```

## 2. Pull The Demo Database

The public demo bundle is a pre-built 1000 Genomes BRCA1 cohort.

```bash
vcfclick db pull demo \
  https://github.com/nuin/vcfclick/releases/download/v0.1.0/1000g-brca1-demo.tar.gz
```

Confirm it exists:

```bash
vcfclick db list
vcfclick db info demo
```

## 3. Run Your First Queries

Count variants:

```bash
vcfclick db query demo \
  "SELECT count(*) AS variants FROM variants"
```

Count samples with a non-reference call in the BRCA1 interval:

```bash
vcfclick db query demo \
  "SELECT count(DISTINCT (ingest_id, sample_id)) AS samples
   FROM genotypes
   WHERE chrom = 'chr17'
     AND pos BETWEEN 43044295 AND 43170245"
```

Show high allele-frequency variants from typed VCF INFO fields:

```bash
vcfclick db query demo \
  "SELECT chrom, pos, ref, alt, info_AF
   FROM variants
   WHERE info_AF IS NOT NULL
   ORDER BY info_AF DESC
   LIMIT 10"
```

## 4. Open The TUI

Install the TUI extra if you did not already:

```bash
uv tool install --force "vcfclick[tui]"
```

Open the demo database:

```bash
vcfclick tui --db demo
```

Try:

- `BRCA1`
- `chr17:43044295-43125483`

See [Terminal UI](TUI.md) for pane details and keyboard expectations.

## 5. Ingest Your Own VCF

vcfclick expects biallelic rows. Decompose multi-allelic records first:

```bash
bcftools norm -m - input.vcf.gz | bgzip > normalised.vcf.gz
tabix -p vcf normalised.vcf.gz
```

Preview how fields will be routed:

```bash
vcfclick discover normalised.vcf.gz
```

Create and ingest:

```bash
vcfclick db create my-cohort
vcfclick db ingest my-cohort normalised.vcf.gz \
  --cohort study1 \
  --ingest-id batch_a
```

Inspect:

```bash
vcfclick db info my-cohort
vcfclick db stats my-cohort
```

`db stats` is currently chDB-only. If you selected
`VCFCLICK_BACKEND=duckdb`, use `db info` for basic counts.

The full workflow is in [User Guide](USER_GUIDE.md).

## 6. Where Data Lives

By default:

```text
~/.vcfclick/dbs/<database-name>/
```

Use another home directory:

```bash
VCFCLICK_HOME=/data/vcfclick vcfclick db list
```

## 7. Next Steps

- [User Guide](USER_GUIDE.md) for ingestion, querying, export, bundles.
- [Backends](BACKENDS.md) for chDB vs DuckDB.
- [Schema Reference](SCHEMA.md) for table definitions and query patterns.
- [MCP and annotations](MCP.md) for natural-language SQL workflows.
