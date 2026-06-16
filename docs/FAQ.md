# FAQ

## Is vcfclick For Clinical Reporting?

No. vcfclick is a research-preview tool for exploratory cohort queries,
local demos, and bioinformatics workflows. It is not a clinical
reporting system.

## How Do I Install It? Is There A Homebrew Or Conda Package?

Install with `uv tool install vcfclick` or `pipx install vcfclick`;
both pull the native wheels (`cyvcf2`, `chdb`, `duckdb`, `pyarrow`)
automatically and work the same on macOS and Linux.

There is no Homebrew or conda-forge package. `chdb` (the embedded
ClickHouse engine) ships only as a binary wheel with no source build,
so a from-source Homebrew formula isn't possible, and homebrew-core
does not accept binary-only dependency stacks. `uv` / `pipx` are the
supported paths. See [Install](../README.md#install).

## What Kind Of VCF Should I Ingest?

Use normalized, bgzipped VCFs. Multi-allelic sites must be decomposed:

```bash
bcftools norm -m - input.vcf.gz -Oz -o input.norm.vcf.gz
tabix -p vcf input.norm.vcf.gz
```

Then ingest:

```bash
vcfclick db create study
vcfclick db ingest study input.norm.vcf.gz --cohort cases
```

## What's The Difference Between `merge` And `combine`?

They merge different things. `vcfclick merge` joins VCFs with **disjoint
samples** into one multi-sample VCF (it wraps `bcftools merge`) — use it
to assemble a trio or cohort from separate per-sample files. `vcfclick
combine` unions **call sets that may share samples** — e.g. the same
cohort called by two callers — annotating each record with `set=`
provenance and resolving a shared sample by input priority. `combine` is
a reimplementation of GATK3's `CombineVariants` (removed in GATK4). See
[Combining call sets](COMBINE.md).

## Can vcfclick Do Trio / De-Novo Analysis?

Yes — load a pedigree with `vcfclick db ped` and run `vcfclick db trio`
for de-novo, recessive, and dominant candidates. It is candidate
*filtering*, not variant calling, and defensible de-novo needs
`vcfclick db ingest --keep-reference` so a parent's hom-reference is
stored rather than inferred from absence. See [Trio](TRIO.md).

## Why Are There Fewer Genotype Rows Than Samples Times Variants?

`genotypes` is sparse. vcfclick stores non-reference calls only. A
sample missing from `genotypes` at a variant is treated as `0/0` by
convention.

This is intentional: it keeps cohort databases much smaller and makes
non-reference scans fast.

## How Do I Count Homozygous-Reference Samples?

Count total samples from `samples`, count non-reference samples from
`genotypes`, then subtract. Do not search for null rows with a left
join.

See [Common query patterns](SCHEMA.md#common-query-patterns).

## Why Did A `gq >= 20 AND dp >= 10` Query Return Zero Rows?

Many public joint-call releases include genotypes but not per-sample
`GQ` or `DP`. In SQL, comparisons against `NULL` fail, so a quality
filter can remove every row.

Check whether the fields are populated:

```bash
vcfclick db query demo \
  "SELECT count() AS rows,
          sum(CASE WHEN gq IS NOT NULL THEN 1 ELSE 0 END) AS with_gq,
          sum(CASE WHEN dp IS NOT NULL THEN 1 ELSE 0 END) AS with_dp
   FROM genotypes"
```

If `with_gq` and `with_dp` are zero, use raw genotype counts or choose
a quality field that exists in your VCF.

## Should I Query `variants.info_AF` Or Compute AF From Genotypes?

Use `variants.info_AF` for broad ranking and quick scans when the VCF
already carries population-level AF. It avoids aggregating the full
sparse genotype table in the browser or on small machines.

Compute AF from `genotypes` when you need cohort-specific AF for your
loaded samples. Use the denominator pattern in
[Schema reference](SCHEMA.md#cohort-allele-frequency).

## Where Are Databases Stored?

By default:

```text
~/.vcfclick/dbs/<name>/
```

Set `VCFCLICK_HOME` to move them:

```bash
export VCFCLICK_HOME=/data/vcfclick
```

## Can I Share A Database?

Yes. Use backend-neutral bundles:

```bash
vcfclick db push study study.tar.gz
vcfclick db pull restored-study study.tar.gz
```

The bundle contains Parquet files, not backend-specific database files.

## Can I Use DuckDB Instead Of chDB?

Yes:

```bash
VCFCLICK_BACKEND=duckdb vcfclick db create study
```

DuckDB is useful for lightweight installs, and it is **required** for
conda / Bioconda installs because chDB is not packaged for conda —
there a `conda install` gets DuckDB only. chDB remains the default on
PyPI installs and currently supports more operations, including
`vcfclick db stats`.

## How Do I Run The TUI?

Install the optional extra, then run:

```bash
uv tool install "vcfclick[tui]"
vcfclick tui --db demo
```

From source:

```bash
uv sync --extra tui --group dev
uv run vcfclick tui --db demo
```

See [Terminal UI](TUI.md).

## How Do I Use vcfclick With An LLM?

Use the MCP server:

```bash
VCFCLICK_DB_NAME=demo python -m vcfclick_mcp.server
```

Configure your MCP client to launch that command. See
[MCP and annotations](MCP.md).

## Does The Browser Demo Use The Same Backend?

No. The browser demo runs DuckDB-Wasm over Parquet files inside the
browser. The installable CLI runs locally and can use chDB or DuckDB.

The shared concept is the schema and the visible SQL workflow.

## Why Does Ingest Use So Much Disk Temporarily?

Ingest stages typed Parquet before committing rows to the database. This
keeps failed ingests from leaving partial data behind, but it means you
need working disk space during import.

For very large cohorts, ingest in batches, keep `VCFCLICK_HOME` on a
disk with enough free space, and use Parquet bundles for transfer.
