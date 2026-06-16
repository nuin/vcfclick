# vcfclick

[![test](https://github.com/nuin/vcfclick/actions/workflows/test.yml/badge.svg)](https://github.com/nuin/vcfclick/actions/workflows/test.yml)
[![PyPI](https://img.shields.io/pypi/v/vcfclick.svg)](https://pypi.org/project/vcfclick/)
[![PyPI downloads](https://img.shields.io/pypi/dm/vcfclick.svg)](https://pypistats.org/packages/vcfclick)
[![GitHub stars](https://img.shields.io/github/stars/nuin/vcfclick.svg?style=flat)](https://github.com/nuin/vcfclick/stargazers)
[![live demo](https://img.shields.io/badge/demo-try%20in%20browser-b45309)](https://nuin.github.io/vcfclick-demo/)

vcfclick turns VCF cohorts into local, queryable SQL databases for research
labs and bioinformatics teams.

- Ingest joint VCFs or batches of per-sample VCFs.
- Query variants, genotypes, samples, and ingestions with SQL.
- Explore cohorts in an optional terminal UI.
- Share databases as portable Parquet bundles.
- Use MCP to let an LLM write visible, auditable SQL.

Status: **research preview**. vcfclick is intended for exploratory
research workflows, not clinical reporting.

## Try It First

**[Try the live browser demo](https://nuin.github.io/vcfclick-demo/)**.

The browser demo runs DuckDB-Wasm over a public 1000 Genomes Parquet
cohort. It does not require installing vcfclick. It is a quick way to
see the core interaction: ask a genomics question, inspect the generated
SQL, and view the result.

The installable CLI is different: it creates local databases under your
`VCFCLICK_HOME` (default `~/.vcfclick`) and can use embedded chDB
(ClickHouse engine) or DuckDB as the storage backend.

## Install

Recommended for CLI use:

```bash
uv tool install vcfclick
vcfclick --help
```

Alternative:

```bash
pipx install vcfclick
vcfclick --help
```

With the optional terminal UI:

```bash
uv tool install "vcfclick[tui]"
vcfclick tui
```

From a source checkout:

```bash
git clone https://github.com/nuin/vcfclick.git
cd vcfclick
uv sync --extra tui --group dev
uv run vcfclick --help
```

vcfclick depends on native Python wheels (`cyvcf2`, `chdb`, `duckdb`,
`pyarrow`). Wheels are available for common macOS arm64 and Linux x86_64
Python versions. If your platform builds `cyvcf2` from source, install
`htslib` development headers first.

## 30-Second CLI Demo

Pull the public BRCA1 demo bundle, then run a SQL query:

```bash
vcfclick db pull demo \
  https://github.com/nuin/vcfclick/releases/download/v0.1.0/1000g-brca1-demo.tar.gz

vcfclick db query demo \
  "SELECT count(DISTINCT (ingest_id, sample_id)) AS samples
   FROM genotypes
   WHERE chrom = 'chr17'
     AND pos BETWEEN 43044295 AND 43170245"
```

Open the same database in the TUI:

```bash
vcfclick tui --db demo
```

## Documentation

Start here:

- [Getting started](docs/GETTING_STARTED.md) - install, pull the demo,
  run first SQL queries, launch the TUI.
- [User guide](docs/USER_GUIDE.md) - create databases, ingest VCFs,
  query, inspect, compare, export, bundle, and restore.
- [Backends](docs/BACKENDS.md) - chDB vs DuckDB, install paths, conda,
  and moving data between backends.
- [Terminal UI](docs/TUI.md) - install the Textual extra and use the
  Locus, Operations, and SQL panes.
- [MCP and annotations](docs/MCP.md) - configure an MCP client, load
  gene/ClinVar references, and use visible LLM-generated SQL.
- [Trio / family analysis](docs/TRIO.md) - merge per-sample VCFs, load
  a pedigree, and report de-novo / recessive / dominant candidates.
- [Schema reference](docs/SCHEMA.md) - table definitions, query
  conventions, sparse genotype rules, and common SQL patterns.
- [FAQ](docs/FAQ.md) - common install, memory, query, backend, and data
  interpretation questions.

Project and contributor docs:

- [Examples](examples/README.md) - worked BRCA1 natural-language SQL
  session.
- [Benchmarks](bench/BENCHMARK.md) - ingest performance measurements.
- [Contributing](CONTRIBUTING.md) - development setup, tests, releases.
- [Citation](docs/CITATION.md) - DOI and BibTeX.
- [License rationale](LICENSING.md) - Apache 2.0 and what it means.

## Core Concepts

### One Database Per Cohort Or Project

The CLI manages named databases under:

```text
~/.vcfclick/dbs/<name>/
```

Set `VCFCLICK_HOME=/path/to/home` if you want databases somewhere else.

### Four Cohort Tables

Every database has the same logical tables:

| Table | Meaning |
|---|---|
| `variants` | one row per `(ingest_id, chrom, pos, ref, alt)` |
| `genotypes` | sparse non-reference sample calls only |
| `samples` | one row per `(ingest_id, sample_id)` |
| `ingestions` | one row per uploaded VCF or imported dump |

The most important rule: **`genotypes` is sparse**. Homozygous-reference
calls (`0/0`) are not stored. See [schema query patterns](docs/SCHEMA.md#common-query-patterns)
before writing allele-frequency or hom-ref queries by hand.

### Backend Choice

vcfclick can run on either backend:

- **chDB**: embedded ClickHouse engine, default when installed, best fit
  for cohort-scale local databases.
- **DuckDB**: embedded single-file backend, useful for conda/Bioconda
  packaging and lightweight environments.

Choose with:

```bash
VCFCLICK_BACKEND=chdb vcfclick db list
VCFCLICK_BACKEND=duckdb vcfclick db list
```

Backends use different on-disk formats. Move data between them with
`vcfclick db dump` and `vcfclick db ingest-parquet`; details are in
[Backends](docs/BACKENDS.md).

## Architecture

vcfclick separates sample/cohort data from reference annotations.

```text
VCF / Parquet input
        |
        v
named vcfclick database
  - variants
  - genotypes
  - samples
  - ingestions
        |
        +-- SQL CLI / TUI
        +-- MCP tools for visible generated SQL

shared annotation store
  - gene coordinates
  - ClinVar significance table
```

Sample data lives in the selected backend for each named database.
Annotation data lives in an embedded DuckDB reference store shared by
the MCP tools.

## Current Limits

- Multi-allelic sites must be decomposed before ingest:
  `bcftools norm -m - input.vcf.gz`.
- vcfclick does not currently auto-load pedigree/sex metadata.
- DuckDB backend support is useful but not identical to chDB support;
  some operations may be chDB-first.
- The natural-language layer is meant to produce visible SQL, not to
  hide SQL from the user.

## License

Apache License 2.0. See [LICENSE](LICENSE) and [LICENSING.md](LICENSING.md).
