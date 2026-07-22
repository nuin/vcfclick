# CLI reference

Every `vcfclick` command, with a one-line summary, an example, and its full option list (verbatim from `--help`). Commands are grouped by task; deeper how-to guides are linked where they exist.

Run `vcfclick --help` or `vcfclick <command> --help` any time for the same information at the terminal. The terminal UI, web UI, and benchmark commands need their optional extras (`vcfclick[tui]`, `[web]`, `[benchmark]`).

## Commands

- **Databases** — [`db create`](#db-create) · [`db list`](#db-list) · [`db info`](#db-info) · [`db path`](#db-path) · [`db rm`](#db-rm) · [`db query`](#db-query) · [`db stats`](#db-stats) · [`db diff`](#db-diff)
- **Ingesting variants** — [`db ingest`](#db-ingest) · [`db ingest-batch`](#db-ingest-batch) · [`merge`](#merge) · [`combine`](#combine) · [`discover`](#discover)
- **Family / trio analysis** — [`db ped`](#db-ped) · [`db trio`](#db-trio) · [`db qc`](#db-qc)
- **Benchmarking** — [`benchmark`](#benchmark) · [`benchmark-cohort`](#benchmark-cohort)
- **Annotations** — [`annotations load`](#annotations-load) · [`annotations load-clinvar`](#annotations-load-clinvar) · [`annotations load-gnomad`](#annotations-load-gnomad)
- **Export & sharing** — [`db dump`](#db-dump) · [`db ingest-parquet`](#db-ingest-parquet) · [`db push`](#db-push) · [`db pull`](#db-pull)
- **Interactive UIs** — [`tui`](#tui) · [`web`](#web)

## Databases

One named database per cohort, under `~/.vcfclick` (or `$VCFCLICK_HOME`). See the [User guide](USER_GUIDE.md).

### db create

Create an empty database with the vcfclick schema.

```bash
vcfclick db create trio
```

<details><summary>options</summary>

```
Usage: vcfclick db create [OPTIONS] NAME

  Create a new empty database with the vcfclick schema applied.

Options:
  --help  Show this message and exit.
```

</details>

### db list

List all named databases.

```bash
vcfclick db list
```

<details><summary>options</summary>

```
Usage: vcfclick db list [OPTIONS]

  List all named databases.

Options:
  --help  Show this message and exit.
```

</details>

### db info

Row counts, ingestions, and metadata for a database.

```bash
vcfclick db info trio
```

<details><summary>options</summary>

```
Usage: vcfclick db info [OPTIONS] NAME

  Show metadata about a database (row counts, ingestions, size).

Options:
  --help  Show this message and exit.
```

</details>

### db path

Print a database's on-disk path (no existence check).

```bash
vcfclick db path trio
```

<details><summary>options</summary>

```
Usage: vcfclick db path [OPTIONS] NAME

  Print the on-disk path of a named database (no checks).

Options:
  --help  Show this message and exit.
```

</details>

### db rm

Permanently delete a database and all its data.

```bash
vcfclick db rm --yes trio
```

<details><summary>options</summary>

```
Usage: vcfclick db rm [OPTIONS] NAME

  Permanently delete a named database and all its data.

Options:
  --yes   Confirm the action without prompting.
  --help  Show this message and exit.
```

</details>

### db query

Run SQL against a database and print the result. See [SCHEMA.md](SCHEMA.md).

```bash
vcfclick db query trio "SELECT count() FROM variants"
```

<details><summary>options</summary>

```
Usage: vcfclick db query [OPTIONS] NAME SQL

  Run a SQL query against a named database and print the result.

Options:
  --format TEXT  chDB output format (PrettyCompact, JSON, CSV, TSV, Vertical,
                 ...).  [default: PrettyCompact]
  --help         Show this message and exit.
```

</details>

### db stats

Schema-population stats for an ingested cohort (chDB backend).

```bash
vcfclick db stats trio
```

<details><summary>options</summary>

```
Usage: vcfclick db stats [OPTIONS] NAME

  Schema-population stats for an ingested cohort.

Options:
  --top INTEGER  Show at most TOP overflow-Map keys per table.  [default: 20]
  --help         Show this message and exit.
```

</details>

### db diff

Per-variant allele-frequency comparison across two cohorts.

```bash
vcfclick db diff trio --cohort-a A --cohort-b B
```

<details><summary>options</summary>

```
Usage: vcfclick db diff [OPTIONS] NAME

  Per-variant allele-frequency comparison across two cohorts.

Options:
  --cohort-a TEXT  First cohort name (as stored on samples.cohort).
                   [required]
  --cohort-b TEXT  Second cohort name.  [required]
  --top INTEGER    Limit to the top N variants by absolute AF difference. 0 =
                   no limit.  [default: 50]
  --format TEXT    chDB output format (PrettyCompact, JSON, CSV, TSV, ...).
                   [default: PrettyCompact]
  --help           Show this message and exit.
```

</details>

## Ingesting variants

Get VCFs into a database. Multi-allelics must be decomposed first (`bcftools norm -m -`).

### db ingest

Ingest one normalized VCF into a database.

```bash
vcfclick db ingest trio joint.vcf.gz --cohort trio
```

<details><summary>options</summary>

```
Usage: vcfclick db ingest [OPTIONS] NAME VCF_PATH

  Ingest a (normalised) VCF into a named database.

Options:
  --cohort TEXT      Cohort label this VCF belongs to.  [default: default]
  --ingest-id TEXT   Stable upload identifier (UUID4 if omitted). Reuse to
                     replace prior data.
  --workers INTEGER  Parallel worker processes for ingestion.  [default: 4]
  --serial           Use the single-process serial ingester instead of
                     parallel.
  --keep-reference   Also store confident hom-reference (0/0) genotype calls,
                     not just non-reference. Needed for defensible trio de-
                     novo analysis (a parent must be provably 0/0, not merely
                     absent). Stores more rows — use for trios/families from a
                     joint-called VCF, not large cohorts. No-calls (./.) are
                     still dropped.
  --help             Show this message and exit.
```

</details>

### db ingest-batch

Ingest many per-sample VCFs as one cohort (from a dir or manifest).

```bash
vcfclick db ingest-batch trio --from-dir vcfs/
```

<details><summary>options</summary>

```
Usage: vcfclick db ingest-batch [OPTIONS] NAME

  Ingest many per-sample VCFs into NAME as one cohort.

Options:
  --from-dir DIRECTORY  Ingest every *.vcf.gz under DIR. ingest_id is the
                        filename stem.
  --manifest FILE       TSV with required `vcf_path` column; optional
                        `sample_id`/`ingest_id` and `cohort` columns override
                        the defaults.
  --cohort TEXT         Default cohort label. Required with --from-dir; used
                        as the fallback for manifest rows that don't carry
                        their own `cohort`.
  --help                Show this message and exit.
```

</details>

### merge

Merge per-sample VCFs into one joint VCF (needs bcftools).

```bash
vcfclick merge a.vcf.gz b.vcf.gz -o joint.vcf.gz
```

### combine

Combine callers' call sets with `set=` provenance (GATK3 CombineVariants). See [COMBINE.md](COMBINE.md).

```bash
vcfclick combine gatk.vcf dv.vcf -o out.vcf --name gatk --name dv
```

### discover

Report which VCF fields land in typed columns vs the overflow map.

```bash
vcfclick discover calls.vcf.gz
```

## Family / trio analysis

Pedigree-aware inheritance analysis. See [Trio](TRIO.md).

### db ped

Load family relationships from a PED/FAM file. See [TRIO.md](TRIO.md).

```bash
vcfclick db ped trio fam.ped
```

<details><summary>options</summary>

```
Usage: vcfclick db ped [OPTIONS] NAME PED_PATH

  Load family relationships from a PED/FAM file into NAME.

  The PED's individual ids must match sample ids already ingested under the
  target ingest_id (v1 assumes a joint-called trio, so all members share one
  ingest_id). Re-loading replaces the prior pedigree for that ingest_id.

Options:
  --ingest-id TEXT  Ingest_id the pedigree's samples belong to. Inferred when
                    the database has exactly one ingestion.
  --help            Show this message and exit.
```

</details>

### db trio

Candidate variants under Mendelian models (de-novo/recessive/dominant/comp-het). See [TRIO.md](TRIO.md).

```bash
vcfclick db trio trio --proband HG002 --category denovo
```

<details><summary>options</summary>

```
Usage: vcfclick db trio [OPTIONS] NAME

  Report candidate variants under Mendelian inheritance models for a trio,
  with genotype quality gates and an AF rarity filter.

Options:
  --proband TEXT                  Sample id of the affected child.  [required]
  --category [denovo|recessive|dominant|comphet|all]
                                  Inheritance model. 'all' prints per-model
                                  candidate counts.  [default: all]
  --min-gq INTEGER                [default: 20]
  --min-dp INTEGER                [default: 10]
  --max-af FLOAT                  Keep variants with population info_AF <=
                                  this (rarity filter).  [default: 0.01]
  --min-ab FLOAT                  [default: 0.25]
  --max-ab FLOAT                  [default: 0.75]
  --gnomad-max-af FLOAT           Additionally drop candidates whose gnomAD
                                  popmax AF exceeds this (needs `vcfclick
                                  annotations load-gnomad`). Variants absent
                                  from the loaded gnomAD slice are kept as
                                  rare.
  --help                          Show this message and exit.
```

</details>

### db qc

Per-sample QC: het/hom, Ti/Tv, chrX-het sex check. See [QC.md](QC.md).

```bash
vcfclick db qc trio
```

<details><summary>options</summary>

```
Usage: vcfclick db qc [OPTIONS] NAME

  Per-sample QC: het/hom ratio, Ti/Tv, and a chrX-het sex check.

Options:
  --format [table|json]  [default: table]
  --help                 Show this message and exit.
```

</details>

## Benchmarking

Concordance benchmarking (hap.py-style) + analytics. See [Benchmarking](BENCHMARK.md).

### benchmark

Benchmark a query VCF against a truth VCF; recall/precision/F1. See [BENCHMARK.md](BENCHMARK.md).

```bash
vcfclick benchmark --truth t.vcf.gz --query q.vcf.gz --ref g.fa --regions c.bed -o report/
```

<details><summary>options</summary>

```
Usage: vcfclick benchmark [OPTIONS]

  Benchmark a query VCF against a truth VCF (normalized genotype concordance).

Options:
  --truth FILE                    Truth (gold-standard) VCF.  [required]
  --query FILE                    Query (call set under test) VCF.  [required]
  --ref FILE                      Reference FASTA (indexed .fai).  [required]
  --regions FILE                  Confident-region BED. Omit to treat every
                                  call as confident.
  -o, --output DIRECTORY          Output directory for reports.  [required]
  --engine [normalized|haplotype|exact]
                                  Reconciliation engine. 'haplotype' adds
                                  local-haplotype rescue; 'exact' is a
                                  diagnostic that skips reference left-
                                  alignment.  [default: normalized]
  --report-formats TEXT           Comma-separated subset of
                                  csv,json,parquet,html (or 'all').  [default:
                                  all]
  --on-ref-mismatch [error|skip]  Behaviour when a REF allele disagrees with
                                  the reference.  [default: error]
  --conf-containment [start|full]
                                  Confident-region membership: variant start,
                                  or whole ref span.  [default: start]
  --decompose-mnp / --no-decompose-mnp
                                  Atomize MNPs into per-position SNPs (loses
                                  phase; off by default).  [default: no-
                                  decompose-mnp]
  --strict                        Promote warnings (e.g. missing --regions) to
                                  hard errors.
  --pass-only / --all             Show only the PASS (or only the ALL) filter
                                  view in the headline.
  --stratify TEXT                 Comma-separated concordance stratification
                                  axes (gnomad,clinvar,gene) joined against
                                  the annotation store; writes
                                  stratified_<axis>.csv.
  --roc                           Write roc.tsv — a query-quality threshold
                                  sweep (recall/precision trade-off per
                                  variant type).
  --strat-region NAME=REGIONS.bed
                                  Stratify concordance by a genome-region set
                                  (e.g. lowcomplex=lc.bed); repeatable. Writes
                                  stratified_regions.csv.
  --audit                         Write fn_annotated.csv / fp_annotated.csv —
                                  each FN/FP joined to its gene, ClinVar
                                  significance, and gnomAD AF.
  --help                          Show this message and exit.
```

</details>

### benchmark-cohort

Benchmark many callers against one truth; per-caller + blind spots. See [BENCHMARK.md](BENCHMARK.md).

```bash
vcfclick benchmark-cohort --truth t.vcf.gz --ref g.fa --caller gatk=g.vcf.gz -o cohort/
```

<details><summary>options</summary>

```
Usage: vcfclick benchmark-cohort [OPTIONS]

  Benchmark several callers against one truth; report per-caller concordance
  and each caller's relative blind spots (variants others catch).

Options:
  --truth FILE                    Truth VCF.  [required]
  --ref FILE                      Reference FASTA.  [required]
  --regions FILE                  Confident-region BED.
  --caller NAME=QUERY.vcf         A named query call set; repeat for each
                                  caller.  [required]
  -o, --output DIRECTORY          [required]
  --engine [normalized|haplotype|exact]
                                  [default: haplotype]
  --on-ref-mismatch [error|skip]
  --history FILE                  Append per-caller metrics to this history
                                  CSV.
  --label TEXT                    Label for the history row (e.g. a pipeline
                                  version).  [default: run]
  --help                          Show this message and exit.
```

</details>

## Annotations

The shared DuckDB reference store (genes, ClinVar, gnomAD). See [MCP and annotations](MCP.md).

### annotations load

Populate gene coordinates from a GENCODE GFF3.

```bash
vcfclick annotations load
```

<details><summary>options</summary>

```
Usage: vcfclick annotations load [OPTIONS]

  Populate the gene-coordinates table from a GENCODE GFF3.

  Required once after `pip install vcfclick` for the MCP server's
  `position_for_gene` tool to resolve symbols → coordinates. Re-run when
  GENCODE releases a new annotation version (yearly-ish).

Options:
  --gff FILE       Local GENCODE GFF3 (.gff3 or .gff3.gz). Downloads v45 from
                   EBI if omitted.
  --keep-existing  Don't truncate refseq_genes before loading (default:
                   replace).
  --help           Show this message and exit.
```

</details>

### annotations load-clinvar

Populate the ClinVar significance table.

```bash
vcfclick annotations load-clinvar
```

<details><summary>options</summary>

```
Usage: vcfclick annotations load-clinvar [OPTIONS]

  Populate the ClinVar significance table from the NCBI ClinVar VCF.

  Required for the MCP server's `clinvar_lookup` tool to return real
  significance calls. The NCBI VCF refreshes weekly; re-run monthly (or before
  any clinically-adjacent demo) to stay current. Bare numeric contigs are
  normalised to `chr`-prefixed during load so lookups against sample data
  (which uses chr-style) compose.

Options:
  --vcf FILE       Local ClinVar VCF (.vcf.gz). Downloads the current NCBI
                   weekly release if omitted.
  --keep-existing  Don't truncate clinvar_variants before loading (default:
                   replace).
  --help           Show this message and exit.
```

</details>

### annotations load-gnomad

Load gnomAD allele frequencies from a sites VCF.

```bash
vcfclick annotations load-gnomad gnomad.vcf.gz
```

<details><summary>options</summary>

```
Usage: vcfclick annotations load-gnomad [OPTIONS] VCF

  Load gnomAD allele frequencies from a gnomAD sites VCF.

  gnomAD is too large to bundle, so pass a VCF you supply — a region slice, an
  af-only file, or a per-chromosome sites VCF. A small region can be pulled
  with tabix-over-HTTPS from the public gnomAD bucket; see docs/MCP.md. Powers
  the `gnomad_lookup` MCP tool and the `db trio --gnomad-max-af` rarity
  filter.

Options:
  --replace  Truncate gnomad_af before loading (default: append, so several
             per-chromosome slices can be loaded incrementally).
  --help     Show this message and exit.
```

</details>

## Export & sharing

Move databases between machines and backends.

### db dump

Export all tables to Parquet files.

```bash
vcfclick db dump trio -o dump/
```

<details><summary>options</summary>

```
Usage: vcfclick db dump [OPTIONS] NAME

  Export all tables from a named database to Parquet files.

Options:
  --out DIRECTORY  Output directory (default: ./<name>-dump/).
  --help           Show this message and exit.
```

</details>

### db ingest-parquet

Ingest a `db dump` Parquet set into a database.

```bash
vcfclick db ingest-parquet trio dump/
```

<details><summary>options</summary>

```
Usage: vcfclick db ingest-parquet [OPTIONS] NAME DUMP_DIR

  Ingest a Parquet dump (produced by `db dump`) into NAME.

Options:
  --cohort TEXT     Cohort label to assign to the imported data.  [default:
                    default]
  --ingest-id TEXT  Stable upload identifier (UUID4 if omitted). Reuse to
                    replace prior data.
  --help            Show this message and exit.
```

</details>

### db push

Bundle a database as a portable tar.gz.

```bash
vcfclick db push trio -o trio.tar.gz
```

<details><summary>options</summary>

```
Usage: vcfclick db push [OPTIONS] NAME OUT_PATH

  Dump a database and bundle it as a portable tar.gz file.

Options:
  --help  Show this message and exit.
```

</details>

### db pull

Restore a database from a tar.gz bundle (file or HTTPS URL).

```bash
vcfclick db pull demo https://.../demo.tar.gz
```

<details><summary>options</summary>

```
Usage: vcfclick db pull [OPTIONS] NAME SOURCE

  Restore a database from a tar.gz bundle (HTTPS URL or local file).

Options:
  --help  Show this message and exit.
```

</details>

## Interactive UIs

Optional extras (`vcfclick[tui]`, `vcfclick[web]`).

### tui

Launch the terminal UI (Locus / Operations / SQL panes). See [TUI.md](TUI.md).

```bash
vcfclick tui --db trio
```

<details><summary>options</summary>

```
Usage: vcfclick tui [OPTIONS]

  Launch the optional terminal UI.

Options:
  --db TEXT  Database to open initially.
  --help     Show this message and exit.
```

</details>

### web

Launch the local browser UI (SQL, NL→SQL, trio, combine). See [WEB.md](WEB.md).

```bash
vcfclick web trio
```

<details><summary>options</summary>

```
Usage: vcfclick web [OPTIONS] NAME

  Launch the optional web UI for database NAME.

Options:
  --port INTEGER  Port to bind.  [default: 8765]
  --host TEXT     Interface to bind. Defaults to localhost; this is a single-
                  user local tool with no authentication — only change this if
                  you understand the exposure.  [default: 127.0.0.1]
  --no-browser    Do not open a browser automatically.
  --help          Show this message and exit.
```

</details>
