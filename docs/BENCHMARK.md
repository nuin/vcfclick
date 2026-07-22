# Benchmarking VCFs (`vcfclick benchmark`)

`vcfclick benchmark` compares a **query** VCF against a **truth** VCF over a
reference and a confident-region BED, and reports true/false positives and
negatives with precision, recall, and F1 — the same task as
[Illumina `hap.py`](https://github.com/Illumina/hap.py), reimplemented natively
(no hap.py or `vcfeval` at runtime).

It goes beyond a hap.py-style report by making the concordance **queryable**:
stratify by gnomAD frequency / ClinVar significance / gene / genome region,
sweep quality thresholds, audit every error with its annotation context,
benchmark many callers at once, and track metrics across pipeline versions.

> **Status: research preview.** For research QC and pipeline regression testing,
> not clinical reporting. SNP concordance tracks hap.py to ~0.1–0.2%; INDELs to
> ~0.2% with `--engine haplotype`. See [Validation](VALIDATION.md) for the real
> hap.py conformance runs, and the honest caveats at the end of this page.

## Contents

- [Install](#install)
- [Inputs](#inputs)
- [Quick start](#quick-start)
- [Engines](#engines-which-to-use)
- [Output files](#output-files)
- [Reading the metrics](#reading-the-metrics)
- [Core options](#core-options)
- [Concordance analytics](#concordance-analytics-beyond-happy)
- [Multi-caller cohorts](#multi-caller-cohorts)
- [MCP: ask an LLM about the errors](#mcp-ask-an-llm-about-the-errors)
- [Comparing to hap.py](#comparing-to-happy)
- [Troubleshooting](#troubleshooting)

## Install

The benchmark commands need the optional `benchmark` extra (`pyfaidx`, `numpy`):

```bash
uv tool install "vcfclick[benchmark]"      # or: pipx install "vcfclick[benchmark]"
```

The annotation-driven analytics (`--stratify`, `--audit`) additionally need the
annotation store loaded (genes / ClinVar / gnomAD) — see [MCP and
annotations](MCP.md).

## Inputs

| Input | What | Notes |
|---|---|---|
| `--truth` | Gold-standard VCF | e.g. a GIAB benchmark VCF. Bgzipped or plain. |
| `--query` | The call set under test | The caller output you're scoring. |
| `--ref` | Reference FASTA | Must be **indexed** (`.fai`); `samtools faidx genome.fa`. Contig names are matched with `chr`/no-`chr` and `MT`/`chrM` aliasing. |
| `--regions` | Confident-region BED | **Strongly recommended.** Where truth is complete; calls outside are UNK (not scored), matching hap.py. Omitting it scores everything and is flagged. |

Multi-allelic records are decomposed automatically; indels are left-aligned
against the reference before matching.

## Quick start

```bash
vcfclick benchmark \
  --truth HG002_benchmark.vcf.gz \
  --query my_calls.vcf.gz \
  --ref GRCh38.fa \
  --regions HG002_confident.bed \
  --engine haplotype \
  -o report/
```

Console output is the headline per variant type and filter view:

```
SNP   PASS  recall=0.9985 precision=0.9990 f1=0.9987
SNP   ALL   recall=0.9985 precision=0.9990 f1=0.9987
INDEL PASS  recall=0.9951 precision=0.9948 f1=0.9949
INDEL ALL   recall=0.9951 precision=0.9948 f1=0.9949
reports → report/
```

## Engines — which to use

| `--engine` | What it does | Use when |
|---|---|---|
| `normalized` (default) | Reference-normalize both sides, then a genotype-aware keyed match. Fast. | SNP-grade benchmarking; a conservative INDEL lower bound. |
| `haplotype` | Adds a local-haplotype pass that resolves representation-different calls (an MNP vs two SNPs; a `1/2` het-alt spelled as two `0/1`). | **hap.py-comparable numbers, especially INDELs.** Recommended for most benchmarking. |
| `exact` | Diagnostic: keys on the trimmed representation *without* left-alignment. | Measuring how much reference normalization is buying you. |

## Output files

Written to the `-o` directory:

| File | Contents |
|---|---|
| `summary.csv` | Strict GA4GH-shaped summary (SNP/INDEL × PASS/ALL): TRUTH.TP/FN, QUERY.TP/FP/UNK, Recall/Precision/F1. Drop-in comparable to hap.py's. |
| `vcfclick.summary.csv` | Same rows plus an `Engine` provenance column. |
| `metrics.json` | The summary as JSON. |
| `run_meta.json` | Provenance: engine, inputs, excluded counts, which extras ran. |
| `index.html` | Self-contained report (headline cards + tables), stamped with the engine and an INDEL caveat. |
| `benchmark.parquet` | The full per-variant classified frame (with `--report-formats parquet`) — the substrate for the analytics below. |

Analytics add more files (`stratified_*.csv`, `roc.tsv`, `fn_annotated.csv`, …),
covered below.

## Reading the metrics

- **TP** — a query call that matches a truth call (genotype-aware).
- **FN** — a truth call the query missed (recall = TP / (TP + FN)).
- **FP** — a query call with no truth match (precision = TP / (TP + FP)).
- **UNK** — a query call **outside** the confident region: truth is unknown
  there, so it is *not scored* (as in hap.py). This is why a confident BED
  matters — without one everything is forced in-region.
- **F1** — harmonic mean of precision and recall.
- **PASS vs ALL** — the `PASS` view scores only query calls that pass FILTER
  (a filtered query call becomes an FN there); `ALL` scores every query call.
  Truth is always scored (never filtered), matching hap.py.
- **SNP vs INDEL** — reported separately. Complex/MNP variants that aren't
  SNP/INDEL are surfaced in `run_meta.json` (`unsummarized_types`), never
  silently dropped.

## Core options

| Option | Default | Meaning |
|---|---|---|
| `--engine {normalized,haplotype,exact}` | `normalized` | See [Engines](#engines-which-to-use). |
| `--regions BED` | none | Confident regions. Omit → everything in-region (warned). |
| `--conf-containment {start,full}` | `start` | Is a variant "in-region" by its start position, or must its whole ref span be inside? |
| `--decompose-mnp / --no-decompose-mnp` | off | Atomize MNPs into per-position SNPs (loses phase; leave off for the haplotype engine). |
| `--on-ref-mismatch {error,skip}` | `error` | If a query/truth REF disagrees with the reference: hard error, or drop the record (count surfaced). |
| `--strict` | off | Promote warnings (e.g. missing `--regions`) to errors. |
| `--pass-only / --all` | both | Restrict the console headline to one filter view. |
| `--report-formats` | all | Subset of `csv,json,parquet,html`. Use `parquet` to keep the per-variant frame. |

## Concordance analytics (beyond hap.py)

These make the concordance queryable against annotations and quality. They read
the per-variant frame, so they work on any run.

### Stratify by annotation — `--stratify`

Recall/precision per **gnomAD AF bin**, **ClinVar significance**, or **gene**,
joined from the annotation store. Requires annotations loaded (see [MCP](MCP.md)).

```bash
vcfclick benchmark ... --engine haplotype --stratify gnomad,clinvar,gene -o report/
```

Writes `stratified_gnomad.csv`, `stratified_clinvar.csv`, `stratified_gene.csv`,
each with `stratum, truth_tp, truth_fn, query_tp, query_fp, recall, precision`.
Answers *"what's my recall on rare variants? on ClinVar-pathogenic sites? in
BRCA1?"* — questions hap.py's design can't.

- **gnomad** bins: `novel` (not in gnomAD), `rare` (<0.001), `low` (<0.05), `common`.
- **clinvar**: one row per clinical-significance value.
- **gene**: one row per gene the variants fall in.

### Stratify by genome region — `--strat-region`

Recall/precision per user-supplied region set (low-complexity, segmental
duplications, GA4GH genome stratifications, …). Repeatable.

```bash
vcfclick benchmark ... \
  --strat-region lowcomplex=low_complexity.bed \
  --strat-region segdup=segmental_dups.bed \
  -o report/
```

Writes `stratified_regions.csv`. A variant in overlapping strata counts in each;
variants in no region fall under `none`.

### Quality ROC / PR curve — `--roc`

Sweep a minimum query-quality (QUAL) threshold: as it rises, low-quality query
calls drop out (recall falls, precision rises), tracing the precision-recall
trade-off.

```bash
vcfclick benchmark ... --roc -o report/
```

Writes `roc.tsv` (`Type, Threshold, TP, FP, Recall, Precision`) — feed it to any
plotter (or `happyR`).

### Audit the errors — `--audit`

Every false negative and false positive, joined to its **gene, ClinVar
significance, and gnomAD AF** — so you can see not just *how many* errors but
*which ones and whether they matter*.

```bash
vcfclick benchmark ... --audit -o report/
```

Writes `fn_annotated.csv` and `fp_annotated.csv`
(`chrom, pos, ref, alt, vtype, gene, clin_sig, af`). *"Did my caller miss any
ClinVar-pathogenic variants?"* is now one file away.

## Multi-caller cohorts

`benchmark-cohort` scores **several callers against one truth in a single pass**
— what would be N separate hap.py runs plus a pile of stitched CSVs.

```bash
vcfclick benchmark-cohort \
  --truth HG002_benchmark.vcf.gz \
  --ref GRCh38.fa \
  --regions HG002_confident.bed \
  --caller gatk=gatk.vcf.gz \
  --caller deepvariant=dv.vcf.gz \
  --caller dragen=dragen.vcf.gz \
  --engine haplotype \
  -o cohort/ \
  --history history.csv --label pipeline-v2.3
```

- Prints per-caller recall/precision by type, and each caller's **relative blind
  spots** — variants it misses that other callers recover.
- Writes `per_caller_metrics.csv`.
- With `--history`/`--label`, appends the metrics to a CSV keyed by label so you
  can **track concordance across pipeline versions** (regression testing).

## MCP: ask an LLM about the errors

With an MCP client configured (see [MCP](MCP.md)) and a concordance parquet
(`--report-formats parquet`), the `benchmark_errors` tool returns the annotated
FN/FP rows, so a language model can answer *"which variants did this caller miss,
and do any of them matter?"* directly:

```
benchmark_errors(concordance_parquet="report/benchmark.parquet", kind="FN")
```

## Comparing to hap.py

vcfclick's `haplotype` engine was validated directly against real hap.py
(v0.3.12, `xcmp`) on synthetic representation cases, GIAB HG002 chr20/chr1, and
an **independent** assembly-based caller (GIAB T2T-Q100). SNP and INDEL
concordance track hap.py to **~0.2%**. Full numbers, method, and the honest
residual are in [Validation](VALIDATION.md).

**Honest caveats:**

- Use `--engine haplotype` for hap.py-comparable INDEL numbers; `normalized` is
  a conservative lower bound.
- Always pass `--regions` — the UNK gating is what makes precision meaningful.
- For caller-vs-caller comparison without a gold truth, you are measuring
  *concordance relative to the baseline you name `--truth`*, not accuracy.
- Validated on GRCh38/HG002; not a blanket guarantee for every
  organism/reference. It's a research preview — for regulatory-grade INDEL
  benchmarking, cross-check against hap.py itself.

## Troubleshooting

- **`contig 'chr1' is not in the reference`** — the VCF and FASTA use
  incompatible contig naming and aliasing couldn't reconcile them. Use a
  reference matching your VCFs' build/naming.
- **`REF ... != reference`** — the VCF's REF allele doesn't match the given
  reference (wrong build, or a different reference). Fix the reference, or pass
  `--on-ref-mismatch skip` to drop those records (the skipped count is surfaced).
- **`No module named 'pyfaidx'`** — install the extra: `vcfclick[benchmark]`.
- **Empty `stratified_*.csv` / audit files** — the annotation store isn't loaded
  or isn't reachable; see [MCP and annotations](MCP.md).
- **Precision looks implausibly high without `--regions`** — without a confident
  BED there is no UNK bucket; supply one.
