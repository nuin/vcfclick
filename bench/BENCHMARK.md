# vcfclick ingest benchmark

*Date: 2026-05-31 (Mac M-series arm64, 16 GB RAM, default settings unless noted).*

## Workload

A joint VCF from the **1000 Genomes 30x phased release**, chr17:40M-50M:

- VCF: `CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz`,
  tabix-sliced to `chr17:40,000,000-50,000,000`.
- Source URL:
  `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20201028_3202_phased/`
- 235,768 variants × 3,202 samples = 44,986,737 sparse non-reference calls.
- ~114 MB bgzip-compressed.
- Includes BRCA1, BRCA2 partial, NF1 partial, and other clinically
  recognised genes — representative of real-cohort gene-panel ingestion.

## vcfclick ingest

Native arm64, Python 3.14 + chDB 4.x + cyvcf2 + pyarrow. The parallel
splitter reads the tabix `.tbi` index directly (~1 ms) to balance
work across workers by variant count. Three configurations re-run
from cold (working DB truncated each time):

```
Serial (1 process):       247.7 s   (952 v/s)
Parallel 4 workers:        91.6 s  (2,575 v/s)   2.7×
Parallel 8 workers:        69.4 s  (3,396 v/s)   3.6×
```

Counts match across all runs: **235,768 variants / 44,986,737 calls
stored**. Sparse-table compression: 6.2% of the dense theoretical
maximum (3,202 samples × 235,768 variants ≈ 754M dense cells vs 45M
stored).

Reproduce:

```bash
git clone https://github.com/nuin/vcfclick.git
cd vcfclick
uv sync
./bench/run.sh                          # all three configs
./bench/run.sh serial                   # one config only
./bench/run.sh parallel-4 parallel-8    # subset
```

The script caches the chr17:40-50M slice (~110 MB) under `bench/data/`
on first run, then ingests it into an isolated `bench/.vcfclick/` so
your real `~/.vcfclick/dbs/` is never touched. The ingester's own
`[ingest] done. N variants in X.Xs (Y/s)` line is the measurement.

## Open work

1. **Query-workload benchmark** — "samples with non-ref call in BRCA1"
   on the same data, comparing per-query latency across vcfclick
   configurations.
2. **Larger workload** — chr17 full chromosome or whole genome — to
   see how vcfclick's per-configuration speedup curve scales.
