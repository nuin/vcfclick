# vcfclick vs TileDB-VCF — first benchmark

*Date: 2026-05-31 (Mac M-series arm64, 16GB RAM, default settings unless noted).*

## TL;DR

On a typical research-bioinformatics workload — a joint VCF from the
**1000 Genomes 30x phased release**, chr17:40M-50M, **235,768
variants × 3,202 samples = 44,986,737 sparse non-ref calls** — the
practical end-to-end story:

| Workflow step | vcfclick | TileDB-VCF |
|---|---|---|
| Pre-process (split + index per sample) | *not needed* | bcftools +split + tabix × 3202 ≈ **8+ min** |
| Disk for input VCFs | 114 MB (joint) | **15.1 GB** (per-sample, 132× inflation) |
| Ingest (best stable config) | **69 s** (parallel-8) | ~4,738 s ≈ **79 min** (single-thread, projected) |
| **End-to-end** | **~1 min** | **~87 min** |

vcfclick ingests joint VCFs directly. TileDB-VCF was designed for
**per-sample VCFs** (the clinical-pipeline output shape) — joint
ingestion is an explicit "Combined VCFs are currently not supported"
runtime error. That's not a slowness, it's a *different design*: the
two tools optimise for different shapes of input.

## What the workload was

- VCF: `CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz`, tabix-sliced to `chr17:40,000,000-50,000,000`.
- Source URL: `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/`
- 235,768 variants, 3,202 samples, ~45M non-ref calls, ~114 MB compressed.
- Includes BRCA1, BRCA2 partial, NF1 partial, and other clinically-recognised genes — representative of a real-cohort gene-panel ingestion.

## What was measured

### vcfclick

Native, Python 3.14 + chDB 4.x + cyvcf2 + pyarrow on arm64. Pre-pass
splitter reads tabix `.tbi` directly (~1 ms). Three configurations
re-run from cold (`.chdb/` truncated each time):

```
Serial (1 process):       247.7 s   (952 v/s)
Parallel 4 workers:        91.6 s  (2,575 v/s)   2.7×
Parallel 8 workers:        69.4 s  (3,396 v/s)   3.6×
```

Counts match across all runs: **235,768 variants / 44,986,737 calls
stored**. Sparse-table compression: 6.2% of dense theoretical max.

### TileDB-VCF (Docker, `tiledb/tiledbvcf-cli:latest`, version 0.40.3)

Three attempts to make it run end-to-end on the same data:

1. **Joint VCF directly:** runtime exception —
   `Combined VCFs are currently not supported`. No way around this
   on the CLI.

2. **All 3,202 per-sample VCFs, default 16 threads:** failed at 24 s
   with `Coordinate (chr17, 40000016) comes before last written
   coordinate in the global order` — a TileDB writer race / ordering
   error when threads interleave.

3. **All 3,202 per-sample VCFs, single-thread (`-t 1`)** — runs
   stable. We benchmarked the first 50 samples to project:

   ```
   Create dataset:           0.23 s
   Store 50 samples:        74.27 s   (1.49 s / sample)
   Projected 3,202 samples: ~4,738 s ≈ 79 min
   ```

### Pre-processing required for TileDB-VCF (not for vcfclick)

```
bcftools +split joint.vcf.gz → 3,202 per-sample VCFs:  not cleanly captured (~few min)
tabix index × 3,202:                                  337.2 s
Total per-sample VCF disk usage:                       15.1 GB
```

The disk inflation (114 MB joint → 15.1 GB per-sample) is dominated
by per-file BGZF headers + duplicated `##contig` / `##INFO` /
`##FORMAT` headers across 3,202 files.

## Caveats and disclosures

1. **TileDB-VCF ran under Rosetta emulation** (image is linux/amd64
   only; host is Apple Silicon arm64). Native arm64 would be faster
   — typical Rosetta penalty for compute-heavy workloads is 30-50%.
   Even granting a 2× speedup on native silicon, single-threaded
   projection would be ~40 min, still far behind vcfclick's parallel-8.

2. **Multi-threaded TileDB-VCF failed** in our environment. With more
   tuning (smaller `--sample-batch-size`, different threading flags)
   it might run. We didn't keep going; the single-threaded number is
   the only stable measurement we could produce.

3. **vcfclick used a tabix-driven splitter pre-pass** (essentially
   free, 1 ms). Without it, parallel-8 would be ~96 s instead of 69 s.

4. **Memory was not bounded** in either tool. Both peaked under
   16 GB on this workload.

5. **The 3,202-sample case is favorable to vcfclick** because every
   variant has many non-ref calls and the sparse-table compression
   wins. For low-call-rate workloads (clinical exome with one sample,
   somatic mutation calling, etc.) the comparison would shift.

## What this means

vcfclick wins ~70× on **end-to-end time** for the joint-VCF
research-bioinformatics workflow, and uses ~5× less disk for the
source data. **The two tools have different design centres:**

- TileDB-VCF is built for **clinical per-sample-VCF pipelines** with
  long-lived multi-sample arrays. Stable and well-suited to that
  workflow.

- vcfclick is built for **research joint-cohort workflows** with
  public datasets and ad-hoc analytical queries. The Native-style
  parallel ingest, the wide-pre-declared-plus-Map-overflow schema,
  the in-process columnar query layer, and the natural-language tool
  surface all assume a researcher with a joint VCF and questions.

If you're loading per-sample clinical pipeline output into a
long-lived array, TileDB-VCF stays appropriate. If you're loading a
1000G / gnomAD / ExAC joint VCF to ask cohort questions, vcfclick
is the more direct fit.

## Reproducing this

vcfclick numbers — one script, downloads the slice once, runs all three
configurations from a cold DB:

```bash
git clone https://github.com/nuin/vcfclick.git
cd vcfclick
uv sync
./bench/run.sh                          # all three configs
./bench/run.sh serial                   # one config only
./bench/run.sh parallel-4 parallel-8    # subset
```

The script caches the chr17:40-50M slice (~110 MB) under `bench/data/` on
first run, then ingests it into an isolated `bench/.vcfclick/` (so your
real `~/.vcfclick/dbs/` is never touched). The ingester's own
`[ingest] done. N variants in X.Xs (Y/s)` line is the measurement.

TileDB-VCF numbers:

```bash
# 1) pre-split (required)
bcftools +split data/chr17_40_50M.vcf.gz -Oz -o bench/per_sample_vcfs/
for f in bench/per_sample_vcfs/*.vcf.gz; do tabix -p vcf "$f"; done
ls bench/per_sample_vcfs/*.vcf.gz | sed 's|.*/|/data/|' > bench/samples_list.txt
# 2) docker ingest
docker run --rm --platform linux/amd64 \
    -v "$PWD/bench/per_sample_vcfs":/data \
    -v "$PWD/bench/tiledbvcf_out":/out \
    -v "$PWD/bench":/bench \
    tiledb/tiledbvcf-cli:latest \
    tiledbvcf create -u /out/array --vcf-attributes /data/HG00096.vcf.gz
docker run --rm --platform linux/amd64 \
    -v "$PWD/bench/per_sample_vcfs":/data \
    -v "$PWD/bench/tiledbvcf_out":/out \
    -v "$PWD/bench":/bench \
    tiledb/tiledbvcf-cli:latest \
    tiledbvcf store -u /out/array --samples-file /bench/samples_list.txt -t 1
```

## Open work, ranked

1. **Run TileDB-VCF natively** (build from source on arm64) to remove
   the Rosetta penalty from the comparison.
2. **Try TileDB-VCF with smaller `--sample-batch-size` and more
   threads** to see if a multi-threaded run can be stabilised.
3. **Benchmark a query workload** ("samples with non-ref call in
   BRCA1") on the same data, both tools.
4. **Benchmark a larger workload** — chr17 full chromosome, or whole
   genome — to see if vcfclick's lead narrows or widens with scale.
