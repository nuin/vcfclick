"""Parallel VCF ingestion with index-aware region splitting.

Architecture:

    [worker 1] parse region 1 → variants_*.parquet
                              → genotypes_*.parquet  ┐
    [worker 2] parse region 2 → variants_*.parquet   │  main:
                              → genotypes_*.parquet  │  INSERT INTO t
    [worker N] parse region N → variants_*.parquet   │  SELECT * FROM
                              → genotypes_*.parquet  ┘  file('staging/*.parquet')

Workers have NO chDB dependency — they only need cyvcf2 + pyarrow.
Concurrent writers to one chDB session aren't safe; concurrent writers
to N independent Parquet files are trivial.

Two design choices the earlier version got wrong, fixed here:

1. Region splitting is variant-count-aware, not uniform-position.
   The uniform splitter divided each contig into N equal position
   ranges; when data lives in a dense subregion (BRCA1 panel, exome,
   etc.) all variants land in one range and N-1 workers do nothing.
   `split_by_variant_count()` does a single pre-pass counting variants
   per 100Kb bucket, then greedy-splits each contig into ranges of
   approximately equal variant count. The pre-pass cost (~10s for the
   chr17 test slice) amortizes immediately against the parallelism
   that follows.

2. Workers flush Parquet batches during parsing, not all at the end.
   Buffering 45M genotype rows in a single worker before the final
   write costs ~10GB of Python memory and hits paging. The new worker
   flushes every BATCH_SIZE variants — bounded memory, multiple small
   Parquet files per worker, glob-import the lot.

Bonus: staging Parquet files ARE valid exports. `--keep-staging`
retains them for downstream DuckDB / Snowflake / Spark / Iceberg use.
"""

from __future__ import annotations

import argparse
import shutil
import time
import uuid
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from cyvcf2 import VCF

from ingest._arrow import (
    GENOTYPES_ARROW_SCHEMA,
    INGESTIONS_ARROW_SCHEMA,
    SAMPLES_ARROW_SCHEMA,
    VARIANTS_ARROW_SCHEMA,
    write_parquet,
)
from ingest.vcf_load import (
    BATCH_SIZE,
    build_genotype_rows,
    build_variant_row,
    classify_header,
)
from storage import apply_schema, db_path, get_session, insert_via_parquet

DEFAULT_BUCKET_SIZE = 100_000  # 100Kb position buckets for the splitter
SPARSE_CONTIG_THRESHOLD = 1_000  # skip splitting if a contig has fewer


def split_by_variant_count(
    vcf_path: str,
    n_workers: int,
    bucket_size: int = DEFAULT_BUCKET_SIZE,
) -> list[str]:
    """Single-pass count of variants per position bucket → balanced ranges.

    Cost: one VCF iteration without genotype access (`for variant in vcf`
    only touches the binary record header, not the per-sample arrays).
    For the chr17 10Mb / 235k-variant test slice this is ~10 seconds.
    """
    vcf = VCF(vcf_path)
    counts: dict[str, dict[int, int]] = {}
    for variant in vcf:
        contig = variant.CHROM
        bucket = variant.POS // bucket_size
        c = counts.setdefault(contig, {})
        c[bucket] = c.get(bucket, 0) + 1

    regions: list[str] = []
    for contig in vcf.seqnames:
        bucket_counts = counts.get(contig)
        if not bucket_counts:
            continue
        regions.extend(
            _split_contig_balanced(contig, bucket_counts, n_workers, bucket_size)
        )
    return regions


def _split_contig_balanced(
    contig: str,
    bucket_counts: dict[int, int],
    n_workers: int,
    bucket_size: int,
) -> list[str]:
    sorted_buckets = sorted(bucket_counts.keys())
    if not sorted_buckets:
        return []

    total = sum(bucket_counts.values())
    if total < SPARSE_CONTIG_THRESHOLD or n_workers <= 1:
        start = sorted_buckets[0] * bucket_size + 1
        end = (sorted_buckets[-1] + 1) * bucket_size
        return [f"{contig}:{start}-{end}"]

    target = total / n_workers
    regions: list[str] = []
    cur_start = sorted_buckets[0]
    cur_count = 0

    for i, b in enumerate(sorted_buckets):
        cur_count += bucket_counts[b]
        is_last_region = len(regions) == n_workers - 1
        is_last_bucket = i == len(sorted_buckets) - 1

        # Cut here if we hit the target AND we have remaining workers,
        # OR if this is the last bucket (must emit the final region).
        if (cur_count >= target and not is_last_region) or is_last_bucket:
            start = cur_start * bucket_size + 1
            end = (b + 1) * bucket_size
            regions.append(f"{contig}:{start}-{end}")
            cur_start = b + 1
            cur_count = 0

    return regions


def _worker(args: tuple) -> tuple[str, int, int]:
    """Parse one region, emit Parquet files in BATCH_SIZE batches.

    Returns (region, n_variants, n_batches) — n_batches is informational,
    used to log per-worker progress in the main process log.
    """
    region, vcf_path, ingest_id, staging_dir, extra_format_fields, batch_size = args
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)

    safe_region = region.replace(":", "_").replace("-", "_")
    staging = Path(staging_dir)

    variants_batch: list[list] = []
    genotypes_batch: list[list] = []
    total_variants = 0
    batch_idx = 0

    def flush() -> None:
        nonlocal batch_idx
        if not variants_batch:
            return
        v_path = staging / f"variants_{safe_region}_{batch_idx:04d}.parquet"
        g_path = staging / f"genotypes_{safe_region}_{batch_idx:04d}.parquet"
        write_parquet(variants_batch, VARIANTS_ARROW_SCHEMA, v_path)
        write_parquet(genotypes_batch, GENOTYPES_ARROW_SCHEMA, g_path)
        variants_batch.clear()
        genotypes_batch.clear()
        batch_idx += 1

    for variant in vcf(region):
        if len(variant.ALT) != 1:
            raise ValueError(
                f"Multi-allelic at {variant.CHROM}:{variant.POS}. "
                f"Normalise with bcftools norm -m -."
            )
        variants_batch.append(build_variant_row(variant, ingest_id))
        genotypes_batch.extend(
            build_genotype_rows(variant, samples, extra_format_fields, ingest_id)
        )
        total_variants += 1
        if len(variants_batch) >= batch_size:
            flush()

    flush()
    return region, total_variants, batch_idx


def _ensure_schema() -> None:
    sess = get_session()
    n = (
        sess.query(
            "SELECT count() FROM system.tables "
            "WHERE database = currentDatabase() AND name = 'variants'",
            "CSV",
        )
        .bytes()
        .decode()
        .strip()
    )
    if n == "0":
        apply_schema(Path(__file__).parent.parent / "schema")


def ingest_parallel(
    vcf_path: str,
    cohort: str,
    ingest_id: str | None = None,
    workers: int = 4,
    keep_staging: bool = False,
    staging_dir: str | None = None,
    bucket_size: int = DEFAULT_BUCKET_SIZE,
    batch_size: int = BATCH_SIZE,
) -> str:
    if ingest_id is None:
        ingest_id = str(uuid.uuid4())

    _ensure_schema()
    sess = get_session()

    vcf = VCF(vcf_path)
    classification = classify_header(vcf)
    extra_format_fields = classification["extra_format"]
    samples = list(vcf.samples)

    staging = Path(staging_dir) if staging_dir else db_path() / "staging" / ingest_id
    staging.mkdir(parents=True, exist_ok=True)

    print(f"[parallel-ingest] {vcf_path}")
    print(f"[parallel-ingest] ingest_id: {ingest_id}")
    print(f"[parallel-ingest] cohort:    {cohort}")
    print(f"[parallel-ingest] samples:   {len(samples)}")
    print(f"[parallel-ingest] workers:   {workers}")
    print(f"[parallel-ingest] staging:   {staging}")
    print(
        f"[parallel-ingest] INFO typed={len(classification['typed_info'])} "
        f"→ info_extra={len(classification['extra_info'])}"
    )
    print(
        f"[parallel-ingest] FORMAT typed={len(classification['typed_format'])} "
        f"→ format_extra={len(classification['extra_format'])}"
    )

    # Plan balanced regions. Prefer the tabix index (essentially free)
    # and fall back to a cyvcf2 pre-pass if the .tbi is missing.
    from ingest._tabix import split_via_tbi

    started_split = time.time()
    regions = split_via_tbi(vcf_path, workers, bucket_size)
    if regions is None:
        regions = split_by_variant_count(vcf_path, workers, bucket_size)
        split_source = "cyvcf2 pre-pass"
    else:
        split_source = "tabix .tbi"
    split_elapsed = time.time() - started_split
    print(
        f"[parallel-ingest] split:  {split_elapsed:.2f}s via {split_source} → "
        f"{len(regions)} regions across {len(set(r.split(':')[0] for r in regions))} contigs"
    )

    # Samples table — safe Parquet-staged insert (no string interpolation
    # of VCF-supplied sample IDs).
    insert_via_parquet(
        "samples",
        SAMPLES_ARROW_SCHEMA,
        [
            {"ingest_id": ingest_id, "sample_id": s, "cohort": cohort, "sex": None}
            for s in samples
        ],
    )

    args_list = [
        (r, vcf_path, ingest_id, str(staging), extra_format_fields, batch_size)
        for r in regions
    ]

    started_parse = time.time()
    total = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for region, n, n_batches in pool.map(_worker, args_list):
            total += n
            if n > 0:
                print(
                    f"[parallel-ingest]   {region}: {n:,} variants "
                    f"({n_batches} batches)"
                )
    parse_elapsed = time.time() - started_parse

    started_import = time.time()
    sess.query(
        f"INSERT INTO variants "
        f"SELECT * FROM file('{staging}/variants_*.parquet', 'Parquet')"
    )
    sess.query(
        f"INSERT INTO genotypes "
        f"SELECT * FROM file('{staging}/genotypes_*.parquet', 'Parquet')"
    )
    import_elapsed = time.time() - started_import

    insert_via_parquet(
        "ingestions",
        INGESTIONS_ARROW_SCHEMA,
        [
            {
                "ingest_id": ingest_id,
                "cohort": cohort,
                "vcf_path": vcf_path,
                "n_variants": total,
                "n_samples": len(samples),
            }
        ],
    )

    if not keep_staging:
        shutil.rmtree(staging)
    else:
        print(f"[parallel-ingest] staging Parquet files kept at {staging}")

    total_elapsed = split_elapsed + parse_elapsed + import_elapsed
    print(
        f"[parallel-ingest] split:  {split_elapsed:.1f}s  "
        f"parse: {parse_elapsed:.1f}s ({total / max(parse_elapsed, 0.001):,.0f}/s)  "
        f"import: {import_elapsed:.1f}s"
    )
    print(
        f"[parallel-ingest] total:  {total_elapsed:.1f}s, {total:,} variants "
        f"({total / max(total_elapsed, 0.001):,.0f}/s overall with {workers} workers)"
    )
    return ingest_id


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf_path")
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--ingest-id", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument(
        "--bucket-size",
        type=int,
        default=DEFAULT_BUCKET_SIZE,
        help=f"Splitter position-bucket size (default {DEFAULT_BUCKET_SIZE:,} bp). "
        f"Smaller = finer load balance, larger = cheaper pre-pass.",
    )
    ap.add_argument(
        "--batch-size",
        type=int,
        default=BATCH_SIZE,
        help=f"Variants per worker Parquet batch (default {BATCH_SIZE:,}). "
        f"Bounds per-worker memory.",
    )
    ap.add_argument(
        "--keep-staging",
        action="store_true",
        help="Keep the worker Parquet files for downstream use (they "
        "ARE the export format).",
    )
    ap.add_argument(
        "--staging-dir",
        default=None,
        help="Override staging directory (default: <DB>/staging/<ingest_id>/).",
    )
    args = ap.parse_args()
    ingest_parallel(
        args.vcf_path,
        args.cohort,
        args.ingest_id,
        args.workers,
        args.keep_staging,
        args.staging_dir,
        args.bucket_size,
        args.batch_size,
    )


if __name__ == "__main__":
    main()
