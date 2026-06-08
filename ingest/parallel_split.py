"""Balanced region splitting for parallel VCF ingestion."""

from __future__ import annotations

from cyvcf2 import VCF

DEFAULT_BUCKET_SIZE = 100_000
SPARSE_CONTIG_THRESHOLD = 1_000


def split_by_variant_count(
    vcf_path: str,
    n_workers: int,
    bucket_size: int = DEFAULT_BUCKET_SIZE,
) -> list[str]:
    """Single-pass count of variants per position bucket."""
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
        if bucket_counts:
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

        if (cur_count >= target and not is_last_region) or is_last_bucket:
            start = cur_start * bucket_size + 1
            end = (b + 1) * bucket_size
            regions.append(f"{contig}:{start}-{end}")
            cur_start = b + 1
            cur_count = 0

    return regions
