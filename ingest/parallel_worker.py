"""Parse one indexed VCF region into independent Parquet batches."""

from pathlib import Path

from cyvcf2 import VCF

from ingest._arrow import (
    GENOTYPES_ARROW_SCHEMA,
    MISSING_GENOTYPES_ARROW_SCHEMA,
    VARIANTS_ARROW_SCHEMA,
    write_parquet,
)
from ingest.vcf_rows import MAX_BATCH_ROWS, build_record_rows


def _worker(args: tuple) -> tuple[str, int, int]:
    """Parse one region, emit Parquet files in BATCH_SIZE batches.

    Returns (region, n_variants, n_batches) — n_batches is informational,
    used to log per-worker progress in the main process log.
    """
    (
        region,
        vcf_path,
        ingest_id,
        staging_dir,
        extra_format_fields,
        batch_size,
        keep_reference,
        record_missing,
    ) = args
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)

    safe_region = region.replace(":", "_").replace("-", "_")
    staging = Path(staging_dir)

    variants_batch: list[list] = []
    genotypes_batch: list[list] = []
    missing_batch: list[list] = []
    total_variants = 0
    batch_idx = 0

    def flush() -> None:
        nonlocal batch_idx
        if not variants_batch:
            return
        v_path = staging / f"variants_{safe_region}_{batch_idx:04d}.parquet"
        g_path = staging / f"genotypes_{safe_region}_{batch_idx:04d}.parquet"
        m_path = staging / f"missing_{safe_region}_{batch_idx:04d}.parquet"
        write_parquet(variants_batch, VARIANTS_ARROW_SCHEMA, v_path)
        write_parquet(genotypes_batch, GENOTYPES_ARROW_SCHEMA, g_path)
        write_parquet(missing_batch, MISSING_GENOTYPES_ARROW_SCHEMA, m_path)
        variants_batch.clear()
        genotypes_batch.clear()
        missing_batch.clear()
        batch_idx += 1

    for variant in vcf(region):
        if len(variant.ALT) != 1:
            raise ValueError(
                f"Multi-allelic at {variant.CHROM}:{variant.POS}. "
                f"Normalise with bcftools norm -m -."
            )
        variant_row, genotype_rows, missing_rows = build_record_rows(
            variant,
            samples,
            extra_format_fields,
            ingest_id,
            keep_reference,
            record_missing,
        )
        variants_batch.append(variant_row)
        genotypes_batch.extend(genotype_rows)
        missing_batch.extend(missing_rows)
        total_variants += 1
        if (
            len(variants_batch) >= batch_size
            or len(genotypes_batch) + len(missing_batch) >= MAX_BATCH_ROWS
        ):
            flush()

    flush()
    return region, total_variants, batch_idx
