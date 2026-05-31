"""Parallel VCF ingestion.

Architecture, simpler now that chDB is in-process:

    [worker 1] parse region 1 → Parquet file
    [worker 2] parse region 2 → Parquet file       →  main:
    [worker 3] parse region 3 → Parquet file          INSERT INTO t
    [worker N] parse region N → Parquet file          SELECT * FROM
                                                      file('staging/*.parquet')

Workers have NO chDB dependency — they only need cyvcf2 + pyarrow.
That keeps the parallel-friendly part (CPU-bound parsing) cleanly
separated from the serial part (chDB write to MergeTree storage).
Concurrent writers to a single chDB session aren't safe; concurrent
writers to N independent Parquet files are trivial.

Bonus: the staging Parquet files ARE valid exports. Add `--keep-staging`
to retain them after import — they'll show up in `<staging>/variants_*.parquet`
and `<staging>/genotypes_*.parquet` for downstream use in DuckDB / Snowflake / Spark.
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
    build_genotype_rows,
    build_variant_row,
    classify_header,
)
from storage import DB_PATH, apply_schema, get_session, insert_via_parquet


DEFAULT_CONTIG_LENGTH = 300_000_000  # > longest human chromosome (~250Mb)


def auto_split_regions(vcf_path: str, n_chunks: int) -> list[str]:
    """Split each contig into n_chunks position ranges. Uses vcf.seqnames
    as the source of truth (always populated); falls back to a wide
    default range when ##contig headers lack `length=` (which the 1000G
    phased files notably do)."""
    vcf = VCF(vcf_path)
    contig_lengths: dict[str, int] = {}
    for h in vcf.header_iter():
        try:
            d = h.info(extra=True)
        except Exception:
            continue
        if str(d.get("HeaderType", "")).upper() != "CONTIG":
            continue
        cid = d.get("ID")
        if not cid:
            continue
        length = d.get("length") or d.get("Length")
        contig_lengths[cid] = int(length) if length else DEFAULT_CONTIG_LENGTH

    regions: list[str] = []
    for name in vcf.seqnames:
        length = contig_lengths.get(name, DEFAULT_CONTIG_LENGTH)
        chunk_size = max(length // n_chunks, 1_000_000)
        start = 1
        while start <= length:
            end = min(start + chunk_size - 1, length)
            regions.append(f"{name}:{start}-{end}")
            start = end + 1
    return regions


def _worker(args: tuple) -> tuple[str, int]:
    """Parse one region and emit two Parquet files (variants, genotypes).

    Workers have no chDB dependency by design — keeps process-pool
    semantics simple and lets the parser parallelise without
    coordinating writers.
    """
    region, vcf_path, ingest_id, staging_dir, extra_format_fields = args
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)

    variants_rows: list[list] = []
    genotypes_rows: list[list] = []

    for variant in vcf(region):
        if len(variant.ALT) != 1:
            raise ValueError(
                f"Multi-allelic at {variant.CHROM}:{variant.POS}. "
                f"Normalise with bcftools norm -m -."
            )
        variants_rows.append(build_variant_row(variant, ingest_id))
        genotypes_rows.extend(
            build_genotype_rows(variant, samples, extra_format_fields, ingest_id)
        )

    safe_region = region.replace(":", "_").replace("-", "_")
    v_path = Path(staging_dir) / f"variants_{safe_region}.parquet"
    g_path = Path(staging_dir) / f"genotypes_{safe_region}.parquet"
    write_parquet(variants_rows, VARIANTS_ARROW_SCHEMA, v_path)
    write_parquet(genotypes_rows, GENOTYPES_ARROW_SCHEMA, g_path)

    return region, len(variants_rows)


def _ensure_schema() -> None:
    sess = get_session()
    n = sess.query(
        "SELECT count() FROM system.tables "
        "WHERE database = currentDatabase() AND name = 'variants'",
        "CSV",
    ).bytes().decode().strip()
    if n == "0":
        apply_schema(Path(__file__).parent.parent / "schema")


def ingest_parallel(
    vcf_path: str,
    cohort: str,
    ingest_id: str | None = None,
    workers: int = 4,
    keep_staging: bool = False,
    staging_dir: str | None = None,
) -> str:
    if ingest_id is None:
        ingest_id = str(uuid.uuid4())

    _ensure_schema()
    sess = get_session()

    vcf = VCF(vcf_path)
    classification = classify_header(vcf)
    extra_format_fields = classification["extra_format"]
    samples = list(vcf.samples)

    regions = auto_split_regions(vcf_path, workers)

    staging = Path(staging_dir) if staging_dir else DB_PATH / "staging" / ingest_id
    staging.mkdir(parents=True, exist_ok=True)

    print(f"[parallel-ingest] {vcf_path}")
    print(f"[parallel-ingest] ingest_id: {ingest_id}")
    print(f"[parallel-ingest] cohort:    {cohort}")
    print(f"[parallel-ingest] samples:   {len(samples)}")
    print(f"[parallel-ingest] workers:   {workers}")
    print(f"[parallel-ingest] regions:   {len(regions)}")
    print(f"[parallel-ingest] staging:   {staging}")
    print(
        f"[parallel-ingest] INFO typed={len(classification['typed_info'])} "
        f"→ info_extra={len(classification['extra_info'])}"
    )
    print(
        f"[parallel-ingest] FORMAT typed={len(classification['typed_format'])} "
        f"→ format_extra={len(classification['extra_format'])}"
    )

    # Samples go through the safe Parquet-staged path (no string
    # interpolation of VCF-supplied sample IDs into SQL).
    insert_via_parquet(
        "samples",
        SAMPLES_ARROW_SCHEMA,
        [
            {"ingest_id": ingest_id, "sample_id": s,
             "cohort": cohort, "sex": None}
            for s in samples
        ],
    )

    args_list = [
        (r, vcf_path, ingest_id, str(staging), extra_format_fields)
        for r in regions
    ]

    # Parse phase — parallel.
    started_parse = time.time()
    total = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for region, n in pool.map(_worker, args_list):
            total += n
            if n > 0:
                print(f"[parallel-ingest]   {region}: {n:,} variants")
    parse_elapsed = time.time() - started_parse

    # Import phase — serial bulk import via glob.
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
        [{
            "ingest_id": ingest_id,
            "cohort": cohort,
            "vcf_path": vcf_path,
            "n_variants": total,
            "n_samples": len(samples),
        }],
    )

    if not keep_staging:
        shutil.rmtree(staging)
    else:
        print(f"[parallel-ingest] staging Parquet files kept at {staging}")

    total_elapsed = parse_elapsed + import_elapsed
    print(
        f"[parallel-ingest] parse:  {parse_elapsed:.1f}s "
        f"({total / max(parse_elapsed, 0.001):,.0f}/s with {workers} workers)"
    )
    print(f"[parallel-ingest] import: {import_elapsed:.1f}s (bulk Parquet → chDB)")
    print(
        f"[parallel-ingest] total:  {total_elapsed:.1f}s, {total:,} variants "
        f"({total / max(total_elapsed, 0.001):,.0f}/s overall)"
    )
    return ingest_id


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf_path")
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--ingest-id", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument(
        "--keep-staging", action="store_true",
        help="Keep the worker Parquet files for downstream use (they "
             "ARE the export format).",
    )
    ap.add_argument(
        "--staging-dir", default=None,
        help="Override staging directory (default: <DB>/staging/<ingest_id>/).",
    )
    args = ap.parse_args()
    ingest_parallel(
        args.vcf_path, args.cohort, args.ingest_id,
        args.workers, args.keep_staging, args.staging_dir,
    )


if __name__ == "__main__":
    main()
