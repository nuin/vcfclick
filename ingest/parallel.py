"""Parallel VCF ingestion with index-aware region splitting.

Architecture:

    [worker 1] parse region 1 → variants_*.parquet
                              → genotypes_*.parquet  ┐
    [worker 2] parse region 2 → variants_*.parquet   │  main:
                              → genotypes_*.parquet  │  INSERT INTO t
    [worker N] parse region N → variants_*.parquet   │  SELECT * FROM
                              → genotypes_*.parquet  ┘  parquet_file_expr(
                                                          'staging/*.parquet')

Workers have NO storage-engine dependency — they only need cyvcf2 and
pyarrow. The only step that touches the active backend is the main
process Phase 2 `INSERT … SELECT FROM <parquet_file_expr>` against the
staged glob, which is rendered per-backend by storage.parquet_file_expr
(`file('p', 'Parquet')` for chDB, `read_parquet('p')` for DuckDB).
The path works identically on either backend; concurrent writers to
one engine session aren't safe, but concurrent writers to N
independent Parquet files are trivial.

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

import logging
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
from ingest.parallel_split import DEFAULT_BUCKET_SIZE, split_by_variant_count
from ingest.routing import classify_header
from ingest.vcf_load import BATCH_SIZE
from ingest.vcf_rows import build_genotype_rows, build_variant_row
from storage import (
    apply_schema,
    db_path,
    get_session,
    ingest_id_lock,
    insert_via_parquet,
    parquet_file_expr,
    rollback_ingest,
    validate_ingest_id,
)

log = logging.getLogger(__name__)


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
    from storage import table_exists

    if not table_exists("variants"):
        apply_schema()


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
    validate_ingest_id(ingest_id)

    # Serialise concurrent ingests under the same (DB, ingest_id) —
    # see ingest.vcf_load.ingest. Without this two parallel runs sharing
    # an ingest_id would race on the default staging dir, rollback_ingest,
    # and the bulk-import glob, producing a mixed corrupt state.
    with ingest_id_lock(ingest_id):
        return _ingest_parallel_locked(
            vcf_path,
            cohort,
            ingest_id=ingest_id,
            workers=workers,
            keep_staging=keep_staging,
            staging_dir=staging_dir,
            bucket_size=bucket_size,
            batch_size=batch_size,
        )


def _ingest_parallel_locked(
    vcf_path: str,
    cohort: str,
    *,
    ingest_id: str,
    workers: int,
    keep_staging: bool,
    staging_dir: str | None,
    bucket_size: int,
    batch_size: int,
) -> str:
    """Real parallel-ingest body — already holds the per-ingest_id
    file lock. See ingest_parallel() for the public API."""

    _ensure_schema()
    sess = get_session()

    # Open + classify BEFORE touching the engine. A corrupt header / bad
    # classification raises here, before the rollback later runs, so a
    # re-ingest under an existing id whose new VCF fails to read leaves
    # the prior rows intact. See ingest.vcf_load.ingest docstring for
    # the full replacement-semantics contract.
    vcf = VCF(vcf_path)
    classification = classify_header(vcf)
    extra_format_fields = classification["extra_format"]
    samples = list(vcf.samples)

    staging = Path(staging_dir) if staging_dir else db_path() / "staging" / ingest_id
    # Wipe a stale staging directory from a previous failed run under
    # the same ingest_id. Without this the bulk-import globs at the end
    # of Phase 2 would pick up the old Parquets alongside the new ones —
    # documented corruption path. Only clear when we chose the staging
    # location ourselves; a user-supplied --staging-dir is theirs.
    if not staging_dir and staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True, exist_ok=True)

    log.info("[parallel-ingest] %s", vcf_path)
    log.info("[parallel-ingest] ingest_id: %s", ingest_id)
    log.info("[parallel-ingest] cohort:    %s", cohort)
    log.info("[parallel-ingest] samples:   %d", len(samples))
    log.info("[parallel-ingest] workers:   %d", workers)
    log.info("[parallel-ingest] staging:   %s", staging)
    log.info(
        "[parallel-ingest] INFO typed=%d → info_extra=%d",
        len(classification["typed_info"]),
        len(classification["extra_info"]),
    )
    log.info(
        "[parallel-ingest] FORMAT typed=%d → format_extra=%d",
        len(classification["typed_format"]),
        len(classification["extra_format"]),
    )

    # Plan balanced regions. Prefer the tabix index (essentially free)
    # and fall back to a cyvcf2 pre-pass if the .tbi is missing.
    from ingest._tabix import split_via_tbi

    started_split = time.time()
    regions = split_via_tbi(vcf_path, workers, bucket_size)
    # Tabix-derived split returns None when no .tbi is present and an
    # empty list when the index exists but the linear-index buckets are
    # too sparse to balance (typical for tiny VCFs — a 5-record fixture
    # may fit in a single 16Kb bucket with no measurable density to
    # split on). Both cases fall back to the cyvcf2 count pre-pass,
    # which can split any non-empty VCF into at least one region.
    if not regions:
        regions = split_by_variant_count(vcf_path, workers, bucket_size)
        split_source = "cyvcf2 pre-pass"
    else:
        split_source = "tabix .tbi"
    split_elapsed = time.time() - started_split
    log.info(
        "[parallel-ingest] split:  %.2fs via %s → %d regions across %d contigs",
        split_elapsed,
        split_source,
        len(regions),
        len({r.split(":")[0] for r in regions}),
    )

    # Tracks whether Phase 2 (engine writes) have started — see the same
    # flag in ingest.vcf_load.ingest. Determines whether the except arm
    # rolls back at all.
    commit_started = False

    try:
        # ── Phase 1: workers stage Parquets to disk. NO engine writes. ──
        # If any worker raises (multi-allelic record, malformed row,
        # subprocess crash), pool.map surfaces the exception and we
        # exit the try block without ever touching the engine — prior data
        # under this ingest_id stays intact. Matches the stage-then-
        # commit contract documented in ingest.vcf_load.ingest.
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
                    log.info(
                        "[parallel-ingest]   %s: %s variants (%d batches)",
                        region,
                        f"{n:,}",
                        n_batches,
                    )
        parse_elapsed = time.time() - started_parse

        # ── Phase 2: workers succeeded. Commit atomically. ──
        # Wipe prior rows under this ingest_id, then bulk-import the
        # staged Parquets (variants + genotypes), then write samples
        # and the ingestions catalog. An engine-side failure during any
        # of these still wipes prior data but the window is narrow.
        commit_started = True
        rollback_ingest(ingest_id)

        insert_via_parquet(
            "samples",
            SAMPLES_ARROW_SCHEMA,
            [
                {"ingest_id": ingest_id, "sample_id": s, "cohort": cohort, "sex": None}
                for s in samples
            ],
        )

        started_import = time.time()
        # Explicit column lists on both sides — same rationale as
        # ingest.vcf_load._import_parquet: immune to either engine ever shifting
        # Parquet imports from name-based to positional mapping.
        from ingest._arrow import (
            GENOTYPES_COLUMNS,
            VARIANTS_COLUMNS,
            column_list_sql,
        )

        v_cols = column_list_sql(VARIANTS_COLUMNS)
        g_cols = column_list_sql(GENOTYPES_COLUMNS)
        # SQL-quote the glob paths — `staging_dir` can come from user
        # input (the ingest_parallel public API exposes it), so a
        # single quote in the path would otherwise close the file()
        # literal and let an attacker inject arbitrary SQL.
        v_expr = parquet_file_expr(f"{staging}/variants_*.parquet")
        g_expr = parquet_file_expr(f"{staging}/genotypes_*.parquet")
        sess.query(f"INSERT INTO variants ({v_cols}) SELECT {v_cols} FROM {v_expr}")
        sess.query(f"INSERT INTO genotypes ({g_cols}) SELECT {g_cols} FROM {g_expr}")
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
    except BaseException:
        # Only roll back if Phase 2 (engine writes) actually started. If
        # a worker raised during Phase 1, nothing was written to the engine
        # and rolling back would destroy prior data under the same
        # ingest_id and silently turn a "failed re-ingest" into a wipe.
        if commit_started:
            log.warning(
                "[parallel-ingest] failed mid-commit — rolling back %s",
                ingest_id,
            )
            try:
                rollback_ingest(ingest_id)
            except Exception as rb_err:
                log.error("[parallel-ingest] rollback FAILED: %s", rb_err)
        else:
            log.warning(
                "[parallel-ingest] failed during parse — no engine writes "
                "occurred, prior data under ingest_id=%s preserved",
                ingest_id,
            )
        raise

    if not keep_staging:
        shutil.rmtree(staging)
    else:
        log.info("[parallel-ingest] staging Parquet files kept at %s", staging)

    total_elapsed = split_elapsed + parse_elapsed + import_elapsed
    log.info(
        "[parallel-ingest] split:  %.1fs  parse: %.1fs (%s/s)  import: %.1fs",
        split_elapsed,
        parse_elapsed,
        f"{total / max(parse_elapsed, 0.001):,.0f}",
        import_elapsed,
    )
    log.info(
        "[parallel-ingest] total:  %.1fs, %s variants (%s/s overall with %d workers)",
        total_elapsed,
        f"{total:,}",
        f"{total / max(total_elapsed, 0.001):,.0f}",
        workers,
    )
    return ingest_id


# Library module — invoke via `vcfclick db ingest <name> <vcf>`.
# The public CLI lives in cli/main.py.
