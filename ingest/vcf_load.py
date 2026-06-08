"""VCF → chDB ingestion via Parquet staging.

Single-process serial loader. The data path is:

    cyvcf2 record  →  row tuple  →  Parquet batch  →  chDB import

Each batch (BATCH_SIZE variants) is written to a temporary Parquet file,
then bulk-imported with `INSERT INTO t SELECT * FROM file('batch.parquet')`.
That's the fastest way into chDB MergeTree storage — and it shares the
exact code path used by the parallel loader (workers write Parquet
files; main imports the glob).

Schema policy: VCF 4.3 reserved INFO/FORMAT fields and common GATK
fields land in typed columns. Anything else lands in info_extra /
format_extra Maps. See ingest.routing for the routing tables — the
single source of truth for which fields are typed.

Pre-requisite: multi-allelic sites decomposed via
    bcftools norm -m - input.vcf.gz | bgzip > normalised.vcf.gz
"""

from __future__ import annotations

import logging
import tempfile
import time
import uuid
from pathlib import Path

from cyvcf2 import VCF

from ingest._arrow import (
    GENOTYPES_ARROW_SCHEMA,
    INGESTIONS_ARROW_SCHEMA,
    SAMPLES_ARROW_SCHEMA,
    VARIANTS_ARROW_SCHEMA,
    write_parquet,
)
from ingest.routing import classify_header
from ingest.vcf_rows import build_genotype_rows, build_variant_row
from storage import (
    apply_schema,
    db_path,
    get_session,
    ingest_id_lock,
    insert_via_parquet,
    rollback_ingest,
    validate_ingest_id,
)

log = logging.getLogger(__name__)


def _import_parquet(table: str, parquet_path: Path) -> None:
    """Bulk-import a single Parquet file into the active backend's table."""
    from ingest._arrow import TABLE_COLUMNS, column_list_sql
    from storage import parquet_file_expr

    sess = get_session()
    cols = column_list_sql(TABLE_COLUMNS[table])
    # Explicit column lists on both sides so the import is immune to the
    # engine ever shifting Parquet handling from name-based to positional
    # mapping. `parquet_file_expr` renders the backend-specific reader
    # (`file(..., 'Parquet')` for chDB; `read_parquet(...)` for DuckDB)
    # and quotes the path through sql_quote_str.
    sess.query(
        f"INSERT INTO {table} ({cols}) "
        f"SELECT {cols} FROM {parquet_file_expr(str(parquet_path))}"
    )


def _ensure_schema() -> None:
    """Apply the schema if the variants table isn't already there."""
    from storage import table_exists

    if not table_exists("variants"):
        apply_schema()


BATCH_SIZE = 10_000


def ingest(
    vcf_path: str,
    cohort: str,
    ingest_id: str | None = None,
) -> str:
    """Load a normalised VCF into the embedded chDB store.

    Stage-then-commit: the variant loop only writes Parquet files to
    a temp directory. chDB writes (delete prior rows under this
    ingest_id, insert samples, bulk-import the staged Parquets,
    insert into the ingestions catalog) all happen at the end, after
    the full VCF has parsed successfully.

    Replacement semantics: re-running under an existing `ingest_id`
    truly replaces the prior data. The rollback runs only after Phase 1
    (parsing) succeeds, so failures during parsing — bad header,
    multi-allelic record, malformed row, KeyboardInterrupt — leave the
    prior data intact. The narrow remaining window is a chDB-side
    failure DURING the Phase 2 imports themselves (disk full mid-
    import, chDB session crash); a failure there leaves the database
    in the same state as any partial-commit DB system, and the
    `except` arm's rollback still cleans up whatever new rows did
    land. For belt-and-suspenders atomicity across that narrow
    window, ingest into a fresh `ingest_id` and remove the old one
    once the new ingest succeeds.
    """
    if ingest_id is None:
        ingest_id = str(uuid.uuid4())
    validate_ingest_id(ingest_id)

    _ensure_schema()

    # Serialise concurrent ingests under the same (DB, ingest_id) —
    # the file lock blocks a second invocation until the first
    # finishes, so workflow runners firing parallel `vcfclick db
    # ingest` calls can't race on the staging dir, rollback, or
    # bulk-import. See storage.db.ingest_id_lock docstring.
    with ingest_id_lock(ingest_id):
        return _ingest_locked(vcf_path, cohort, ingest_id)


def _prepare_vcf(vcf_path: str):
    """Open and classify before touching chDB."""
    vcf = VCF(vcf_path)
    classification = classify_header(vcf)
    return vcf, classification, list(vcf.samples)


def _log_ingest_plan(
    vcf_path: str,
    cohort: str,
    ingest_id: str,
    samples: list[str],
    classification: dict,
) -> None:
    log.info("[ingest] %s", vcf_path)
    log.info("[ingest] ingest_id: %s", ingest_id)
    log.info("[ingest] cohort:    %s", cohort)
    log.info("[ingest] samples:   %d", len(samples))
    log.info(
        "[ingest] INFO typed=%d → info_extra=%d",
        len(classification["typed_info"]),
        len(classification["extra_info"]),
    )
    log.info(
        "[ingest] FORMAT typed=%d → format_extra=%d",
        len(classification["typed_format"]),
        len(classification["extra_format"]),
    )
    if classification["extra_info"]:
        log.info("[ingest]   info_extra keys: %s", classification["extra_info"])
    if classification["extra_format"]:
        log.info("[ingest]   format_extra keys: %s", classification["extra_format"])


def _write_stage_batch(
    staging_path: Path,
    n_variants: int,
    variants_batch: list[list],
    genotypes_batch: list[list],
) -> None:
    if not variants_batch:
        return
    v_path = staging_path / f"v_{n_variants:09d}.parquet"
    g_path = staging_path / f"g_{n_variants:09d}.parquet"
    write_parquet(variants_batch, VARIANTS_ARROW_SCHEMA, v_path)
    write_parquet(genotypes_batch, GENOTYPES_ARROW_SCHEMA, g_path)
    variants_batch.clear()
    genotypes_batch.clear()


def _stage_vcf(
    vcf,
    vcf_path: str,
    samples: list[str],
    extra_format_fields: list[str],
    ingest_id: str,
    staging_path: Path,
    started: float,
) -> int:
    variants_batch: list[list] = []
    genotypes_batch: list[list] = []
    n_variants = 0

    for variant in vcf:
        if len(variant.ALT) != 1:
            raise ValueError(
                f"Multi-allelic site at {variant.CHROM}:{variant.POS} "
                f"({len(variant.ALT)} ALTs). Re-normalise with: "
                f"bcftools norm -m - {vcf_path} | bgzip > out.vcf.gz"
            )
        variants_batch.append(build_variant_row(variant, ingest_id))
        genotypes_batch.extend(
            build_genotype_rows(variant, samples, extra_format_fields, ingest_id)
        )
        n_variants += 1

        if len(variants_batch) >= BATCH_SIZE:
            _write_stage_batch(
                staging_path, n_variants, variants_batch, genotypes_batch
            )
            elapsed = time.time() - started
            log.info(
                "[ingest] %s variants (%s/s)",
                f"{n_variants:>10,}",
                f"{n_variants / elapsed:>8,.0f}",
            )

    _write_stage_batch(staging_path, n_variants, variants_batch, genotypes_batch)
    return n_variants


def _commit_staged_ingest(
    staging_path: Path,
    vcf_path: str,
    cohort: str,
    ingest_id: str,
    samples: list[str],
    n_variants: int,
) -> None:
    rollback_ingest(ingest_id)
    insert_via_parquet(
        "samples",
        SAMPLES_ARROW_SCHEMA,
        [
            {"ingest_id": ingest_id, "sample_id": s, "cohort": cohort, "sex": None}
            for s in samples
        ],
    )

    for v_path in sorted(staging_path.glob("v_*.parquet")):
        _import_parquet("variants", v_path)
        g_path = staging_path / v_path.name.replace("v_", "g_", 1)
        if g_path.exists() and g_path.stat().st_size > 0:
            _import_parquet("genotypes", g_path)

    insert_via_parquet(
        "ingestions",
        INGESTIONS_ARROW_SCHEMA,
        [
            {
                "ingest_id": ingest_id,
                "cohort": cohort,
                "vcf_path": vcf_path,
                "n_variants": n_variants,
                "n_samples": len(samples),
            }
        ],
    )


def _handle_ingest_failure(commit_started: bool, ingest_id: str) -> None:
    if commit_started:
        log.warning("[ingest] failed mid-commit — rolling back %s", ingest_id)
        try:
            rollback_ingest(ingest_id)
        except Exception as rb_err:
            log.error("[ingest] rollback FAILED: %s", rb_err)
    else:
        log.warning(
            "[ingest] failed during parse — no chDB writes occurred, "
            "prior data under ingest_id=%s preserved",
            ingest_id,
        )


def _ingest_locked(vcf_path: str, cohort: str, ingest_id: str) -> str:
    """Real ingest body — already holds the per-ingest_id file lock."""
    vcf, classification, samples = _prepare_vcf(vcf_path)
    _log_ingest_plan(vcf_path, cohort, ingest_id, samples, classification)

    started = time.time()
    n_variants = 0
    commit_started = False

    try:
        db_dir = db_path()
        with tempfile.TemporaryDirectory(
            prefix="vcfclick_ingest_", dir=str(db_dir.parent)
        ) as staging:
            staging_path = Path(staging)
            n_variants = _stage_vcf(
                vcf,
                vcf_path,
                samples,
                classification["extra_format"],
                ingest_id,
                staging_path,
                started,
            )
            commit_started = True
            _commit_staged_ingest(
                staging_path, vcf_path, cohort, ingest_id, samples, n_variants
            )
    except BaseException:
        _handle_ingest_failure(commit_started, ingest_id)
        raise

    elapsed = time.time() - started
    log.info(
        "[ingest] done. %s variants in %.1fs (%s/s)",
        f"{n_variants:,}",
        elapsed,
        f"{n_variants / max(elapsed, 0.001):,.0f}",
    )
    return ingest_id
