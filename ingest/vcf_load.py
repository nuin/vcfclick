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
    GENOTYPES_COLUMNS,
    INGESTIONS_ARROW_SCHEMA,
    SAMPLES_ARROW_SCHEMA,
    VARIANTS_ARROW_SCHEMA,
    VARIANTS_COLUMNS,
    write_parquet,
)
from ingest.routing import FORMAT_PAIR, FORMAT_SCALAR, FORMAT_TRIPLE
from ingest.routing import INFO_FLAG, INFO_PAIR, INFO_SCALAR, classify_header
from storage import (
    apply_schema,
    db_path,
    get_session,
    insert_via_parquet,
    rollback_ingest,
    validate_ingest_id,
)

log = logging.getLogger(__name__)


def _stringify(value) -> str:
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    return str(value)


def _format_arr(variant, field):
    """Read a FORMAT field as a per-sample array via cyvcf2, or None."""
    try:
        arr = variant.format(field)
    except KeyError:
        return None
    return arr


def _cell(arr, i):
    """Extract sample `i`'s value from a FORMAT array, unwrapping (n,1) shapes."""
    v = arr[i]
    if hasattr(v, "__len__") and not isinstance(v, (str, bytes)) and len(v) == 1:
        v = v[0]
    return v


def _scalar_int(arr, i):
    if arr is None:
        return None
    try:
        v = int(_cell(arr, i))
    except (IndexError, TypeError, ValueError):
        return None
    return v if v >= 0 else None


def _pair_int(arr, i):
    if arr is None:
        return (None, None)
    try:
        v = arr[i]
        if len(v) < 2:
            return (None, None)
        a, b = int(v[0]), int(v[1])
    except (IndexError, TypeError, ValueError):
        return (None, None)
    return (a if a >= 0 else None, b if b >= 0 else None)


def _triple_int(arr, i):
    if arr is None:
        return (None, None, None)
    try:
        v = arr[i]
        if len(v) < 3:
            return (None, None, None)
        a, b, c = int(v[0]), int(v[1]), int(v[2])
    except (IndexError, TypeError, ValueError):
        return (None, None, None)
    return (
        a if a >= 0 else None,
        b if b >= 0 else None,
        c if c >= 0 else None,
    )


def _triple_float(arr, i):
    if arr is None:
        return (None, None, None)
    try:
        v = arr[i]
        if len(v) < 3:
            return (None, None, None)
        a, b, c = float(v[0]), float(v[1]), float(v[2])
    except (IndexError, TypeError, ValueError):
        return (None, None, None)
    # cyvcf2 uses NaN as the float-missing sentinel.
    return (
        None if a != a else a,
        None if b != b else b,
        None if c != c else c,
    )


def build_variant_row(variant, ingest_id: str) -> list:
    """One row for the variants table. Caller has verified len(ALT) == 1."""
    info = dict(variant.INFO)

    row: dict = {col: None for col in VARIANTS_COLUMNS}
    row["ingest_id"] = ingest_id
    row["chrom"] = variant.CHROM
    row["pos"] = variant.POS
    row["ref"] = variant.REF
    row["alt"] = variant.ALT[0]
    row["vcf_id"] = variant.ID
    row["qual"] = variant.QUAL
    # cyvcf2 conflates PASS and '.' as None for FILTER. Pass through —
    # the filter column is Nullable.
    row["filter"] = variant.FILTER

    extra: dict[str, str] = {}
    for field, value in info.items():
        if field in INFO_SCALAR:
            row[INFO_SCALAR[field]] = value
        elif field in INFO_PAIR:
            ref_col, alt_col = INFO_PAIR[field]
            if isinstance(value, (list, tuple)) and len(value) >= 2:
                row[ref_col] = value[0]
                row[alt_col] = value[1]
        elif field in INFO_FLAG:
            row[INFO_FLAG[field]] = 1 if value else 0
        else:
            if value is not None:
                extra[field] = _stringify(value)

    # Flags are non-Nullable in the schema; default any unset to 0.
    for flag_col in INFO_FLAG.values():
        if row[flag_col] is None:
            row[flag_col] = 0

    row["info_extra"] = extra
    return [row[c] for c in VARIANTS_COLUMNS]


# cyvcf2 gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT.
# Sparse table convention: skip HOM_REF and UNKNOWN.
GT_ENCODE = {1: 1, 3: 2}


def build_genotype_rows(
    variant, samples: list[str], extra_format_fields: list[str], ingest_id: str
) -> list[list]:
    gt_types = variant.gt_types
    gt_phases = variant.gt_phases

    # Read every typed FORMAT field through variant.format(). The cyvcf2
    # shortcuts (gt_quals/gt_depths/gt_ref_depths/gt_alt_depths) silently
    # return -1 when the FORMAT column ordering varies between records in
    # the same VCF — observed on the routing-test fixture, where a later
    # record with `GT:GQ:DP` came after one with `GT:GQ:DP:AD:PL:...`.
    scalar_arrs = {f: _format_arr(variant, f) for f in FORMAT_SCALAR}
    pair_arrs = {f: _format_arr(variant, f) for f in FORMAT_PAIR}
    triple_arrs = {f: _format_arr(variant, f) for f in FORMAT_TRIPLE}

    extra_arrays: dict[str, object] = {}
    for f in extra_format_fields:
        arr = _format_arr(variant, f)
        if arr is not None:
            extra_arrays[f] = arr

    rows = []
    for i, sample_id in enumerate(samples):
        encoded = GT_ENCODE.get(int(gt_types[i]))
        if encoded is None:
            continue

        row: dict = {col: None for col in GENOTYPES_COLUMNS}
        row["ingest_id"] = ingest_id
        row["chrom"] = variant.CHROM
        row["pos"] = variant.POS
        row["ref"] = variant.REF
        row["alt"] = variant.ALT[0]
        row["sample_id"] = sample_id
        row["gt"] = encoded
        row["phased"] = 1 if gt_phases[i] else 0

        for src, col in FORMAT_SCALAR.items():
            if col == "ft":
                # ft is a string; handled below.
                continue
            row[col] = _scalar_int(scalar_arrs[src], i)

        # ft: per-sample FILTER string. cyvcf2 returns it as ndarray of
        # bytes or list of str; treat "." and "" as missing.
        ft_arr = scalar_arrs.get("FT")
        if ft_arr is not None:
            try:
                v = _cell(ft_arr, i)
                if isinstance(v, bytes):
                    v = v.decode()
                s = str(v).strip() if v is not None else ""
                if s and s != ".":
                    row["ft"] = s
            except (IndexError, TypeError):
                pass

        for src, (ref_col, alt_col) in FORMAT_PAIR.items():
            row[ref_col], row[alt_col] = _pair_int(pair_arrs[src], i)

        for src, (c1, c2, c3) in FORMAT_TRIPLE.items():
            if src == "GL":
                row[c1], row[c2], row[c3] = _triple_float(triple_arrs[src], i)
            else:
                row[c1], row[c2], row[c3] = _triple_int(triple_arrs[src], i)

        extra: dict[str, str] = {}
        for f, arr in extra_arrays.items():
            try:
                extra[f] = _stringify(arr[i])
            except (IndexError, TypeError):
                continue
        row["format_extra"] = extra

        rows.append([row[c] for c in GENOTYPES_COLUMNS])

    return rows


def _import_parquet(table: str, parquet_path: Path) -> None:
    """Bulk-import a single Parquet file into a chDB table."""
    from ingest._arrow import TABLE_COLUMNS, column_list_sql

    sess = get_session()
    cols = column_list_sql(TABLE_COLUMNS[table])
    # Explicit column lists on both sides so the import is immune to
    # chDB ever shifting Parquet handling from name-based to positional
    # mapping. `file()` resolves relative to chDB's user_files directory
    # by default; an absolute path bypasses that. Single-quote the path
    # for SQL safety (Parquet paths aren't expected to contain quotes).
    sess.query(
        f"INSERT INTO {table} ({cols}) "
        f"SELECT {cols} FROM file('{parquet_path}', 'Parquet')"
    )


def _ensure_schema() -> None:
    """Apply the schema if the variants table isn't already there."""
    sess = get_session()
    result = (
        sess.query(
            "SELECT count() FROM system.tables "
            "WHERE database = currentDatabase() AND name = 'variants'",
            "CSV",
        )
        .bytes()
        .decode()
        .strip()
    )
    if result == "0":
        apply_schema(Path(__file__).parent.parent / "schema")


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

    # Open + classify BEFORE touching chDB. Bad headers / unreadable
    # files raise here, before the rollback below runs, so a re-ingest
    # under an existing id whose new VCF is corrupt leaves the prior
    # rows queryable.
    vcf = VCF(vcf_path)
    classification = classify_header(vcf)
    extra_format_fields = classification["extra_format"]
    samples = list(vcf.samples)

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

    # Tracks whether Phase 2 (chDB writes) has started. Determines
    # whether the except arm runs rollback at all: if Phase 1 raises
    # we've written nothing and rollback would destroy prior data
    # belonging to this ingest_id.
    commit_started = False

    try:
        # Two-phase stage-then-commit. Phase 1 writes Parquet files to
        # a tempdir without touching chDB — if the variant loop raises
        # (multi-allelic, malformed body, KeyboardInterrupt), prior data
        # under this ingest_id is still untouched. Phase 2 runs only if
        # the full VCF parsed successfully: rollback prior rows, then
        # bulk-import the staged Parquets in one go.
        variants_batch: list[list] = []
        genotypes_batch: list[list] = []
        n_variants = 0
        started = time.time()

        # Stage on the same volume as the destination DB rather than
        # the system /tmp. With the stage-then-commit flow the whole
        # VCF's Parquet is held on disk before Phase 2 imports — on
        # systems where /tmp is on a small ramdisk or a separate
        # partition that's smaller than the DB volume, a fresh large
        # ingest would otherwise fail before commit purely from
        # staging-disk exhaustion. The DB volume necessarily has room
        # for the destination data, so it has room for the staging
        # equivalent.
        db_dir = db_path()
        with tempfile.TemporaryDirectory(
            prefix="vcfclick_ingest_", dir=str(db_dir.parent)
        ) as staging:
            staging_path = Path(staging)

            def flush_to_disk() -> None:
                if not variants_batch:
                    return
                v_path = staging_path / f"v_{n_variants:09d}.parquet"
                g_path = staging_path / f"g_{n_variants:09d}.parquet"
                write_parquet(variants_batch, VARIANTS_ARROW_SCHEMA, v_path)
                write_parquet(genotypes_batch, GENOTYPES_ARROW_SCHEMA, g_path)
                variants_batch.clear()
                genotypes_batch.clear()

            # ── Phase 1: parse VCF, stage Parquets, NO chDB writes. ──
            for variant in vcf:
                if len(variant.ALT) != 1:
                    raise ValueError(
                        f"Multi-allelic site at {variant.CHROM}:{variant.POS} "
                        f"({len(variant.ALT)} ALTs). Re-normalise with: "
                        f"bcftools norm -m - {vcf_path} | bgzip > out.vcf.gz"
                    )
                variants_batch.append(build_variant_row(variant, ingest_id))
                genotypes_batch.extend(
                    build_genotype_rows(
                        variant, samples, extra_format_fields, ingest_id
                    )
                )
                n_variants += 1

                if len(variants_batch) >= BATCH_SIZE:
                    flush_to_disk()
                    elapsed = time.time() - started
                    log.info(
                        "[ingest] %s variants (%s/s)",
                        f"{n_variants:>10,}",
                        f"{n_variants / elapsed:>8,.0f}",
                    )

            flush_to_disk()

            # ── Phase 2: parse succeeded. Commit atomically. ──
            # Prior rows under this ingest_id wiped; samples + staged
            # variants + staged genotypes + ingestions catalog all
            # written in one go. Any chDB failure here still leaves
            # prior data gone, but the window is narrow vs. the parse.
            commit_started = True
            rollback_ingest(ingest_id)

            insert_via_parquet(
                "samples",
                SAMPLES_ARROW_SCHEMA,
                [
                    {
                        "ingest_id": ingest_id,
                        "sample_id": s,
                        "cohort": cohort,
                        "sex": None,
                    }
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
    except BaseException:
        # Only roll back if Phase 2 (chDB writes) actually started. If
        # Phase 1 raised (multi-allelic, malformed body, KeyboardInterrupt
        # mid-parse) we wrote nothing to chDB — rolling back here would
        # destroy prior data under the same ingest_id and silently turn
        # a "failed re-ingest" into a wipe. Catches BaseException so
        # Ctrl-C during Phase 2 still cleans up partial commits.
        if commit_started:
            log.warning("[ingest] failed mid-commit — rolling back %s", ingest_id)
            try:
                rollback_ingest(ingest_id)
            except Exception as rb_err:
                # Surface rollback failures separately — the user needs
                # to know the DB is dirty so they can manually `db rm`.
                log.error("[ingest] rollback FAILED: %s", rb_err)
        else:
            log.warning(
                "[ingest] failed during parse — no chDB writes occurred, "
                "prior data under ingest_id=%s preserved",
                ingest_id,
            )
        raise

    elapsed = time.time() - started
    log.info(
        "[ingest] done. %s variants in %.1fs (%s/s)",
        f"{n_variants:,}",
        elapsed,
        f"{n_variants / max(elapsed, 0.001):,.0f}",
    )
    return ingest_id
