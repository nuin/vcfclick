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
format_extra Maps. The routing tables below are the single source of
truth — extend them when promoting fields from overflow to typed.

Pre-requisite: multi-allelic sites decomposed via
    bcftools norm -m - input.vcf.gz | bgzip > normalised.vcf.gz
"""

from __future__ import annotations

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
from storage import apply_schema, get_session, insert_via_parquet

INFO_SCALAR = {
    "AC": "info_AC",
    "AF": "info_AF",
    "AN": "info_AN",
    "DP": "info_DP",
    "MQ": "info_MQ",
    "MQ0": "info_MQ0",
    "NS": "info_NS",
    "BQ": "info_BQ",
    "SB": "info_SB",
    "END": "info_END",
    "CIGAR": "info_CIGAR",
    "AA": "info_AA",
    "QD": "info_QD",
    "FS": "info_FS",
    "SOR": "info_SOR",
    "MQRankSum": "info_MQRankSum",
    "ReadPosRankSum": "info_ReadPosRankSum",
    "ExcessHet": "info_ExcessHet",
    "InbreedingCoeff": "info_InbreedingCoeff",
    "MLEAC": "info_MLEAC",
    "MLEAF": "info_MLEAF",
    "BaseQRankSum": "info_BaseQRankSum",
    "ClippingRankSum": "info_ClippingRankSum",
}

INFO_PAIR = {
    "AD": ("info_AD_ref", "info_AD_alt"),
}

INFO_FLAG = {
    "SOMATIC": "info_SOMATIC",
    "VALIDATED": "info_VALIDATED",
    "DB": "info_DB",
    "H2": "info_H2",
    "H3": "info_H3",
    "1000G": "info_1000G",
}

FORMAT_SCALAR = {
    "GQ": "gq",
    "DP": "dp",
    "MQ": "mq",
    "FT": "ft",
    "PS": "ps",
    "PQ": "pq",
}

FORMAT_PAIR = {
    "AD": ("ad_ref", "ad_alt"),
    "ADF": ("adf_ref", "adf_alt"),
    "ADR": ("adr_ref", "adr_alt"),
}

FORMAT_TRIPLE = {
    "PL": ("pl_ref_ref", "pl_ref_alt", "pl_alt_alt"),
    "GL": ("gl_ref_ref", "gl_ref_alt", "gl_alt_alt"),
}


def _typed_format_fields() -> set[str]:
    return {"GT"} | FORMAT_SCALAR.keys() | FORMAT_PAIR.keys() | FORMAT_TRIPLE.keys()


def classify_header(vcf: VCF) -> dict[str, list[str]]:
    typed_info, extra_info = [], []
    typed_format, extra_format = [], []
    typed_format_ids = _typed_format_fields()

    for h in vcf.header_iter():
        try:
            d = h.info(extra=True)
        except Exception:
            continue
        kind = str(d.get("HeaderType", "")).upper()
        field = d.get("ID")
        if not field:
            continue
        if kind == "INFO":
            if field in INFO_SCALAR or field in INFO_PAIR or field in INFO_FLAG:
                typed_info.append(field)
            else:
                extra_info.append(field)
        elif kind == "FORMAT":
            if field in typed_format_ids:
                typed_format.append(field)
            else:
                extra_format.append(field)

    return {
        "typed_info": sorted(typed_info),
        "extra_info": sorted(extra_info),
        "typed_format": sorted(typed_format),
        "extra_format": sorted(extra_format),
    }


def _stringify(value) -> str:
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    return str(value)


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
    gqs = variant.gt_quals
    dps = variant.gt_depths
    ad_refs = variant.gt_ref_depths
    ad_alts = variant.gt_alt_depths

    extra_arrays: dict[str, object] = {}
    for f in extra_format_fields:
        try:
            arr = variant.format(f)
        except KeyError:
            continue
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

        if gqs is not None:
            v = gqs[i]
            if v is not None and v >= 0:
                row["gq"] = int(v)
        if dps is not None:
            v = dps[i]
            if v is not None and v >= 0:
                row["dp"] = int(v)
        if ad_refs is not None:
            v = ad_refs[i]
            if v is not None and v >= 0:
                row["ad_ref"] = int(v)
        if ad_alts is not None:
            v = ad_alts[i]
            if v is not None and v >= 0:
                row["ad_alt"] = int(v)

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
    sess = get_session()
    # `file()` resolves relative to the chDB user_files directory by
    # default; passing an absolute path bypasses that. Wrap in single
    # quotes for SQL safety (paths shouldn't contain quotes in practice).
    sess.query(f"INSERT INTO {table} SELECT * FROM file('{parquet_path}', 'Parquet')")


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
    """Load a normalised VCF into the embedded chDB store."""
    if ingest_id is None:
        ingest_id = str(uuid.uuid4())

    _ensure_schema()

    vcf = VCF(vcf_path)
    classification = classify_header(vcf)
    extra_format_fields = classification["extra_format"]
    samples = list(vcf.samples)

    print(f"[ingest] {vcf_path}")
    print(f"[ingest] ingest_id: {ingest_id}")
    print(f"[ingest] cohort:    {cohort}")
    print(f"[ingest] samples:   {len(samples)}")
    print(
        f"[ingest] INFO typed={len(classification['typed_info'])} "
        f"→ info_extra={len(classification['extra_info'])}"
    )
    print(
        f"[ingest] FORMAT typed={len(classification['typed_format'])} "
        f"→ format_extra={len(classification['extra_format'])}"
    )
    if classification["extra_info"]:
        print(f"[ingest]   info_extra keys: {classification['extra_info']}")
    if classification["extra_format"]:
        print(f"[ingest]   format_extra keys: {classification['extra_format']}")

    # Samples + ingestions catalog go through the same Parquet-staged
    # bulk-insert as variants/genotypes. Avoids string-interpolating
    # sample IDs (which come from the VCF header) into SQL.
    insert_via_parquet(
        "samples",
        SAMPLES_ARROW_SCHEMA,
        [
            {"ingest_id": ingest_id, "sample_id": s, "cohort": cohort, "sex": None}
            for s in samples
        ],
    )

    variants_batch: list[list] = []
    genotypes_batch: list[list] = []
    n_variants = 0
    started = time.time()

    with tempfile.TemporaryDirectory(prefix="vcfclick_ingest_") as staging:
        staging_path = Path(staging)

        def flush() -> None:
            if not variants_batch:
                return
            v_path = staging_path / f"v_{n_variants}.parquet"
            g_path = staging_path / f"g_{n_variants}.parquet"
            write_parquet(variants_batch, VARIANTS_ARROW_SCHEMA, v_path)
            write_parquet(genotypes_batch, GENOTYPES_ARROW_SCHEMA, g_path)
            _import_parquet("variants", v_path)
            if g_path.stat().st_size > 0:
                _import_parquet("genotypes", g_path)
            variants_batch.clear()
            genotypes_batch.clear()

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
                flush()
                elapsed = time.time() - started
                print(
                    f"[ingest] {n_variants:>10,} variants "
                    f"({n_variants / elapsed:>8,.0f}/s)"
                )

        flush()

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

    elapsed = time.time() - started
    print(
        f"[ingest] done. {n_variants:,} variants in {elapsed:.1f}s "
        f"({n_variants / max(elapsed, 0.001):,.0f}/s)"
    )
    return ingest_id


# Library module — invoke via `vcfclick db ingest <name> <vcf> --serial`.
# The public CLI lives in cli/main.py.
