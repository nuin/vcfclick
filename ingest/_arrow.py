"""pyarrow schemas + row-batch → Parquet conversion.

The ingestion data path is:
    cyvcf2 record → row tuple → batched into Arrow Table → Parquet file
    → chDB `INSERT INTO ... SELECT * FROM file('*.parquet', 'Parquet')`

Parquet as the intermediate format means workers in the parallel path
have no chDB dependency (they only need cyvcf2 + pyarrow), AND the
same Parquet files double as the user-facing export format.

Nullability matches the ClickHouse schema exactly. Non-nullable
columns (the UInt8 flags, gt, phased, ingest_id, coords) are declared
`nullable=False` so a programming bug that emits None gets caught at
Arrow conversion time rather than producing a corrupt insert.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pyarrow as pa
import pyarrow.parquet as pq


def _f(name: str, type_, nullable: bool = True) -> pa.Field:
    return pa.field(name, type_, nullable=nullable)


# Order MUST match VARIANTS_COLUMNS in ingest.vcf_load.
VARIANTS_ARROW_SCHEMA = pa.schema(
    [
        _f("ingest_id", pa.string(), False),
        _f("chrom", pa.string(), False),
        _f("pos", pa.uint32(), False),
        _f("ref", pa.string(), False),
        _f("alt", pa.string(), False),
        _f("vcf_id", pa.string()),
        _f("qual", pa.float32()),
        _f("filter", pa.string()),
        _f("info_AC", pa.uint32()),
        _f("info_AF", pa.float32()),
        _f("info_AN", pa.uint32()),
        _f("info_DP", pa.uint32()),
        _f("info_MQ", pa.float32()),
        _f("info_MQ0", pa.uint32()),
        _f("info_NS", pa.uint32()),
        _f("info_BQ", pa.float32()),
        _f("info_SB", pa.float32()),
        _f("info_END", pa.uint32()),
        _f("info_CIGAR", pa.string()),
        _f("info_AA", pa.string()),
        _f("info_QD", pa.float32()),
        _f("info_FS", pa.float32()),
        _f("info_SOR", pa.float32()),
        _f("info_MQRankSum", pa.float32()),
        _f("info_ReadPosRankSum", pa.float32()),
        _f("info_ExcessHet", pa.float32()),
        _f("info_InbreedingCoeff", pa.float32()),
        _f("info_MLEAC", pa.uint32()),
        _f("info_MLEAF", pa.float32()),
        _f("info_BaseQRankSum", pa.float32()),
        _f("info_ClippingRankSum", pa.float32()),
        # DRAGEN-specific record-level INFO scalars.
        _f("info_FractionInformativeReads", pa.float32()),
        _f("info_HAPCOMP", pa.uint32()),
        _f("info_HAPDOM", pa.float32()),
        _f("info_DragenSnvHardQUAL", pa.float32()),
        _f("info_DragenIndelHardQUAL", pa.float32()),
        _f("info_AD_ref", pa.uint32()),
        _f("info_AD_alt", pa.uint32()),
        _f("info_SOMATIC", pa.uint8(), False),
        _f("info_VALIDATED", pa.uint8(), False),
        _f("info_DB", pa.uint8(), False),
        _f("info_H2", pa.uint8(), False),
        _f("info_H3", pa.uint8(), False),
        _f("info_1000G", pa.uint8(), False),
        _f("info_extra", pa.map_(pa.string(), pa.string()), False),
    ]
)


# Order MUST match GENOTYPES_COLUMNS in ingest.vcf_load.
GENOTYPES_ARROW_SCHEMA = pa.schema(
    [
        _f("ingest_id", pa.string(), False),
        _f("chrom", pa.string(), False),
        _f("pos", pa.uint32(), False),
        _f("ref", pa.string(), False),
        _f("alt", pa.string(), False),
        _f("sample_id", pa.string(), False),
        _f("gt", pa.int8(), False),
        _f("phased", pa.uint8(), False),
        _f("gq", pa.uint16()),
        _f("dp", pa.uint16()),
        _f("mq", pa.uint16()),
        _f("ft", pa.string()),
        _f("ps", pa.uint32()),
        _f("pq", pa.uint16()),
        _f("ad_ref", pa.uint16()),
        _f("ad_alt", pa.uint16()),
        _f("adf_ref", pa.uint16()),
        _f("adf_alt", pa.uint16()),
        _f("adr_ref", pa.uint16()),
        _f("adr_alt", pa.uint16()),
        _f("pl_ref_ref", pa.uint16()),
        _f("pl_ref_alt", pa.uint16()),
        _f("pl_alt_alt", pa.uint16()),
        _f("gl_ref_ref", pa.float32()),
        _f("gl_ref_alt", pa.float32()),
        _f("gl_alt_alt", pa.float32()),
        _f("format_extra", pa.map_(pa.string(), pa.string()), False),
    ]
)


def column_names(schema: pa.Schema) -> list[str]:
    return [f.name for f in schema]


SAMPLES_ARROW_SCHEMA = pa.schema(
    [
        _f("ingest_id", pa.string(), False),
        _f("sample_id", pa.string(), False),
        _f("cohort", pa.string(), False),
        _f("sex", pa.string()),
    ]
)

INGESTIONS_ARROW_SCHEMA = pa.schema(
    [
        _f("ingest_id", pa.string(), False),
        _f("cohort", pa.string(), False),
        _f("vcf_path", pa.string(), False),
        _f("n_variants", pa.uint64(), False),
        _f("n_samples", pa.uint32(), False),
    ]
)


VARIANTS_COLUMNS = column_names(VARIANTS_ARROW_SCHEMA)
GENOTYPES_COLUMNS = column_names(GENOTYPES_ARROW_SCHEMA)
SAMPLES_COLUMNS = column_names(SAMPLES_ARROW_SCHEMA)
INGESTIONS_COLUMNS = column_names(INGESTIONS_ARROW_SCHEMA)


def write_parquet(rows: Iterable[list], schema: pa.Schema, path: Path) -> int:
    """Convert row tuples (in `schema` column order) to a Parquet file.

    Returns the number of rows written. Empty `rows` writes a valid
    empty Parquet file with the schema, which the chDB importer
    silently no-ops on.
    """
    rows = list(rows)
    if not rows:
        # An empty Parquet is fine — chDB's file() function returns 0 rows.
        pq.write_table(
            pa.Table.from_arrays(
                [pa.array([], type=f.type) for f in schema],
                schema=schema,
            ),
            path,
        )
        return 0

    columns = list(zip(*rows))
    arrays = [
        pa.array(list(col), type=field.type) for col, field in zip(columns, schema)
    ]
    table = pa.Table.from_arrays(arrays, schema=schema)
    pq.write_table(table, path)
    return len(rows)
