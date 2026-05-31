-- Variants table: one row per unique (ingest_id, chrom, pos, ref, alt).
--
-- Crucially: rows are NOT merged across ingestions. The same coordinate
-- + ref + alt observed in two different VCFs is two rows, because their
-- QUAL, FILTER, INFO annotations, and QC origin can all differ — merging
-- them would silently lose information. Cross-ingestion correlation is a
-- query-time concern, not a storage-time one.
--
-- Dedup on same-ID re-upload is automatic via ReplacingMergeTree on
-- (ingest_id, chrom, pos, ref, alt). Re-running an ingestion with the
-- same --ingest-id silently replaces rows; using a new --ingest-id
-- adds new rows alongside the old.

CREATE TABLE variants (
    ingest_id    LowCardinality(String),           -- one per VCF upload

    chrom        LowCardinality(String),
    pos          UInt32,
    ref          String,
    alt          String,

    vcf_id       Nullable(String),
    qual         Nullable(Float32),
    filter       LowCardinality(Nullable(String)),   -- NULL = unspecified ('.'), distinct from 'PASS'

    -- VCF 4.3 reserved INFO. Per-ALT fields are scalars (decomposed).
    info_AC      Nullable(UInt32),
    info_AF      Nullable(Float32),
    info_AN      Nullable(UInt32),
    info_AD_ref  Nullable(UInt32),
    info_AD_alt  Nullable(UInt32),
    info_DP      Nullable(UInt32),
    info_MQ      Nullable(Float32),
    info_MQ0     Nullable(UInt32),
    info_NS      Nullable(UInt32),
    info_BQ      Nullable(Float32),
    info_SB      Nullable(Float32),
    info_END     Nullable(UInt32),
    info_CIGAR   Nullable(String),
    info_AA      Nullable(String),

    info_SOMATIC   UInt8 DEFAULT 0,
    info_VALIDATED UInt8 DEFAULT 0,
    info_DB        UInt8 DEFAULT 0,
    info_H2        UInt8 DEFAULT 0,
    info_H3        UInt8 DEFAULT 0,
    info_1000G     UInt8 DEFAULT 0,

    -- Common GATK metrics.
    info_QD               Nullable(Float32),
    info_FS               Nullable(Float32),
    info_SOR              Nullable(Float32),
    info_MQRankSum        Nullable(Float32),
    info_ReadPosRankSum   Nullable(Float32),
    info_ExcessHet        Nullable(Float32),
    info_InbreedingCoeff  Nullable(Float32),
    info_MLEAC            Nullable(UInt32),
    info_MLEAF            Nullable(Float32),
    info_BaseQRankSum     Nullable(Float32),
    info_ClippingRankSum  Nullable(Float32),

    -- Overflow for non-reserved INFO fields.
    info_extra   Map(String, String),

    ingested_at  DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(ingested_at)
ORDER BY (ingest_id, chrom, pos, ref, alt)
SETTINGS deduplicate_merge_projection_mode = 'rebuild';
-- ^ Required as of ClickHouse 24.x to allow ADD PROJECTION on a
-- ReplacingMergeTree. 'rebuild' regenerates the projection during merge
-- dedup so projection reads stay consistent with the main table.

-- Cross-ingestion range scan: "all variants at this gene region across
-- all uploads." Drops ingest_id from the sort to let scans coalesce.
ALTER TABLE variants ADD PROJECTION p_range_scan (
    SELECT chrom, pos, ref, alt, ingest_id, filter, info_AF, info_DP
    ORDER BY (chrom, pos)
);
