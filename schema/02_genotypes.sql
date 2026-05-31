-- Sparse non-reference genotype calls. Scoped by ingest_id.
--
-- Sample identity is (ingest_id, sample_id). A sample named "001" in
-- ingestion A is NOT the same as "001" in ingestion B — they're whatever
-- the uploader said they were, and we don't try to reconcile. Cross-
-- ingestion matching (e.g., the same patient sequenced twice) is a
-- higher-level concern that a user-defined alias table would handle.
--
-- Re-ingestion with the same ingest_id is idempotent: ReplacingMergeTree
-- on the full sorting key (ingest_id, chrom, pos, ref, alt, sample_id)
-- dedupes silently. Re-ingestion with a new ingest_id appends.

CREATE TABLE genotypes (
    ingest_id    LowCardinality(String),

    chrom        LowCardinality(String),
    pos          UInt32,
    ref          String,
    alt          String,
    sample_id    LowCardinality(String),

    -- Genotype encoding: alt-allele count (Int8).
    --   1  = heterozygous (0/1, 1/0)
    --   2  = homozygous alt (1/1)
    --   -1 = hemizygous or mixed-missing
    -- 0/0 and ./. are NOT stored — absence is the signal.
    gt           Int8,
    phased       UInt8 DEFAULT 0,

    -- Reserved FORMAT fields.
    gq           Nullable(UInt16),
    dp           Nullable(UInt16),
    ad_ref       Nullable(UInt16),
    ad_alt       Nullable(UInt16),
    adf_ref      Nullable(UInt16),
    adf_alt      Nullable(UInt16),
    adr_ref      Nullable(UInt16),
    adr_alt      Nullable(UInt16),
    pl_ref_ref   Nullable(UInt16),
    pl_ref_alt   Nullable(UInt16),
    pl_alt_alt   Nullable(UInt16),
    gl_ref_ref   Nullable(Float32),
    gl_ref_alt   Nullable(Float32),
    gl_alt_alt   Nullable(Float32),
    mq           Nullable(UInt16),
    ft           LowCardinality(Nullable(String)),   -- per-sample FILTER; NULL when unspecified
    ps           Nullable(UInt32),
    pq           Nullable(UInt16),

    -- Overflow for non-reserved FORMAT fields.
    format_extra Map(String, String),

    ingested_at  DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(ingested_at)
ORDER BY (ingest_id, chrom, pos, ref, alt, sample_id)
SETTINGS deduplicate_merge_projection_mode = 'rebuild';
-- Required as of ClickHouse 24.x to allow ADD PROJECTION on a
-- ReplacingMergeTree (same rationale as variants).

-- Per-sample queries within an ingestion. Sample-major sort order.
ALTER TABLE genotypes ADD INDEX idx_sample sample_id TYPE bloom_filter GRANULARITY 4;

ALTER TABLE genotypes ADD PROJECTION p_sample_major (
    SELECT *
    ORDER BY (ingest_id, sample_id, chrom, pos, ref, alt)
);
