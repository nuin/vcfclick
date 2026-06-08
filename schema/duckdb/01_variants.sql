-- Variants table — DuckDB schema, byte-for-byte equivalent to the chDB
-- schema/01_variants.sql in terms of which columns exist and what their
-- Arrow round-trip types are. The differences are all engine-level:
--
--   * No `ReplacingMergeTree` — DuckDB has no merge-tree engine. The
--     same idempotent-replace semantics on `--ingest-id` are enforced
--     in Python by `rollback_ingest()` (DELETE WHERE ingest_id = X
--     followed by INSERT), independent of backend.
--
--   * No `ORDER BY (...)` — DuckDB has no clustering/storage-order
--     concept of the kind ClickHouse uses. Range scans run on
--     min/max statistics over column row-groups instead.
--
--   * `LowCardinality(String)` becomes `VARCHAR` — DuckDB applies its
--     own dictionary encoding internally. The shape of the data on
--     disk is similar; the SQL surface is plain VARCHAR.
--
--   * `Nullable(X)` collapses to `X` — DuckDB columns are nullable by
--     default unless `NOT NULL` is declared. Columns that were
--     `LowCardinality(String)` (NOT NULL in chDB) get an explicit
--     `NOT NULL` here to match.
--
--   * No projections — DuckDB does not yet support materialised
--     projection equivalents. Range scans over (chrom, pos) rely on
--     min/max row-group pruning, which has been good enough at
--     vcfclick's per-cohort scale on the bench.

CREATE TABLE variants (
    ingest_id    VARCHAR NOT NULL,

    chrom        VARCHAR NOT NULL,
    pos          UINTEGER NOT NULL,
    ref          VARCHAR NOT NULL,
    alt          VARCHAR NOT NULL,

    vcf_id       VARCHAR,
    qual         REAL,
    filter       VARCHAR,                       -- NULL = unspecified ('.'), distinct from 'PASS'

    -- VCF 4.3 reserved INFO. Per-ALT fields are scalars (decomposed).
    info_AC      UINTEGER,
    info_AF      REAL,
    info_AN      UINTEGER,
    info_AD_ref  UINTEGER,
    info_AD_alt  UINTEGER,
    info_DP      UINTEGER,
    info_MQ      REAL,
    info_MQ0     UINTEGER,
    info_NS      UINTEGER,
    info_BQ      REAL,
    info_SB      REAL,
    info_END     UINTEGER,
    info_CIGAR   VARCHAR,
    info_AA      VARCHAR,

    info_SOMATIC   UTINYINT DEFAULT 0,
    info_VALIDATED UTINYINT DEFAULT 0,
    info_DB        UTINYINT DEFAULT 0,
    info_H2        UTINYINT DEFAULT 0,
    info_H3        UTINYINT DEFAULT 0,
    info_1000G     UTINYINT DEFAULT 0,

    -- Common GATK metrics.
    info_QD               REAL,
    info_FS               REAL,
    info_SOR              REAL,
    info_MQRankSum        REAL,
    info_ReadPosRankSum   REAL,
    info_ExcessHet        REAL,
    info_InbreedingCoeff  REAL,
    info_MLEAC            UINTEGER,
    info_MLEAF            REAL,
    info_BaseQRankSum     REAL,
    info_ClippingRankSum  REAL,

    -- DRAGEN-specific.
    info_FractionInformativeReads REAL,
    info_HAPCOMP                  UINTEGER,
    info_HAPDOM                   REAL,
    info_DragenSnvHardQUAL        REAL,
    info_DragenIndelHardQUAL      REAL,

    -- Overflow for non-reserved INFO fields.
    info_extra   MAP(VARCHAR, VARCHAR),

    ingested_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
