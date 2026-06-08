-- Samples + ingestions — DuckDB equivalent of schema/03_samples.sql.
-- The cohort-size-MV rationale from the chDB schema still applies
-- (see vcfclick_mcp/server.py SCHEMA_DESCRIPTION and docs/SCHEMA.md);
-- there's no materialised view here either.

CREATE TABLE samples (
    ingest_id    VARCHAR NOT NULL,
    sample_id    VARCHAR NOT NULL,
    cohort       VARCHAR NOT NULL,
    sex          VARCHAR,
    ingested_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE ingestions (
    ingest_id    VARCHAR NOT NULL,
    cohort       VARCHAR NOT NULL,
    vcf_path     VARCHAR NOT NULL,
    n_variants   UBIGINT DEFAULT 0,
    n_samples    UINTEGER DEFAULT 0,
    ingested_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
