-- Pedigree / family relationships — DuckDB equivalent of
-- schema/04_pedigree.sql. Same columns and order; LowCardinality(X)
-- becomes VARCHAR, Nullable(X) collapses to X (nullable by default).

CREATE TABLE pedigree (
    ingest_id    VARCHAR NOT NULL,
    sample_id    VARCHAR NOT NULL,
    family_id    VARCHAR NOT NULL,
    father_id    VARCHAR NOT NULL,
    mother_id    VARCHAR NOT NULL,
    sex          VARCHAR,
    affected     VARCHAR,
    ingested_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
