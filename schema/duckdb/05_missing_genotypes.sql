-- Missing genotype calls — DuckDB equivalent of
-- schema/05_missing_genotypes.sql. Same columns and order.

CREATE TABLE missing_genotypes (
    ingest_id    VARCHAR NOT NULL,
    chrom        VARCHAR NOT NULL,
    pos          UINTEGER NOT NULL,
    ref          VARCHAR NOT NULL,
    alt          VARCHAR NOT NULL,
    sample_id    VARCHAR NOT NULL,
    ingested_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
