-- Population panel — DuckDB equivalent of schema/06_populations.sql.
-- Same columns and order.

CREATE TABLE populations (
    ingest_id        VARCHAR NOT NULL,
    sample_id        VARCHAR NOT NULL,
    population       VARCHAR NOT NULL,
    super_population VARCHAR,
    sex              VARCHAR,
    ingested_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
