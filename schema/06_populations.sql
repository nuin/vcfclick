-- Sample -> population panel, loaded with `vcfclick db panel <name>
-- <file>` (1000 Genomes panel format or a generic TSV/CSV). Like
-- `pedigree`, it is loaded separately from VCF ingest and is NOT wiped
-- when a VCF is re-ingested under the same ingest_id.
--
-- Sample identity is (ingest_id, sample_id), matching `samples`.
-- Re-loading a panel replaces the rows for the affected samples.

CREATE TABLE populations (
    ingest_id        LowCardinality(String),
    sample_id        LowCardinality(String),
    population       LowCardinality(String),
    super_population LowCardinality(Nullable(String)),
    sex              LowCardinality(Nullable(String)),  -- 'male' / 'female' / NULL

    ingested_at      DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(ingested_at)
ORDER BY (ingest_id, sample_id);
