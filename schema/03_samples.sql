-- Samples table. Sample identity is (ingest_id, sample_id).
-- Cohort is an orthogonal grouping that can span ingestions.

CREATE TABLE samples (
    ingest_id    LowCardinality(String),
    sample_id    LowCardinality(String),
    cohort       LowCardinality(String),
    sex          LowCardinality(Nullable(String)),
    ingested_at  DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(ingested_at)
ORDER BY (ingest_id, sample_id);

-- Cohort sizes view: cohort → n_samples. Counts distinct (ingest_id,
-- sample_id) pairs so the same cohort name across two ingestions sums
-- correctly. Use this as the denominator for cohort allele frequencies.
CREATE MATERIALIZED VIEW cohort_sizes_mv
ENGINE = SummingMergeTree()
ORDER BY cohort
AS SELECT cohort, count() AS n_samples
   FROM samples
   GROUP BY cohort;


-- Ingestion catalog. One row per VCF upload. Useful for management
-- queries ("what's loaded?") and for the MCP server's get_ingestions tool.
CREATE TABLE ingestions (
    ingest_id    String,
    cohort       LowCardinality(String),
    vcf_path     String,
    n_variants   UInt64 DEFAULT 0,
    n_samples    UInt32 DEFAULT 0,
    ingested_at  DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(ingested_at)
ORDER BY ingest_id;
