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

-- NOTE on cohort-size denominators: there is intentionally no
-- materialized cohort_sizes view. An earlier version had a
-- `cohort_sizes_mv` SummingMergeTree, but SummingMergeTree never
-- decrements on DELETE — rolling back a samples insert (failed
-- ingest, --ingest-id replacement) silently left the view's count
-- inflated, breaking any AF query that used it as the denominator.
-- Compute cohort size directly when you need it:
--     SELECT cohort, count(DISTINCT (ingest_id, sample_id)) AS n
--     FROM samples WHERE cohort = ... GROUP BY cohort
-- At realistic cohort scale (~10^4 samples) this is microseconds.

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
