-- Per-sample fully missing genotype calls (./. or .), one row per
-- (ingest_id, chrom, pos, ref, alt, sample_id). Written at ingest by
-- default; `db ingest --no-record-missing` skips it.
--
-- `genotypes` keeps its meaning (a stored row is a called non-reference
-- genotype, or 0/0 under --keep-reference), so a sample absent from
-- `genotypes` is 0/0 unless it has a row here. That lets per-group
-- called counts be derived without storing every 0/0:
--     called_in_group = group_size - missing_in_group
--
-- Partially missing calls (./1, 0/.) are NOT rows here: they are not
-- fully missing. The per-site variants.an_called / ac_called counts stay
-- exact for them.

CREATE TABLE missing_genotypes (
    ingest_id    LowCardinality(String),

    chrom        LowCardinality(String),
    pos          UInt32,
    ref          String,
    alt          String,
    sample_id    LowCardinality(String),

    ingested_at  DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(ingested_at)
ORDER BY (ingest_id, chrom, pos, ref, alt, sample_id);
