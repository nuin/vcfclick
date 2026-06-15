-- Pedigree / family relationships, loaded from a PED/FAM file via
-- `vcfclick db ped <name> <file.ped>`. Separate from VCF ingest: a PED
-- describes who-is-whose-parent among already-ingested samples, so it's
-- loaded independently and is NOT wiped when a VCF is re-ingested under
-- the same ingest_id.
--
-- Sample identity is (ingest_id, sample_id), matching `samples`.
-- father_id / mother_id are the sample_id of the parent within the same
-- ingest_id, or '0' for a founder (PED convention). v1 supports trios
-- from a JOINT-called VCF, so all three members share one ingest_id.

CREATE TABLE pedigree (
    ingest_id    LowCardinality(String),
    sample_id    LowCardinality(String),
    family_id    LowCardinality(String),
    father_id    LowCardinality(String),           -- sample_id or '0'
    mother_id    LowCardinality(String),           -- sample_id or '0'
    sex          LowCardinality(Nullable(String)), -- 'male' / 'female' / NULL
    affected     LowCardinality(Nullable(String)), -- 'affected' / 'unaffected' / NULL

    ingested_at  DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(ingested_at)
ORDER BY (ingest_id, sample_id);
