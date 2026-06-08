-- Genotypes table — DuckDB equivalent of schema/02_genotypes.sql.
-- Same sparse encoding: 0/0 and ./. are NOT stored; absence is the
-- signal. The bloom-filter index and sample-major projection from the
-- chDB schema don't have direct DuckDB equivalents; per-sample queries
-- rely on DuckDB's row-group statistics for pruning at this scale.

CREATE TABLE genotypes (
    ingest_id    VARCHAR NOT NULL,

    chrom        VARCHAR NOT NULL,
    pos          UINTEGER NOT NULL,
    ref          VARCHAR NOT NULL,
    alt          VARCHAR NOT NULL,
    sample_id    VARCHAR NOT NULL,

    -- Genotype encoding: alt-allele count (TINYINT, signed).
    --   1  = heterozygous (0/1, 1/0)
    --   2  = homozygous alt (1/1)
    --   -1 = hemizygous or mixed-missing
    -- 0/0 and ./. are NOT stored — absence is the signal.
    gt           TINYINT NOT NULL,
    phased       UTINYINT DEFAULT 0,

    -- Reserved FORMAT fields.
    gq           USMALLINT,
    dp           USMALLINT,
    ad_ref       USMALLINT,
    ad_alt       USMALLINT,
    adf_ref      USMALLINT,
    adf_alt      USMALLINT,
    adr_ref      USMALLINT,
    adr_alt      USMALLINT,
    pl_ref_ref   USMALLINT,
    pl_ref_alt   USMALLINT,
    pl_alt_alt   USMALLINT,
    gl_ref_ref   REAL,
    gl_ref_alt   REAL,
    gl_alt_alt   REAL,
    mq           USMALLINT,
    ft           VARCHAR,
    ps           UINTEGER,
    pq           USMALLINT,

    -- Overflow for non-reserved FORMAT fields.
    format_extra MAP(VARCHAR, VARCHAR),

    ingested_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
