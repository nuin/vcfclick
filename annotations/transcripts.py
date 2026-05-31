"""Transcript / exon / CDS / UTR annotation hierarchy.

Phase 2 work. Same license as the engine. No tier boundary.

Ships open alongside the rest. The stubs below capture the planned API
so the design is visible in the repo, and so anyone reading the code
can see what's coming without having to dig through a roadmap doc.

Why transcript-level matters:
  - "non-ref in BRCA1" includes deep intronic + UTR calls that are
    usually noise for downstream interpretation.
  - "non-ref in BRCA1 CDS, AF < 0.01" is the question a research
    bioinformatician actually wants.
  - Canonical transcript (MANE Select) avoids spurious hits on rare
    isoforms.
  - Splice-site distance is needed for any consequence prediction.

Data source: GTF/GFF (RefSeq GFF3 + Ensembl GTF, both GRCh38).
NCBI's GFF3 is public domain; Ensembl GTF is Apache 2 / open data.
The integration — curated, GRCh38-pinned, canonical-tagged DuckDB
artefact — is shipped under the same OSS license as the engine.
"""

from __future__ import annotations


# ─────────────────────────────────────────────────────────────────────
# Planned DuckDB schema. Lives in the same annotations.duckdb file as
# the gene + ClinVar tables.
# ─────────────────────────────────────────────────────────────────────

TRANSCRIPTS_SCHEMA_DDL = """
CREATE TABLE IF NOT EXISTS transcripts (
    transcript_id    VARCHAR PRIMARY KEY,    -- e.g. 'NM_007294.4'
    gene_symbol      VARCHAR NOT NULL,
    chrom            VARCHAR NOT NULL,
    start_pos        UINTEGER NOT NULL,
    end_pos          UINTEGER NOT NULL,
    strand           VARCHAR,
    biotype          VARCHAR,                -- 'protein_coding', 'lncRNA', ...
    is_canonical     BOOLEAN DEFAULT FALSE,  -- MANE Select tag
    is_mane_plus     BOOLEAN DEFAULT FALSE   -- MANE Plus Clinical
);

CREATE INDEX IF NOT EXISTS idx_transcripts_gene
    ON transcripts (gene_symbol);

CREATE TABLE IF NOT EXISTS exons (
    transcript_id    VARCHAR NOT NULL,
    exon_number      USMALLINT NOT NULL,     -- 1-indexed in transcription order
    chrom            VARCHAR NOT NULL,
    start_pos        UINTEGER NOT NULL,
    end_pos          UINTEGER NOT NULL,
    PRIMARY KEY (transcript_id, exon_number)
);

CREATE INDEX IF NOT EXISTS idx_exons_range
    ON exons (chrom, start_pos, end_pos);

CREATE TABLE IF NOT EXISTS cds (
    transcript_id    VARCHAR NOT NULL,
    exon_number      USMALLINT NOT NULL,
    chrom            VARCHAR NOT NULL,
    start_pos        UINTEGER NOT NULL,      -- coding-only, excludes UTR
    end_pos          UINTEGER NOT NULL,
    phase            USMALLINT,              -- 0, 1, or 2
    PRIMARY KEY (transcript_id, exon_number)
);

CREATE INDEX IF NOT EXISTS idx_cds_range
    ON cds (chrom, start_pos, end_pos);
"""


# ─────────────────────────────────────────────────────────────────────
# Planned public API. Phase 2 implementation.
# ─────────────────────────────────────────────────────────────────────

def transcripts_for_gene(symbol: str) -> list:
    """All transcripts of a gene."""
    raise NotImplementedError("Phase 2.")


def canonical_transcript(symbol: str):
    """The MANE Select transcript for a gene, if defined."""
    raise NotImplementedError("Phase 2.")


def cds_regions_for_gene(symbol: str) -> list:
    """Disjoint CDS ranges for a gene, suitable for a SQL range filter.
    The clinically meaningful version of position_for_gene()."""
    raise NotImplementedError("Phase 2.")


def exon_at(chrom: str, pos: int) -> list:
    """All (transcript_id, exon_number) pairs containing a position."""
    raise NotImplementedError("Phase 2.")


def splice_site_distance(chrom: str, pos: int) -> int | None:
    """Distance in bp to the nearest exon/intron boundary (signed)."""
    raise NotImplementedError("Phase 2.")
