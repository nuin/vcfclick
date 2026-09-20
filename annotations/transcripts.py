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

from dataclasses import dataclass


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


@dataclass(frozen=True)
class Transcript:
    """One transcript row."""

    transcript_id: str
    gene_symbol: str
    chrom: str
    start_pos: int
    end_pos: int
    strand: str | None
    biotype: str | None
    is_canonical: bool
    is_mane_plus: bool


_TX_COLS = (
    "transcript_id, gene_symbol, chrom, start_pos, end_pos, strand, biotype, "
    "is_canonical, is_mane_plus"
)


def transcripts_for_gene(symbol: str) -> list[Transcript]:
    """All transcripts of a gene, canonical first."""
    from annotations.db import get_connection

    rows = (
        get_connection()
        .execute(
            f"SELECT {_TX_COLS} FROM transcripts WHERE gene_symbol = ? "
            "ORDER BY is_canonical DESC, transcript_id",
            [symbol],
        )
        .fetchall()
    )
    return [Transcript(*r) for r in rows]


def canonical_transcript(symbol: str) -> Transcript | None:
    """The MANE Select transcript for a gene, if defined.

    Using it avoids spurious hits on rare isoforms.
    """
    from annotations.db import get_connection

    row = (
        get_connection()
        .execute(
            f"SELECT {_TX_COLS} FROM transcripts "
            "WHERE gene_symbol = ? AND is_canonical ORDER BY transcript_id LIMIT 1",
            [symbol],
        )
        .fetchone()
    )
    return Transcript(*row) if row else None


def cds_regions_for_gene(symbol: str) -> list[tuple[str, int, int]]:
    """Disjoint CDS ranges for a gene, suitable for a SQL range filter.

    The clinically meaningful version of `position_for_gene()`: excludes
    introns and UTRs, so "non-ref in BRCA1 CDS" is expressible. Ranges from
    the canonical transcript when one is defined, else every transcript;
    overlapping ranges are merged.
    """
    from annotations.db import get_connection

    conn = get_connection()
    canon = canonical_transcript(symbol)
    if canon is not None:
        rows = conn.execute(
            "SELECT chrom, start_pos, end_pos FROM cds WHERE transcript_id = ? "
            "ORDER BY chrom, start_pos",
            [canon.transcript_id],
        ).fetchall()
    else:
        rows = conn.execute(
            "SELECT c.chrom, c.start_pos, c.end_pos FROM cds c "
            "JOIN transcripts t ON c.transcript_id = t.transcript_id "
            "WHERE t.gene_symbol = ? ORDER BY c.chrom, c.start_pos",
            [symbol],
        ).fetchall()

    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in rows:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2] + 1:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, int(start), int(end)))
    return merged


def exon_at(chrom: str, pos: int) -> list[tuple[str, int]]:
    """All (transcript_id, exon_number) pairs whose exon contains `pos`."""
    from annotations.db import get_connection

    rows = (
        get_connection()
        .execute(
            "SELECT transcript_id, exon_number FROM exons "
            "WHERE chrom = ? AND ? BETWEEN start_pos AND end_pos "
            "ORDER BY transcript_id, exon_number",
            [chrom, pos],
        )
        .fetchall()
    )
    return [(t, int(n)) for t, n in rows]


def splice_site_distance(chrom: str, pos: int) -> int | None:
    """Distance in bp to the nearest exon/intron boundary, or None if no exon
    is near. 0 means the position sits on a boundary.

    Needed for any consequence prediction: a variant a couple of bases from a
    boundary is a candidate splice variant even when it looks intronic.
    """
    from annotations.db import get_connection

    row = (
        get_connection()
        .execute(
            "SELECT min(least(abs(? - start_pos), abs(? - end_pos))) FROM exons "
            "WHERE chrom = ?",
            [pos, pos, chrom],
        )
        .fetchone()
    )
    return None if row is None or row[0] is None else int(row[0])
