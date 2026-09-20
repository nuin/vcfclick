"""Transcript / exon / CDS annotations (the Phase-2 hierarchy).

Why it matters: "non-ref in BRCA1" spans deep intronic and UTR calls; the
question a bioinformatician actually asks is "non-ref in BRCA1 *CDS*". These
tests use a small synthetic GENCODE-shaped GFF3 so they need no download.

Synthetic gene GENEX on chr1, + strand, two transcripts:
  T1 (MANE Select, protein_coding): exons 100-200, 300-400; CDS 150-200, 300-350
  T2 (alt isoform):                 exons 100-250
"""

from __future__ import annotations

import gzip

import pytest

GFF = """\
##gff-version 3
chr1\tHAVANA\tgene\t100\t400\t.\t+\t.\tID=g1;gene_id=G1;gene_name=GENEX;gene_type=protein_coding
chr1\tHAVANA\ttranscript\t100\t400\t.\t+\t.\tID=T1;Parent=g1;transcript_id=T1;gene_name=GENEX;transcript_type=protein_coding;tag=basic,Ensembl_canonical,MANE_Select
chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tParent=T1;transcript_id=T1;gene_name=GENEX;exon_number=1
chr1\tHAVANA\texon\t300\t400\t.\t+\t.\tParent=T1;transcript_id=T1;gene_name=GENEX;exon_number=2
chr1\tHAVANA\tCDS\t150\t200\t.\t+\t0\tParent=T1;transcript_id=T1;gene_name=GENEX;exon_number=1
chr1\tHAVANA\tCDS\t300\t350\t.\t+\t2\tParent=T1;transcript_id=T1;gene_name=GENEX;exon_number=2
chr1\tHAVANA\ttranscript\t100\t250\t.\t+\t.\tID=T2;Parent=g1;transcript_id=T2;gene_name=GENEX;transcript_type=protein_coding;tag=basic
chr1\tHAVANA\texon\t100\t250\t.\t+\t.\tParent=T2;transcript_id=T2;gene_name=GENEX;exon_number=1
"""


@pytest.fixture
def loaded(tmp_path, monkeypatch):
    """Point the annotation store at a temp DuckDB and load the fixture GFF3."""
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    gff = tmp_path / "t.gff3.gz"
    with gzip.open(gff, "wt") as fh:
        fh.write(GFF)
    from annotations.loaders.gencode_transcripts import load

    counts = load(gff, replace=True)
    return counts


def test_load_reports_feature_counts(loaded):
    assert loaded["transcripts"] == 2
    assert loaded["exons"] == 3
    assert loaded["cds"] == 2


def test_transcripts_for_gene(loaded):
    from annotations.transcripts import transcripts_for_gene

    rows = transcripts_for_gene("GENEX")
    assert {r.transcript_id for r in rows} == {"T1", "T2"}
    t1 = next(r for r in rows if r.transcript_id == "T1")
    assert t1.is_canonical is True and t1.biotype == "protein_coding"


def test_canonical_transcript_is_mane_select(loaded):
    from annotations.transcripts import canonical_transcript

    t = canonical_transcript("GENEX")
    assert t is not None and t.transcript_id == "T1"


def test_canonical_transcript_none_when_absent(loaded):
    from annotations.transcripts import canonical_transcript

    assert canonical_transcript("NOSUCHGENE") is None


def test_cds_regions_for_gene_excludes_utr_and_introns(loaded):
    """The clinically meaningful version of position_for_gene(): CDS only,
    so the 100-149 UTR and the 201-299 intron are excluded."""
    from annotations.transcripts import cds_regions_for_gene

    regions = cds_regions_for_gene("GENEX")
    assert regions == [("chr1", 150, 200), ("chr1", 300, 350)]


def test_exon_at_returns_transcript_and_number(loaded):
    from annotations.transcripts import exon_at

    hits = set(exon_at("chr1", 150))  # inside T1 exon 1 and T2 exon 1
    assert ("T1", 1) in hits and ("T2", 1) in hits
    assert exon_at("chr1", 260) == []  # intron of T1, past T2


def test_splice_site_distance(loaded):
    from annotations.transcripts import splice_site_distance

    # Exon 1 of T1 spans 100-200: position 198 is 2bp from the 200 boundary.
    assert splice_site_distance("chr1", 198) == 2
    # Right on a boundary -> 0.
    assert splice_site_distance("chr1", 200) == 0
    # No transcript anywhere near -> None.
    assert splice_site_distance("chr9", 12345) is None
