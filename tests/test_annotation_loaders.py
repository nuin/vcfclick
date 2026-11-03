"""Annotation loaders replace a table only once the new file has parsed.

A loader parses its input first and swaps the rows in inside one
transaction, so a truncated download or a file with no usable rows raises
and leaves the store as it was. Small synthetic inputs; no network.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

GENE_GFF = (
    "##gff-version 3\n"
    "chr1\tHAVANA\tgene\t100\t400\t.\t+\t.\tID=g1;gene_id=G1;gene_name=GENEX;"
    "gene_type=protein_coding\n"
    "chr1\tHAVANA\ttranscript\t100\t400\t.\t+\t.\tID=T1;Parent=g1;transcript_id=T1;"
    "gene_name=GENEX;transcript_type=protein_coding;tag=basic,MANE_Select\n"
    "chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tParent=T1;transcript_id=T1;"
    "gene_name=GENEX;exon_number=1\n"
)
NO_FEATURES_GFF = "##gff-version 3\nchr1\tHAVANA\tregion\t1\t10\t.\t+\t.\tID=r\n"


def _gff(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    with gzip.open(p, "wt") as fh:
        fh.write(text)
    return p


def _truncated(tmp_path: Path, src: Path) -> Path:
    data = src.read_bytes()
    p = tmp_path / f"cut_{src.name}"
    p.write_bytes(data[: len(data) // 2])
    return p


def _count(table: str) -> int:
    from annotations.db import read_connection

    conn = read_connection()
    try:
        return conn.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
    finally:
        conn.close()


def _big_gff(tmp_path: Path) -> Path:
    lines = ["##gff-version 3"]
    for i in range(3000):
        lines.append(
            f"chr1\tHAVANA\tgene\t{i * 10 + 1}\t{i * 10 + 5}\t.\t+\t.\t"
            f"ID=g{i};gene_id=G{i};gene_name=GENE{i};gene_type=protein_coding"
        )
    return _gff(tmp_path, "big.gff3.gz", "\n".join(lines) + "\n")


@pytest.mark.parametrize("bad", ["empty", "truncated"])
def test_gene_load_keeps_the_table_when_the_file_fails(
    tmp_path, isolated_annotation_db, bad
):
    from annotations.loaders.gencode_genes import load

    good = _big_gff(tmp_path)
    assert load(good) == 3000
    path = (
        _gff(tmp_path, "none.gff3.gz", NO_FEATURES_GFF)
        if bad == "empty"
        else _truncated(tmp_path, good)
    )
    with pytest.raises((RuntimeError, EOFError, OSError)):
        load(path)
    assert _count("refseq_genes") == 3000


@pytest.mark.parametrize("bad", ["empty", "truncated"])
def test_transcript_load_keeps_the_tables_when_the_file_fails(
    tmp_path, isolated_annotation_db, bad
):
    from annotations.loaders.gencode_transcripts import load

    good = _gff(tmp_path, "good.gff3.gz", GENE_GFF * 1)
    assert load(good)["transcripts"] == 1
    if bad == "empty":
        path = _gff(tmp_path, "none.gff3.gz", NO_FEATURES_GFF)
    else:
        path = _truncated(tmp_path, _big_gff(tmp_path))
    with pytest.raises((RuntimeError, EOFError, OSError)):
        load(path)
    assert _count("transcripts") == 1
    assert _count("exons") == 1


def test_a_failed_insert_rolls_the_delete_back(
    tmp_path, isolated_annotation_db, monkeypatch
):
    """The DELETE and the INSERT are one transaction: an INSERT that fails
    after the DELETE leaves the old rows."""
    import annotations.loaders.gencode_genes as gg

    assert gg.load(_big_gff(tmp_path)) == 3000
    # two rows with one primary key: the INSERT fails after the DELETE ran
    dup = _gff(
        tmp_path,
        "dup.gff3.gz",
        "chr1\tH\tgene\t1\t5\t.\t+\t.\tgene_id=A;gene_name=DUP\n"
        "chr1\tH\tgene\t9\t15\t.\t+\t.\tgene_id=B;gene_name=OTHER\n",
    )
    real = gg.iter_genes

    def twice(path):
        rows = list(real(path))
        return rows + rows[:1]

    monkeypatch.setattr(gg, "iter_genes", twice)
    with pytest.raises(Exception):  # noqa: B017 - the constraint error
        gg.load(dup)
    assert _count("refseq_genes") == 3000


def test_hgnc_load_keeps_the_table_when_the_file_fails(
    tmp_path, isolated_annotation_db
):
    from annotations.loaders.hgnc_aliases import load

    fixture = Path(__file__).parent / "fixtures" / "hgnc_aliases.tsv"
    load(fixture)
    before = _count("gene_aliases")
    assert before > 0
    empty = tmp_path / "empty.tsv"
    empty.write_text(fixture.read_text().splitlines()[0] + "\n")
    with pytest.raises(ValueError):
        load(empty)
    assert _count("gene_aliases") == before


def test_loaders_release_the_store(tmp_path, isolated_annotation_db, run_python):
    """After a load returns, another process can read the store: the loader
    closed its writer."""
    from annotations.loaders.gencode_genes import load

    load(_big_gff(tmp_path))
    out = run_python(
        tmp_path,
        "import annotations.db as a; print(a.position_for_gene('GENE7').start_pos)",
    )
    assert out.strip() == "71"


PAR_GFF = (
    "##gff-version 3\n"
    "chrX\tHAVANA\tgene\t624344\t659411\t.\t+\t.\tgene_id=ENSG00000185960.14;"
    "gene_name=SHOX;gene_type=protein_coding\n"
    "chrY\tHAVANA\tgene\t624344\t659411\t.\t+\t.\tgene_id=ENSG00000185960.14_PAR_Y;"
    "gene_name=SHOX;gene_type=protein_coding\n"
    "chr1\tHAVANA\tgene\t100\t200\t.\t+\t.\tgene_id=ENSG1.1;gene_name=TWICE;"
    "gene_type=lncRNA\n"
    "chr1\tHAVANA\tgene\t500\t600\t.\t+\t.\tgene_id=ENSG2.1;gene_name=TWICE;"
    "gene_type=lncRNA\n"
)


def test_a_par_gene_keeps_both_loci(tmp_path, isolated_annotation_db, caplog):
    """SHOX is listed on chrX and on chrY (its `_PAR_Y` copy): both loci are
    stored. A symbol repeated on one contig keeps its first locus, and the
    loader says how many it dropped."""
    from annotations.db import gene_at
    from annotations.loaders.gencode_genes import load

    with caplog.at_level("INFO", logger="annotations.loaders.gencode_genes"):
        assert load(_gff(tmp_path, "par.gff3.gz", PAR_GFF)) == 3
    assert [g.chrom for g in gene_at("chrY", 630_000)] == ["chrY"]
    assert [g.chrom for g in gene_at("chrX", 630_000)] == ["chrX"]
    assert [g.start_pos for g in gene_at("chr1", 150)] == [100]
    assert gene_at("chr1", 550) == []
    assert "1 gene row" in caplog.text and "TWICE" in caplog.text


def test_a_store_keyed_by_symbol_alone_is_rebuilt(tmp_path, isolated_annotation_db):
    """A store made before the key included the contig: a replacing load
    rebuilds the table, so the chrY copy fits."""
    import duckdb

    from annotations.loaders.gencode_genes import load

    conn = duckdb.connect(str(isolated_annotation_db))
    conn.execute(
        "CREATE TABLE refseq_genes (gene_symbol VARCHAR PRIMARY KEY, chrom VARCHAR "
        "NOT NULL, start_pos UINTEGER NOT NULL, end_pos UINTEGER NOT NULL, strand "
        "VARCHAR, refseq_id VARCHAR, description VARCHAR)"
    )
    conn.execute("INSERT INTO refseq_genes VALUES ('OLD','chr2',1,5,'+','E','pc')")
    conn.close()
    assert load(_gff(tmp_path, "par.gff3.gz", PAR_GFF)) == 3
    assert _count("refseq_genes") == 3


def test_genes_scope_takes_every_locus_of_a_symbol(tmp_path, isolated_annotation_db):
    from annotations.db import resolve_gene_tokens
    from annotations.loaders.gencode_genes import load
    from cli import scope as sc

    load(_gff(tmp_path, "par.gff3.gz", PAR_GFF))
    found = resolve_gene_tokens(["SHOX", "ENSG00000185960"])
    assert sorted(g.chrom for _, g in found["SHOX"]) == ["chrX", "chrY"]
    assert sorted(g.chrom for _, g in found["ENSG00000185960"]) == ["chrX", "chrY"]
    scope = sc.build_scope(
        sc.ScopeRequest.from_values({"genes": ["SHOX"]}), build="GRCh38"
    )
    assert scope.contains("chrX", 630_000) and scope.contains("chrY", 630_000)
    assert scope.echo["genes"]["genes"] == 1
    assert not any("chrX only" in n for n in scope.notes)
