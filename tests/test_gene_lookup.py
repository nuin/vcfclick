"""Gene lookups that ignore how a symbol is cased and how a contig is spelt.

Two bugs in shared code, fixed for every caller (db trio, db gene, the MCP
tools, the TUI):

  * `position_for_gene` upper-cased the symbol it was given and matched the
    stored one exactly, but the GENCODE loader stores a symbol as HGNC spells
    it: `C9orf72`, `C19orf12` and `C12orf57` never resolved.
  * `GeneIndex.at` (and `gene_at`) needed the store's contig name exactly
    (`chr7`), so a VCF with bare contigs (`7`) got no genes: no gene column in
    `db trio` and no compound het at all.

The annotation store is DuckDB and never chDB, so the lookups run in process
against an isolated store; the `db trio` check runs the CLI.
"""

from __future__ import annotations

import json

import pytest

from storage.contigs import canonical_chrom, chrom_aliases
from tests.conftest import run_cli
from tests.fixtures._vcfgen import TRIO_PED, TRIO_SAMPLES, Site, write_vcf

GENES = [
    ("C9orf72", "chr9", 27_546_545, 27_573_866, "-", "ENSG00000147894.16", "pc"),
    ("C19orf12", "chr19", 29_698_937, 29_715_789, "-", "ENSG00000131943.19", "pc"),
    ("C12orf57", "chr12", 6_942_978, 6_946_003, "+", "ENSG00000111678.11", "pc"),
    ("BRCA1", "chr17", 43_044_295, 43_125_483, "-", "ENSG00000012048.24", "pc"),
    ("G7", "chr7", 1_000, 2_000, "+", "ENSG00000000007.1", "pc"),
    ("XG", "chrX", 2_800_000, 2_900_000, "+", "ENSG00000124343.14", "pc"),
    ("MT-ND1", "chrM", 3_307, 4_262, "+", "ENSG00000198888.2", "pc"),
]


def _seed(genes=GENES) -> None:
    from annotations.db import get_connection

    conn = get_connection()
    conn.executemany("INSERT INTO refseq_genes VALUES (?, ?, ?, ?, ?, ?, ?)", genes)
    conn.close()


@pytest.mark.parametrize(
    "name, canonical",
    [
        ("7", "chr7"),
        ("chr7", "chr7"),
        ("22", "chr22"),
        ("X", "chrX"),
        ("chr23", "chrX"),
        ("23", "chrX"),
        ("Y", "chrY"),
        ("MT", "chrM"),
        ("M", "chrM"),
        ("chrUn_KI270302v1", "chrUn_KI270302v1"),
        ("07", "07"),
        ("PAR", "PAR"),
    ],
)
def test_canonical_chrom(name, canonical):
    assert canonical_chrom(name) == canonical


def test_chrom_aliases_cover_every_spelling():
    assert chrom_aliases("7") == ("7", "chr7")
    assert chrom_aliases("chr7") == ("7", "chr7")
    assert "23" in chrom_aliases("X") and "chrX" in chrom_aliases("23")
    assert set(chrom_aliases("MT")) >= {"chrM", "M", "MT"}
    assert chrom_aliases("chrUn_x") == ("chrUn_x",)


@pytest.mark.parametrize("asked", ["C9orf72", "c9orf72", "C9ORF72"])
def test_position_for_gene_ignores_case(isolated_annotation_db, asked):
    from annotations.db import position_for_gene

    _seed()
    g = position_for_gene(asked)
    assert g is not None
    assert (g.gene_symbol, g.chrom, g.start_pos) == ("C9orf72", "chr9", 27_546_545)


@pytest.mark.parametrize("asked", ["C19orf12", "C12ORF57", "brca1"])
def test_position_for_gene_resolves_mixed_case_symbols(isolated_annotation_db, asked):
    from annotations.db import position_for_gene

    _seed()
    assert position_for_gene(asked).gene_symbol.upper() == asked.upper()


def test_position_for_gene_prefers_the_exact_spelling(isolated_annotation_db):
    from annotations.db import position_for_gene

    _seed([("Abc1", "chr1", 10, 20, "+", None, None), ("ABC1", "chr2", 30, 40, "+", None, None)])  # fmt: skip
    assert position_for_gene("ABC1").chrom == "chr2"
    assert position_for_gene("Abc1").chrom == "chr1"
    assert position_for_gene("abc1").gene_symbol in ("ABC1", "Abc1")


def test_position_for_gene_closes_its_connection(isolated_annotation_db, monkeypatch):
    import annotations.db as adb

    _seed()
    opened = []
    real = adb.get_connection

    def tracking():
        conn = real()
        opened.append(conn)
        return conn

    monkeypatch.setattr(adb, "get_connection", tracking)
    adb.position_for_gene("BRCA1")
    adb.position_for_gene("NOPE")
    assert len(opened) == 2
    for conn in opened:
        with pytest.raises(Exception):
            conn.execute("SELECT 1")


@pytest.mark.parametrize("chrom", ["7", "chr7"])
def test_gene_index_and_gene_at_match_any_contig_spelling(
    isolated_annotation_db, chrom
):
    from annotations.db import GeneIndex, gene_at

    _seed()
    index = GeneIndex()
    assert [g.gene_symbol for g in index.at(chrom, 1_500)] == ["G7"]
    assert [g.gene_symbol for g in gene_at(chrom, 1_500)] == ["G7"]
    assert index.at(chrom, 2_001) == []


@pytest.mark.parametrize(
    "chrom, pos, symbol",
    [("X", 2_850_000, "XG"), ("23", 2_850_000, "XG"), ("MT", 4_000, "MT-ND1")],
)
def test_gene_index_maps_sex_and_mito_aliases(
    isolated_annotation_db, chrom, pos, symbol
):
    from annotations.db import GeneIndex, gene_at

    _seed()
    assert [g.gene_symbol for g in GeneIndex().at(chrom, pos)] == [symbol]
    assert [g.gene_symbol for g in gene_at(chrom, pos)] == [symbol]


@pytest.mark.parametrize("backend", ["chdb", "duckdb"])
def test_db_trio_finds_genes_on_a_bare_contig_vcf(
    backend, vcfclick_home, tmp_path, monkeypatch
):
    """A compound-het pair on contig `7` (not `chr7`), one het from each
    parent, inside G7: before the fix GeneIndex.at('7', ...) found nothing, so
    the gene column was empty and comphet had no gene."""
    store = tmp_path / "ann.duckdb"
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(store))
    import annotations.db as adb

    monkeypatch.setattr(adb, "DUCKDB_PATH", store)
    _seed()
    het, ref = "0/1:15,15:99:30", "0/0:30,0:99:30"
    sites = [
        Site("7", 1_200, "A", "G", (het, het, ref), info="AF=0.002", fmt="GT:AD:GQ:DP"),
        Site("7", 1_800, "C", "T", (het, ref, het), info="AF=0.002", fmt="GT:AD:GQ:DP"),
    ]
    vcf = write_vcf(tmp_path / "bare.vcf", TRIO_SAMPLES, sites, contigs=("7",))
    ped = tmp_path / "trio.ped"
    ped.write_text(TRIO_PED)
    run_cli(vcfclick_home, backend, "db", "create", "t")
    run_cli(
        vcfclick_home, backend, "db", "ingest", "t", str(vcf),
        "--ingest-id", "g1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(vcfclick_home, backend, "db", "ped", "t", str(ped))
    out = run_cli(
        vcfclick_home,
        backend,
        "db",
        "trio",
        "t",
        "--proband",
        "CHILD",
        "--format",
        "json",
    ).stdout
    doc = json.loads(out)
    assert [g["gene"] for g in doc["models"]["comphet"]["genes"]] == ["G7"]
    dominant = doc["models"]["dominant"]["candidates"]
    assert {c["pos"]: c["genes"] for c in dominant} == {1_200: ["G7"], 1_800: ["G7"]}


def test_db_gene_uses_the_same_alias_groups():
    """storage.gene_query (db gene, the MCP variants_in_gene tool) takes its
    spellings from storage.contigs: a chrX gene also matches chr23 and 23."""
    from storage.gene_query import chrom_aliases as gene_query_aliases

    assert gene_query_aliases("chrX") == ["chrX", "X", "chr23", "23"]
    assert gene_query_aliases("chr17") == ["17", "chr17"]
    assert set(gene_query_aliases("chrM")) >= {"chrM", "M", "MT"}
