"""`annotations load-hgnc` and the batched gene resolver.

In process against an isolated annotation store (DuckDB, never chDB). The HGNC
file is the committed extract tests/fixtures/hgnc_aliases.tsv; nothing is
downloaded.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from tests.conftest import FIXTURES

GFF = FIXTURES / "scope_genes.gff3"
HGNC = FIXTURES / "hgnc_aliases.tsv"


@pytest.fixture
def genes_loaded(isolated_annotation_db):
    from annotations.loaders.gencode_genes import load

    load(GFF)
    return isolated_annotation_db


def _resolve(tokens):
    from annotations.db import resolve_gene_tokens

    return resolve_gene_tokens(tokens)


def _symbols(found, token):
    from annotations.db import gene_token_key

    return [(kind, g.gene_symbol) for kind, g in found.get(gene_token_key(token), [])]


def test_loader_reads_approved_rows_only(genes_loaded):
    from annotations.db import get_connection
    from annotations.loaders.hgnc_aliases import load

    counts = load(HGNC)
    assert counts == {"hgnc_id": 3, "entrez": 3, "previous": 3, "alias": 5}
    conn = get_connection()
    symbols = {
        r[0]
        for r in conn.execute("SELECT DISTINCT symbol FROM gene_aliases").fetchall()
    }
    conn.close()
    assert "WITHDRAWN1" not in symbols  # an Entry Withdrawn row is not read


def test_loader_splits_quoted_pipe_lists(genes_loaded):
    from annotations.loaders.hgnc_aliases import iter_aliases

    rows = set(iter_aliases(HGNC))
    assert ("NBIA4", "previous", "C19orf12") in rows
    assert ("SPG43", "previous", "C19orf12") in rows


def test_loader_refuses_a_file_without_symbols(tmp_path, isolated_annotation_db):
    from annotations.loaders.hgnc_aliases import load

    bad = tmp_path / "x.tsv"
    bad.write_text("a\tb\n1\t2\n")
    with pytest.raises(ValueError, match="not an HGNC TSV"):
        load(bad)


def test_current_symbols_and_ensembl_ids_resolve_without_the_table(genes_loaded):
    found = _resolve(["c9orf72", "C19ORF12", "ENSG00000111678.99", "ENSG00000147894"])
    assert _symbols(found, "c9orf72") == [("symbol", "C9orf72")]
    assert _symbols(found, "C19ORF12") == [("symbol", "C19orf12")]
    assert _symbols(found, "ENSG00000111678.99") == [("ensembl", "C12orf57")]
    assert _symbols(found, "ENSG00000147894") == [("ensembl", "C9orf72")]


@pytest.mark.parametrize(
    "token", ["FTDALS1", "ALSFTD", "HGNC:28337", "203228", "NBIA4"]
)
def test_aliases_resolve_only_with_the_table(genes_loaded, token):
    from annotations.loaders.hgnc_aliases import load

    assert _resolve([token]) == {}
    load(HGNC)
    kinds = _symbols(_resolve([token]), token)
    assert len(kinds) == 1 and kinds[0][1] in ("C9orf72", "C19orf12")


def test_an_alias_of_two_genes_lists_both(genes_loaded):
    from annotations.loaders.hgnc_aliases import load

    load(HGNC)
    assert sorted(s for _, s in _symbols(_resolve(["SHARED1"]), "SHARED1")) == [
        "C12orf57",
        "C19orf12",
    ]


def test_a_current_symbol_outranks_a_previous_one(genes_loaded):
    """WITHDRAWN1's row would make C9orf72 a previous symbol of another gene;
    it is not approved, and a current symbol wins anyway."""
    from annotations.loaders.hgnc_aliases import load

    load(HGNC)
    assert _symbols(_resolve(["C9ORF72"]), "C9ORF72")[0] == ("symbol", "C9orf72")


def test_two_thousand_tokens_are_one_query_on_one_connection(genes_loaded, monkeypatch):
    import annotations.db as adb

    statements = []
    real = adb.get_connection
    opened = []

    class Spy:
        def __init__(self, conn):
            self._conn = conn

        def execute(self, sql, *a):
            statements.append(sql)
            return self._conn.execute(sql, *a)

        def __getattr__(self, name):
            return getattr(self._conn, name)

    def spying():
        opened.append(1)
        return Spy(real())

    monkeypatch.setattr(adb, "get_connection", spying)
    tokens = [f"NOPE{i}" for i in range(1990)] + ["C9orf72"] * 10
    found = adb.resolve_gene_tokens(tokens)
    assert len(opened) == 1
    lookups = [s for s in statements if "refseq_genes" in s]
    assert len(lookups) == 1
    assert list(found) == ["C9ORF72"]


def test_load_hgnc_cli(genes_loaded):
    from cli.main import cli

    r = CliRunner().invoke(cli, ["annotations", "load-hgnc", str(HGNC)])
    assert r.exit_code == 0, r.output
    assert "loaded   14 gene aliases" in r.output


def test_par_gene_keeps_the_first_copy(genes_loaded):
    found = _resolve(["XG"])
    [(kind, gene)] = found["XG"]
    assert (kind, gene.chrom) == ("symbol", "chrX")
