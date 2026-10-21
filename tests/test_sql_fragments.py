"""storage.sql fragment builders, as text. Their behaviour on the engines is
tested in tests/test_family_checks.py (which runs them on chDB and DuckDB)."""

from __future__ import annotations

import pytest

from storage import sql


def test_str_list_quotes_every_value():
    assert sql.str_list(["a", "b"]) == "'a', 'b'"
    assert sql.str_list([]) == ""
    # quotes and backslashes go through sql_quote_str, never concatenated
    assert sql.str_list(["it's", "a\\b"]) == "'it''s', 'a\\\\b'"
    assert sql.str_list(iter(("x",))) == "'x'"


def test_autosome_filter_lists_every_spelling():
    f = sql.autosome_filter("g.chrom")
    assert f.startswith("g.chrom IN ('1', '2', ")
    assert "'chr22'" in f and "'chrX'" not in f and "'X'" not in f
    with pytest.raises(ValueError):
        sql.autosome_filter("chrom; DROP TABLE x")


def test_sample_join_is_an_inner_join_keyed_on_the_site():
    j = sql.sample_join("f", "DAD")
    assert j == (
        "INNER JOIN genotypes f ON f.ingest_id = g.ingest_id AND f.chrom = g.chrom "
        "AND f.pos = g.pos AND f.ref = g.ref AND f.alt = g.alt "
        "AND f.sample_id = 'DAD'"
    )
    assert "LEFT" not in j
    assert sql.sample_join("m", "o'brien", base="p").endswith(
        "m.sample_id = 'o''brien'"
    )
    for bad in ("f; --", "f x", "1f", ""):
        with pytest.raises(ValueError):
            sql.sample_join(bad, "S")


@pytest.mark.parametrize(
    "backend,expected",
    [
        (
            "duckdb",
            "(CAST(CAST(g.ad_alt AS BIGINT) AS DOUBLE) / "
            "CAST(CAST(g.ad_alt AS BIGINT) + CAST(g.ad_ref AS BIGINT) AS DOUBLE))",
        ),
        (
            "chdb",
            "(toFloat64(toInt64(g.ad_alt)) / "
            "toFloat64(toInt64(g.ad_alt) + toInt64(g.ad_ref)))",
        ),
    ],
)
def test_allele_fraction_widens_before_adding(backend, expected, monkeypatch):
    """Both operands are widened to 64 bits before the add: DuckDB adds two
    USMALLINT columns in USMALLINT and raises on 65535 + 1."""
    monkeypatch.setenv("VCFCLICK_BACKEND", backend)
    assert sql.allele_fraction_sql("g.ad_alt", "g.ad_ref") == expected


def test_genotype_pivot_sql_shape():
    q = sql.genotype_pivot_sql("i1", {"p": "CHILD", "f": "DAD"}, ("gt", "gq"))
    assert "max(CASE WHEN sample_id = 'CHILD' THEN gt END) AS p_gt" in q
    assert "max(CASE WHEN sample_id = 'DAD' THEN gq END) AS f_gq" in q
    assert "WHERE ingest_id = 'i1' AND sample_id IN ('CHILD', 'DAD')" in q
    assert q.endswith("GROUP BY ingest_id, chrom, pos, ref, alt")
    assert "LEFT" not in q
    with pytest.raises(ValueError):
        sql.genotype_pivot_sql("i1", {})
    with pytest.raises(ValueError):
        sql.genotype_pivot_sql("i1", {"p x": "C"})


def test_db_qc_uses_the_shared_x_names():
    """db qc's chrX test is built from contigs.X_NAMES, not its own list."""
    from cli import db_qc

    assert not hasattr(db_qc, "_CHRX")
    sql = db_qc._qc_sql()
    assert "chrom IN ('chrX', 'X', 'chr23', '23')" in sql
