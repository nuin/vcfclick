"""storage.sql fragment builders, as text. Their behaviour on the engines is
tested in tests/test_family_checks.py (which runs them on chDB and DuckDB)."""

from __future__ import annotations

import pytest

from storage import sql


@pytest.mark.parametrize(
    "backend,backslash,backslash_quote",
    [
        # ClickHouse reads `\` as an escape inside a literal: double it.
        ("chdb", "'a\\\\b'", "'\\\\'''"),
        # DuckDB reads `\` literally: doubling it changes the value.
        ("duckdb", "'a\\b'", "'\\'''"),
    ],
)
def test_sql_quote_str_escapes_for_the_backend(
    backend, backslash, backslash_quote, monkeypatch
):
    monkeypatch.setenv("VCFCLICK_BACKEND", backend)
    assert sql.sql_quote_str("plain") == "'plain'"
    assert sql.sql_quote_str("o'brien") == "'o''brien'"
    assert sql.sql_quote_str("a\\b") == backslash
    # A quote is always doubled, never backslash-escaped (`\'` is a literal
    # backslash plus the string terminator on DuckDB).
    assert sql.sql_quote_str("\\'") == backslash_quote


@pytest.mark.parametrize("backend", ["chdb", "duckdb"])
def test_str_list_quotes_every_value(backend, monkeypatch):
    monkeypatch.setenv("VCFCLICK_BACKEND", backend)
    assert sql.str_list(["a", "b"]) == "'a', 'b'"
    assert sql.str_list([]) == ""
    # quotes and backslashes go through sql_quote_str, never concatenated
    assert sql.str_list(["it's", "a\\b"]) == ("'it''s', " + sql.sql_quote_str("a\\b"))
    assert sql.str_list(iter(("x",))) == "'x'"


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
