from __future__ import annotations

import pytest

from storage import apply_schema, get_session
from tui.services import (
    DatabaseError,
    LocusInputError,
    ParsedLocus,
    TuiServiceError,
    database_summary,
    execute_sql,
    list_database_names,
    parse_locus_input,
    validate_database,
)


@pytest.fixture(autouse=True)
def force_duckdb_backend(monkeypatch):
    monkeypatch.setenv("VCFCLICK_BACKEND", "duckdb")


def test_parse_range_with_commas():
    locus = parse_locus_input("chr17:43,044,295-43,125,483")
    assert locus == ParsedLocus(
        kind="range",
        label="chr17:43044295-43125483",
        chrom="chr17",
        start_pos=43044295,
        end_pos=43125483,
        gene_symbol=None,
    )


def test_parse_gene_symbol():
    locus = parse_locus_input("BRCA1")
    assert locus == ParsedLocus(
        kind="gene",
        label="BRCA1",
        chrom=None,
        start_pos=None,
        end_pos=None,
        gene_symbol="BRCA1",
    )


@pytest.mark.parametrize(
    "raw",
    ["", "chr17", "chr17:431-430", "chr17:start-end", "BRCA1 BRCA2"],
)
def test_parse_invalid_locus(raw):
    with pytest.raises(LocusInputError):
        parse_locus_input(raw)


def test_list_database_names_uses_vcfclick_home(vcfclick_home):
    (vcfclick_home / "dbs" / "alpha").mkdir(parents=True)
    (vcfclick_home / "dbs" / "beta").mkdir(parents=True)
    assert list_database_names() == ["alpha", "beta"]


def test_validate_database_rejects_missing(vcfclick_home):
    with pytest.raises(DatabaseError):
        validate_database("missing")


def test_validate_database_wraps_unsafe_name(vcfclick_home):
    with pytest.raises(DatabaseError):
        validate_database("../unsafe")


def test_execute_sql_returns_structured_rows(vcfclick_home, monkeypatch):
    name = "smoke-sql"
    validate_home = vcfclick_home / "dbs" / name
    validate_home.mkdir(parents=True)
    monkeypatch.setenv("VCFCLICK_DB_NAME", name)
    get_session(name)
    apply_schema()

    result = execute_sql(name, "SELECT count() AS n FROM variants")
    assert result.sql == "SELECT count() AS n FROM variants"
    assert result.columns == ["n"]
    assert result.rows == [[0]]
    assert result.row_count == 1


def test_execute_sql_wraps_malformed_sql(vcfclick_home, monkeypatch):
    name = "bad-sql"
    validate_home = vcfclick_home / "dbs" / name
    validate_home.mkdir(parents=True)
    monkeypatch.setenv("VCFCLICK_DB_NAME", name)
    get_session(name)
    apply_schema()

    with pytest.raises(TuiServiceError):
        execute_sql(name, "SELECT FROM")


def test_database_summary_counts_missing_schema_as_none(vcfclick_home):
    name = "schema-missing"
    (vcfclick_home / "dbs" / name).mkdir(parents=True)

    summary = database_summary(name)
    assert summary.name == name
    assert summary.variants is None
    assert summary.genotypes is None
    assert summary.samples is None
    assert summary.ingestions is None


def test_database_summary_counts_empty_schema(vcfclick_home, monkeypatch):
    name = "smoke-summary"
    (vcfclick_home / "dbs" / name).mkdir(parents=True)
    monkeypatch.setenv("VCFCLICK_DB_NAME", name)
    get_session(name)
    apply_schema()

    summary = database_summary(name)
    assert summary.name == name
    assert summary.variants == 0
    assert summary.genotypes == 0
    assert summary.samples == 0
    assert summary.ingestions == 0
