"""One list of inheritance models, four consumers.

`cli.trio_analysis.CATEGORIES` is the definition. The `db trio --category`
choices, the model registry, the web UI's category select and the MCP
briefing all derive from it; this file fails when one of them is hand-edited
out of step (the web page once listed three of the four models and nothing
noticed).

In-process only and no chDB: nothing here opens a database.
"""

from __future__ import annotations

import re

import pytest

from cli.trio_analysis import (
    CATEGORIES,
    CATEGORY_CHOICES,
    DEFAULT_MAX_AF,
    GENE_CATEGORIES,
    ROW_CATEGORIES,
    X_CATEGORIES,
    _quality_gate,
)


def test_the_constant_is_the_row_models_then_the_gene_model():
    # then the sex- and region-aware X model (step X)
    assert CATEGORIES == ROW_CATEGORIES + GENE_CATEGORIES + X_CATEGORIES
    assert X_CATEGORIES == ("xlinked",)
    assert CATEGORY_CHOICES == CATEGORIES + ("all",)
    assert len(set(CATEGORY_CHOICES)) == len(CATEGORY_CHOICES)


def test_cli_choices_are_the_constant():
    from cli.db_trio import db_trio

    option = next(p for p in db_trio.params if p.name == "category")
    assert tuple(option.type.choices) == CATEGORY_CHOICES
    assert option.default == "all"


def test_the_model_registry_is_the_constant():
    from cli.trio_output import MODELS

    assert tuple(MODELS) == CATEGORIES


def test_the_web_select_lists_every_category_once_in_order():
    pytest.importorskip("fastapi")
    from vcfclick_web.page import INDEX_HTML

    select = re.search(r'<select id="t-cat">(.*?)</select>', INDEX_HTML, re.S)
    assert select, "the Trio tab has a category select"
    assert re.findall(r"<option>([^<]*)</option>", select.group(1)) == list(CATEGORIES)
    assert "__TRIO_CATEGORY_OPTIONS__" not in INDEX_HTML


def test_the_web_endpoint_accepts_exactly_the_categories(monkeypatch):
    """The whitelist is checked before any database is touched. The store
    lookup is stubbed to "no single ingestion", so a valid name gets past
    the whitelist and stops there with that message, and an unknown one is
    refused by name; nothing opens a database."""
    pytest.importorskip("fastapi")
    pytest.importorskip("httpx")
    from fastapi.testclient import TestClient

    from vcfclick_web.app import app

    monkeypatch.setattr("cli.db_trio._sole_ingest_id", lambda name: None)
    client = TestClient(app)
    refused = client.get("/api/trio", params={"proband": "C", "category": "nope"})
    assert refused.json() == {"error": "unknown category 'nope'"}
    for category in CATEGORIES:
        out = client.get("/api/trio", params={"proband": "C", "category": category})
        assert "single-ingestion" in out.json()["error"], category
    # "all" is a CLI summary mode, not a web listing
    out = client.get("/api/trio", params={"proband": "C", "category": "all"})
    assert out.json() == {"error": "unknown category 'all'"}


def test_the_mcp_briefing_names_every_category_the_cli_takes():
    from vcfclick_mcp.server import SCHEMA_DESCRIPTION

    assert (
        f"`vcfclick db trio --category` takes: {', '.join(CATEGORY_CHOICES)}."
        in SCHEMA_DESCRIPTION
    )


def test_the_mcp_briefing_states_the_gates_the_cli_applies():
    """The briefing once told the model that `db trio` applies GQ >= 20
    exactly, but the CLI is NULL-tolerant (a missing GQ/DP passes) and the
    generic convention drops those rows. Both are now stated, and the form
    quoted is the CLI's own."""
    from vcfclick_mcp.server import SCHEMA_DESCRIPTION
    from storage.sql import AF_EPSILON

    cli_form = _quality_gate("g", 20, 10).replace("g.", "")
    assert cli_form == "(gq IS NULL OR gq >= 20) AND (dp IS NULL OR dp >= 10)"
    assert cli_form in SCHEMA_DESCRIPTION
    assert "NULL-TOLERANT" in SCHEMA_DESCRIPTION
    assert "exactly these" not in SCHEMA_DESCRIPTION
    assert f"info_AF <= {DEFAULT_MAX_AF:g} + {AF_EPSILON!r}" in SCHEMA_DESCRIPTION
    # the cohort-derived AF behaviour and its enum
    from cli.trio_analysis import AF_GATE_VALUES

    for state in AF_GATE_VALUES:
        assert state in SCHEMA_DESCRIPTION, state
    # convention 2 (the strict form) is still there for generic queries
    assert "AND gq >= 20 AND dp >= 10" in SCHEMA_DESCRIPTION
