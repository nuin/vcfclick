"""MCP server integration tests.

Two tiers of coverage:

1. **In-process** — invoke `FastMCP.call_tool()` directly. Catches
   tool-registration drift, argument-schema mismatches, return-shape
   bugs, and SCHEMA_DESCRIPTION drift. Fast, no subprocess.

2. **Subprocess** — spawn the actual MCP server over stdio (the
   transport Claude Desktop uses) and list tools via the real
   protocol. Catches server startup / transport / tool-registration
   wire issues that the in-process path can't see.

Notes:
  - The DuckDB-backed tool tests use the conftest `isolated_annotation_db`
    fixture so they don't touch the user's real annotation store.
  - The `run_sql` tool is NOT exercised here: chDB holds a single
    EmbeddedServer per process, so running it in-process clashes with
    the rest of the test suite. The subprocess startup test below
    proves run_sql is registered and reachable; deeper coverage of
    run_sql against real data is in tests/test_cli.py via the CLI.
"""

from __future__ import annotations

import asyncio
import os
import shutil
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
CLINVAR_FIXTURE = Path(__file__).parent / "fixtures" / "clinvar_mini.vcf.gz"

EXPECTED_TOOLS = {
    "get_schema",
    "run_sql",
    "position_for_gene",
    "gene_at",
    "clinvar_lookup",
}


def _run(coro):
    """Run an async coroutine inside a fresh loop — keeps the tests
    sync-callable regardless of pytest's anyio plugin state."""
    return asyncio.new_event_loop().run_until_complete(coro)


def _unwrap(structured):
    """FastMCP wraps every structured tool return in {"result": <value>}.
    This unwraps that consistently. None and list returns both come back
    as {"result": None} / {"result": [...]}."""
    if isinstance(structured, dict) and set(structured.keys()) == {"result"}:
        return structured["result"]
    return structured


# ─────────────────────── tool registry / schemas ───────────────────────


def test_all_advertised_tools_are_registered():
    from vcfclick_mcp.server import mcp

    tools = _run(mcp.list_tools())
    assert {t.name for t in tools} == EXPECTED_TOOLS


def test_tool_argument_schemas_match_signatures():
    """The JSON-Schema FastMCP generates from each tool's Python
    signature is what the LLM sees — drift here means the LLM calls
    tools with the wrong arg shapes."""
    from vcfclick_mcp.server import mcp

    schemas = {t.name: t.inputSchema for t in _run(mcp.list_tools())}

    assert schemas["get_schema"]["properties"] == {}

    assert schemas["run_sql"]["properties"]["query"]["type"] == "string"
    assert "query" in schemas["run_sql"]["required"]

    assert schemas["position_for_gene"]["properties"]["symbol"]["type"] == "string"

    ga = schemas["gene_at"]["properties"]
    assert ga["chrom"]["type"] == "string"
    assert ga["pos"]["type"] == "integer"

    cv = schemas["clinvar_lookup"]["properties"]
    assert cv["chrom"]["type"] == "string"
    assert cv["pos"]["type"] == "integer"
    assert cv["ref"]["type"] == "string"
    assert cv["alt"]["type"] == "string"


# ─────────────────────── get_schema (briefing) ───────────────────────


def test_get_schema_returns_briefing_text():
    from vcfclick_mcp.server import SCHEMA_DESCRIPTION, mcp

    _, structured = _run(mcp.call_tool("get_schema", {}))
    assert _unwrap(structured) == SCHEMA_DESCRIPTION


def test_briefing_includes_non_obvious_invariants():
    """The pieces of the schema briefing the LLM most commonly gets
    wrong without explicit guidance — sparse-table semantics and the
    GQ/DP-NULL silent-failure trap added after the 1000G demo."""
    from vcfclick_mcp.server import mcp

    _, structured = _run(mcp.call_tool("get_schema", {}))
    text = str(_unwrap(structured))

    # Sparse-table convention — most non-obvious thing about the schema.
    assert "SPARSE TABLE" in text
    assert "0/0" in text
    # GQ/DP-NULL trap — added after the BRCA1 demo silently returned 0.
    assert "gq >= 20" in text
    assert "NULL silently fails" in text
    # Cross-ingestion: rows are NOT merged across uploads.
    assert "NOT merged across ingestions" in text
    # Cohort AF denominator: must teach the LLM to compute cohort size
    # against `samples` ALONE, not through the join to sparse genotypes.
    # Counting via the join only sees non-reference samples → inflated
    # AF. Codex round 6 caught this exact failure mode in an earlier
    # iteration of the briefing.
    assert "cohort_size" in text or "FROM samples\n" in text
    assert "CROSS JOIN" in text


# ─────────────────────── DuckDB-backed tools ───────────────────────


def _seed_gene(symbol="BRCA1", chrom="chr17", start=43044295, end=43170245):
    from annotations.db import get_connection

    get_connection().execute(
        "INSERT INTO refseq_genes VALUES (?, ?, ?, ?, '-', 'ENSG_test', 'protein_coding')",
        [symbol, chrom, start, end],
    )


def test_position_for_gene_returns_coords(isolated_annotation_db):
    from vcfclick_mcp.server import mcp

    _seed_gene()

    _, structured = _run(mcp.call_tool("position_for_gene", {"symbol": "BRCA1"}))
    assert _unwrap(structured) == {
        "gene_symbol": "BRCA1",
        "chrom": "chr17",
        "start_pos": 43044295,
        "end_pos": 43170245,
        "strand": "-",
    }


def test_position_for_gene_returns_none_for_unknown(isolated_annotation_db):
    from vcfclick_mcp.server import mcp

    _seed_gene()  # one gene, but not the one we'll query

    _, structured = _run(mcp.call_tool("position_for_gene", {"symbol": "NOSUCHGENE"}))
    assert _unwrap(structured) is None


def test_gene_at_returns_overlapping_genes(isolated_annotation_db):
    from vcfclick_mcp.server import mcp

    _seed_gene()

    _, structured = _run(mcp.call_tool("gene_at", {"chrom": "chr17", "pos": 43100000}))
    rows = _unwrap(structured)
    assert len(rows) == 1
    assert rows[0]["gene_symbol"] == "BRCA1"


def test_gene_at_empty_outside_any_gene(isolated_annotation_db):
    from vcfclick_mcp.server import mcp

    _seed_gene()

    _, structured = _run(mcp.call_tool("gene_at", {"chrom": "chr17", "pos": 1}))
    assert _unwrap(structured) == []


def test_clinvar_lookup_hits_known_variant(isolated_annotation_db):
    """Full chain: ClinVar loader → DuckDB → MCP tool."""
    from annotations.loaders.clinvar import load
    from vcfclick_mcp.server import mcp

    load(CLINVAR_FIXTURE)

    _, structured = _run(
        mcp.call_tool(
            "clinvar_lookup",
            {"chrom": "chr17", "pos": 43044295, "ref": "C", "alt": "T"},
        )
    )
    result = _unwrap(structured)
    assert result["clin_sig"] == "Pathogenic"
    assert "Hereditary" in result["condition"]
    assert result["clinvar_id"] == "1001"


def test_clinvar_lookup_miss_returns_none(isolated_annotation_db):
    from annotations.loaders.clinvar import load
    from vcfclick_mcp.server import mcp

    load(CLINVAR_FIXTURE)

    _, structured = _run(
        mcp.call_tool(
            "clinvar_lookup",
            {"chrom": "chr1", "pos": 99999, "ref": "A", "alt": "T"},
        )
    )
    assert _unwrap(structured) is None


# ─────────────── subprocess + stdio transport smoke test ───────────────


def test_server_starts_and_lists_tools_over_stdio(tmp_path):
    """Spawn the actual MCP server, talk JSON-RPC over stdio, and
    confirm the tool registry comes through the wire. This is the only
    place we exercise the real transport — everything else is in-process.

    Covers: server bootstrap, stdio framing, tool registration over the
    protocol. Does NOT cover: tool execution correctness (the in-process
    tests above do that) or chDB-backed tools (one-server-per-process
    constraint plus this test's own process state would conflict)."""
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    # Isolated VCFCLICK_HOME so the server doesn't try to open any
    # real cohort DB during startup.
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(tmp_path / "vcfclick")
    env.pop("VCFCLICK_DB_NAME", None)
    env["PYTHONPATH"] = str(REPO)

    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "vcfclick_mcp.server"],
        env=env,
    )

    async def _list_tools_via_stdio():
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                resp = await session.list_tools()
                return {t.name for t in resp.tools}

    names = _run(_list_tools_via_stdio())
    assert names == EXPECTED_TOOLS
