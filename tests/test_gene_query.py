"""Tests for gene-level queries: `vcfclick db gene` and the MCP tool
`variants_in_gene`.

Fixture: tests/fixtures/tiny.vcf.gz (chr1:100, 250, 500, 750, 900; samples
S1-S3). A throwaway annotation store defines GENEX = chr1:200-800, which
covers 250, 500 and 750.
"""

from __future__ import annotations

import asyncio
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")


@pytest.fixture
def gene_store(tmp_path, monkeypatch) -> Path:
    """An annotation store with one gene, reachable by CLI subprocesses."""
    path = tmp_path / "ann.duckdb"
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(path))
    import annotations.db as adb

    conn = adb.get_connection()
    conn.execute(
        "INSERT INTO refseq_genes VALUES ('GENEX', 'chr1', 200, 800, '+', '1', 'test gene')"
    )
    conn.close()  # DuckDB allows one writer; release it for the subprocesses
    return path


def _vc(home: Path, *args: str, ok: bool = True):
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    if ok:
        assert r.returncode == 0, f"{' '.join(args)} failed:\n{r.stderr}"
    return r


def _db(home: Path, vcf: Path, name: str = "g") -> None:
    _vc(home, "db", "create", name)
    _vc(
        home,
        "db",
        "ingest",
        name,
        str(vcf),
        "--cohort",
        "c",
        "--ingest-id",
        "i1",
        "--serial",
    )


def _gene(home: Path, *args: str) -> dict:
    return json.loads(_vc(home, "db", "gene", "g", *args, "--format", "json").stdout)


def test_variants_in_gene_with_carrier_counts(vcfclick_home, tiny_vcf, gene_store):
    _db(vcfclick_home, tiny_vcf)
    res = _gene(vcfclick_home, "GENEX")
    assert res["gene"]["gene_symbol"] == "GENEX"
    assert res["region"] == {"chrom": "chr1", "start": 200, "end": 800}
    by_pos = {r["pos"]: r for r in res["variants"]}
    assert sorted(by_pos) == [250, 500, 750]
    assert by_pos[250]["carriers"] == 1 and by_pos[250]["hom_alt"] == 0
    assert by_pos[500]["carriers"] == 3 and by_pos[500]["hom_alt"] == 0
    assert by_pos[750]["carriers"] == 3 and by_pos[750]["hom_alt"] == 2
    assert by_pos[500]["ref"] == "G" and by_pos[500]["alt"] == "A"
    assert res["row_count"] == 3 and res["truncated"] is False
    assert "FROM variants" in res["sql"]  # the SQL is returned for auditability


def test_flank_limit_and_case(vcfclick_home, tiny_vcf, gene_store):
    _db(vcfclick_home, tiny_vcf)
    assert [
        r["pos"] for r in _gene(vcfclick_home, "genex", "--flank", "200")["variants"]
    ] == [100, 250, 500, 750, 900]
    limited = _gene(vcfclick_home, "GENEX", "--limit", "2")
    assert limited["row_count"] == 2 and limited["truncated"] is True


def test_unknown_gene_says_how_to_load_genes(vcfclick_home, tiny_vcf, gene_store):
    _db(vcfclick_home, tiny_vcf)
    r = _vc(vcfclick_home, "db", "gene", "g", "NOPE1", ok=False)
    assert r.returncode != 0
    assert "NOPE1" in r.stderr and "annotations load" in r.stderr


def test_chromosome_naming_mismatch_is_bridged(vcfclick_home, tmp_path, gene_store):
    """Genes are stored as chr1; a VCF may use 1. Both must match."""
    vcf = tmp_path / "nochr.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"
        "1\t300\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "1\t5000\t.\tC\tT\t50\tPASS\t.\tGT\t1/1\n"
    )
    _db(vcfclick_home, vcf)
    res = _gene(vcfclick_home, "GENEX")
    assert [(r["chrom"], r["pos"]) for r in res["variants"]] == [("1", 300)]


def test_table_and_tsv_output(vcfclick_home, tiny_vcf, gene_store):
    _db(vcfclick_home, tiny_vcf)
    table = _vc(vcfclick_home, "db", "gene", "g", "GENEX").stdout
    assert "GENEX" in table and "chr1:200-800" in table and "carriers" in table
    tsv = _vc(
        vcfclick_home, "db", "gene", "g", "GENEX", "--format", "tsv"
    ).stdout.splitlines()
    assert tsv[0].split("\t")[:4] == ["chrom", "pos", "ref", "alt"]
    assert len(tsv) == 4


# --- MCP tool ----------------------------------------------------------------


def test_mcp_tool_is_registered_with_a_symbol_argument():
    from vcfclick_mcp.server import mcp

    tools = {
        t.name: t for t in asyncio.new_event_loop().run_until_complete(mcp.list_tools())
    }
    assert "variants_in_gene" in tools
    props = tools["variants_in_gene"].inputSchema["properties"]
    assert props["symbol"]["type"] == "string"
    assert props["flank"]["type"] == "integer" and props["limit"]["type"] == "integer"


def test_mcp_tool_end_to_end_over_stdio(vcfclick_home, tiny_vcf, gene_store):
    """The real server process, over the MCP protocol, against a real database."""
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    _db(vcfclick_home, tiny_vcf)
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(vcfclick_home)
    env["VCFCLICK_DB_NAME"] = "g"
    env["PYTHONPATH"] = str(REPO)
    params = StdioServerParameters(
        command=sys.executable, args=["-m", "vcfclick_mcp.server"], env=env
    )

    async def call():
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                found = await session.call_tool("variants_in_gene", {"symbol": "GENEX"})
                missing = await session.call_tool(
                    "variants_in_gene", {"symbol": "NOPE1"}
                )
                return found.structuredContent, missing.structuredContent

    found, missing = asyncio.new_event_loop().run_until_complete(call())
    found = found.get("result", found)
    assert [v["pos"] for v in found["variants"]] == [250, 500, 750]
    assert found["variants"][2]["hom_alt"] == 2
    assert missing.get("result", missing) is None
