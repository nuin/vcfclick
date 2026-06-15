# Examples

Walk-throughs of vcfclick's natural-language-over-SQL workflow. Each
example uses the same dataset (the 1000 Genomes Phase 3 BRCA1 cohort
demo bundle the README points at) so a reader can reproduce the
queries verbatim.

- [`brca1-cohort.md`](brca1-cohort.md) — five canonical cohort
  questions, the MCP tools the LLM calls for each, the SQL it
  generates, and the actual chDB output.

Every SQL block and every result block in these files is **real** —
queries were run against the demo bundle on a vanilla install, not
mocked or projected. The LLM's English framing and tool-call sequence
are *illustrative* of what Claude (or any MCP client) would produce
against the schema in
[`vcfclick_mcp/server.py`](../vcfclick_mcp/server.py); a different
client will phrase things slightly differently but should arrive at
the same SQL.

For setup and reference material, see:

- [Getting started](../docs/GETTING_STARTED.md)
- [MCP and annotations](../docs/MCP.md)
- [Schema reference](../docs/SCHEMA.md)
