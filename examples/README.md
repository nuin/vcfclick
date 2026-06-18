# Examples

- [`vcfclick-demo.ipynb`](vcfclick-demo.ipynb) — a **self-contained,
  runnable demo notebook** ([open in Colab](https://colab.research.google.com/github/nuin/vcfclick/blob/main/examples/vcfclick-demo.ipynb)).
  It installs vcfclick, downloads its own example VCFs, and walks the
  headline flows end to end — ingest → SQL, GIAB-validated trio de novo,
  compound-het + gnomAD rarity filtering, sample QC, and `combine`. Every
  cell carries the real output of a v0.7.0 run, so it also doubles as a
  release smoke-test.

The markdown walk-throughs below cover the natural-language-over-SQL
workflow against the 1000 Genomes BRCA1 demo bundle:

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
