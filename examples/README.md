# Examples

A **self-contained, runnable demo** — it downloads its own example VCFs
and walks the headline flows end to end: ingest → SQL, GIAB-validated
trio de novo, compound-het + gnomAD rarity filtering, sample QC, and
`combine`. Two formats:

- [`vcfclick-demo.ipynb`](vcfclick-demo.ipynb) — Jupyter
  ([open in Colab](https://colab.research.google.com/github/nuin/vcfclick/blob/main/examples/vcfclick-demo.ipynb)).
  Zero-setup one-click; every cell carries the real output of a v0.7.0
  run, so it renders with results on GitHub and doubles as a smoke-test.
- [`vcfclick_demo.py`](vcfclick_demo.py) — an **interactive**
  [marimo](https://marimo.io) version: plain-Python, git-friendly, and
  reactive. A live SQL editor, a trio model dropdown + gnomAD-rarity
  slider, a `--keep-reference` toggle that flips the de-novo count
  between 4 and 0, and a `--min-callsets` slider — change a control and
  the dependent cells re-run. Run it (deps pulled via inline script
  metadata): `uvx marimo edit --sandbox examples/vcfclick_demo.py`. Needs
  a real Python kernel — vcfclick's native engines can't run in
  browser-WASM.

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
