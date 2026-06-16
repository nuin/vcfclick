# Web UI (`vcfclick web`)

`vcfclick web` starts a small local server and opens a browser UI over a
cohort database — a SQL explorer, a natural-language→SQL box, and trio /
combine panels. It is an optional extra; everything runs on your machine
against the embedded backend, with no hosted service.

## Install and run

```bash
uv tool install "vcfclick[web]"     # or: pipx install "vcfclick[web]"

vcfclick web epilepsy_2026
# vcfclick web → http://127.0.0.1:8765  (opens your browser)
```

Options:

| Flag | Default | Meaning |
|---|---|---|
| `--port` | `8765` | Port to bind. |
| `--host` | `127.0.0.1` | Interface to bind. Localhost only by default — there is **no authentication**, so only change this if you understand the exposure. |
| `--no-browser` | off | Don't open a browser automatically. |

Stop the server with `Ctrl-C`.

## What's in the UI

- **Schema sidebar** — the cohort's tables and columns (from the locked
  Arrow schemas), click to expand.
- **SQL** — write a query, run it, see the result table. Read-only:
  `INSERT`/`UPDATE`/`DELETE`/`DROP`/`ALTER`/`CREATE` are rejected.
- **Ask** — type a question in English; the server uses the same schema
  briefing the MCP server gives an LLM to generate SQL, shows you the
  SQL, and runs it. Bring your own API key (Gemini or Anthropic) — it
  stays in your browser and is sent only to the provider you pick.
- **Trio** — pick a proband, an inheritance model (de-novo / recessive /
  dominant) and quality thresholds; reuses `db trio`. De-novo / dominant
  show a note if the cohort lacks stored hom-reference rows (re-ingest
  with `--keep-reference`).
- **Combine** — paste two call sets and run the GATK3-style
  `CombineVariants` merge with `set=` provenance.

## How it works

The server is a thin HTTP layer (FastAPI) over existing vcfclick
internals — it adds no new query logic. The SQL path is the same
`get_session(...).query(...)` used by `vcfclick db query` and the MCP
`run_sql` tool; the trio panel calls the `db trio` SQL builders; the
combine panel calls the same `combine_vcfs` engine as the CLI. The UI is
a single self-contained page served by the server.

It binds to `127.0.0.1` and is single-user with no auth — it is meant to
run on the same machine as your databases, like the TUI, not to be
exposed to a network.

## Relation to the other interfaces

| Interface | Use it for |
|---|---|
| CLI (`vcfclick db query`) | scripting, pipelines, headless servers |
| [TUI](TUI.md) (`vcfclick tui`) | a terminal cohort browser |
| Web (`vcfclick web`) | a browser UI: SQL, NL→SQL, trio, combine |
| [MCP](MCP.md) | letting an LLM client write visible SQL for you |
