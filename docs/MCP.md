# MCP And Annotations

vcfclick includes an MCP server so an MCP client can ask genomics
questions, receive visible SQL, and run that SQL against a local cohort
database.

The MCP layer is intentionally transparent:

1. annotation tools resolve gene symbols and ClinVar alleles from a
   local DuckDB reference store;
2. generated SQL runs against the selected vcfclick database;
3. the client can show the plan, SQL, columns, and rows.

MCP is optional. The CLI and TUI work without it.

## Tools

The server exposes:

| Tool | Purpose |
|---|---|
| `get_schema` | Return the schema and query rules given to the client. |
| `run_sql` | Execute SQL against the active vcfclick database. |
| `position_for_gene` | Resolve an HGNC gene symbol to GRCh38 coordinates. |
| `gene_at` | Return genes overlapping one position. |
| `clinvar_lookup` | Return ClinVar significance for one allele. |
| `gnomad_lookup` | Return the gnomAD allele frequency (and popmax) for one allele. |

## Load Annotation Data

Gene lookup requires a gene annotation table:

```bash
vcfclick annotations load
```

If `--gff` is omitted, vcfclick downloads the default GENCODE GFF3.
Use a local file when working offline or when you need a specific
release:

```bash
vcfclick annotations load --gff gencode.v45.annotation.gff3.gz
```

ClinVar lookup requires the ClinVar table:

```bash
vcfclick annotations load-clinvar
```

If `--vcf` is omitted, vcfclick downloads the current NCBI ClinVar VCF.
For reproducible demos, keep a local copy and pass it explicitly:

```bash
vcfclick annotations load-clinvar --vcf clinvar.vcf.gz
```

gnomAD population allele frequencies power `gnomad_lookup` and the
`db trio --gnomad-max-af` rarity filter. gnomAD is too large to bundle,
so you supply a sites VCF — a region slice pulls in seconds via
tabix-over-HTTPS from the public gnomAD bucket:

```bash
tabix https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr7.vcf.bgz \
    chr7:117480000-117670000 | bgzip > cftr.gnomad.vcf.gz
vcfclick annotations load-gnomad cftr.gnomad.vcf.gz
```

`load-gnomad` appends by default, so per-chromosome slices load
incrementally; pass `--replace` to start fresh.

The annotation store is separate from cohort databases. Loading
annotations does not modify any `variants` or `genotypes` table.

## Start The MCP Server

Run the server with a Python executable that can import vcfclick:

```bash
python -m vcfclick_mcp.server
```

From a source checkout:

```bash
uv run python -m vcfclick_mcp.server
```

Select the active cohort database with `VCFCLICK_DB_NAME`:

```bash
VCFCLICK_DB_NAME=demo python -m vcfclick_mcp.server
```

Use `VCFCLICK_HOME` and `VCFCLICK_BACKEND` the same way you do for the
CLI:

```bash
VCFCLICK_HOME=/data/vcfclick \
VCFCLICK_BACKEND=chdb \
VCFCLICK_DB_NAME=study \
python -m vcfclick_mcp.server
```

## Example MCP Client Configuration

The exact JSON location depends on the client. The important pieces are
the Python executable, the module name, and the environment.

```json
{
  "mcpServers": {
    "vcfclick": {
      "command": "/path/to/python",
      "args": ["-m", "vcfclick_mcp.server"],
      "env": {
        "VCFCLICK_HOME": "/Users/me/.vcfclick",
        "VCFCLICK_BACKEND": "chdb",
        "VCFCLICK_DB_NAME": "demo"
      }
    }
  }
}
```

For a source checkout managed by `uv`, use the checkout environment's
Python executable. One way to find it is:

```bash
uv run python -c "import sys; print(sys.executable)"
```

## Query Rules The MCP Server Teaches

The server's schema briefing emphasizes several rules that matter for
correct cohort queries:

- `genotypes` is sparse: only non-reference calls are stored.
- A missing genotype row means `0/0` by convention.
- Count samples with `COUNT(DISTINCT (ingest_id, sample_id))`.
- Compute allele-frequency denominators from `samples`, not from a
  join through sparse `genotypes`.
- Treat `gq` and `dp` filters carefully because public joint-call VCFs
  often leave them null.
- Use annotation tools before writing SQL when the user asks about a
  gene symbol or ClinVar significance.

See [Schema reference](SCHEMA.md) for the SQL patterns behind those
rules.

## Example Prompt

After the server is configured and `VCFCLICK_DB_NAME=demo` points to a
database, ask the MCP client:

```text
How many samples have a non-reference call in BRCA1? Show the SQL.
```

The expected plan is:

1. call `position_for_gene("BRCA1")`;
2. query `genotypes` in that coordinate range;
3. return both the SQL and the result.

## Troubleshooting

If the server starts but `run_sql` cannot find data, check:

```bash
vcfclick db list
vcfclick db info demo
echo "$VCFCLICK_HOME"
echo "$VCFCLICK_DB_NAME"
```

If gene lookup returns nothing, load annotations or verify the gene
symbol is present in the chosen GENCODE release.

If ClinVar lookup returns nothing, load ClinVar and verify the allele
uses the same normalized representation as the cohort (`chr` prefix,
1-based position, decomposed `ref`/`alt`).
