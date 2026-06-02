"""MCP server exposing the vcfclick database to LLM clients.

Backed by the embedded chDB session in `storage.db` — no server, no
port, no Docker. Composes two stores:

  - chDB:   sample data (variants, genotypes, samples, ingestions)
  - DuckDB: reference data (RefSeq genes, ClinVar)

The LLM is taught (via SCHEMA_DESCRIPTION) to translate bioinformatics
questions into a two-step plan: annotation lookups first (DuckDB),
then a parameterised SQL query against chDB. The UI shows each step.
"""

from __future__ import annotations

import json

from mcp.server.fastmcp import FastMCP

from annotations import (
    clinvar_lookup as _clinvar_lookup,
    gene_at as _gene_at,
    position_for_gene as _position_for_gene,
)
from storage import get_session


mcp = FastMCP("vcfclick")


SCHEMA_DESCRIPTION = """\
You translate bioinformatics questions into a two-step plan:
  1. Annotation lookups (gene → coordinates, variant → ClinVar) against
     the local DuckDB reference store, via the provided tools.
  2. A SQL query against the chDB sample-data store, parameterised by
     the results of step 1.

Annotations are NOT joined into chDB at storage time. If you need gene
coordinates or ClinVar significance, call the annotation tools.

──────────────────── CHDB TABLES (sample data) ────────────

EVERY sample-data row carries `ingest_id` — the identifier of the VCF
upload it came from. Rows are NOT merged across ingestions. The same
(chrom, pos, ref, alt) observed in two different VCFs is two rows
because their QUAL, FILTER, INFO, and QC origin can all differ.
Cross-ingestion correlation is a query-time concern, not storage.

variants — one row per (ingest_id, chrom, pos, ref, alt).
  Raw VCF INFO fields as typed columns (info_AC, info_AF, info_AN,
  info_DP, info_MQ, info_QD, info_FS, info_SOR, etc.) plus VCF
  mandatory fields (vcf_id, qual, filter). Multi-allelic sites are
  pre-decomposed: per-ALT INFO fields are scalars, not arrays.
  Non-reserved INFO fields live in `info_extra Map(String, String)`.

genotypes — sparse: ONLY non-reference calls are stored.
  One row per (ingest_id, chrom, pos, ref, alt, sample_id).
  Per-sample FORMAT fields: gt (Int8 alt-allele count), gq, dp,
  ad_ref/ad_alt, adf_*/adr_*, pl_ref_ref/pl_ref_alt/pl_alt_alt,
  mq, ft, ps, pq, phased. Non-reserved FORMAT fields in
  `format_extra Map(String, String)`.

samples — (ingest_id, sample_id, cohort, sex).
  Sample identity is (ingest_id, sample_id).

ingestions — catalog: (ingest_id, cohort, vcf_path, n_variants,
  n_samples, ingested_at). Query when the user asks "what's loaded?".

cohort_sizes_mv — materialised view: cohort → n_samples. Denominator
for cohort allele frequencies.

──────────────── ANNOTATION TOOLS (DuckDB, reference) ────────────

position_for_gene(symbol)   → (chrom, start_pos, end_pos)
gene_at(chrom, pos)         → overlapping gene symbols
clinvar_lookup(chrom, pos, ref, alt) → ClinVar significance

Use these BEFORE writing SQL whenever the user mentions a gene symbol,
gene region, or clinical significance. Do not invent gene_symbol or
clin_sig columns in chDB — they do not exist there.

──────────────────── CRITICAL CONVENTIONS ────────────────────

1. SPARSE TABLE: a sample absent from `genotypes` at a given
   (chrom, pos, ref, alt) is 0/0 by convention. NEVER write
   LEFT JOIN ... IS NULL. Do NOT add `AND gt != 0`.

2. DEFAULT QUALITY FILTER: every query against `genotypes` includes
       AND gq >= 20 AND dp >= 10
   unless the user explicitly overrides. Sensible defaults for germline
   calls.

   IMPORTANT: NULL silently fails the comparison, so if `gq` or `dp` is
   NULL for every row in the result set (common in phased / joint-call
   public releases like 1000 Genomes Phase 3, which ship genotype-only),
   the filter silently drops everything and returns 0.

   Always validate a "0" or suspiciously low filtered result by running
   the same query without the quality filter. If the raw count is high
   and the filtered count is 0, the dataset does not carry per-sample
   GQ/DP and you should report both numbers, explain why, and ask the
   user whether they want raw counts or a different quality gate (e.g.
   `ft = 'PASS'` if FT is populated, or `ad_alt / (ad_ref + ad_alt)` for
   allele-balance filtering).

3. COUNTING SAMPLES: always COUNT(DISTINCT (ingest_id, sample_id)).

4. INGESTION SCOPE: if the user doesn't name a specific ingest_id or
   cohort, the query spans ALL ingestions and the result should be
   labelled accordingly ("across N ingestions"). If the user names a
   cohort, JOIN samples to filter by cohort. If they name an ingest_id,
   filter directly.

5. ALLELE FREQUENCY: sum(gt) / (2 * cohort_sizes_mv.n_samples). AF is
   meaningful PER COHORT — always scope it.

6. COORDINATES: GRCh38, UCSC-style ('chr' prefix).

7. NON-RESERVED FIELDS: query Map column with
       info_extra['SOMETHING']     or     format_extra['SOMETHING']
   Returns String; CAST if numeric comparison.

──────────────────────────── EXAMPLE ────────────────────────────

User: "How many samples in cohort 'demo_brca' have a non-ref call in BRCA1?"

Plan:
  1. position_for_gene("BRCA1") → ("chr17", 43044295, 43125483)
  2. Scoped SQL: JOIN samples to filter by cohort.

SQL:
    SELECT count(DISTINCT (g.ingest_id, g.sample_id)) AS n_samples
    FROM genotypes g
    INNER JOIN samples s
        ON s.ingest_id = g.ingest_id AND s.sample_id = g.sample_id
    WHERE s.cohort = 'demo_brca'
      AND g.chrom = 'chr17'
      AND g.pos BETWEEN 43044295 AND 43125483
      AND g.gq >= 20 AND g.dp >= 10;

Always return BOTH the plan (annotation lookups performed) AND the SQL
string. The UI shows them side-by-side — the chain of reasoning is
part of the answer, not a debug trace.
"""


@mcp.tool()
def get_schema() -> str:
    """Return the briefing the LLM uses to plan queries."""
    return SCHEMA_DESCRIPTION


@mcp.tool()
def run_sql(query: str) -> dict:
    """Execute a chDB SQL query and return rows + the SQL that ran."""
    sess = get_session()
    raw = sess.query(query, "JSONCompact").bytes().decode()
    parsed = json.loads(raw)
    return {
        "sql": query,
        "columns": [m["name"] for m in parsed.get("meta", [])],
        "rows": parsed.get("data", []),
        "row_count": len(parsed.get("data", [])),
    }


@mcp.tool()
def position_for_gene(symbol: str) -> dict | None:
    """Translate an HGNC gene symbol to GRCh38 coordinates."""
    g = _position_for_gene(symbol)
    if g is None:
        return None
    return {
        "gene_symbol": g.gene_symbol,
        "chrom": g.chrom,
        "start_pos": g.start_pos,
        "end_pos": g.end_pos,
        "strand": g.strand,
    }


@mcp.tool()
def gene_at(chrom: str, pos: int) -> list[dict]:
    """Return all genes overlapping a single GRCh38 position."""
    return [
        {
            "gene_symbol": g.gene_symbol,
            "chrom": g.chrom,
            "start_pos": g.start_pos,
            "end_pos": g.end_pos,
            "strand": g.strand,
        }
        for g in _gene_at(chrom, pos)
    ]


@mcp.tool()
def clinvar_lookup(chrom: str, pos: int, ref: str, alt: str) -> dict | None:
    """ClinVar significance + review status for a specific allele."""
    r = _clinvar_lookup(chrom, pos, ref, alt)
    if r is None:
        return None
    return {
        "chrom": r.chrom,
        "pos": r.pos,
        "ref": r.ref,
        "alt": r.alt,
        "clin_sig": r.clin_sig,
        "review_status": r.review_status,
        "clinvar_id": r.clinvar_id,
        "condition": r.condition,
    }


if __name__ == "__main__":
    mcp.run()
