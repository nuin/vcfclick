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
    gnomad_af as _gnomad_af,
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

──────────────── ANNOTATION TOOLS (DuckDB, reference) ────────────

position_for_gene(symbol)   → (chrom, start_pos, end_pos)
gene_at(chrom, pos)         → overlapping gene symbols
clinvar_lookup(chrom, pos, ref, alt) → ClinVar significance
gnomad_lookup(chrom, pos, ref, alt)  → gnomAD allele frequency (popmax)

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

5. ALLELE FREQUENCY:
       ac = sum(gt) over `genotypes` JOINed to cohort samples.
       an = 2 * count(DISTINCT (ingest_id, sample_id)) computed
            FROM `samples` ALONE, filtered by cohort — NOT from
            the join with genotypes. Because genotypes is sparse
            (0/0 calls are absent), counting samples through the
            join only sees non-reference samples and gives a too-
            small denominator → inflated AF. Compute the cohort
            size as its own CTE / subquery against `samples` and
            bring it in via CROSS JOIN.

   Canonical pattern (this is what `vcfclick db diff` does):
       WITH cohort_size AS (
           SELECT 2 * count(DISTINCT (ingest_id, sample_id)) AS an
           FROM samples WHERE cohort = 'study1'
       )
       SELECT
           sum(g.gt) AS ac,
           cs.an     AS an,
           sum(g.gt) / cs.an AS af
       FROM genotypes g
       INNER JOIN samples s
           ON s.ingest_id = g.ingest_id AND s.sample_id = g.sample_id
       CROSS JOIN cohort_size cs
       WHERE s.cohort = 'study1'
         AND g.chrom = 'chr17' AND g.pos = 43044295
       GROUP BY cs.an;

   AF is meaningful PER COHORT — always scope it. There is
   intentionally no materialized cohort_sizes view — SummingMergeTree
   wouldn't decrement on rollback or replacement.

6. COORDINATES: GRCh38, UCSC-style ('chr' prefix).

7. NON-RESERVED FIELDS: query Map column with
       info_extra['SOMETHING']     or     format_extra['SOMETHING']
   Returns String; CAST if numeric comparison.

8. FAMILY / TRIO ANALYSIS: the `pedigree` table maps relationships:
       (ingest_id, sample_id, family_id, father_id, mother_id, sex, affected)
   father_id/mother_id are sample_ids ('0' = founder). Resolve a
   proband's parents from pedigree, then self-join genotypes on
   (ingest_id, chrom, pos, ref, alt) aliasing the proband (g), father
   (f), mother (m). Genotype encoding: gt 0=hom-ref, 1=het, 2=hom-alt.

   Inheritance models (gt comes from `genotypes`):
     - recessive:  g.gt=2 AND f.gt=1 AND m.gt=1
     - dominant:   g.gt=1 AND exactly one parent carries (gt>0), other gt=0
     - de novo:    g.gt>0 AND f.gt=0 AND m.gt=0

   CRITICAL de-novo caveat: genotypes is sparse, so a parent ABSENT at a
   site is 0/0 OR ./. (no-call) — indistinguishable. Confident de novo
   needs the parents to have a STORED gt=0 row (requires the cohort was
   ingested with --keep-reference). An INNER JOIN to both parents
   requiring f.gt=0 AND m.gt=0 naturally excludes no-call parents (no
   row → no join). If the DB has no gt=0 rows, de novo cannot be proven;
   say so rather than reporting absence-based de novo.

   Standard candidate filters (slivar-style): GQ>=20, DP>=10, het
   allele-balance ad_alt/(ad_ref+ad_alt) near 0.5, and population
   rarity info_AF<=0.01. The `vcfclick db trio` CLI command applies
   exactly these.

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


@mcp.tool()
def gnomad_lookup(chrom: str, pos: int, ref: str, alt: str) -> dict | None:
    """gnomAD allele frequency for a specific allele. `popmax` is the
    rarity-relevant value (highest genetic-ancestry-group AF). None means
    the allele is absent from the loaded gnomAD slice — treat as rare, not
    AF 0."""
    r = _gnomad_af(chrom, pos, ref, alt)
    if r is None:
        return None
    return {
        "chrom": r.chrom,
        "pos": r.pos,
        "ref": r.ref,
        "alt": r.alt,
        "af": r.af,
        "af_grpmax": r.af_grpmax,
        "popmax": r.popmax,
    }


if __name__ == "__main__":
    mcp.run()
