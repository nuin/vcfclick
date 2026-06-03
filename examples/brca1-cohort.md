# BRCA1 cohort questions — vcfclick MCP session

Five canonical questions a researcher might ask of the 1000 Genomes
Phase 3 BRCA1 cohort, with the MCP tool calls, the SQL the LLM
generates, and the actual result.

**Dataset:** 3,014 variants × 3,202 samples × 588,584 sparse non-ref
calls, scoped to BRCA1 (chr17:43,044,295-43,170,245).

## Setup

```bash
# Install + load reference annotations (one-time, ~5 min total).
pip install vcfclick
vcfclick annotations load            # GENCODE genes (~60 MB)
vcfclick annotations load-clinvar    # ClinVar weekly VCF (~80 MB)

# Pull the cohort demo bundle.
vcfclick db pull demo \
    https://github.com/nuin/vcfclick/releases/download/v0.1.0/1000g-brca1-demo.tar.gz

# Point Claude Desktop at it. In ~/Library/Application Support/Claude/claude_desktop_config.json:
#
#   "vcfclick": {
#     "command": "/path/to/vcfclick/.venv/bin/python",
#     "args": ["-m", "vcfclick_mcp.server"],
#     "cwd": "/path/to/vcfclick",
#     "env": {
#       "PYTHONPATH": "/path/to/vcfclick",
#       "VCFCLICK_DB_NAME": "demo"
#     }
#   }
```

The five questions below are then literally typed into Claude Desktop
(or any MCP-capable client). The SQL and outputs shown are what the
client renders.

---

## Q1: How many samples have a non-reference call anywhere in BRCA1?

**LLM workflow:**

1. `position_for_gene("BRCA1")` → `chr17:43,044,295-43,170,245`
2. `run_sql(...)` with the range.

**Generated SQL:**

```sql
SELECT count(DISTINCT (ingest_id, sample_id)) AS samples_with_brca1_call
FROM genotypes
WHERE chrom = 'chr17'
  AND pos BETWEEN 43044295 AND 43170245
```

**Result:**

```
   ┌─samples_with_brca1_call─┐
1. │                    3202 │
   └─────────────────────────┘
```

Every sample in the cohort has at least one non-reference call in the
BRCA1 region — unsurprising for an 84 Kb gene at 30× joint-called
depth, but the question gets answered in milliseconds.

---

## Q2: What are the five most common BRCA1 variants in this cohort?

**LLM workflow:**

1. `position_for_gene("BRCA1")` (cached from Q1).
2. `run_sql(...)` ordering by `info_AF` — note that `info_AF` is a
   typed column on `variants`, not a Map key, because AF is a VCF 4.3
   reserved field.

**Generated SQL:**

```sql
SELECT pos, ref, alt, info_AC, info_AF, info_AN
FROM variants
WHERE chrom = 'chr17'
  AND pos BETWEEN 43044295 AND 43170245
ORDER BY info_AF DESC NULLS LAST
LIMIT 5
```

**Result:**

```
   ┌──────pos─┬─ref─┬─alt─┬─info_AC─┬──info_AF─┬─info_AN─┐
1. │ 43165485 │ G   │ GT  │    6384 │ 0.996877 │    6404 │
2. │ 43161270 │ A   │ G   │    6379 │ 0.996096 │    6404 │
3. │ 43080327 │ G   │ C   │    6273 │ 0.979544 │    6404 │
4. │ 43169893 │ C   │ T   │    5615 │ 0.876796 │    6404 │
5. │ 43159815 │ G   │ C   │    5404 │ 0.843848 │    6404 │
   └──────────┴─────┴─────┴─────────┴──────────┴─────────┘
```

The top variant has AF ≈ 0.997 — almost every chromosome in the cohort
carries it. (`info_AN = 6404` because 3202 diploid samples have 6404
alleles.)

---

## Q3: At that top variant, how many samples are het vs. hom-alt vs. hom-ref?

**Note on the sparse encoding:** the `genotypes` table only stores
non-reference calls (`gt = 1` for het, `gt = 2` for hom-alt). Hom-ref
calls (`0/0`) are *absent rows*, not zero-valued rows, so the
hom-ref count is derived as `(total samples) − (stored rows for this
variant)`. The MCP schema description explicitly tells the LLM this.

**Generated SQL:**

```sql
SELECT
    sum(gt = 1) AS het_samples,
    sum(gt = 2) AS hom_alt_samples,
    3202 - count() AS hom_ref_samples
FROM genotypes
WHERE chrom = 'chr17' AND pos = 43165485 AND ref = 'G' AND alt = 'GT'
```

**Result:**

```
   ┌─het_samples─┬─hom_alt_samples─┬─hom_ref_samples─┐
1. │          20 │            3182 │               0 │
   └─────────────┴─────────────────┴─────────────────┘
```

20 het + 3,182 hom-alt + 0 hom-ref = 3,202 samples — consistent with
AF ≈ 0.997.

---

## Q4: Show me low-quality variants — GATK QD below 2 — in the BRCA1 region

**Generated SQL:** No annotation tools needed; this is a pure
typed-column query. `info_QD` is one of the common GATK fields the
ingester routes to a typed column.

```sql
SELECT pos, ref, alt, info_QD, info_FS
FROM variants
WHERE chrom = 'chr17'
  AND pos BETWEEN 43044295 AND 43170245
  AND info_QD IS NOT NULL
  AND info_QD < 2
ORDER BY info_QD
LIMIT 5
```

**Result:**

```
   ┌──────pos─┬─ref─┬─alt─┬─info_QD─┬─info_FS─┐
1. │ 43098484 │ A   │ G   │    0.25 │   4.194 │
2. │ 43057493 │ A   │ T   │     0.3 │  23.612 │
3. │ 43057489 │ T   │ C   │    0.43 │  24.398 │
4. │ 43128393 │ G   │ C   │    0.64 │   5.542 │
5. │ 43128378 │ GGT │ G   │    1.08 │       0 │
   └──────────┴─────┴─────┴─────────┴─────────┘
```

The `info_QD` column was readable directly because QD is in the
ingester's typed list. Custom fields (CSQ, ClinVar annotations baked
into the VCF, etc.) would live under `info_extra['<KEY>']` instead —
the discover output (`vcfclick discover input.vcf.gz`) tells you which
is which before you query.

---

## Q5: Per-sub-cohort allele count for that top variant

**Generated SQL:** Joins `genotypes` to `samples` on
`(ingest_id, sample_id)` so per-cohort grouping is available — cohort
is a sample-table column the user sets at ingest time.

```sql
SELECT
    s.cohort,
    sum(g.gt) AS allele_count_in_cohort,
    count(DISTINCT s.sample_id) * 2 AS allele_total
FROM genotypes g
JOIN samples s
  ON s.ingest_id = g.ingest_id AND s.sample_id = g.sample_id
WHERE g.chrom = 'chr17'
  AND g.pos = 43165485
  AND g.ref = 'G' AND g.alt = 'GT'
GROUP BY s.cohort
```

**Result:**

```
   ┌─cohort─┬─allele_count_in_cohort─┬─allele_total─┐
1. │ 1000g  │                   6384 │         6404 │
   └────────┴────────────────────────┴──────────────┘
```

The demo bundle is a single cohort (`1000g`), so the grouping
collapses to one row. With multiple ingested cohorts side-by-side
under the same database, this same query yields per-cohort AF and is
the basis for ad-hoc case/control comparisons.

---

## Composing with ClinVar

After `vcfclick annotations load-clinvar`, the MCP server's
`clinvar_lookup` tool returns real pathogenicity calls. A
multi-tool question like

> "Which BRCA1 variants in this cohort match a ClinVar-pathogenic call?"

then becomes: `position_for_gene("BRCA1")` → `run_sql(...)` to pull
the cohort's BRCA1 variants → `clinvar_lookup(chrom, pos, ref, alt)`
for each → filter where `clin_sig LIKE '%athogenic%'`. The exact
count depends on the ClinVar weekly release at the time of the
query.

This composition — chDB sample data plus DuckDB reference data,
joined by the MCP server at question time — is the architectural
point of separating the two stores. Reference data refreshes on a
ClinVar/GENCODE cadence; sample data grows independently as new
cohorts arrive.
