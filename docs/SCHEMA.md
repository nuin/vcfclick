# Schema reference

Every chDB cohort vcfclick creates uses the same four tables plus one
materialised view. This page lists every column. The authoritative
sources are the SQL files in [`schema/`](../schema/); this doc just
flattens them for SQL writers.

  | Table | Cardinality | Sort key |
  |---|---|---|
  | [`variants`](#variants) | one row per `(ingest_id, chrom, pos, ref, alt)` | `(ingest_id, chrom, pos, ref, alt)` |
  | [`genotypes`](#genotypes) | sparse — only non-reference calls | `(ingest_id, chrom, pos, ref, alt, sample_id)` |
  | [`samples`](#samples) | one row per `(ingest_id, sample_id)` | `(ingest_id, sample_id)` |
  | [`ingestions`](#ingestions) | one row per VCF upload | `ingest_id` |
  | [`cohort_sizes_mv`](#cohort_sizes_mv) | materialised view; `cohort → n_samples` | `cohort` |

All four tables use `ReplacingMergeTree` keyed by `ingested_at`. Re-
ingesting under the same `ingest_id` is idempotent — chDB dedupes on
merge. Use `SELECT … FROM <table> FINAL` to force dedup at query time
if you need an immediately-consistent count.

---

## Three conventions you have to internalise

1. **`genotypes` is sparse.** Only non-reference calls are stored.
   A sample absent from `genotypes` at `(chrom, pos, ref, alt)` is
   `0/0` by convention. **Never write `LEFT JOIN … WHERE g.gt IS NULL`**
   to find hom-ref samples — they're just not there. Derive hom-ref
   counts by subtraction:
   `hom_ref = total_samples_in_cohort − count(genotypes)`.

2. **Rows are NOT merged across ingestions.** The same
   `(chrom, pos, ref, alt)` ingested from two different VCFs is
   *two rows* on `variants` and *two rows per sample* on `genotypes`
   because their QUAL / FILTER / INFO and QC origin can all differ.
   Cross-ingestion correlation is your job at query time, not chDB's
   at storage time.

3. **`gq` / `dp` are often `NULL` on public joint-call releases**
   (1000G Phase 3, gnomAD, etc.) — those datasets ship genotype-only.
   `NULL` silently fails any comparison, so a filter like
   `WHERE g.gq >= 20 AND g.dp >= 10` will drop *every* row. If a
   suspiciously low count comes back, re-run without the filter to
   check whether `gq`/`dp` are populated at all on that cohort. The
   MCP `SCHEMA_DESCRIPTION` instructs the LLM about this trap; SQL-by-
   hand users hit it on their own.

---

## `variants`

One row per unique `(ingest_id, chrom, pos, ref, alt)`. Multi-allelic
sites must be decomposed before ingest (`bcftools norm -m -`); the
ingester errors with a helpful message otherwise.

### Identity + VCF mandatory fields

| Column | Type | Meaning |
|---|---|---|
| `ingest_id` | `LowCardinality(String)` | the `--ingest-id` passed at ingest time |
| `chrom` | `LowCardinality(String)` | UCSC-style (`chr` prefix). GRCh38. |
| `pos` | `UInt32` | 1-based |
| `ref` | `String` | |
| `alt` | `String` | single ALT per row (decomposed) |
| `vcf_id` | `Nullable(String)` | VCF `ID` column; usually rsID or `.` |
| `qual` | `Nullable(Float32)` | VCF `QUAL` column |
| `filter` | `LowCardinality(Nullable(String))` | `NULL` means `.` (unspecified), distinct from `'PASS'` |
| `ingested_at` | `DateTime` | wall-clock at ingest |

### Reserved INFO scalars (VCF 4.3)

| Column | Type | VCF field |
|---|---|---|
| `info_AC` | `Nullable(UInt32)` | `INFO/AC` (per-ALT, scalar after decomposition) |
| `info_AF` | `Nullable(Float32)` | `INFO/AF` |
| `info_AN` | `Nullable(UInt32)` | `INFO/AN` |
| `info_AD_ref` | `Nullable(UInt32)` | `INFO/AD[0]` |
| `info_AD_alt` | `Nullable(UInt32)` | `INFO/AD[1]` |
| `info_DP` | `Nullable(UInt32)` | `INFO/DP` |
| `info_MQ` | `Nullable(Float32)` | `INFO/MQ` |
| `info_MQ0` | `Nullable(UInt32)` | `INFO/MQ0` |
| `info_NS` | `Nullable(UInt32)` | `INFO/NS` |
| `info_BQ` | `Nullable(Float32)` | `INFO/BQ` |
| `info_SB` | `Nullable(Float32)` | `INFO/SB` |
| `info_END` | `Nullable(UInt32)` | `INFO/END` |
| `info_CIGAR` | `Nullable(String)` | `INFO/CIGAR` |
| `info_AA` | `Nullable(String)` | `INFO/AA` (ancestral allele) |

### Reserved INFO flags

Flags are `UInt8 DEFAULT 0`. Present → `1`, absent → `0`. Never `NULL`.

| Column | VCF field |
|---|---|
| `info_SOMATIC` | `INFO/SOMATIC` |
| `info_VALIDATED` | `INFO/VALIDATED` |
| `info_DB` | `INFO/DB` (in dbSNP) |
| `info_H2` | `INFO/H2` (HapMap2) |
| `info_H3` | `INFO/H3` (HapMap3) |
| `info_1000G` | `INFO/1000G` |

### Common GATK metrics

| Column | Type | VCF field |
|---|---|---|
| `info_QD` | `Nullable(Float32)` | `INFO/QD` |
| `info_FS` | `Nullable(Float32)` | `INFO/FS` |
| `info_SOR` | `Nullable(Float32)` | `INFO/SOR` |
| `info_MQRankSum` | `Nullable(Float32)` | `INFO/MQRankSum` |
| `info_ReadPosRankSum` | `Nullable(Float32)` | `INFO/ReadPosRankSum` |
| `info_ExcessHet` | `Nullable(Float32)` | `INFO/ExcessHet` |
| `info_InbreedingCoeff` | `Nullable(Float32)` | `INFO/InbreedingCoeff` |
| `info_MLEAC` | `Nullable(UInt32)` | `INFO/MLEAC` |
| `info_MLEAF` | `Nullable(Float32)` | `INFO/MLEAF` |
| `info_BaseQRankSum` | `Nullable(Float32)` | `INFO/BaseQRankSum` |
| `info_ClippingRankSum` | `Nullable(Float32)` | `INFO/ClippingRankSum` |

### DRAGEN-specific (germline + somatic)

| Column | Type | VCF field |
|---|---|---|
| `info_FractionInformativeReads` | `Nullable(Float32)` | `INFO/FractionInformativeReads` |
| `info_HAPCOMP` | `Nullable(UInt32)` | `INFO/HAPCOMP` |
| `info_HAPDOM` | `Nullable(Float32)` | `INFO/HAPDOM` |
| `info_DragenSnvHardQUAL` | `Nullable(Float32)` | `INFO/DragenSnvHardQUAL` |
| `info_DragenIndelHardQUAL` | `Nullable(Float32)` | `INFO/DragenIndelHardQUAL` |

### Overflow

| Column | Type | Meaning |
|---|---|---|
| `info_extra` | `Map(String, String)` | every non-routed `INFO` field, value stringified |

Read a key:

```sql
SELECT info_extra['CSQ'] FROM variants WHERE pos = 43044295;
```

Numeric Map values are stored as strings; cast at query time:

```sql
SELECT toFloat32OrNull(info_extra['AS_VQSLOD']) AS vqslod FROM variants;
```

To find out which keys are in `info_extra` for an ingested cohort,
run `vcfclick db stats <name>` — it lists the most frequent overflow
keys. To preview where a *new* VCF's fields will land, run
`vcfclick discover <vcf>`.

### Projections + indexes

- `p_range_scan` — secondary sort by `(chrom, pos)`. Used implicitly
  when a query has a region predicate (`chrom = 'chr17' AND pos
  BETWEEN ...`) without an `ingest_id` filter.

---

## `genotypes`

Sparse: only non-reference calls. One row per
`(ingest_id, chrom, pos, ref, alt, sample_id)`.

### Identity + genotype encoding

| Column | Type | Meaning |
|---|---|---|
| `ingest_id` | `LowCardinality(String)` | same as on `variants` |
| `chrom`, `pos`, `ref`, `alt` | as in `variants` | |
| `sample_id` | `LowCardinality(String)` | per-VCF sample name |
| `gt` | `Int8` | alt-allele count: `1` = het (0/1 or 1/0), `2` = hom-alt (1/1), `-1` = hemizygous or mixed-missing. **`0` and `./.` are NOT stored.** |
| `phased` | `UInt8 DEFAULT 0` | `1` if `|`-separated in the VCF, else `0` |
| `ingested_at` | `DateTime` | |

### Reserved FORMAT scalars

| Column | Type | VCF field |
|---|---|---|
| `gq` | `Nullable(UInt16)` | `FORMAT/GQ` |
| `dp` | `Nullable(UInt16)` | `FORMAT/DP` |
| `mq` | `Nullable(UInt16)` | `FORMAT/MQ` |
| `ft` | `LowCardinality(Nullable(String))` | `FORMAT/FT` (per-sample filter) |
| `ps` | `Nullable(UInt32)` | `FORMAT/PS` (phase set) |
| `pq` | `Nullable(UInt16)` | `FORMAT/PQ` (phasing quality) |

### Reserved FORMAT pairs (`Number=R`)

| Columns | Type | VCF field |
|---|---|---|
| `ad_ref`, `ad_alt` | `Nullable(UInt16)` | `FORMAT/AD` |
| `adf_ref`, `adf_alt` | `Nullable(UInt16)` | `FORMAT/ADF` |
| `adr_ref`, `adr_alt` | `Nullable(UInt16)` | `FORMAT/ADR` |

### Reserved FORMAT triples (`Number=G`, biallelic ploidy-2)

| Columns | Type | VCF field |
|---|---|---|
| `pl_ref_ref`, `pl_ref_alt`, `pl_alt_alt` | `Nullable(UInt16)` | `FORMAT/PL` |
| `gl_ref_ref`, `gl_ref_alt`, `gl_alt_alt` | `Nullable(Float32)` | `FORMAT/GL` |

### Overflow

| Column | Type | Meaning |
|---|---|---|
| `format_extra` | `Map(String, String)` | every non-routed `FORMAT` field, value stringified per sample |

### Projections + indexes

- `idx_sample` — bloom filter on `sample_id`, granularity 4.
  Cheap pre-filter when you JOIN by sample.
- `p_sample_major` — secondary sort by `(ingest_id, sample_id, chrom, pos, ref, alt)`.
  Used implicitly when a query is sample-major (`WHERE sample_id = ...`).

---

## `samples`

One row per `(ingest_id, sample_id)`. `cohort` is the orthogonal
grouping you set at ingest time; the same cohort name can span
multiple ingestions.

| Column | Type | Meaning |
|---|---|---|
| `ingest_id` | `LowCardinality(String)` | |
| `sample_id` | `LowCardinality(String)` | as read from the VCF header |
| `cohort` | `LowCardinality(String)` | the `--cohort` you passed at ingest |
| `sex` | `LowCardinality(Nullable(String))` | not auto-populated by the ingester; reserved for a future PED loader |
| `ingested_at` | `DateTime` | |

Sample identity is `(ingest_id, sample_id)`. A sample named `S1` in
ingestion `A` is NOT the same as `S1` in ingestion `B`. Cross-
ingestion alias resolution (same patient sequenced twice) is a
user-defined alias table; vcfclick does not auto-merge.

---

## `ingestions`

One row per VCF upload. Useful for management queries ("what's
loaded?") and for the MCP server.

| Column | Type | Meaning |
|---|---|---|
| `ingest_id` | `String` | |
| `cohort` | `LowCardinality(String)` | |
| `vcf_path` | `String` | absolute path passed to `db ingest` |
| `n_variants` | `UInt64 DEFAULT 0` | populated at end-of-ingest |
| `n_samples` | `UInt32 DEFAULT 0` | populated at end-of-ingest |
| `ingested_at` | `DateTime` | start-of-ingest wall-clock |

---

## `cohort_sizes_mv`

Materialised view (`SummingMergeTree`) giving each cohort's sample
count. Use it as the denominator for cohort allele frequency:

```sql
SELECT
    sum(g.gt) AS ac,
    (SELECT 2 * n_samples FROM cohort_sizes_mv WHERE cohort = 'study1') AS an,
    sum(g.gt) / (SELECT 2 * n_samples FROM cohort_sizes_mv WHERE cohort = 'study1') AS af
FROM genotypes g
INNER JOIN samples s
    ON s.ingest_id = g.ingest_id AND s.sample_id = g.sample_id
WHERE s.cohort = 'study1'
  AND g.chrom = 'chr17' AND g.pos BETWEEN 43044295 AND 43170245;
```

`vcfclick db diff <db> --cohort-a A --cohort-b B` does this for you
across two cohorts.

| Column | Type | Meaning |
|---|---|---|
| `cohort` | `LowCardinality(String)` | |
| `n_samples` | `UInt64` | distinct `(ingest_id, sample_id)` count |

---

## Promoting an overflow field to typed

If a field your VCFs always carry shows up under `info_extra` or
`format_extra` and you'd rather query it as a typed column:

1. Add the column to the relevant SQL file under [`schema/`](../schema/).
2. Add a routing entry in [`ingest/routing.py`](../ingest/routing.py)
   — `INFO_SCALAR`, `INFO_PAIR`, `INFO_FLAG`, `FORMAT_SCALAR`,
   `FORMAT_PAIR`, or `FORMAT_TRIPLE`.
3. Add the matching `pyarrow` type entry in
   [`ingest/_arrow.py`](../ingest/_arrow.py).
4. `vcfclick db rm <name>` and `vcfclick db create <name>` to apply
   the new schema, then re-ingest. The DRAGEN columns shipped in
   v0.1.2 followed this exact recipe.

`vcfclick discover <vcf>` emits a draft DDL fragment for any
overflow field it sees, which you can paste into step 1 unchanged.

---

## See also

- [`schema/01_variants.sql`](../schema/01_variants.sql),
  [`schema/02_genotypes.sql`](../schema/02_genotypes.sql),
  [`schema/03_samples.sql`](../schema/03_samples.sql) — canonical DDL.
- [`ingest/routing.py`](../ingest/routing.py) — the routing tables
  that determine which VCF fields land where.
- [`vcfclick_mcp/server.py`](../vcfclick_mcp/server.py) — the
  `SCHEMA_DESCRIPTION` briefing the MCP server hands to an LLM
  client; states the conventions above in prompt form.
- [`examples/brca1-cohort.md`](../examples/brca1-cohort.md) — five
  worked queries against the demo bundle.
