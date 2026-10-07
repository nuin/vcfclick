# `vcfnorm`: in-process multi-allelic splitting and consequence decoding — design

Status: draft for review (branch `feature/vcfnorm`).

## Goal

Add a top-level `vcfnorm/` package with two opt-in features:

1. `--split-multiallelic` at ingest: decompose multi-allelic records in
   process, producing exactly the rows that `bcftools norm -m -` followed by
   a normal ingest would produce.
2. A new `consequences` table holding VEP `CSQ` / SnpEff `ANN` annotations
   decoded into typed columns. It is filled at ingest (`--decode-consequences`)
   or afterwards (`vcfclick db consequences`).

`vcfnorm` imports only cyvcf2, numpy, pyarrow and (feature 2) polars. It never
imports `storage`, `cli`, `ingest`, chDB or DuckDB, so it can later be
extracted into its own package. The glue that stages rows and talks to the
backend lives in `ingest/` and `cli/`.

Non-goals:

- Left-alignment or allele trimming. Both need a reference FASTA. Users who
  need them still run `bcftools norm -f ref.fa`. The docs say so.
- The other bcftools splitting modes: `--multi-overlaps .`, `--keep-sum`,
  `--old-rec-tag`, and `-m -snps` / `-m -indels`. Only plain `-m -` (that is,
  `-m -both`) is mirrored.
- Joining biallelic records (`-m +`).
- bcftools `BCSQ`, symbolic structural-variant alleles (`<DEL>`, …), and
  `Number=G` for ploidy above 2. A record that needs any of these is rejected
  with a clear error.
- Changing existing table semantics. `genotypes` stays sparse, and rows are
  never merged across ingestions. `info_extra` is left exactly as it is
  today: the raw `CSQ` / `ANN` string stays there.
- A new query engine, polars-bio, or Rust.

## 1. Multi-allelic splitting

### Reference behaviour (bcftools 1.24, verified)

Given a record with ALT alleles `a1..an`, `bcftools norm -m -` emits one
record per ALT, in ALT order. For the record of allele `k`:

| field | rule |
|---|---|
| CHROM, POS, ID, REF, QUAL, FILTER | copied (the ID too: `rs1` appears on every split record) |
| ALT | `ak` only. No trimming: `ATG>C` stays `ATG>C`, and `ATG>ATGTG` is not shortened |
| `Number=A` (INFO, FORMAT) | value `k` |
| `Number=R` | values `0` and `k` |
| `Number=G` (diploid) | values at genotype indices `(0,0)`, `(0,k)`, `(k,k)`, using the VCF index `j(j+1)/2 + i` |
| `Number=G` (haploid sample) | values `0` and `k` |
| `Number=1`, fixed `n`, `.`, Flag | copied unchanged. This includes `CSQ` / `ANN`, which get the whole original string on every split record |
| GT | allele `k` → `1`; any other ALT → `0` (`--multi-overlaps 0`, the default); missing `.` and phase are kept. `1/2` → `1/0` and `0/1`; `./2` → `./0` and `./1`; haploid `2` → `0` and `1` |
| missing per-sample value | stays missing. A `.` for an `R` field becomes `.,.`, the way bcftools prints it |

The fields this covers include the typed columns `info_AC`, `info_AF`,
`info_MLEAC`, `info_MLEAF`, `info_AD_ref/alt`, `ad_*`, `adf_*`, `adr_*`,
`pl_*` and `gl_*`, and every `info_extra` / `format_extra` value, all driven
by each field's `Number` in the header. A value whose count does not match its
`Number` raises a `ValueError` naming the site, as other parse errors already
do. That happens before any database write, so earlier data under the same
`ingest_id` is kept.

`gt`, `n_called`, `an_called` and `ac_called` are then computed from the
recoded allele array by the existing code, so they are exact per allele.

### Mechanism: a record adapter, not a temporary VCF

`build_record_rows` (shared by the serial and parallel loaders) reads a small
part of the cyvcf2 `Variant` API: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`,
`FILTER`, `INFO` (mapping), `format(field)`, `gt_types`, `gt_phases` and
`genotype.array()`. `vcfnorm.split` provides `SplitRecord`, which implements
exactly that surface. Its values are sliced out of the parent record's numpy
arrays and INFO values using the rules above. `gt_types` is re-derived with
cyvcf2's rule (missing alleles ignored; none called → UNKNOWN; all called
alleles `0` → HOM_REF; all called alleles ALT → HOM_ALT; a mix → HET, so `./1`
counts as HET). The loaders then call the unchanged `build_record_rows` on it.

Approaches considered:

- **Re-parse split text through cyvcf2** (`Writer.variant_from_string`).
  Rejected: records made that way have no reader attached, so
  `Variant.format()` raises. Verified with cyvcf2 0.33.0.
- **Write a normalised temporary VCF and ingest that.** Rejected: it doubles
  the parse work and disk I/O, and the parallel path would need the temporary
  file bgzipped and tabix-indexed. cyvcf2 cannot write that index.
- **Record adapter (chosen).** Biallelic records pass through as the real
  cyvcf2 objects, so the default path is unchanged and no slower. Only
  multi-allelic records pay for slicing. The parallel workers split inside
  their own region.

`SplitRecord` raises `AttributeError` for anything outside that API surface.
A contract test runs `build_record_rows` on it, so if `vcf_rows` starts using
a new attribute, the test fails instead of the code silently misbehaving.

### Wiring

- `ingest/records.py` (new) holds one iterator, `decomposed(vcf_iter, split,
  fields, vcf_path)`. With `split` off, it raises today's error for a
  multi-allelic record; the message now also names `--split-multiallelic`.
  With `split` on, it yields `vcfnorm` splits. Both `vcf_load._stage_vcf` and
  `parallel_worker._worker` use this iterator in place of their own
  `len(variant.ALT) != 1` checks, so both loaders behave identically.
- The header's `Number` map is built once per VCF by
  `vcfnorm.header.field_numbers(vcf)`.
- `n_variants` in `ingestions` counts split records, matching what the
  bcftools-normalised file reports.

### CLI

`--split-multiallelic` (default off) on `db ingest` (serial and parallel
paths) and `db ingest-batch`. It is threaded through `IngestOptions`,
`ingest()`, `ingest_parallel()` and the batch path, the same way
`keep_reference` / `record_missing` are today.

## 2. Consequence decoding

### Table `consequences`

One row per annotation entry. It joins to `variants` on `(ingest_id, chrom,
pos, ref, alt)`. ClickHouse uses `ReplacingMergeTree(ingested_at)`; DuckDB
uses the same columns in the same order.

| column | ClickHouse type | source (VEP / SnpEff) |
|---|---|---|
| `ingest_id` | LowCardinality(String) | |
| `chrom` | LowCardinality(String) | |
| `pos` | UInt32 | |
| `ref` | String | |
| `alt` | String | the variant's ALT this entry was matched to |
| `source` | LowCardinality(String) | `'vep'` / `'snpeff'` |
| `entry` | UInt16 | 0-based position of the entry in the record's annotation string |
| `allele` | String | `Allele` / `Allele`, raw |
| `consequence` | Array(String) | `Consequence` / `Annotation`, split on `&` |
| `impact` | LowCardinality(Nullable(String)) | `IMPACT` / `Annotation_Impact` |
| `gene_symbol` | Nullable(String) | `SYMBOL` / `Gene_Name` |
| `gene_id` | Nullable(String) | `Gene` / `Gene_ID` |
| `feature_type` | LowCardinality(Nullable(String)) | `Feature_type` / `Feature_Type` |
| `transcript_id` | Nullable(String) | `Feature` / `Feature_ID` |
| `biotype` | LowCardinality(Nullable(String)) | `BIOTYPE` / `Transcript_BioType` |
| `hgvsc` | Nullable(String) | `HGVSc` / `HGVS.c` |
| `hgvsp` | Nullable(String) | `HGVSp` / `HGVS.p` |
| `exon` | Nullable(String) | `EXON` (`"3/12"`) / `Rank`, when no consequence term is `intron_variant` |
| `intron` | Nullable(String) | `INTRON` / `Rank`, when a consequence term is `intron_variant` |
| `canonical` | Nullable(UInt8) | `CANONICAL` (`YES` → 1, otherwise 0) / NULL |
| `mane_select` | Nullable(String) | `MANE_SELECT`, or the older `MANE` field / NULL |
| `csq_extra` | Map(String, String) | every other non-empty subfield, keyed by its header name |
| `ingested_at` | DateTime DEFAULT now() | housekeeping |

`ORDER BY (ingest_id, chrom, pos, ref, alt, source, entry)`.

This differs from the brief, which keyed rows on `transcript`.
`transcript_id` is a column, not part of the key, for two reasons. Intergenic
entries have no feature ID, and ClickHouse sort keys cannot be Nullable.
Also, VEP can emit two entries for the same feature (for example `--merged`
runs, or several alleles at one site). In a ReplacingMergeTree, rows with the
same key collapse into one, which would make chDB and DuckDB disagree.
`source` + `entry` is unique per record, so both backends store the same rows.

Value handling:

- An empty subfield is NULL, and is left out of `csq_extra`.
- VEP subfields are percent-decoded (any `%XX`; for example the `%3D` that VEP
  writes for `=` in synonymous HGVSp). SnpEff subfields are kept as written.
- Only `consequence` is split into a list. Other `&`-joined values (such as
  `Existing_variation`) stay as raw strings in `csq_extra`.

Array queries differ between the backends: `has(consequence, 'x')` in
ClickHouse and `list_contains(consequence, 'x')` in DuckDB. SCHEMA.md and the
MCP `SCHEMA_DESCRIPTION` give both. Parquet import of `list<string>` into
`Array(String)` / `VARCHAR[]` was verified on both backends.

### Field layout

- **VEP:** field names and order come from the `Format: A|B|…` part of the
  `##INFO=<ID=CSQ…>` description. If there is no `Format:`, decoding that
  source raises an error naming the header line.
- **SnpEff `ANN`:** the fixed 16-field spec (`Allele | Annotation |
  Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID |
  Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length |
  CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS /
  INFO`). Any trailing extra subfields go to `csq_extra` as `field_N`.
- If a record carries both `CSQ` and `ANN`, both are decoded, and the `source`
  column tells them apart.

### Allele matching

A split record (whether split by bcftools or by `vcfnorm`) carries the
original record's whole annotation string. Each entry must therefore be
assigned to the right ALT, or dropped.

The **allele group** is the set of variants in one ingestion that share
`(chrom, pos, ref)` and have an identical raw annotation string. Those are
exactly the siblings split from one original record; a record that was never
split forms a group of one. The group's ALT set rebuilds the original
record's alleles. This works the same at ingest and in the backfill, and for
input split by either bcftools or vcfnorm. The ALT order is lost, but VEP's
allele representation only depends on the set of alleles.

Each entry's `Allele` is matched against these candidates, in order:

1. **VEP default representation over the group.** If REF and every ALT in the
   group share their first base, that base is trimmed from every allele, and
   an empty result becomes `-`. Otherwise all alleles are left as they are.
   Within a group this mapping is one-to-one, so a match is never ambiguous.
2. **Per-allele minimal representation** (VEP `--minimal`): shared prefix,
   then shared suffix, trimmed from REF/ALT separately; empty becomes `-`.
3. **The literal ALT**, which is what SnpEff uses.

For SnpEff only step 3 applies. An entry that matches no ALT in the group
(for example an allele that was filtered out before ingest) is dropped and
counted. The command reports `entries decoded / dropped`. `ALLELE_NUM`, when
present, is kept in `csq_extra` but not used for matching, so that ingest-time
decoding and the backfill always give the same result.

### Polars

`vcfnorm.consequences.decode(table: pa.Table, fmt) -> pa.Table`. The input
has the variant key columns plus the raw `CSQ` / `ANN` string; the output
matches `CONSEQUENCES_ARROW_SCHEMA`. Inside, polars does: split on `,` →
explode, keeping the entry index → split on `|` → struct, using the header
field names → join to candidate allele representations computed per group →
build the typed columns and `csq_extra`. This is the only place polars is
used.

`vcfnorm` defines the annotation part of the schema.
`ingest/_arrow.py` combines it with the key columns to form
`CONSEQUENCES_ARROW_SCHEMA`, registers it in `TABLE_COLUMNS`, and the
schema-agreement test covers it.

**Dependency: an optional extra, `vcfclick[csq]` (`polars>=1.20`).** Only
feature 2 needs polars, and the polars wheel is about 40 MB, so the base
install stays lean. Running `--decode-consequences` or `db consequences`
without polars fails before parsing, with `uv tool install
"vcfclick[csq]"` in the message. CI installs the extra in both the chDB and
DuckDB jobs. Bioconda: polars is on conda-forge (2.0.0 as of today), so a
recipe can add it as a run dependency. The in-repo `packaging/bioconda`
recipe is not changed by this work.

### Entry points

- **At ingest:** `--decode-consequences` (default off) on `db ingest`
  (serial and parallel) and `db ingest-batch`. Decoded rows are staged as
  `c_*.parquet` next to the variant batches and imported in the same commit,
  so a failure rolls them back too. Siblings of one record are consecutive in
  the file, so the loaders put off a batch flush until the current allele
  group is closed. That way a group is never split across two Parquet batches.
  If the VCF has no `CSQ` / `ANN` header, the loader logs once and writes no
  rows.
- **Backfill:** `vcfclick db consequences <db> [--ingest-id ID]
  [--csq-format 'Allele|Consequence|…']`. It reads `info_extra['CSQ']` /
  `info_extra['ANN']` for each ingestion (or just one), groups and decodes
  them, then deletes and re-inserts that ingestion's `consequences` rows. Run
  twice, it gives the same result. The VEP field order is taken from the
  header of `ingestions.vcf_path` if that file still exists; otherwise from
  `--csq-format`; otherwise the command errors, saying which ingestion needs
  it. It reports variants scanned, entries decoded and entries dropped per
  ingestion.

### Plumbing for the new table

- `schema/07_consequences.sql` and `schema/duckdb/07_consequences.sql`.
- `storage.db`: add the table to `_ADDED_TABLE_FILES` (`upgrade_schema`
  creates it in older databases) and to `_INGESTION_SCOPED_TABLES` plus
  `_OPTIONAL_INGESTION_TABLES` (rollback clears it, and skips it if absent).
- `export.parquet.TABLES`, `cli.db_bundle.BUNDLE_TABLES`,
  `ingest.parquet_load` (`_EXPECTED_SCHEMAS`, and an optional
  `consequences.parquet` in `DumpFiles`), and the web UI `_TABLES` list.
  Older bundles without the file still load.
- MCP `SCHEMA_DESCRIPTION`: the table, the join key, the array-query
  difference between backends, and "query `consequences`, not
  `info_extra['CSQ']`".

## Edge cases

- Multi-allelic records with `--decode-consequences` but without
  `--split-multiallelic` are still rejected, exactly as today.
- A record that is both split and annotated: each split record gets only the
  entries for its own allele (tested).
- Indels with a shared first base (`-` deletions, trimmed insertions), MNVs
  (`AT>GC`, nothing trimmed), mixed groups (`ATG>A,ATGTG,C`, where VEP does
  not trim because the first bases differ), and a `*` ALT, which matches a
  `*` entry and is otherwise dropped.
- Haploid and mixed-ploidy samples (`Number=G` uses the haploid rule per
  sample), phased genotypes, `./.`, and partial calls.
- A `CSQ` value that is `.` or empty gives no rows. A malformed entry (wrong
  subfield count) is dropped and counted, never fatal.
- Ingesting the same `ingest_id` again replaces its consequences (rollback
  covers the table).
- Older databases gain the empty table on first use. Older bundles load
  without it.

## Test plan

Fixtures (small plain-text VCFs, bgzipped at test time with
`conftest.bgzip_vcf`):

- `multiallelic_fields.vcf`: A/R/G INFO and FORMAT fields, typed and
  overflow, haploid and diploid, phased, partial and missing calls,
  `1/2`-style het-alt calls, an indel and an MNV group, and a `*` allele.
  The existing `multiallelic.vcf.gz` stays as it is, because
  `test_ingest_atomic` relies on its exact 3-record shape for the rejection
  test.
- `multiallelic_fields.split.vcf`: the committed output of `bcftools 1.24
  norm -m -` on the file above. This is the oracle, so CI does not need
  bcftools.
- `vep.vcf` (with a `Format:` header, a canonical and a MANE transcript,
  intergenic and regulatory entries, `%3D` escapes, indels, a deletion, an
  MNV and a multi-allelic record) and `snpeff.vcf` (the same shapes, ANN
  format).

Tests:

- `vcfnorm` unit tests, with no storage: `Number` slicing for A, R and G
  (diploid and haploid); GT recoding and `gt_types` derivation; VEP, minimal
  and SnpEff allele representations; percent-decoding; the `CSQ` header
  parser; and the `SplitRecord` contract.
- **Oracle tests, on both backends and both serial and parallel:** ingest
  `multiallelic_fields.vcf` with `--split-multiallelic`, and ingest the
  committed `.split.vcf` normally. `variants`, `genotypes` and
  `missing_genotypes` must be identical row for row (excluding `ingest_id`
  and `ingested_at`). When bcftools is on PATH, a further test regenerates the
  split file and checks it still matches the committed one, which catches
  changes in bcftools behaviour.
- The default ingest of a multi-allelic file still fails, and the message
  names the flag.
- Consequences: expected rows for the VEP and SnpEff fixtures; ingest-time
  decoding equals backfill decoding; vcfnorm-split input equals
  bcftools-split input; the backfill is idempotent; `--csq-format` is used
  when the VCF is gone; a missing extra gives the install-hint error.
- Round trips: `db dump` → `ingest-parquet`, `db push` → `db pull`, an older
  bundle without `consequences.parquet`, `upgrade_schema` on an older
  database, and rollback on a failed ingest.
- The schema-agreement test covers the new table; the MCP description test
  mentions it.
- ruff, isort and the full suite pass with `VCFCLICK_BACKEND=chdb` and
  `=duckdb`.

## Docs

- `docs/VCFNORM.md` (new): splitting (the bcftools-equivalence rules, the GT
  convention, and what it doesn't do: left-alignment needs a FASTA) and
  consequences (allele matching, backfill, `--csq-format`).
- `docs/CLI.md`, `docs/SCHEMA.md` (the new table, query patterns, and
  `has` / `list_contains`), the README feature list and docs index, the MCP
  `SCHEMA_DESCRIPTION`, and CHANGELOG (Unreleased). No version bump.

## Delivery

Feature 1 lands first (vcfnorm split, the flag, oracle tests, docs), then
feature 2 (table, decoding, backfill, docs). Each has its own tests and docs,
on `feature/vcfnorm`.
