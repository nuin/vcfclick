# Population-genetics support — design

Status: approved for implementation (branch `feature/popgen`).

## Goal

Make vcfclick databases usable for population-genetics analysis, both from
its own CLI and as Parquet bundles read in the browser by PopGeneJS
(DuckDB-Wasm). Four pieces:

1. Called-genotype accounting (missing data is currently indistinguishable
   from 0/0 in the default sparse layout).
2. A sample → population panel.
3. Ancestral-allele handling for polarised statistics (`info_AA` already
   exists as a typed column; it needs normalisation).
4. A `vcfclick db popgen` command group computing standard statistics.

Non-goals: changing the meaning of the existing `genotypes` table (every
row stays a called non-reference genotype, or a 0/0 under
`--keep-reference`); new trio/QC behaviour; a migration framework.

## 1. Called-genotype accounting

### 1a. Per-site cohort counts on `variants`

Add three columns (both backends, Arrow schema, schema-agreement test):

| column | type | meaning |
|---|---|---|
| `n_called` | Nullable(UInt32) | samples with a fully called GT at this site |
| `an_called` | Nullable(UInt32) | called alleles at this site (ploidy-aware: haploid call = 1, diploid = 2; partially missing `./1` counts its called allele) |
| `ac_called` | Nullable(UInt32) | ALT alleles among called alleles, counted from GT (a haploid `1` contributes 1, not 2) |

Computed at ingest from cyvcf2 per-sample allele arrays (`variant.genotypes`,
allele index −1 = missing), not from `gt_types`, so ploidy and partial
missingness are exact. Decomposed (single-ALT) records only, as today.
NULL means "not recorded" (databases or bundles created before this change).

### 1b. Per-sample missing calls: new table `missing_genotypes`

```
ingest_id, chrom, pos, ref, alt, sample_id   (+ ingested_at housekeeping)
```

One row per sample whose GT is fully missing (`./.`, `.`) at a site.
Written by default; `db ingest --no-record-missing` skips it. Same engine /
ORDER BY conventions as `genotypes`. This keeps `genotypes` semantics intact
while letting per-population called counts be derived:
`called_in_pop = pop_size − missing_in_pop`.

Partially missing diploid calls (`./1`): store the called allele in
`genotypes` as today's encoding does (whatever cyvcf2 gt_types yields) and
document the edge case; per-site `an_called` / `ac_called` stay exact.

All ingest paths must populate both: serial `vcf_load`, parallel ingest,
`ingest-batch`, `merge`/`combine` outputs where they write tables, and
Parquet bundle dump/load (`dump`, `ingest-parquet`, `db push`/pull bundles).
Loading an older bundle without these columns/table must still work (NULL /
absent).

## 2. Population panel

New table `populations` (both backends):

```
ingest_id, sample_id, population, super_population (nullable), sex (nullable)
```

Keyed like `pedigree` (`(ingest_id, sample_id)`), ReplacingMergeTree.

CLI: `vcfclick db panel <db> <file> [--ingest-id ID] [--sample-col] [--pop-col]
[--super-pop-col] [--sex-col]`

- Accepts the 1000 Genomes panel format (`sample pop super_pop gender`,
  tab-separated, header) out of the box, plus generic TSV/CSV with column
  options. `gender`/`sex` values normalised to male/female/NULL.
- Without `--ingest-id`, applies to every ingestion that contains the sample.
- Reports panel samples not in the DB and DB samples missing from the panel.
- Idempotent: re-loading replaces rows for the affected (ingest_id, sample_id).
- Included in dumps/bundles; documented in SCHEMA.md; described in the MCP
  `SCHEMA_DESCRIPTION`.

## 3. Ancestral allele

`info_AA` stays as stored. Add a pure helper that normalises it:
take the text before the first `|`; uppercase base = high confidence,
lowercase = low confidence; `.`, `-`, `N`, empty = unknown. A site is
polarisable when the normalised base equals REF (ancestral = REF) or ALT
(derived = REF; flip counts); otherwise unpolarised.

`db popgen` option `--ancestral {aa,aa-high,ref,none}` (default `aa` when any
`info_AA` is present, else `none`):
`aa` uses high+low confidence, `aa-high` only uppercase, `ref` assumes REF
is ancestral, `none` → folded statistics only. Report how many sites were
polarised / dropped.

## 4. `vcfclick db popgen`

Group with subcommands, shared options:

- Scope: `--region chr:start-end` (repeatable), `--gene SYMBOL` (via the
  existing annotation store, if available), `--ingest-id`.
- Site filters (defaults): biallelic SNVs only (`--include-indels` to widen),
  `--pass-only` (FILTER is PASS or NULL), `--min-call-rate 0.9` (per
  population, from called counts), `--maf` (default 0, applied cohort-wide).
- Autosomes only by default (exclude X/Y/MT and non-PAR sex chromosomes);
  `--include-sex-chroms` uses panel sex for ploidy on X/Y (or document as
  unsupported in v1 and refuse).
- Grouping: by `population` (default), `--by super_population`, or whole
  cohort when no panel is loaded (single group "all").
- `--format table|tsv|json`.

Subcommands:

- `summary` — per group: n samples, sites considered, segregating sites S,
  θ_W (per site and total), π, Tajima's D, Fay & Wu's H (normalised, Zeng
  2006; only if polarised), mean observed vs expected heterozygosity, F.
- `sfs` — per group unfolded (if polarised) and folded SFS. Missing data
  handled by hypergeometric projection down to `--project N` haplotypes
  (default: min called haplotypes across retained sites in that group);
  sites with fewer called haplotypes than N are dropped (reported).
- `fst` — pairwise Hudson F_ST (Bhatia et al. 2013, ratio of averages)
  between groups, genome/region-wide, plus per-window with `--window`.
- `windows` — sliding windows (`--window`, `--step`, bp) along each region:
  S, θ_W, π, Tajima's D per group, and pairwise F_ST if `--fst`. TSV-first,
  ready for plotting.

Computation: SQL computes per-site, per-group allele counts
(alt alleles from `genotypes` dosage joined to `populations`; called
haplotypes = 2 × (group size − missing in group) for autosomes) in both
backend dialects; statistics are computed in Python/numpy from those counts.
θ estimators and D use per-site called sample sizes correctly (π from
per-site 2·k·(n−k)/(n(n−1)) with each site's own n; θ_W/D on the projected
SFS or with per-site a_n — document the choice).

Older databases without called counts or `missing_genotypes`: warn once,
treat absent as 0/0, and include `"missing_data_tracked": false` in JSON
output.

## Validation

- Unit tests for every estimator against hand-computed values and textbook
  examples (θ_W constants, Tajima's D, Hudson F_ST, projection).
- End-to-end fixtures: small VCFs (bgzip/tabix) with known population labels,
  missing calls, polarised and unpolarised sites, haploid calls; run each
  subcommand on chDB and DuckDB and assert identical results.
- Cross-check numbers once against an independent implementation if
  available locally (e.g. scikit-allel or vcftools) and record expected
  values as literals in tests (no new runtime dependency).

## Docs

`docs/POPGEN.md` (new), `docs/CLI.md`, `docs/SCHEMA.md`, README feature list
and docs index, MCP `SCHEMA_DESCRIPTION`, CHANGELOG (Unreleased).
