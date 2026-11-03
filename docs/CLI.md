# CLI reference

Every `vcfclick` command, with a one-line summary, an example, and its full option list (verbatim from `--help`). Commands are grouped by task; deeper how-to guides are linked where they exist.

Run `vcfclick --help` or `vcfclick <command> --help` any time for the same information at the terminal. The terminal UI, web UI, and benchmark commands need their optional extras (`vcfclick[tui]`, `[web]`, `[benchmark]`).

## Commands

- **Databases** — [`db create`](#db-create) · [`db list`](#db-list) · [`db info`](#db-info) · [`db path`](#db-path) · [`db rm`](#db-rm) · [`db query`](#db-query) · [`db stats`](#db-stats) · [`db diff`](#db-diff) · [`db gene`](#db-gene)
- **Ingesting variants** — [`db ingest`](#db-ingest) · [`db ingest-batch`](#db-ingest-batch) · [`merge`](#merge) · [`combine`](#combine) · [`discover`](#discover)
- **Population genetics** — [`db panel`](#db-panel) · [`db popgen`](#db-popgen) ([`summary`](#db-popgen-summary) · [`sfs`](#db-popgen-sfs) · [`fst`](#db-popgen-fst) · [`windows`](#db-popgen-windows))
- **Family / trio analysis** — [`db ped`](#db-ped) · [`db families`](#db-families) · [`db trio`](#db-trio) · [`db qc`](#db-qc) · [`db relatedness`](#db-relatedness)
- **Filter presets** — [`preset list`](#preset-list) · [`preset show`](#preset-show) · [`preset save`](#preset-save) · [`preset set` / `unset` / `edit` / `rename` / `rm`](#preset-set-unset-edit-rename-rm) · [`preset export` / `import`](#preset-export-import) · [`preset validate`](#preset-validate) · [`preset path`](#preset-path)
- **Benchmarking** — [`benchmark`](#benchmark) · [`benchmark-cohort`](#benchmark-cohort)
- **Annotations** — [`annotations load`](#annotations-load) · [`annotations load-clinvar`](#annotations-load-clinvar) · [`annotations load-transcripts`](#annotations-load-transcripts) · [`annotations load-gnomad`](#annotations-load-gnomad) · [`annotations load-hgnc`](#annotations-load-hgnc)
- **Export & sharing** — [`db export`](#db-export) · [`db dump`](#db-dump) · [`db ingest-parquet`](#db-ingest-parquet) · [`db push`](#db-push) · [`db pull`](#db-pull)
- **Interactive UIs** — [`tui`](#tui) · [`web`](#web)

## Databases

One named database per cohort, under `~/.vcfclick` (or `$VCFCLICK_HOME`). See the [User guide](USER_GUIDE.md).

### db create

Create an empty database with the vcfclick schema.

```bash
vcfclick db create trio
```

<details><summary>options</summary>

```
Usage: vcfclick db create [OPTIONS] NAME

  Create a new empty database with the vcfclick schema applied.

Options:
  --help  Show this message and exit.
```

</details>

### db list

List all named databases.

```bash
vcfclick db list
```

<details><summary>options</summary>

```
Usage: vcfclick db list [OPTIONS]

  List all named databases.

Options:
  --help  Show this message and exit.
```

</details>

### db info

Row counts, ingestions, and metadata for a database.

```bash
vcfclick db info trio
```

<details><summary>options</summary>

```
Usage: vcfclick db info [OPTIONS] NAME

  Show metadata about a database (row counts, ingestions, size).

Options:
  --help  Show this message and exit.
```

</details>

### db path

Print a database's on-disk path (no existence check).

```bash
vcfclick db path trio
```

<details><summary>options</summary>

```
Usage: vcfclick db path [OPTIONS] NAME

  Print the on-disk path of a named database (no checks).

Options:
  --help  Show this message and exit.
```

</details>

### db rm

Permanently delete a database and all its data.

```bash
vcfclick db rm --yes trio
```

<details><summary>options</summary>

```
Usage: vcfclick db rm [OPTIONS] NAME

  Permanently delete a named database and all its data.

Options:
  --yes   Confirm the action without prompting.
  --help  Show this message and exit.
```

</details>

### db query

Run SQL against a database and print the result. See [SCHEMA.md](SCHEMA.md).

```bash
vcfclick db query trio "SELECT count() FROM variants"
```

<details><summary>options</summary>

```
Usage: vcfclick db query [OPTIONS] NAME SQL

  Run a SQL query against a named database and print the result.

Options:
  --format TEXT  chDB output format (PrettyCompact, JSON, CSV, TSV, Vertical,
                 ...).  [default: PrettyCompact]
  --help         Show this message and exit.
```

</details>

### db stats

Schema-population stats for an ingested cohort (chDB backend).

```bash
vcfclick db stats trio
```

<details><summary>options</summary>

```
Usage: vcfclick db stats [OPTIONS] NAME

  Schema-population stats for an ingested cohort.

Options:
  --top INTEGER  Show at most TOP overflow-Map keys per table.  [default: 20]
  --help         Show this message and exit.
```

</details>

### db diff

Per-variant allele-frequency comparison across two cohorts.

```bash
vcfclick db diff trio --cohort-a A --cohort-b B
```

<details><summary>options</summary>

```
Usage: vcfclick db diff [OPTIONS] NAME

  Per-variant allele-frequency comparison across two cohorts.

Options:
  --cohort-a TEXT  First cohort name (as stored on samples.cohort).
                   [required]
  --cohort-b TEXT  Second cohort name.  [required]
  --top INTEGER    Limit to the top N variants by absolute AF difference. 0 =
                   no limit.  [default: 50]
  --format TEXT    chDB output format (PrettyCompact, JSON, CSV, TSV, ...).
                   [default: PrettyCompact]
  --help           Show this message and exit.
```

</details>

### db gene

Every variant in a gene (HGNC symbol), with how many samples carry it and how many are homozygous. Needs gene coordinates loaded once with [`annotations load`](#annotations-load); `chr17` and `17` naming both match.

```bash
vcfclick db gene brca_cohort BRCA1
vcfclick db gene brca_cohort BRCA1 --flank 5000 --format tsv > brca1_variants.tsv
```

<details><summary>options</summary>

```
Usage: vcfclick db gene [OPTIONS] NAME SYMBOL

  Variants in a gene (by HGNC symbol), with carrier and hom-alt counts.

Options:
  --flank INTEGER            Extend the gene by this many bp on each side.
                             [default: 0]
  --limit INTEGER            At most this many variants.  [default: 500]
  --format [table|json|tsv]  [default: table]
  --help                     Show this message and exit.
```

</details>

## Ingesting variants

Get VCFs into a database. Multi-allelics must be decomposed first (`bcftools norm -m -`).

### db ingest

Ingest one normalized VCF into a database.

```bash
vcfclick db ingest trio joint.vcf.gz --cohort trio
```

<details><summary>options</summary>

```
Usage: vcfclick db ingest [OPTIONS] NAME VCF_PATH

  Ingest a (normalised) VCF into a named database.

Options:
  --cohort TEXT                   Cohort label this VCF belongs to.  [default:
                                  default]
  --ingest-id TEXT                Stable upload identifier (UUID4 if omitted).
                                  Reuse to replace prior data.
  --workers INTEGER               Parallel worker processes for ingestion.
                                  [default: 4]
  --serial                        Use the single-process serial ingester
                                  instead of parallel.
  --keep-reference                Also store confident hom-reference (0/0)
                                  genotype calls, not just non-reference.
                                  Needed for defensible trio de-novo analysis
                                  (a parent must be provably 0/0, not merely
                                  absent). Stores more rows — use for
                                  trios/families from a joint-called VCF, not
                                  large cohorts. No-calls (./.) are still not
                                  stored as genotypes (see --record-missing).
  --record-missing / --no-record-missing
                                  Record each fully missing call (./.) in the
                                  missing_genotypes table, so per-population
                                  called counts can be derived (needed by `db
                                  popgen` on cohorts with missing data). Per-
                                  site called counts on `variants` (n_called,
                                  an_called, ac_called) are always recorded.
                                  Adds one row per sample per site with a no-
                                  call (about samples x sites x missing rate),
                                  so on large cohorts with much missing data
                                  it can outgrow the genotypes table; skip it
                                  if you will not run `db popgen`.  [default:
                                  record-missing]
  --help                          Show this message and exit.
```

</details>

### db ingest-batch

Ingest many per-sample VCFs as one cohort (from a dir or manifest).

```bash
vcfclick db ingest-batch trio --from-dir vcfs/
```

<details><summary>options</summary>

```
Usage: vcfclick db ingest-batch [OPTIONS] NAME

  Ingest many per-sample VCFs into NAME as one cohort.

Options:
  --from-dir DIRECTORY            Ingest every *.vcf.gz under DIR. ingest_id
                                  is the filename stem.
  --manifest FILE                 TSV with required `vcf_path` column;
                                  optional `sample_id`/`ingest_id` and
                                  `cohort` columns override the defaults.
  --cohort TEXT                   Default cohort label. Required with --from-
                                  dir; used as the fallback for manifest rows
                                  that don't carry their own `cohort`.
  --record-missing / --no-record-missing
                                  Record fully missing calls (./.) in
                                  missing_genotypes (see `db ingest --help`).
                                  [default: record-missing]
  --help                          Show this message and exit.
```

</details>

### merge

Merge per-sample VCFs into one joint VCF (needs bcftools).

```bash
vcfclick merge a.vcf.gz b.vcf.gz -o joint.vcf.gz
```

### combine

Combine callers' call sets with `set=` provenance (GATK3 CombineVariants). See [COMBINE.md](COMBINE.md).

```bash
vcfclick combine gatk.vcf dv.vcf -o out.vcf --name gatk --name dv
```

<details><summary>options</summary>

```
Usage: vcfclick combine [OPTIONS] VCFS...

  Combine two or more VCF call sets into one, with set= provenance.

Options:
  -o, --out FILE            Output combined VCF path (.vcf or .vcf.gz). [required]
  --name TEXT               Set name for an input (repeat once per input, in order).
  --min-callsets INTEGER    Keep only sites present in at least this many inputs.
                            [default: 1]
  --pass-only               Count only PASS calls toward --min-callsets; a filtered
                            input is named filterIn<name> in set=. FILTER '.' = PASS.
  --count-by [allele|site]  Count inputs by exact allele (default) or by position
                            ('site' keeps every allele there; GATK --minimumN).
                            [default: allele]
  --reference FILE          Reference FASTA. Split multi-allelics and left-align/trim
                            internally (bcftools norm -m - -f equivalent) instead of
                            refusing. Needs pyfaidx (vcfclick[benchmark]).
  --carry-info              Carry QUAL/FILTER/INFO (plus header lines) from the
                            highest-priority input that called each allele.
  --help                    Show this message and exit.
```

</details>

### discover

Report which VCF fields land in typed columns vs the overflow map.

```bash
vcfclick discover calls.vcf.gz
```

## Population genetics

Per-population diversity, neutrality tests, site-frequency spectra and F_ST, computed from a database. See [POPGEN.md](POPGEN.md) for definitions, filters and the choices made.

### db panel

Load a sample -> population panel (1000 Genomes `.panel` format, or any TSV/CSV with a header).

```bash
vcfclick db panel kg integrated_call_samples_v3.20130502.ALL.panel
```

<details><summary>options</summary>

```
Usage: vcfclick db panel [OPTIONS] NAME PANEL_PATH

  Load a sample -> population panel into NAME.

  Reads the 1000 Genomes panel format (sample, pop, super_pop, gender; tab-
  separated with a header) as is, or any TSV/CSV with a header — columns are
  matched by name, or named with the --*-col options. A panel replaces the
  whole labelling of each ingestion it applies to, so samples it does not list
  end up unlabelled.

Options:
  --ingest-id TEXT      Apply the panel to this ingestion only (default: every
                        ingestion that contains each sample).
  --sample-col TEXT     Sample-id column name.
  --pop-col TEXT        Population column name.
  --super-pop-col TEXT  Super-population column name (optional).
  --sex-col TEXT        Sex/gender column name (optional).
  --help                Show this message and exit.
```

</details>

### db popgen

Group of population-genetics subcommands. They share the scope, site-filter, grouping, polarisation and output options; each analyses one ingestion.

<details><summary>options</summary>

```
Usage: vcfclick db popgen [OPTIONS] COMMAND [ARGS]...

  Population-genetics statistics (θ, π, Tajima's D, SFS, F_ST).

Options:
  --help  Show this message and exit.

Commands:
  fst      Pairwise Hudson F_ST between groups (Bhatia et al.
  sfs      Per group site-frequency spectrum, projected for missing data.
  summary  Per group: S, θ_W, π, Tajima's D, Fay & Wu's H, Ho/He/F.
  windows  Sliding-window S, θ_W, π, Tajima's D per group (+ F_ST).
```

</details>

### db popgen summary

Per group: S, θ_W, π, Tajima's D, normalised Fay & Wu's H, Ho/He/F.

```bash
vcfclick db popgen summary kg --by super_population
```

<details><summary>options</summary>

```
Usage: vcfclick db popgen summary [OPTIONS] NAME

  Per group: S, θ_W, π, Tajima's D, Fay & Wu's H, Ho/He/F.

Options:
  --project INTEGER RANGE         SFS projection size in haplotypes (default:
                                  the smallest number of called haplotypes at
                                  any retained site, per group).  [x>=2]
  --ingest-id TEXT                Ingestion to analyse (required when NAME has
                                  several).
  --region TEXT                   chr, chr:pos or chr:start-end (repeatable).
  --gene TEXT                     Restrict to a gene (needs `vcfclick
                                  annotations load`).
  --by [population|super_population|all]
                                  Grouping. Falls back to one group 'all' when
                                  no panel is loaded for the ingestion.
                                  [default: population]
  --include-indels                Also use biallelic indels (default:
                                  biallelic SNVs only).
  --pass-only / --all-filters     Keep only sites whose FILTER is PASS or
                                  missing.  [default: pass-only]
  --min-call-rate FLOAT RANGE     Minimum call rate a site needs in EVERY
                                  group.  [default: 0.9; 0<=x<=1]
  --maf FLOAT RANGE               Minimum cohort-wide minor-allele frequency
                                  (0 keeps monomorphic sites, which θ per site
                                  needs).  [default: 0.0; 0<=x<=0.5]
  --ancestral [aa|aa-high|ref|none]
                                  Polarisation: INFO/AA at any confidence
                                  (aa), high-confidence only (aa-high), REF as
                                  ancestral (ref), or folded statistics only
                                  (none). Default: aa when any site carries
                                  INFO/AA, else none.
  --include-sex-chroms            Not supported yet: X/Y need per-sample
                                  ploidy.
  --allow-untracked               Compute even when the ingestion has no
                                  missing-call tracking (ingested before it
                                  existed, or with --no-record-missing).
                                  Output then says missing_data_tracked:
                                  false.
  --format [table|tsv|json]       [default: table]
  --help                          Show this message and exit.
```

</details>

### db popgen sfs

Per group folded (and, when polarised, unfolded) SFS, projected to a common sample size for missing data.

```bash
vcfclick db popgen sfs kg --project 20 --format json
```

<details><summary>options</summary>

```
Usage: vcfclick db popgen sfs [OPTIONS] NAME

  Per group site-frequency spectrum, projected for missing data.

Options:
  --project INTEGER RANGE         SFS projection size in haplotypes (default:
                                  the smallest number of called haplotypes at
                                  any retained site, per group).  [x>=2]
  --ingest-id TEXT                Ingestion to analyse (required when NAME has
                                  several).
  --region TEXT                   chr, chr:pos or chr:start-end (repeatable).
  --gene TEXT                     Restrict to a gene (needs `vcfclick
                                  annotations load`).
  --by [population|super_population|all]
                                  Grouping. Falls back to one group 'all' when
                                  no panel is loaded for the ingestion.
                                  [default: population]
  --include-indels                Also use biallelic indels (default:
                                  biallelic SNVs only).
  --pass-only / --all-filters     Keep only sites whose FILTER is PASS or
                                  missing.  [default: pass-only]
  --min-call-rate FLOAT RANGE     Minimum call rate a site needs in EVERY
                                  group.  [default: 0.9; 0<=x<=1]
  --maf FLOAT RANGE               Minimum cohort-wide minor-allele frequency
                                  (0 keeps monomorphic sites, which θ per site
                                  needs).  [default: 0.0; 0<=x<=0.5]
  --ancestral [aa|aa-high|ref|none]
                                  Polarisation: INFO/AA at any confidence
                                  (aa), high-confidence only (aa-high), REF as
                                  ancestral (ref), or folded statistics only
                                  (none). Default: aa when any site carries
                                  INFO/AA, else none.
  --include-sex-chroms            Not supported yet: X/Y need per-sample
                                  ploidy.
  --allow-untracked               Compute even when the ingestion has no
                                  missing-call tracking (ingested before it
                                  existed, or with --no-record-missing).
                                  Output then says missing_data_tracked:
                                  false.
  --format [table|tsv|json]       [default: table]
  --help                          Show this message and exit.
```

</details>

### db popgen fst

Pairwise Hudson F_ST (ratio of averages), overall and optionally per window.

```bash
vcfclick db popgen fst kg --window 1000000
```

<details><summary>options</summary>

```
Usage: vcfclick db popgen fst [OPTIONS] NAME

  Pairwise Hudson F_ST between groups (Bhatia et al. 2013).

Options:
  --ingest-id TEXT                Ingestion to analyse (required when NAME has
                                  several).
  --region TEXT                   chr, chr:pos or chr:start-end (repeatable).
  --gene TEXT                     Restrict to a gene (needs `vcfclick
                                  annotations load`).
  --by [population|super_population|all]
                                  Grouping. Falls back to one group 'all' when
                                  no panel is loaded for the ingestion.
                                  [default: population]
  --include-indels                Also use biallelic indels (default:
                                  biallelic SNVs only).
  --pass-only / --all-filters     Keep only sites whose FILTER is PASS or
                                  missing.  [default: pass-only]
  --min-call-rate FLOAT RANGE     Minimum call rate a site needs in EVERY
                                  group.  [default: 0.9; 0<=x<=1]
  --maf FLOAT RANGE               Minimum cohort-wide minor-allele frequency
                                  (0 keeps monomorphic sites, which θ per site
                                  needs).  [default: 0.0; 0<=x<=0.5]
  --ancestral [aa|aa-high|ref|none]
                                  Polarisation: INFO/AA at any confidence
                                  (aa), high-confidence only (aa-high), REF as
                                  ancestral (ref), or folded statistics only
                                  (none). Default: aa when any site carries
                                  INFO/AA, else none.
  --include-sex-chroms            Not supported yet: X/Y need per-sample
                                  ploidy.
  --allow-untracked               Compute even when the ingestion has no
                                  missing-call tracking (ingested before it
                                  existed, or with --no-record-missing).
                                  Output then says missing_data_tracked:
                                  false.
  --format [table|tsv|json]       [default: table]
  --window INTEGER RANGE          Window size (bp).  [x>=1]
  --step INTEGER RANGE            Window step (bp; default = window).  [x>=1]
  --help                          Show this message and exit.
```

</details>

### db popgen windows

Sliding-window S, θ_W, π, Tajima's D per group, plus pairwise F_ST with `--fst`. TSV by default, one row per window.

```bash
vcfclick db popgen windows kg --region 22:16000000-17000000 --window 100000 --fst > win.tsv
```

<details><summary>options</summary>

```
Usage: vcfclick db popgen windows [OPTIONS] NAME

  Sliding-window S, θ_W, π, Tajima's D per group (+ F_ST).

Options:
  --project INTEGER RANGE         SFS projection size in haplotypes (default:
                                  the smallest number of called haplotypes at
                                  any retained site, per group).  [x>=2]
  --ingest-id TEXT                Ingestion to analyse (required when NAME has
                                  several).
  --region TEXT                   chr, chr:pos or chr:start-end (repeatable).
  --gene TEXT                     Restrict to a gene (needs `vcfclick
                                  annotations load`).
  --by [population|super_population|all]
                                  Grouping. Falls back to one group 'all' when
                                  no panel is loaded for the ingestion.
                                  [default: population]
  --include-indels                Also use biallelic indels (default:
                                  biallelic SNVs only).
  --pass-only / --all-filters     Keep only sites whose FILTER is PASS or
                                  missing.  [default: pass-only]
  --min-call-rate FLOAT RANGE     Minimum call rate a site needs in EVERY
                                  group.  [default: 0.9; 0<=x<=1]
  --maf FLOAT RANGE               Minimum cohort-wide minor-allele frequency
                                  (0 keeps monomorphic sites, which θ per site
                                  needs).  [default: 0.0; 0<=x<=0.5]
  --ancestral [aa|aa-high|ref|none]
                                  Polarisation: INFO/AA at any confidence
                                  (aa), high-confidence only (aa-high), REF as
                                  ancestral (ref), or folded statistics only
                                  (none). Default: aa when any site carries
                                  INFO/AA, else none.
  --include-sex-chroms            Not supported yet: X/Y need per-sample
                                  ploidy.
  --allow-untracked               Compute even when the ingestion has no
                                  missing-call tracking (ingested before it
                                  existed, or with --no-record-missing).
                                  Output then says missing_data_tracked:
                                  false.
  --format [table|tsv|json]       [default: tsv]
  --window INTEGER RANGE          Window size (bp).  [x>=1; required]
  --step INTEGER RANGE            Window step (bp; default = window).  [x>=1]
  --fst                           Add pairwise Hudson F_ST columns.
  --help                          Show this message and exit.
```

</details>

## Family / trio analysis

Pedigree-aware inheritance analysis. See [Trio](TRIO.md).

### db ped

Load family relationships from a PED/FAM file, or declare a family with
flags: several families, grandparents, half siblings and duos are supported.
A PED file replaces the stored pedigree. The flags describe a few people of one
family, and a call that would delete individuals it does not name fails until
you add `--merge`, which rewrites only the people the flags name and leaves
everybody else as stored (with a PED file, `--merge` replaces the families the
file lists and says who it removed). Flag-declared parents are stored with
unknown phenotype. See
[TRIO.md](TRIO.md#families-generations-and-probands).

```bash
vcfclick db ped trio fam.ped
vcfclick db ped trio --proband HG002 --father HG003 --mother HG004
vcfclick db ped fams --family famB --merge --proband K1 --proband K2 \
    --father FB --mother MB --parents FB=GF,GM
```

<details><summary>options</summary>

```
Usage: vcfclick db ped [OPTIONS] NAME [PED_PATH]

  Load family relationships into NAME, from a PED/FAM file or as a family
  given by sample ids (--proband, repeatable for siblings, with --father and
  --mother; --parents adds grandparents and half siblings).

  The individual ids must match sample ids already ingested under the target
  ingest_id, and all members share that one ingest_id. A PED file replaces the
  stored pedigree for that ingest_id. Flags describe some people of one
  family; add --merge to change just them and keep everybody else as stored.
  Everything is checked before anything is deleted.

Options:
  --ingest-id TEXT                Ingest_id the pedigree's samples belong to.
                                  Inferred when the database has exactly one
                                  ingestion.
  --proband TEXT                  Without a PED file: an affected child's
                                  sample id. Repeat it for siblings who share
                                  the same parents.
  --father TEXT                   Without a PED file: the proband's father
                                  (sample id; 0 when only the mother is
                                  known).
  --mother TEXT                   Without a PED file: the proband's mother
                                  (sample id; 0 when only the father is
                                  known).
  --family TEXT                   Without a PED file: the family id of every
                                  listed sample. [default: fam1]
  --parents CHILD=FATHER,MOTHER   Without a PED file: the parents of another
                                  member, to declare grandparents or half
                                  siblings (e.g. F=PGF,PGM). Repeatable;
                                  either side may be 0.
  --merge                         Keep what this call does not name. With
                                  flags: rewrite only the people the flags
                                  name, and only what the flags say about
                                  them, and leave everybody else as stored.
                                  With a PED file: replace the families it
                                  lists as wholes (members it omits are
                                  removed, and listed) and keep the other
                                  families. Without it a PED file replaces the
                                  whole pedigree, and a flag call that would
                                  delete individuals it does not name fails.
  --affected TEXT                 Without a PED file: mark a listed sample
                                  affected. Repeatable.
  --unaffected TEXT               Without a PED file: mark a listed sample
                                  unaffected. Repeatable.
  --unknown TEXT                  Without a PED file: mark a listed sample's
                                  phenotype unknown. Repeatable. Probands
                                  default to affected, everyone else to
                                  unknown.
  --sex ID=male|female|unknown    Without a PED file: a listed sample's sex.
                                  Repeatable. A father is male and a mother
                                  female unless this says otherwise.
  --proband-sex [male|female|unknown]
                                  Without a PED file: the proband's sex (a
                                  single --proband only; use --sex for each
                                  sibling, or a PED file).  [default: unknown]
  --help                          Show this message and exit.
```

</details>

### db families

List the families in a database's pedigree: members by generation, the probands `db trio` analyses, every trio or duo with its siblings, and problems in the pedigree. Read-only; exits 0 even when it finds problems. See [TRIO.md](TRIO.md#probands-and-db-families).

```bash
vcfclick db families trio
vcfclick db families trio --format json
```

<details><summary>options</summary>

```
Usage: vcfclick db families [OPTIONS] NAME

  List the families in NAME's pedigree: members by generation, the probands
  `db trio` analyses, each trio or duo with its siblings, and problems in the
  pedigree. Read-only; exits 0 even when it finds problems.

Options:
  --ingest-id TEXT      Ingest_id whose pedigree to list. Inferred when the
                        database has exactly one ingestion.
  --format [text|json]  json: {ingest_id, families: [{family_id, members,
                        probands, trios}], problems}.  [default: text]
  --help                Show this message and exit.
```

</details>

### db trio

Candidate variants under Mendelian models (de-novo/recessive/dominant/comp-het) and the sex-aware X-linked model. See [TRIO.md](TRIO.md).

`--category xlinked` lists the X-linked candidates (hemizygous sons, het and hom-alt daughters, X de novo, female compound het on non-PAR chrX); `--build auto|GRCh37|GRCh38` decides where non-PAR X is; `--min-male-x-dp`, `--max-af-x`, `--min-hemi-ab` and `--max-parent-alt-reads` are its gates. The X guard (on by default, `--no-x-guard` off) drops a het father and a het male proband on non-PAR X from the four classic models. See [TRIO.md](TRIO.md#chromosome-x).

```bash
vcfclick db trio trio --proband HG002 --category denovo
vcfclick db trio fams                      # the probands of every family
vcfclick db trio fams --family famB        # one family
vcfclick db trio fams --all-trios          # every child with both parents
```

Without `--proband`, `db trio` analyses the **probands**: children with both
parents declared who are not marked unaffected and who either have no children
of their own or are affected (`db families` lists them). `--all-trios` is the
older rule, every child with both parents declared. `--format json` adds
`trio.family_id`, `trio.kind` and `trio.siblings` to each document, and a
`families` list to the `trios` document; see
[TRIO.md](TRIO.md#running-db-trio-over-a-family).

`--preset NAME` applies a [filter preset](PRESETS.md): the built-in starting points `broad-first-pass`, `strict-denovo`, `rare-recessive` and `rare-dominant`, or your own. A flag typed next to it wins. A run with a preset says, on stdout, which of its gates cannot bind in this data (a `strict-denovo` run refuses instead). `--category` takes a comma list (`recessive,comphet`). `--strict-missing`, `--pass-only`, `--parent-alt-screen` and `--require-bound-gates` are opt-in gates; with none of them and no preset the output is what it always was. See [PRESETS.md](PRESETS.md).

`--genes`, `--gene-file`, `--region`, `--region-file` and their `--exclude-*`
twins limit a run to a gene panel or regions, for every model: gene symbols in
any case (`C9orf72`), Ensembl ids and, after
[`annotations load-hgnc`](#annotations-load-hgnc), previous symbols, HGNC and
Entrez ids. A short list that loses a gene is an error
(`--allow-unresolved-genes` runs anyway). See
[TRIO.md](TRIO.md#genes-and-regions).

```bash
vcfclick db trio fam1 --proband C1 --genes C9orf72,C19orf12,ENSG00000111678
vcfclick db trio fam1 --proband C1 --gene-file id_panel.txt --gene-flank 500
vcfclick db trio fam1 --proband C1 --region chr7:117100000-117400000
vcfclick db trio fam1 --proband C1 --exclude-region chr6:28510120-33480577
```

`--format json` reports, for every model, `blocked` and `blocked_reason`
(`null`, or `needs_keep_reference` / `needs_genes`), and under `gates` the
`af_gate` outcome (`applied`, `cohort_derived_suspected`,
`skipped_cohort_derived`, `unavailable`). When the VCF's `INFO/AF` is the
cohort's own frequency, the `--max-af` gate is skipped with a note on stderr
unless `--max-af` is given; see
[TRIO.md](TRIO.md#cohort-derived-infoaf). The JSON is additive-only; the
[compatibility policy](TRIO.md#compatibility-policy) is the contract the
desktop apps rely on.

<details><summary>options</summary>

```
Usage: vcfclick db trio [OPTIONS] NAME

  Report candidate variants under Mendelian inheritance models for a trio,
  with genotype quality gates and an AF rarity filter. Without --proband,
  every proband in the pedigree is analysed in turn (--family limits that to
  one family, --all-trios widens it to every child).

Options:
  --proband TEXT                  Sample id of the affected child. Omitted:
                                  every proband in the pedigree, one report
                                  each. A proband is a child with both parents
                                  declared who is not marked unaffected and
                                  who either has no children of their own or
                                  is affected (`db families` lists them).
  --family TEXT                   Analyse every proband of this family (see
                                  `db families`). Not with --proband.
  --all-trios                     Without --proband: analyse every child with
                                  both parents declared, including children
                                  marked unaffected and parents who have
                                  children of their own (the rule before
                                  probands were defined).
  --category [denovo|recessive|dominant|comphet|xlinked|all]
                                  Inheritance model: a name, a comma list
                                  (recessive,comphet) or 'all'. 'all' prints
                                  per-model candidate counts.  [default: all]
  --min-gq INTEGER                Minimum genotype quality (GQ) of the proband
                                  and both parents. A call with no GQ passes
                                  unless --strict-missing.  [default: 20]
  --min-dp INTEGER                Minimum read depth (DP) of the proband and
                                  both parents. A call with no DP passes
                                  unless --strict-missing.  [default: 10]
  --max-af FLOAT                  Keep variants with population info_AF <=
                                  this (rarity filter). If INFO/AF turns out
                                  to be this cohort's own frequency (GATK-
                                  style family VCFs), the filter is skipped
                                  with a note unless this option is given
                                  explicitly.  [default: 0.01]
  --min-ab FLOAT                  Lowest allele balance (alt reads / all
                                  reads) a proband het call may have.
                                  [default: 0.25]
  --max-ab FLOAT                  Highest allele balance a proband het call
                                  may have.  [default: 0.75]
  --gnomad-max-af FLOAT           Additionally drop candidates whose gnomAD
                                  popmax AF exceeds this (needs `vcfclick
                                  annotations load-gnomad`). Variants absent
                                  from the loaded gnomAD slice are kept as
                                  rare.
  --strict-missing / --no-strict-missing
                                  A missing GQ, DP or (for a het call) AD on
                                  the proband or a parent counts as a fail,
                                  not a pass, for each of those gates that is
                                  asked for. A gate whose field is absent from
                                  more than 10% of the rows the selected
                                  models read at their candidate sites cannot
                                  be enforced and is an error (--no-require-
                                  bound-gates runs anyway, with that gate
                                  tolerant for the run).
  --pass-only / --no-pass-only    Keep only variants whose VCF FILTER is PASS
                                  or unset.
  --parent-alt-screen / --no-parent-alt-screen
                                  De novo only: a parent's stored 0/0 is
                                  trusted only when its alt reads stay under
                                  --max-parent-alt-frac (a parent 0/0 with >=
                                  2 alt reads and a high alt fraction is not
                                  proof of absence). Shows the parents' read
                                  counts on de novo candidates.
  --max-parent-alt-frac FLOAT     With --parent-alt-screen: the parent alt-
                                  read fraction at and above which a stored
                                  0/0 is not trusted.  [default: 0.05]
  --require-bound-gates / --no-require-bound-gates
                                  Error, instead of a warning, when a
                                  requested gate cannot bind in this data (the
                                  field is missing on too many rows), a
                                  selected model cannot run on this database,
                                  or no population-frequency filter is in
                                  effect. --strict-missing implies it for the
                                  GQ, DP and AD gates; --no-require-bound-
                                  gates overrides both for this run.
  --min-male-x-dp INTEGER         xlinked: minimum DP of a hemizygous call (a
                                  male proband, and the father of a daughter)
                                  on non-PAR chrX. Default: half of --min-dp,
                                  rounded down (one X copy gives half the
                                  reads).
  --max-af-x FLOAT                xlinked: rarity cap on non-PAR chrX; the
                                  model uses the lower of this and --max-af. A
                                  hemizygous male is affected at the allele
                                  frequency q, not q squared, so X needs a
                                  stricter ceiling.  [default: 0.001]
  --min-hemi-ab FLOAT             xlinked: lowest alt fraction (alt reads /
                                  all reads) of a male proband's hemizygous
                                  ALT call. A call without allele depths
                                  passes.  [default: 0.75]
  --max-parent-alt-reads INTEGER  xlinked de novo: most alt reads a parent's
                                  stored 0/0 may carry (a missing allele depth
                                  passes).  [default: 1]
  --genes SYMBOL[,SYMBOL...]      Only variants inside these genes: symbols
                                  (any case), Ensembl gene ids (ENSG...) and,
                                  after `vcfclick annotations load-hgnc`,
                                  previous symbols, HGNC ids (HGNC:1234) and
                                  Entrez ids. Comma list, repeatable. A gene
                                  is its whole span (introns and UTRs);
                                  promoters need --gene-flank. Needs `vcfclick
                                  annotations load`. --gene-file FILE reads
                                  one per line (# comments); typed with
                                  --genes, the two are joined.
  --gene-file FILE                Read --genes values from a file (see
                                  --genes).
  --exclude-genes SYMBOL[,SYMBOL...]
                                  Leave out variants inside these genes (same
                                  tokens as --genes). --exclude-gene-file FILE
                                  reads one per line.
  --exclude-gene-file FILE        Read --exclude-genes values from a file (see
                                  --exclude-genes).
  --gene-flank INTEGER            Widen every gene of --genes and --exclude-
                                  genes by this many base pairs on each side
                                  (promoters, upstream elements).  [default:
                                  0]
  --region REGION                 Only variants inside this region, chrN:FROM-
                                  TO, 1-based and inclusive (repeatable). With
                                  --genes, a variant in either is kept.
                                  --region-file BED reads 0-based half-open
                                  BED intervals.
  --region-file BED               Read --region values from a file (see
                                  --region).
  --exclude-region REGION         Leave out variants inside this region,
                                  chrN:FROM-TO (repeatable).
  --allow-unresolved-genes        Run even when tokens of --genes / --exclude-
                                  genes resolve to no gene (by default a list
                                  of 20 or fewer may miss none, a longer one
                                  at most 5%); they are left out of the scope
                                  and reported.
  --build [auto|GRCh37|GRCh38]    Assembly of chrX, for its non-PAR region
                                  (xlinked and the X guard). auto: GRCh38 when
                                  a chrX position exceeds the GRCh37 chrX
                                  length (155,270,560), else unknown: the two
                                  bands where GRCh37 and GRCh38 disagree are
                                  left out and counted. CHM13/T2T is not
                                  supported.  [default: auto]
  --x-guard / --no-x-guard        On sure non-PAR chrX, drop a het father, and
                                  a het male proband, from denovo, recessive,
                                  dominant and comphet: a male has one X, so
                                  those calls are artefacts. --no-x-guard
                                  restores the unguarded models.  [default:
                                  x-guard]
  --format [text|json]            json: every model's count and candidates
                                  (with gene, gnomAD and ClinVar when loaded)
                                  in one document; without --proband, a
                                  `trios` list with one document per proband,
                                  and the `families` they are in.  [default:
                                  text]
  --limit INTEGER                 json: at most this many candidates per
                                  model.  [default: 1000]
  --preset NAME|PATH              Apply a filter preset (`vcfclick preset
                                  list`): its values are the defaults, and a
                                  flag typed here wins. NAME@FINGERPRINT fails
                                  the run if the preset's content changed; a
                                  value with a path separator or ending .json
                                  is a file.
  --preset-fingerprint HEX12      Fail unless the preset's fingerprint is this
                                  (also for a preset given as a path).
  --preset-ignore KEY             Switch one preset key off for this run
                                  (repeatable).
  --preset-skip-unavailable       Apply the preset keys this vcfclick knows
                                  and skip the rest (they are NOT applied, and
                                  the output says so); the default is an
                                  error.
  --help                          Show this message and exit.
```

</details>

### db qc

Per-sample QC: het/hom, Ti/Tv, chrX-het sex check. See [QC.md](QC.md).

```bash
vcfclick db qc trio
vcfclick db qc trio --x-nonpar     # chrX het fraction on non-PAR X only
```

<details><summary>options</summary>

```
Usage: vcfclick db qc [OPTIONS] NAME

  Per-sample QC: het/hom ratio, Ti/Tv, and a chrX-het sex check.

Options:
  --format [table|json]         [default: table]
  --x-nonpar                    Take the chrX het fraction on sure non-PAR X
                                only (the region `db trio` infers sex on), so
                                both commands report the same sex.
  --build [auto|GRCh37|GRCh38]  With --x-nonpar: the assembly whose PAR is
                                left out. auto: GRCh38 when a chrX position
                                exceeds the GRCh37 chrX length, else only the
                                region that is non-PAR on both builds.
                                CHM13/T2T is not supported.  [default: auto]
  --help                        Show this message and exit.
```

</details>

### db relatedness

Pairwise kinship between the samples of an ingestion (KING-robust): finds duplicates, parent/child and sibling pairs, and checks a loaded pedigree against the genotypes. Needs genome-wide data: on a single gene or small region it warns instead of classifying, because unrelated people sharing a haplotype would look related. See [QC.md](QC.md#relatedness).

```bash
vcfclick db relatedness cohort
vcfclick db relatedness trio --all --format json
```

<details><summary>options</summary>

```
Usage: vcfclick db relatedness [OPTIONS] NAME

  Pairwise kinship (KING-robust): duplicates, relatives, pedigree errors.

Options:
  --ingest-id TEXT       Only this ingestion (default: each ingestion
                         separately).
  --min-kinship FLOAT    Report pairs at or above this kinship (third degree
                         by default).  [default: 0.0442]
  --all                  Report every pair, including unrelated ones.
  --max-sites INTEGER    Thin to at most this many SNVs (evenly spread) to
                         bound memory.  [default: 200000]
  --force                Classify pairs even when the data covers too small a
                         region.
  --format [table|json]  [default: table]
  --help                 Show this message and exit.
```

</details>

## Filter presets

Named, versioned, shareable sets of `db trio` filters. A preset is a small JSON
file of filter values, with no sample, family, path or database in it. See
[Presets](PRESETS.md) for the guide: precedence, the effective-gate report,
the file format and sharing.

### preset list

List the built-in presets (starting points, not validated), your own in
`$VCFCLICK_HOME/presets`, and any in the directories on `VCFCLICK_PRESET_PATH`.
Always exits 0; a file that cannot be read is listed with its problem.

```bash
vcfclick preset list
vcfclick preset list --source user --format json
```

<details><summary>options</summary>

```
Usage: vcfclick preset list [OPTIONS]

  List presets. Always exits 0: a file that cannot be read is listed with its
  problem, never hidden.

Options:
  --source [builtin|user|shared|all]
                                  Which presets: the built-ins, your own, the
                                  shared directories (VCFCLICK_PRESET_PATH),
                                  or all.  [default: all]
  --command [trio|segregate|case]
                                  Only presets meant for this command (a hint
                                  the author set).
  --format [text|json]            [default: text]
  --help                          Show this message and exit.
```

</details>

### preset show

Show one preset: its filters, what it covers and does not cover, its revision
history with fingerprints, and who approved it.

```bash
vcfclick preset show strict-denovo
vcfclick preset show id-panel --resolve --format json
```

<details><summary>options</summary>

```
Usage: vcfclick preset show [OPTIONS] NAME|PATH

  Show one preset: its values, what it covers and does not, its revision
  history and fingerprint.

Options:
  --command [trio|segregate|case]
                                  Mark the keys this command does not take.
  --resolve                       Show the flag each key maps to.
  --format [text|json]            [default: text]
  --help                          Show this message and exit.
```

</details>

### preset save

Save a preset from filter flags, from another preset, from a saved
`db trio --format json` result, or from JSON.

```bash
vcfclick preset save id-panel --from-preset rare-recessive --max-af 0.003 \
    --description "ID recessive first pass"
vcfclick preset save redo --from-result result.json
```

<details><summary>options</summary>

```
Usage: vcfclick preset save [OPTIONS] NAME

  Save a preset under NAME in your presets directory.

  The filter flags are the ones `db trio` takes: only those you type are saved
  (a value you leave out is inherited from the defaults when the preset is
  used). Built-in names are taken; so is an existing name unless --force.

Options:
  --from-preset BASE              Start from this preset's filters (copied;
                                  there is no inheritance).
  --from-result RESULT.json       Rebuild the filters from a saved `db trio
                                  --format json` result that used a preset:
                                  the values that came from the preset or a
                                  flag.
  --from-json FILE|-              Read {"filters": {...}, "title": ...,
                                  "description": ...} from a file or stdin.
  --title TEXT                    A short title.
  --description TEXT              What the preset is for.
  --force                         Replace a preset of that name.
  --category [denovo|recessive|dominant|comphet|xlinked|all]
                                  Inheritance model: a name, a comma list
                                  (recessive,comphet) or 'all'. 'all' prints
                                  per-model candidate counts.
  --min-gq INTEGER                Minimum genotype quality (GQ) of the proband
                                  and both parents. A call with no GQ passes
                                  unless --strict-missing.
  --min-dp INTEGER                Minimum read depth (DP) of the proband and
                                  both parents. A call with no DP passes
                                  unless --strict-missing.
  --max-af FLOAT                  Keep variants with population info_AF <=
                                  this (rarity filter). If INFO/AF turns out
                                  to be this cohort's own frequency (GATK-
                                  style family VCFs), the filter is skipped
                                  with a note unless this option is given
                                  explicitly.
  --min-ab FLOAT                  Lowest allele balance (alt reads / all
                                  reads) a proband het call may have.
  --max-ab FLOAT                  Highest allele balance a proband het call
                                  may have.
  --gnomad-max-af FLOAT           Additionally drop candidates whose gnomAD
                                  popmax AF exceeds this (needs `vcfclick
                                  annotations load-gnomad`). Variants absent
                                  from the loaded gnomAD slice are kept as
                                  rare.
  --strict-missing / --no-strict-missing
                                  A missing GQ, DP or (for a het call) AD on
                                  the proband or a parent counts as a fail,
                                  not a pass, for each of those gates that is
                                  asked for. A gate whose field is absent from
                                  more than 10% of the rows the selected
                                  models read at their candidate sites cannot
                                  be enforced and is an error (--no-require-
                                  bound-gates runs anyway, with that gate
                                  tolerant for the run).
  --pass-only / --no-pass-only    Keep only variants whose VCF FILTER is PASS
                                  or unset.
  --parent-alt-screen / --no-parent-alt-screen
                                  De novo only: a parent's stored 0/0 is
                                  trusted only when its alt reads stay under
                                  --max-parent-alt-frac (a parent 0/0 with >=
                                  2 alt reads and a high alt fraction is not
                                  proof of absence). Shows the parents' read
                                  counts on de novo candidates.
  --max-parent-alt-frac FLOAT     With --parent-alt-screen: the parent alt-
                                  read fraction at and above which a stored
                                  0/0 is not trusted.
  --require-bound-gates / --no-require-bound-gates
                                  Error, instead of a warning, when a
                                  requested gate cannot bind in this data (the
                                  field is missing on too many rows), a
                                  selected model cannot run on this database,
                                  or no population-frequency filter is in
                                  effect. --strict-missing implies it for the
                                  GQ, DP and AD gates; --no-require-bound-
                                  gates overrides both for this run.
  --min-male-x-dp INTEGER         xlinked: minimum DP of a hemizygous call (a
                                  male proband, and the father of a daughter)
                                  on non-PAR chrX. Default: half of --min-dp,
                                  rounded down (one X copy gives half the
                                  reads).
  --max-af-x FLOAT                xlinked: rarity cap on non-PAR chrX; the
                                  model uses the lower of this and --max-af. A
                                  hemizygous male is affected at the allele
                                  frequency q, not q squared, so X needs a
                                  stricter ceiling.
  --min-hemi-ab FLOAT             xlinked: lowest alt fraction (alt reads /
                                  all reads) of a male proband's hemizygous
                                  ALT call. A call without allele depths
                                  passes.
  --max-parent-alt-reads INTEGER  xlinked de novo: most alt reads a parent's
                                  stored 0/0 may carry (a missing allele depth
                                  passes).
  --genes SYMBOL[,SYMBOL...]      Only variants inside these genes: symbols
                                  (any case), Ensembl gene ids (ENSG...) and,
                                  after `vcfclick annotations load-hgnc`,
                                  previous symbols, HGNC ids (HGNC:1234) and
                                  Entrez ids. Comma list, repeatable. A gene
                                  is its whole span (introns and UTRs);
                                  promoters need --gene-flank. Needs `vcfclick
                                  annotations load`. --gene-file FILE reads
                                  one per line (# comments); typed with
                                  --genes, the two are joined.
  --gene-file FILE                Read --genes values from a file (see
                                  --genes).
  --genes-label TEXT              A label for the gene list, such as the panel
                                  and its version (`PanelApp ID v4.2`).
                                  Descriptive only.
  --exclude-genes SYMBOL[,SYMBOL...]
                                  Leave out variants inside these genes (same
                                  tokens as --genes). --exclude-gene-file FILE
                                  reads one per line.
  --exclude-gene-file FILE        Read --exclude-genes values from a file (see
                                  --exclude-genes).
  --gene-flank INTEGER            Widen every gene of --genes and --exclude-
                                  genes by this many base pairs on each side
                                  (promoters, upstream elements).
  --region REGION                 Only variants inside this region, chrN:FROM-
                                  TO, 1-based and inclusive (repeatable). With
                                  --genes, a variant in either is kept.
                                  --region-file BED reads 0-based half-open
                                  BED intervals.
  --region-file BED               Read --region values from a file (see
                                  --region).
  --exclude-region REGION         Leave out variants inside this region,
                                  chrN:FROM-TO (repeatable).
  --regions-build [GRCh37|GRCh38]
                                  The assembly the regions are written in
                                  (GRCh37 or GRCh38), checked against the
                                  database's build when it is known.
  --help                          Show this message and exit.
```

</details>

### preset set / unset / edit / rename / rm

Change a preset of yours. Built-in and shared presets are read-only, and a
preset with `meta.locked` refuses these without `--force`.

```bash
vcfclick preset set id-panel max_af=0.002 strict_missing=true
vcfclick preset unset id-panel gnomad_max_af
vcfclick preset edit id-panel
vcfclick preset rename id-panel id-panel-v2
vcfclick preset rm id-panel-v2 --yes
```

<details><summary>options</summary>

```
Usage: vcfclick preset set [OPTIONS] NAME KEY=VALUE...

  Set filter values in a preset of yours: `preset set id-panel max_af=0.003`.

Options:
  --force  Change a locked preset.
  --help   Show this message and exit.
```

</details>

<details><summary>options</summary>

```
Usage: vcfclick preset unset [OPTIONS] NAME KEYS...

  Remove filter keys from a preset of yours (they fall back to defaults).

Options:
  --force  Change a locked preset.
  --help   Show this message and exit.
```

</details>

<details><summary>options</summary>

```
Usage: vcfclick preset edit [OPTIONS] NAME

  Open a preset of yours in $EDITOR; it is validated before it is saved.

Options:
  --force  Edit a locked preset.
  --help   Show this message and exit.
```

</details>

<details><summary>options</summary>

```
Usage: vcfclick preset rename [OPTIONS] OLD NEW

  Rename a preset of yours.

Options:
  --force  Rename a locked preset.
  --help   Show this message and exit.
```

</details>

<details><summary>options</summary>

```
Usage: vcfclick preset rm [OPTIONS] NAME

  Delete a preset of yours. `rm` and `delete` are the same command.

Options:
  --yes    Do not ask.
  --force  Delete a locked preset.
  --help   Show this message and exit.
```

</details>

### preset export / import

Hand presets to a colleague. Export leaves out patient-specific keys unless
asked; import is all or nothing and records only the source file's name.

```bash
vcfclick preset export id-panel -o id-panel.vcfclick-preset.json
vcfclick preset import id-panel.vcfclick-preset.json --as colleague-panel
```

<details><summary>options</summary>

```
Usage: vcfclick preset export [OPTIONS] NAMES...

  Write presets to a file to hand to a colleague.

Options:
  -o, --output FILE|DIR  Where to write: a file, or a directory (one file per
                         preset). Without it the preset is printed.
  --bundle               Put several presets in one file.
  --with-patient-data    Keep patient-specific keys (such as HPO terms); they
                         are left out of an export by default.
  --help                 Show this message and exit.
```

</details>

<details><summary>options</summary>

```
Usage: vcfclick preset import [OPTIONS] FILES...

  Import presets from files or bundles. All or nothing.

Options:
  --as NAME             Store it under this name.
  --force               Replace a preset of the same name.
  --allow-unknown-keys  Store a preset with filter keys this vcfclick does not
                        know. Using it is then an error unless --preset-skip-
                        unavailable.
  --help                Show this message and exit.
```

</details>

### preset validate

Check a preset file; exits 1 on an error.

```bash
vcfclick preset validate id-panel.vcfclick-preset.json --command trio
```

<details><summary>options</summary>

```
Usage: vcfclick preset validate [OPTIONS] FILE

  Check a preset file. Exits 1 when it has an error.

Options:
  --command [trio|segregate|case]
                                  Also report the keys this command does not
                                  take.
  --format [text|json]            [default: text]
  --help                          Show this message and exit.
```

</details>

### preset path

Print where presets are looked up: your directory, then each shared one.

```bash
vcfclick preset path
```

<details><summary>options</summary>

```
Usage: vcfclick preset path [OPTIONS]

  Print where presets are looked up: your directory, then each shared one
  (VCFCLICK_PRESET_PATH).

Options:
  --help  Show this message and exit.
```

</details>

## Benchmarking

Concordance benchmarking (hap.py-style) + analytics. See [Benchmarking](BENCHMARK.md).

### benchmark

Benchmark a query VCF against a truth VCF; recall/precision/F1. See [BENCHMARK.md](BENCHMARK.md).

```bash
vcfclick benchmark --truth t.vcf.gz --query q.vcf.gz --ref g.fa --regions c.bed -o report/
```

<details><summary>options</summary>

```
Usage: vcfclick benchmark [OPTIONS]

  Benchmark a query VCF against a truth VCF (normalized genotype concordance).

Options:
  --truth FILE                    Truth (gold-standard) VCF.  [required]
  --query FILE                    Query (call set under test) VCF.  [required]
  --ref FILE                      Reference FASTA (indexed .fai).  [required]
  --regions FILE                  Confident-region BED. Omit to treat every
                                  call as confident.
  -o, --output DIRECTORY          Output directory for reports.  [required]
  --engine [normalized|haplotype|exact]
                                  Reconciliation engine. 'haplotype' adds
                                  local-haplotype rescue; 'exact' is a
                                  diagnostic that skips reference left-
                                  alignment.  [default: normalized]
  --report-formats TEXT           Comma-separated subset of
                                  csv,json,parquet,html (or 'all').  [default:
                                  all]
  --on-ref-mismatch [error|skip]  Behaviour when a REF allele disagrees with
                                  the reference.  [default: error]
  --conf-containment [start|full]
                                  Confident-region membership: variant start,
                                  or whole ref span.  [default: start]
  --decompose-mnp / --no-decompose-mnp
                                  Atomize MNPs into per-position SNPs (loses
                                  phase; off by default).  [default: no-
                                  decompose-mnp]
  --strict                        Promote warnings (e.g. missing --regions) to
                                  hard errors.
  --pass-only / --all             Show only the PASS (or only the ALL) filter
                                  view in the headline.
  --stratify TEXT                 Comma-separated concordance stratification
                                  axes (gnomad,clinvar,gene) joined against
                                  the annotation store; writes
                                  stratified_<axis>.csv.
  --roc                           Write roc.tsv — a query-quality threshold
                                  sweep (recall/precision trade-off per
                                  variant type).
  --strat-region NAME=REGIONS.bed
                                  Stratify concordance by a genome-region set
                                  (e.g. lowcomplex=lc.bed); repeatable. Writes
                                  stratified_regions.csv.
  --audit                         Write fn_annotated.csv / fp_annotated.csv —
                                  each FN/FP joined to its gene, ClinVar
                                  significance, and gnomAD AF.
  --help                          Show this message and exit.
```

</details>

### benchmark-cohort

Benchmark many callers against one truth; per-caller + blind spots. See [BENCHMARK.md](BENCHMARK.md).

```bash
vcfclick benchmark-cohort --truth t.vcf.gz --ref g.fa --caller gatk=g.vcf.gz -o cohort/
```

<details><summary>options</summary>

```
Usage: vcfclick benchmark-cohort [OPTIONS]

  Benchmark several callers against one truth; report per-caller concordance
  and each caller's relative blind spots (variants others catch).

Options:
  --truth FILE                    Truth VCF.  [required]
  --ref FILE                      Reference FASTA.  [required]
  --regions FILE                  Confident-region BED.
  --caller NAME=QUERY.vcf         A named query call set; repeat for each
                                  caller.  [required]
  -o, --output DIRECTORY          [required]
  --engine [normalized|haplotype|exact]
                                  [default: haplotype]
  --on-ref-mismatch [error|skip]
  --history FILE                  Append per-caller metrics to this history
                                  CSV.
  --label TEXT                    Label for the history row (e.g. a pipeline
                                  version).  [default: run]
  --help                          Show this message and exit.
```

</details>

## Annotations

The shared DuckDB reference store (genes, ClinVar, gnomAD). See [MCP and annotations](MCP.md).

### annotations load

Populate gene coordinates from a GENCODE GFF3.

```bash
vcfclick annotations load
```

<details><summary>options</summary>

```
Usage: vcfclick annotations load [OPTIONS]

  Populate the gene-coordinates table from a GENCODE GFF3.

  Required once after `pip install vcfclick` for the MCP server's
  `position_for_gene` tool to resolve symbols → coordinates. Re-run when
  GENCODE releases a new annotation version (yearly-ish).

Options:
  --gff FILE       Local GENCODE GFF3 (.gff3 or .gff3.gz). Downloads v45 from
                   EBI if omitted.
  --keep-existing  Don't truncate refseq_genes before loading (default:
                   replace).
  --help           Show this message and exit.
```

</details>

### annotations load-clinvar

Populate the ClinVar significance table.

```bash
vcfclick annotations load-clinvar
```

<details><summary>options</summary>

```
Usage: vcfclick annotations load-clinvar [OPTIONS]

  Populate the ClinVar significance table from the NCBI ClinVar VCF.

  Required for the MCP server's `clinvar_lookup` tool to return real
  significance calls. The NCBI VCF refreshes weekly; re-run monthly (or before
  any clinically-adjacent demo) to stay current. Bare numeric contigs are
  normalised to `chr`-prefixed during load so lookups against sample data
  (which uses chr-style) compose.

Options:
  --vcf FILE       Local ClinVar VCF (.vcf.gz). Downloads the current NCBI
                   weekly release if omitted.
  --keep-existing  Don't truncate clinvar_variants before loading (default:
                   replace).
  --help           Show this message and exit.
```

</details>

### annotations load-transcripts

Populate the transcript / exon / CDS tables from the same GENCODE GFF3, enabling
CDS-only filtering, MANE Select lookup, and splice-site distance. See [MCP.md](MCP.md).

```bash
vcfclick annotations load-transcripts
```

### annotations load-gnomad

Load gnomAD allele frequencies from a sites VCF.

```bash
vcfclick annotations load-gnomad gnomad.vcf.gz
```

<details><summary>options</summary>

```
Usage: vcfclick annotations load-gnomad [OPTIONS] VCF

  Load gnomAD allele frequencies from a gnomAD sites VCF.

  gnomAD is too large to bundle, so pass a VCF you supply — a region slice, an
  af-only file, or a per-chromosome sites VCF. A small region can be pulled
  with tabix-over-HTTPS from the public gnomAD bucket; see docs/MCP.md. Powers
  the `gnomad_lookup` MCP tool and the `db trio --gnomad-max-af` rarity
  filter.

Options:
  --replace  Truncate gnomad_af before loading (default: append, so several
             per-chromosome slices can be loaded incrementally).
  --help     Show this message and exit.
```

</details>

### annotations load-hgnc

Optional: previous symbols, aliases, HGNC ids and Entrez ids, so a gene list
written with retired symbols resolves (`db trio --genes`). Pass the HGNC
complete-set TSV you downloaded from genenames.org; nothing is fetched. Without
it, current symbols (any case) and Ensembl ids resolve.

```bash
vcfclick annotations load-hgnc hgnc_complete_set.txt
```

<details><summary>options</summary>

```
Usage: vcfclick annotations load-hgnc [OPTIONS] TSV

  Load HGNC previous symbols, aliases, HGNC ids and Entrez ids.

  Optional. Pass the HGNC complete-set TSV (hgnc_complete_set.txt from
  genenames.org; nothing is downloaded). Without it, a `db trio --genes` list
  resolves current symbols (any case) and Ensembl gene ids only; with it, a
  retired symbol, an alias, `HGNC:1234` or an Entrez id resolves to the
  current symbol. Needs `vcfclick annotations load` for the coordinates.

Options:
  --keep-existing  Don't truncate gene_aliases before loading (default:
                   replace).
  --help           Show this message and exit.
```

</details>

## Export & sharing

Move databases between machines and backends.

### db export

Write an ingestion back to VCF, whole or sliced by region, BED, gene, samples or a SQL condition. `.vcf.gz` output is BGZF, so `tabix -p vcf` indexes it. Genotypes are stored sparsely, so samples with no stored call are written `0/0` by default (`--absent-as nocall` for `./.`); ingestions made with `--keep-reference` keep true no-calls. Phase, PL/GL, and the difference between FILTER `PASS` and `.` are not recoverable and the header says so.

```bash
vcfclick db export cohort -o cohort.vcf.gz
vcfclick db export cohort -o brca1.vcf.gz --gene BRCA1 --samples S1,S2 --pass-only
vcfclick db export cohort -o - --where "info_AF < 0.01" --sites-only | bcftools view -H | wc -l
```

<details><summary>options</summary>

```
Usage: vcfclick db export [OPTIONS] NAME

  Export an ingestion to VCF, optionally sliced by region, gene, samples or
  SQL.

Options:
  -o, --out TEXT            Output path: .vcf.gz (BGZF, tabix-indexable),
                            .vcf, or - for stdout.  [required]
  --ingest-id TEXT          Ingestion to export (required when the database
                            has several).
  --samples TEXT            Comma-separated sample IDs to keep.
  --region TEXT             chr1, chr1:1000 or chr1:1000-2000 (repeatable).
  --regions-bed FILE        BED file of regions.
  --gene TEXT               Only variants in this gene (needs `vcfclick
                            annotations load`).
  --flank INTEGER           Extend --gene by this many bp each side.
                            [default: 0]
  --where TEXT              Extra SQL condition on the variants table, e.g.
                            "info_AF < 0.01".
  --pass-only               Only variants whose FILTER is PASS or unset.
  --sites-only              No genotypes: write the first eight columns only.
  --absent-as [ref|nocall]  How to write samples with no stored genotype
                            (ignored for --keep-reference ingestions, where
                            absent always means no-call).  [default: ref]
  --help                    Show this message and exit.
```

</details>

### db dump

Export all tables to Parquet files.

```bash
vcfclick db dump trio -o dump/
```

<details><summary>options</summary>

```
Usage: vcfclick db dump [OPTIONS] NAME

  Export all tables from a named database to Parquet files.

Options:
  --out DIRECTORY  Output directory (default: ./<name>-dump/).
  --help           Show this message and exit.
```

</details>

### db ingest-parquet

Ingest a `db dump` Parquet set into a database.

```bash
vcfclick db ingest-parquet trio dump/
```

<details><summary>options</summary>

```
Usage: vcfclick db ingest-parquet [OPTIONS] NAME DUMP_DIR

  Ingest a Parquet dump (produced by `db dump`) into NAME.

Options:
  --cohort TEXT     Cohort label to assign to the imported data.  [default:
                    default]
  --ingest-id TEXT  Stable upload identifier (UUID4 if omitted). Reuse to
                    replace prior data.
  --help            Show this message and exit.
```

</details>

### db push

Bundle a database as a portable tar.gz.

```bash
vcfclick db push trio -o trio.tar.gz
```

<details><summary>options</summary>

```
Usage: vcfclick db push [OPTIONS] NAME OUT_PATH

  Dump a database and bundle it as a portable tar.gz file.

Options:
  --help  Show this message and exit.
```

</details>

### db pull

Restore a database from a tar.gz bundle (file or HTTPS URL).

```bash
vcfclick db pull demo https://.../demo.tar.gz
```

<details><summary>options</summary>

```
Usage: vcfclick db pull [OPTIONS] NAME SOURCE

  Restore a database from a tar.gz bundle (HTTPS URL or local file).

Options:
  --help  Show this message and exit.
```

</details>

## Interactive UIs

Optional extras (`vcfclick[tui]`, `vcfclick[web]`).

### tui

Launch the terminal UI (Locus / Operations / SQL panes). See [TUI.md](TUI.md).

```bash
vcfclick tui --db trio
```

<details><summary>options</summary>

```
Usage: vcfclick tui [OPTIONS]

  Launch the optional terminal UI.

Options:
  --db TEXT  Database to open initially.
  --help     Show this message and exit.
```

</details>

### web

Launch the local browser UI (SQL, NL→SQL, trio, combine). See [WEB.md](WEB.md).

```bash
vcfclick web trio
```

<details><summary>options</summary>

```
Usage: vcfclick web [OPTIONS] NAME

  Launch the optional web UI for database NAME.

Options:
  --port INTEGER  Port to bind.  [default: 8765]
  --host TEXT     Interface to bind. Defaults to localhost; this is a single-
                  user local tool with no authentication — only change this if
                  you understand the exposure.  [default: 127.0.0.1]
  --no-browser    Do not open a browser automatically.
  --help          Show this message and exit.
```

</details>
