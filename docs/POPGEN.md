# Population genetics

vcfclick can compute the standard population-genetics summaries of a
cohort — Watterson's θ, nucleotide diversity π, Tajima's D, normalised
Fay & Wu's H, observed/expected heterozygosity and F, the site-frequency
spectrum (SFS) and Hudson's F_ST — directly from a database, per
population. The same tables travel in Parquet bundles, so a browser tool
reading a bundle (for example PopGeneJS on DuckDB-Wasm) sees the same
counts.

Three pieces make this possible:

1. **Called-genotype accounting** — at ingest every site records how many
   samples and alleles were actually called (`variants.n_called`,
   `an_called`, `ac_called`), and every fully missing call is recorded in
   `missing_genotypes`. Without this, a no-call (`./.`) is
   indistinguishable from `0/0` in the sparse `genotypes` table.
2. **A population panel** — `vcfclick db panel` loads sample → population
   labels into `populations`.
3. **Ancestral alleles** — `INFO/AA` (already stored as `variants.info_AA`)
   is normalised to polarise sites for unfolded statistics.

## Quick start

```bash
vcfclick db create kg
vcfclick db ingest kg chr22.vcf.gz --ingest-id phase3
vcfclick db panel  kg integrated_call_samples_v3.20130502.ALL.panel

vcfclick db popgen summary kg                        # per population
vcfclick db popgen summary kg --by super_population
vcfclick db popgen sfs     kg --project 20 --format json
vcfclick db popgen fst     kg --by super_population --window 1000000
vcfclick db popgen windows kg --window 100000 --step 50000 --fst > win.tsv
```

`--format table` (default, except `windows` which defaults to `tsv`),
`tsv` (plot-ready, `NA` for undefined values) or `json` (everything,
including the site report and warnings).

## Loading a panel

```bash
vcfclick db panel NAME FILE [--ingest-id ID] [--sample-col C] [--pop-col C]
                            [--super-pop-col C] [--sex-col C]
```

- The 1000 Genomes panel format (`sample pop super_pop gender`,
  tab-separated, header row, trailing empty fields allowed) works as is.
- Any other delimited file with a header works too: tab, comma (a `.csv`
  file or a comma-only header) or whitespace. Columns are found by name
  (`sample`/`sample_id`/`IID`, `pop`/`population`,
  `super_pop`/`super_population`, `gender`/`sex`) or named with the
  `--*-col` options. Only the sample and population columns are required.
- `gender`/`sex` is normalised to `male`/`female`/NULL (`male`, `m`, `1`;
  `female`, `f`, `2`; anything else NULL).
- Without `--ingest-id` the panel applies to every ingestion that
  contains at least one of its samples.
- **A panel replaces the whole labelling of each ingestion it applies
  to.** Re-loading the same file is idempotent; loading a file that lists
  fewer samples leaves the samples it omits unlabelled (no stale labels).
- Re-ingesting a VCF under the same `ingest_id` keeps the labels of the
  samples that are still in it and **removes the labels of samples that
  are not** (a warning lists them). Statistics only ever count samples
  currently in the ingestion, whatever the `populations` table holds.
- The command reports panel samples that are not in the database and
  database samples that are missing from the panel. Unlabelled samples
  are excluded from per-population statistics (and counted in the
  report).

## Scope and site filters

Every subcommand analyses **one ingestion**: sample identity is
`(ingest_id, sample_id)` and two ingestions do not share a site list.
`--ingest-id` is required when the database holds several.

| Option | Default | Meaning |
|---|---|---|
| `--region chr:start-end` | whole ingestion | repeatable; `chr1` and `1` both match |
| `--gene SYMBOL` | — | gene span from the annotation store (`vcfclick annotations load`) |
| `--by` | `population` | `super_population`, or `all` for one cohort-wide group; falls back to `all` (with a warning) when no panel is loaded |
| `--include-indels` | off | biallelic SNVs only by default |
| `--pass-only / --all-filters` | pass-only | FILTER is `PASS` or missing |
| `--min-call-rate` | 0.9 | a site must reach it in **every** group |
| `--maf` | 0 | cohort-wide minor-allele frequency (all samples of the ingestion) |
| `--ancestral` | `aa` if any autosomal site in scope has INFO/AA, else `none` | see below |
| `--allow-untracked` | off | compute on an ingestion without missing-call tracking (see [older databases](#older-databases)) |

Filters are applied in this order, and the number of sites each one
drops is reported (`sites.dropped` in JSON):

1. `non_autosomal` — **autosomes only.** A chromosome is classified by
   the part of its name before the first `_` (so GRCh38 `chrX_..._alt`
   follows X and `chr1_..._random` follows 1), with or without `chr`,
   case-insensitively. Excluded: X, Y, XY, W, Z, M, MT, PAR1/PAR2,
   unplaced and decoy sequence (`chrUn_*`, `chrEBV`, `hs37d5`, `HLA-*`,
   GRCh37 `GL0*`, RefSeq `NT_`/`NW_` scaffolds) and the human RefSeq
   accessions NC_000023 (X), NC_000024 (Y) and NC_012920 (MT);
   NC_000001–NC_000022 are autosomes. **Any other naming is the user's
   responsibility**: an unrecognised name counts as an autosome, so
   restrict with `--region` when a reference uses other names for sex
   chromosomes (e.g. a non-human RefSeq assembly). `--include-sex-chroms`
   is refused in this version: hemizygous calls need per-sample ploidy,
   which is not modelled yet. A `--region` on an excluded chromosome is
   an error rather than an empty result.
2. `not_biallelic` — any position with more than one record (a
   multi-allelic site split by `bcftools norm -m -`).
3. `variant_type` — not a SNV (or, with `--include-indels`, not a plain
   sequence change; symbolic alleles and `*` are never used).
4. `filter` — FILTER is not PASS/missing (with `--pass-only`).
5. `inexact_group_counts` — see [missing data](#missing-data-and-exactness).
6. `call_rate` — called individuals / group size below `--min-call-rate`
   in any group.
7. `maf` — only when `--maf` > 0.

The retained site set is **shared by all groups**, so per-population
statistics and pairwise F_ST are computed over exactly the same sites.
Monomorphic sites are kept (unless `--maf` > 0): they matter for any
per-site normalisation, and they are where most of the genome is.

## Missing data and exactness

`genotypes` stores only called non-reference genotypes. A group's called
haplotypes at a diploid site are derived as

    n = 2 × (group size − samples of the group in missing_genotypes)
    k = Σ gt over the group's rows in genotypes

That derivation is exact for diploid genotypes that are either fully
called or fully missing. It is **not** exact for:

- partially missing calls (`./1`, `0/.`) — `genotypes` stores what
  cyvcf2's `gt_types` reports (`./1` is a het, `0/.` hom-ref), and the
  sample is not "fully missing";
- haploid calls on an autosome (a haploid `1` is stored as `gt = 2`);
- polyploid calls.

The per-site counts on `variants` *are* exact in all of these cases
(they are computed from the allele arrays: a haploid call contributes one
allele, `./1` contributes its one called allele). So every site is
checked: if the cohort-wide derived `n`/`k` disagree with `an_called`/
`ac_called`, the site is dropped as `inexact_group_counts` and a warning
says how many. Statistics are therefore always computed from exact
counts — at the price of dropping those (rare) sites.

### Older databases

Databases ingested before this feature have NULL `n_called` / `an_called`
/ `ac_called` (or lack the columns) and no `missing_genotypes` rows. On
such an ingestion every subcommand **refuses** with an explanation,
unless `--allow-untracked` is given: then it warns once that missing
calls are treated as `0/0`, skips the exactness check, and reports
`"missing_data_tracked": false` in JSON. Re-ingest the VCF to get exact
counts. Ingesting into an older database adds the new columns and tables
in place (`[storage] upgraded schema: ...`, under a per-database lock and
idempotent, so concurrent ingests are safe); older dumps and bundles still
load (the counts are NULL).

An ingestion loaded with `--no-record-missing` has exact per-site counts
but no per-sample missing rows. It also needs `--allow-untracked`; sites
with missing calls then fail the exactness check and are dropped, with a
warning and `missing_data_tracked: false`.

## Ancestral alleles and polarisation

`info_AA` is stored verbatim. It is normalised by taking the text before
the first `|` (the 1000 Genomes form `AA=G|||`); upper case is a
high-confidence call, lower case low confidence (Ensembl EPO
convention); `.`, `-`, `?`, `N` or empty mean unknown.

| `--ancestral` | Polarised when |
|---|---|
| `aa` | AA (either case) equals REF (ancestral = REF) or ALT (derived = REF: counts flipped) |
| `aa-high` | as `aa`, upper-case AA only |
| `ref` | always: REF is assumed ancestral |
| `none` | never: folded statistics only |

Unpolarised sites still count for every folded statistic (S, θ_W, π,
Tajima's D, folded SFS, F_ST, heterozygosity); unfolded statistics
(unfolded SFS, Fay & Wu's H) use the polarised sites only. The report
gives the number of polarised and unpolarised sites.

## Statistics

Notation: at a site, a group has `n` called haplotypes of which `k`
carry ALT (or, once polarised, the derived allele). `a_n = Σ_{i=1}^{n−1}
1/i`, `b_n = Σ_{i=1}^{n−1} 1/i²`.

### summary

Per group:

- `sites`, `segregating_sites` (S: 0 < k < n in the group).
- **π** = Σ_sites 2k(n−k) / (n(n−1)), each site with its **own n**
  (Tajima 1983). `pi_per_site` divides by the retained sites.
- **θ_W** = Σ_segregating 1 / a_{n_i}, each site with its **own n**
  (Watterson 1975 generalised to missing data, as in Ferretti et al.
  2012). With no missing data this is exactly S / a_n.
- **Tajima's D** (Tajima 1989) and **Fay & Wu's H** need one sample size,
  so they are computed on the SFS **projected** to `projection_n`
  haplotypes (below): D = (θ_π − θ_W) / sqrt(e1·S + e2·S(S−1)), with the
  usual a1, a2, b1, b2, c1, c2, e1, e2. Projection leaves π unchanged
  (π is the probability two distinct haplotypes differ, which subsampling
  preserves), so θ_π in D equals the reported π over the sites used.
  With no missing data the projection is the identity and D uses exactly
  the reported π and θ_W.
  The projected S, θ_π, θ_L and θ_H are computed in closed form from the
  hypergeometric moments (S = Σ 1 − P(J=0) − P(J=m), Σ_{j<m} j·P(J=j) =
  mk/n − m·P(J=m), and the second moment for θ_H), which equals building
  the projected spectrum but costs O(1) per site instead of O(m); `sfs`
  builds the full spectrum. D is undefined (`null`/`NA`) for S = 0 or
  n < 4.
- **Fay & Wu's H**, normalised (Zeng, Fu, Shi & Wu 2006, eq. 11), on the
  polarised sites only:
  H = (θ_π − θ_L) / sqrt(Var), θ_L = Σ i ξ_i / (n−1),
  Var = (n−2)/(6(n−1)) θ + [18n²(3n+2)b_{n+1} − (88n³+9n²−13n+6)] /
  (9n(n−1)²) θ², θ = S/a_n, θ² = S(S−1)/(a_n² + b_n). JSON also carries
  the unnormalised `fay_wu_h_raw` = θ_π − θ_H (Fay & Wu 2000).
- **Ho** = mean over sites of heterozygous / called individuals;
  **He** = mean over sites of Nei's (1978) unbiased 2k(n−k)/(n(n−1));
  **F** = 1 − Ho/He (ratio of averages).

### sfs

Per group, the **folded** SFS (always) and the **unfolded** SFS (when
polarised), projected to `projection_n` haplotypes by the hypergeometric
distribution: a site with k of n contributes C(k,j)·C(n−k,m−j)/C(n,m) to
class j (Marth et al. 2004; Gutenkunst et al. 2009 — the dadi
projection). `--project N` sets m; the default is the smallest number of
called haplotypes at any retained site in that group, so no site is
lost. If that minimum is below 4 (for example a nearly
uncalled site kept by `--min-call-rate 0`), a warning says which group
and why D, H (n < 4) or the spectrum (n < 2) cannot be computed. Sites with fewer than m called haplotypes are dropped and counted
(`sites_dropped`). Spectra include the monomorphic classes (j = 0 and,
unfolded, j = m) and are expected counts, so they can be fractional.

### fst

Pairwise **Hudson's F_ST** in the formulation of Bhatia, Patterson,
Sankararaman & Price (2013), per site

    N = (p1 − p2)² − p1(1−p1)/(n1−1) − p2(1−p2)/(n2−1)
    D = p1(1−p2) + p2(1−p1)

combined as a **ratio of averages**, Σ N / Σ D (Bhatia et al. recommend
this over averaging per-site ratios, which is dominated by rare
variants). Small or closely related samples can give slightly negative
values — that is the sampling-bias correction, not an error.
`--window W [--step S]` adds per-window values.

### windows

Sliding windows of `--window` bp every `--step` bp (default: the window
size) along each `--region` (or each chromosome from position 1 to its
last retained site). Per window and group: S, θ_W, π (sums over the
window's retained sites), `*_per_bp` versions, and Tajima's D on the SFS
projected to the group's genome-wide `projection_n`. `--fst` adds
pairwise F_ST columns. Output is one row per window (wide), TSV by
default. Tiling stops at the first window that reaches the end of the
span.

**Per-site and per-bp values.** A VCF normally lists only variable
records, so "per site" means per retained VCF record and "per bp"
divides by the window length assuming every position without a record
is callable and monomorphic. Neither is a per-callable-base estimate:
for that, use an all-sites VCF (with invariant records) or a
callability mask, as in pixy (Korunes & Samuk 2021). The totals (θ_W,
π) and the ratio statistics (D, H, F_ST) do not depend on this.

## Output conventions

Chromosome names in every output are the names **as stored** in the
database (`--region chr1:...` on a database that stores `1` reports `1`).
JSON carries every warning; table/TSV print them on stderr.

## Implementation and scale

Per chromosome, one SQL statement each reads the sites (`variants`), the
ALT dosage and heterozygote counts (`genotypes`, scanned once, grouped by
site and panel label with a LEFT JOIN so unlabelled samples form one
extra bucket; cohort totals are the sum over buckets) and the missing
calls (`missing_genotypes`). Results come back as Arrow. Each chromosome
is filtered before the next is read, and only retained sites are kept,
as 16-bit counts per group (32-bit for groups of more than 32,767
samples), so memory is one chromosome's raw counts
plus about 6 bytes × retained sites × groups. Windows are binary searches
in each chromosome's position-sorted slice.

## Choices made (and why)

| Choice | Reason |
|---|---|
| One ingestion per run | sample identity is `(ingest_id, sample_id)`; ingestions don't share sites |
| Autosomes only; sex chromosomes refused | correct X/Y needs per-sample ploidy; refusing beats a silently wrong diploid assumption |
| Shared site set across groups | groups (and F_ST) are compared over the same sites |
| Drop sites with inexact per-group counts | every statistic comes from exact counts; the drop is reported |
| π and θ_W with per-site n | the standard unbiased treatment of missing data |
| D and H on the projected SFS | their variances assume a single n; projection is the standard (dadi/moments) way to get one without discarding data |
| Default projection = minimum called haplotypes | uses every retained site |
| Hudson F_ST, ratio of averages | Bhatia et al. 2013: no dependence on sample-size ratio, robust to rare variants |
| He unbiased (n/(n−1)) | equals per-site π, so Ho, He and π are on one scale |
| MAF filter cohort-wide | one site set for every group (filtering per group would bias F_ST) |

## Validation

- Every estimator is unit-tested against values worked out by hand with
  exact fractions (`tests/test_popgen_estimators.py`): a_n, the Tajima
  constants for n = 10, D and normalised H worked examples, per-site π
  and θ_W with variable n, the projection, folding, Hudson's N/D and the
  ratio of averages, Ho/He/F.
- Tajima's c1/c2 and Zeng's variance of θ_π − θ_L are re-derived exactly
  from Fu's (1995) covariances of the neutral SFS for several n.
- End to end (`tests/test_popgen_cli.py`), on a fixture with three
  populations, missing (`./.`), partial (`./1`), haploid and phased
  calls, polarised and unpolarised sites: every subcommand runs on both
  chDB and DuckDB and must give byte-identical JSON; values are compared
  with an independent computation straight from the VCF text, with
  scikit-allel 1.3.13 (π, θ_W, Tajima's D, Ho, Hudson F_ST — identical to
  the 12 digits recorded, where the definitions coincide) and with dadi 2.4.4 (projected
  folded and unfolded SFS and Tajima's D with missing data — identical
  to 12 digits). Neither tool is a dependency; the numbers are recorded
  as literals.

## References

- Bhatia G, Patterson N, Sankararaman S, Price AL (2013). Estimating and
  interpreting F_ST: the impact of rare variants. *Genome Res* 23:1514.
- Fay JC, Wu C-I (2000). Hitchhiking under positive Darwinian selection.
  *Genetics* 155:1405.
- Ferretti L, Raineri E, Ramos-Onsins S (2012). Neutrality tests for
  sequences with missing data. *Genetics* 191:1397.
- Fu Y-X (1995). Statistical properties of segregating sites. *Theor Popul
  Biol* 48:172.
- Gutenkunst RN et al. (2009). Inferring the joint demographic history of
  multiple populations from multidimensional SNP frequency data. *PLoS
  Genet* 5:e1000695.
- Hudson RR, Slatkin M, Maddison WP (1992). Estimation of levels of gene
  flow from DNA sequence data. *Genetics* 132:583.
- Korunes KL, Samuk K (2021). pixy: unbiased estimation of nucleotide
  diversity and divergence in the presence of missing data. *Mol Ecol
  Resour* 21:1359.
- Marth GT et al. (2004). The allele frequency spectrum in genome-wide
  human variation data reveals signals of differential demographic
  history in three large world populations. *Genetics* 166:351.
- Nei M (1978). Estimation of average heterozygosity and genetic distance
  from a small number of individuals. *Genetics* 89:583.
- Tajima F (1983). Evolutionary relationship of DNA sequences in finite
  populations. *Genetics* 105:437.
- Tajima F (1989). Statistical method for testing the neutral mutation
  hypothesis by DNA polymorphism. *Genetics* 123:585.
- Watterson GA (1975). On the number of segregating sites in genetical
  models without recombination. *Theor Popul Biol* 7:256.
- Zeng K, Fu Y-X, Shi S, Wu C-I (2006). Statistical tests for detecting
  positive selection by utilizing high-frequency variants. *Genetics*
  174:1431.
