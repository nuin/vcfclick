# Combining call sets (`vcfclick combine`)

`vcfclick combine` merges multiple VCF **call sets** into one, recording
where each variant came from. It is a native reimplementation of GATK3's
`CombineVariants` — a tool GATK4 dropped and never replaced, and which
`bcftools` does not cover.

## `combine` vs `merge`

vcfclick has two ways to put several VCFs together; they solve different
problems.

| | `vcfclick merge` | `vcfclick combine` |
|---|---|---|
| Inputs | per-sample VCFs, **disjoint** samples | call sets that may **share** samples |
| Use case | assemble a trio / cohort joint VCF | two callers, or pre/post-filter, of the same cohort |
| Same sample in two inputs | rejected (must be disjoint) | resolved by **priority** |
| Engine | wraps `bcftools merge` | native (cyvcf2 read + write) |
| Provenance | — | `set=` INFO field per record |

Reach for `merge` to build one multi-sample VCF from separate samples;
reach for `combine` to reconcile multiple *callings* of the same samples.

## Usage

```bash
# Two callers over the same cohort; keep all sites, annotate origin.
vcfclick combine gatk.vcf.gz deepvariant.vcf.gz -o all.vcf.gz

# Consensus: keep only sites called by at least 2 of 3 callers.
vcfclick combine gatk.vcf.gz dv.vcf.gz strelka.vcf.gz \
    -o consensus.vcf.gz --min-callsets 2

# Name the sets explicitly (default names come from the filenames).
vcfclick combine a.vcf.gz b.vcf.gz -o out.vcf.gz --name gatk --name dv
```

Inputs must be on the same reference and **decomposed** (one ALT per
record) — the same precondition as ingest. Decompose first with
`bcftools norm -m -` if needed.

A plain `.vcf` output is fully native (no external tools). A `.vcf.gz`
output is BGZF-compressed and tabix-indexed via htslib's `bgzip` /
`tabix`, so it is region-queryable and usable by parallel ingest like
every other vcfclick VCF — install htslib (`brew install htslib` /
`conda install -c bioconda htslib`) or write a plain `.vcf`.

## What it does

**Union of sites.** Every `(chrom, pos, ref, alt)` across all inputs
appears once in the output.

**`set=` provenance.** Each output record's INFO carries `set=`:

- `set=Intersection` — present in every input;
- `set=<name>` — present in one input only;
- `set=<a>-<b>` — present in those inputs (priority order), but not all.

**PRIORITY genotype resolution.** Input order is priority, highest
first. When a sample appears in more than one input at a site, its
genotype is taken from the **first input with a non-missing call**. A
`./.` in a higher-priority input does not win — the next input fills it.

**Consensus filter.** `--min-callsets N` keeps only sites present in at
least N inputs, for "agreed by ≥N callers" workflows.

## Output and limits

- Output is a fresh VCF with the union of input samples (priority
  order), `GT` per sample, the `set=` INFO field, and the `GQ`/`DP`/`AD`
  FORMAT fields.
- **Quality travels with the genotype.** When a sample's genotype comes
  from a given input, that same input's `GQ`/`DP`/`AD` come with it — so
  the depth and allele balance describe the call that was actually kept,
  and a combined VCF feeds straight into the trio quality gates (`db
  trio --min-gq/--min-dp/--min-ab`). Output FORMAT lists only the
  passthrough fields that some input actually carries (just `GT` if none
  do).
- A sample absent from every input at a site is `./.` with `.` for every
  FORMAT field, never silently `0/0` or a borrowed quality value —
  `combine` does not assert hom-reference where no input has a record.
- The site union is held in memory, which suits the typical "a few call
  sets" use; whole-genome many-caller merges of very large cohorts are
  not the target.

## Comparison with GATK3 CombineVariants

`combine` is a reimplementation of GATK3's `CombineVariants`. The two
produce the same site union, the same `set=` provenance, and the same
PRIORITIZE genotype resolution. The equivalent invocations:

```bash
# GATK3 (removed in GATK4) — needs Java 8 and a reference
java -jar GenomeAnalysisTK.jar -T CombineVariants -R ref.fasta \
    -V:first first.vcf.gz -V:second second.vcf.gz \
    -o combined.vcf \
    -genotypeMergeOptions PRIORITIZE -priority first,second

# vcfclick — no Java, no reference; input order is the priority list
vcfclick combine first.vcf second.vcf -o combined.vcf \
    --name first --name second
```

On two call sets that share sample `S2` (where `S2` is `0/1` in `first`
and `1/1` in `second` at chr1:202), both tools produce the same core —
`set=` and genotypes, with PRIORITIZE keeping `first`'s call for `S2`:

| pos | `set=` | S1 | S2 | S3 |
|---|---|---|---|---|
| 101 | `first` | 0/1 | 0/1 | ./. |
| 202 | `Intersection` | 1/1 | **0/1** | 0/1 |
| 303 | `Intersection` | 0/1 | 1/1 | 0/0 |
| 404 | `second` | ./. | 0/0 | 1/1 |

The `-V:name` binding GATK3 uses for the `set=` label is `--name` in
vcfclick; `-priority a,b` is just the input order. Differences are
representational, not behavioural: GATK3 also recomputes `AC/AF/AN`
INFO, whereas vcfclick emits `set=` plus the `GQ/DP/AD` passthrough.

This equivalence is pinned by a test against **real** GATK 3.8-1 output:
`tests/test_combine_gatk3.py` runs `combine` on the checked-in call sets
and asserts it reproduces the frozen GATK3 golden file
(`tests/fixtures/gatk3/combine_priority.gatk3.vcf`, whose header records
the exact generating command). GATK3 needs Java 8 and is not run in CI;
the golden output is captured once and asserted against.
