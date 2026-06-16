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
  order), `GT` per sample, and the `set=` INFO field.
- A sample absent from every input at a site is `./.`, never silently
  `0/0` — `combine` does not assert hom-reference where no input has a
  record.
- v1 carries **GT only**; per-sample FORMAT fields (`GQ`, `DP`, `AD`)
  from the priority source are not yet propagated. If you need those for
  downstream quality gating, that passthrough is planned.
- The site union is held in memory, which suits the typical "a few call
  sets" use; whole-genome many-caller merges of very large cohorts are
  not the target.
