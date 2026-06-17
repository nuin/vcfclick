# Trio / family analysis

vcfclick can report candidate variants under Mendelian inheritance
models for a trio (affected child + two parents), with the genotype
quality gates and population-frequency rarity filter that rare-disease
analysis uses in practice.

This is **candidate filtering, not variant calling**. It works
downstream of your caller, over an already-genotyped cohort. Quality
gating is only as strong as the FORMAT fields your VCF carries
(`gq`, `dp`, `ad`); where they are absent the gates pass through.

## The pipeline

```bash
# 1. Trios usually arrive as three per-sample VCFs. Combine them into
#    one joint multi-sample VCF (decomposed, ingest-ready).
vcfclick merge proband.vcf.gz father.vcf.gz mother.vcf.gz -o trio.vcf.gz

# 2. Ingest as ONE cohort. Use --keep-reference so a parent's
#    confident 0/0 calls are stored (see "Why --keep-reference" below).
vcfclick db create fam1
vcfclick db ingest fam1 trio.vcf.gz \
    --cohort trio --ingest-id fam1 --keep-reference

# 3. Load the pedigree (standard PED/FAM).
vcfclick db ped fam1 fam1.ped

# 4. Report candidates.
vcfclick db trio fam1 --proband PROBAND
vcfclick db trio fam1 --proband PROBAND --category denovo
```

If your trio was already joint-called into one multi-sample VCF, skip
step 1 and ingest it directly.

## The PED file

Standard six-column PED/FAM, whitespace-delimited:

```
# family  individual  paternal  maternal  sex  phenotype
FAM1  CHILD   FATHER  MOTHER  1  2
FAM1  FATHER  0       0       1  1
FAM1  MOTHER  0       0       2  1
```

- paternal/maternal `0` = founder (no parent in the file)
- sex: 1 = male, 2 = female
- phenotype: 1 = unaffected, 2 = affected

`db ped` validates that every individual and named parent is a real
sample in the cohort, and re-loading replaces the prior pedigree.

## Inheritance models

`db trio` reports, per model, candidate variants where the proband and
parent genotypes fit the pattern. Genotype encoding: `gt` 0 = hom-ref,
1 = het, 2 = hom-alt.

| Model | Pattern | Needs `--keep-reference`? |
|---|---|---|
| `denovo` | proband carries; both parents provably hom-ref | yes |
| `recessive` | proband hom-alt; both parents heterozygous carriers | no |
| `dominant` | proband het; exactly one parent carries, other hom-ref | yes |
| `comphet` | two rare proband hets in the same gene, one from each parent | yes |

`--category all` (default) prints per-model counts; a specific
`--category` prints the variants.

### Compound heterozygous

`--category comphet` finds genes where the proband carries **two** rare
heterozygous variants in *trans* — one inherited from the father, one
from the mother — so both copies of the gene are hit (a recessive
mechanism that per-variant filters miss). Each candidate variant is a
rare proband het with a clear parent-of-origin (one parent carries, the
other is provably hom-ref, hence `--keep-reference`); a gene is reported
only when it has both a paternal and a maternal one. Results are grouped
**per gene**, not per variant.

Because genes live in the annotation store (separate from the cohort),
this needs gene coordinates loaded first:

```bash
vcfclick annotations load          # GENCODE gene coordinates
vcfclick db trio fam1 --proband CHILD --category comphet
```

Parent-of-origin here is genotype-based (no read-backed phasing), so it
is honest candidate filtering: it cannot distinguish true *trans* from
two variants that happen to be on the same parental copy. The annotation
store path can be overridden with `VCFCLICK_ANNOTATIONS_DB`.

## Quality gates

Defaults follow common rare-disease practice and are all tunable:

| Flag | Default | Meaning |
|---|---|---|
| `--min-gq` | 20 | genotype quality (≈99% accurate) |
| `--min-dp` | 10 | read depth |
| `--min-ab` / `--max-ab` | 0.25 / 0.75 | het allele balance `ad_alt/(ad_ref+ad_alt)`, expected ≈0.5 |
| `--max-af` | 0.01 | keep variants with population `info_AF` at or below this |

Gates are lenient on NULL: a VCF lacking `gq`/`dp`/`ad` is not silently
zeroed out — the gate simply doesn't filter on the missing field.

## Why `--keep-reference` (and the de-novo caveat)

vcfclick's `genotypes` table is **sparse** — only non-reference calls
are stored, so a sample absent at a site is `0/0` (reference) **or**
`./.` (no-call), indistinguishable. De novo means "the child carries a
variant neither parent has," which requires the parents to be
*confidently* hom-reference. Absence alone can't prove that: a no-call
parent is not evidence.

`--keep-reference` additionally stores confident `0/0` calls (gt=0) at
variant sites, while still dropping no-calls. Then de-novo's join
(`f.gt = 0 AND m.gt = 0`) matches only when both parents have a stored
reference call — a no-call parent has no row, so that site is correctly
**excluded** from de novo rather than falsely reported.

Without `--keep-reference`, `db trio` still runs but de-novo and
dominant return nothing (no stored parent reference rows) and print a
note to re-ingest with the flag. Recessive works either way.

Rigorous de-novo *calling* (PL-based Bayesian posteriors using genotype
likelihoods) is future work; this is honest candidate filtering.

## Validation against the GIAB benchmark trio

The inheritance models are tested against **real** genotypes from the
Genome in a Bottle Ashkenazi trio — HG002 (son), HG003 (father), HG004
(mother) — taken from the NIST v4.2.1 GRCh38 small-variant benchmark, not
hand-built fixtures. In a CFTR slice the benchmark genotypes reproduce
every model: two recessive sites (proband hom-alt, both parents het), a
paternal-origin and maternal-origin het that together make CFTR a
compound-het candidate gene, and no de novo. The golden fixture
(`tests/fixtures/giab/cftr_trio.vcf.gz`) records its GIAB source URLs in
the header.

## Natural language

The MCP server teaches an LLM the pedigree table, the inheritance
models, and the de-novo sparse caveat, so you can ask trio questions in
English and get the SQL back to inspect. See [MCP](MCP.md).
