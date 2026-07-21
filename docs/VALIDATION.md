# Trio analysis validation

vcfclick's trio inheritance models are validated against **real published
benchmark data** with an **independent ground truth**, not only
hand-built fixtures — and cross-checked against an established external
tool. All of this runs in CI from small checked-in fixtures
(`tests/fixtures/giab/`); this page records the methodology and the
genome-wide numbers behind it.

## Data and ground truth

- **Trio**: the Genome in a Bottle (GIAB) Ashkenazi trio — **HG002**
  (son / proband), **HG003** (father), **HG004** (mother).
- **Calls**: the NIST **v4.2.1** GRCh38 small-variant benchmark VCFs,
  fetched directly from the GIAB FTP (the fixture headers record the
  exact URLs).
- **Ground truth for parent hom-reference**: each sample's v4.2.1
  high-confidence BED (`*_benchmark_noinconsistent.bed`). Inside its BED
  a sample with no variant call is a **confident hom-reference (0/0)**;
  outside it, absence is a **no-call (./.)**. This is the independent
  signal — vcfclick never sees the BED.
- **External tool**: `bcftools +mendelian2` (htslib), the standard
  Mendelian-consistency checker.

## 1. De-novo recovery and the `--keep-reference` defensibility claim

De novo means the child carries a variant **neither parent has** — which
requires the parents to be *provably* hom-reference, not merely absent
from a sparse genotype table. We tested this on real data.

In `chr20:1,000,000–6,000,000` (5 Mb), HG002 carries **50** variants
absent from both parents' benchmark calls — every one a candidate de novo
under the naive "neither parent carries" rule. Classifying each by the
parents' BEDs:

| | sites |
|---|---:|
| HG002-only ("naive de novo") | 50 |
| both parents **in-BED** (confident hom-ref → real de novo) | **4** |
| ≥1 parent **out-of-BED** (no-call → unsupported) | 46 |

The naive rule over-calls de novo by **~12×**. The 50/46 split is from a
one-time scan, reproducible by intersecting the HG002 variant calls with
the parents' BEDs:

```bash
B=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio
# HG002-only SNVs in the region, then keep those inside both parents' BEDs
tabix $B/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz chr20:1000000-6000000
tabix $B/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz chr20:1000000-6000000
# parent confident-hom-ref regions:
curl -s $B/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

The **checked-in CI test** exercises a representative 7-site subset.
`tests/fixtures/giab/denovo_trio.vcf.gz` encodes the 4 confident + 3
no-call sites with real HG002 `GT:GQ:DP:AD`. On that fixture:

- `vcfclick db trio --category denovo` **(ingested `--keep-reference`)**
  returns **exactly the 4** confident sites and **excludes the 3**
  no-call ones.
- The **same data without `--keep-reference`** returns **0** — vcfclick
  refuses to call de novo when it cannot prove a parent is 0/0.

## 2. External cross-checks — `bcftools +mendelian2` and `slivar`

Two independent, established tools were run on the same fixtures and
agree with vcfclick **site-for-site**. The `bcftools +mendelian2` check
is asserted in CI; the `slivar` comparison is a documented one-time run
(slivar ships a Linux binary only), captured below.

`bcftools +mendelian2 -m e` flags as Mendelian-erroneous **exactly the 4
sites** vcfclick calls de novo, and — like vcfclick — does not flag the 3
no-call sites (a missing parent genotype is undecidable). Asserted in CI
(`test_giab_denovo_matches_bcftools_mendelian`, skipped where bcftools is
absent).

[slivar](https://github.com/brentp/slivar), the de-facto standard for
rare-disease trio filtering, was run on the same fixtures:

```bash
# de novo — reports HG002: 4, exactly the confident sites
slivar expr --vcf denovo_trio.vcf.gz --ped denovo_trio.ped --pass-only \
    --trio "denovo:kid.het && dad.hom_ref && mom.hom_ref" -o out.vcf

# recessive — reports HG002: 2 (chr7:117,488,888 and 117,507,446)
slivar expr --vcf cftr_trio.vcf.gz --ped cftr_trio.ped --pass-only \
    --trio "recessive:kid.hom_alt && dad.het && mom.het" -o out.vcf
```

slivar's de-novo expression flags the same **4** confident sites and
excludes the 3 no-call ones (its `dad.hom_ref` / `mom.hom_ref` is false
for a `./.` parent — the same defensibility logic), and its recessive
expression flags the same **2** sites. So vcfclick, bcftools, and slivar
produce identical trio calls on real benchmark data. (slivar ships a
Linux binary only, so this comparison is documented rather than run in
CI; the bcftools cross-check above is the CI-runnable one.)

## 3. Recessive / dominant / compound-het on real genotypes

`tests/fixtures/giab/cftr_trio.vcf.gz` holds real GIAB trio genotypes at
CFTR sites. The models reproduce on real data:

- **recessive** — chr7:117,488,888 and 117,507,446: HG002 hom-alt, both
  parents heterozygous carriers.
- **dominant** — three single-origin hets (proband het, exactly one
  parent carrying).
- **compound het** — **CFTR** is a genuine candidate gene: a
  paternal-origin het (chr7:117,489,437) plus maternal-origin hets, in
  trans.

## Reproducing

```bash
# de-novo cross-check
bcftools +mendelian2 -P tests/fixtures/giab/denovo_trio.ped \
    tests/fixtures/giab/denovo_trio.vcf.gz -m e        # → the 4 de-novo sites

# vcfclick on the same data
vcfclick db create dn
vcfclick db ingest dn tests/fixtures/giab/denovo_trio.vcf.gz \
    --cohort fam --ingest-id g1 --serial --keep-reference
vcfclick db ped  dn tests/fixtures/giab/denovo_trio.ped
vcfclick db trio dn --proband HG002 --category denovo --max-af 1.0
```

The full validation lives in `tests/test_trio.py`
(`test_giab_*`). GIAB data: Zook et al., *Sci. Data* (2016) and the
NIST/GIAB v4.2.1 benchmark; see
<https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/>.

# Benchmarking (`vcfclick benchmark`)

`vcfclick benchmark` compares a query VCF against a truth VCF over a
reference FASTA, restricted to a confident-region BED, and reports
TP/FP/FN with precision/recall/F1 in a GA4GH-shaped `summary.csv`.

Two engines:

- **`normalized`** — reference normalization (left-align + trim + optional
  MNP decomposition) and a genotype-aware, allele/locus-keyed match.
  Resolves representation-shifted indels but not cross-representation
  rewrites (an MNP vs two SNPs), so its INDEL numbers are a conservative
  lower bound.
- **`haplotype`** — adds a vcfeval/xcmp-style local-haplotype pass that
  rescues representation-different but sequence-equivalent calls.

Provenance is explicit: `summary.csv` is the strict canonical shape; the
engine and inputs are stamped in `run_meta.json`, `vcfclick.summary.csv`
(which carries an `Engine` column), and the HTML report.

```bash
vcfclick benchmark \
  --truth truth.vcf.gz --query calls.vcf.gz \
  --ref genome.fa --regions confident.bed -o report/ --engine haplotype
```

## Conformance against hap.py

The `haplotype` engine's **matching** was checked directly against real
hap.py (Illumina `hap.py` v0.3.12, `--engine=xcmp`, run from the
`jmcdani20/hap.py` Docker image) on a purpose-built representation-
equivalence case: a synthetic reference with (a) a truth MNP `AC>GT`
written in the query as two SNPs `A>G` + `C>T`, and (b) a single-base
deletion in an 8-bp homopolymer written right-shifted in the truth and
left-aligned in the query. Both are cases a naive keyed match miscounts
as FP+FN.

| stratum | hap.py (xcmp) | vcfclick `--engine haplotype` |
|---|---|---|
| INDEL recall / precision | 1.0 / 1.0 | 1.0 / 1.0 |
| SNP precision | 1.0 (0 FP) | 1.0 (0 FP) |
| SNP recall | 1.0 | 0.0 default · **1.0 with `--decompose-mnp`** |

Both tools agree the calls are **concordant — zero FP, zero FN** — on
both the MNP-vs-two-SNP rescue and the shifted-homopolymer deletion. The
one default divergence is *quantification, not matching*: hap.py
decomposes the truth MNP into SNPs for BVT counting, so it lands in the
SNP recall denominator; vcfclick keeps it whole (typed `UNK`, surfaced in
`run_meta.json` `unsummarized_types`), so SNP recall reads 0/0. Running
with `--decompose-mnp` makes vcfclick bin it the same way and the numbers
match hap.py exactly (SNP and INDEL 1.0 / 1.0).

**Scope and honesty.** This validates the comparison *logic* on a small
synthetic case, not accuracy over a real genome. A full conformance run
against GIAB HG002 on GRCh38 (whole-genome truth VCF + high-confidence
BED) has not been done; SNP concordance is expected to match hap.py, and
the INDEL/complex quantification gap above is the one known difference to
quantify at scale. hap.py/vcfeval are dev-only oracles, never runtime
dependencies.

The self-benchmark invariant (truth == query ⇒ Recall = Precision =
F1 = 1.0 for every stratum) is asserted in
`tests/benchmark/test_pipeline.py`; the module unit tests live under
`tests/benchmark/`.
