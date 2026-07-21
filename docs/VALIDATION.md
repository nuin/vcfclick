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

### Real GIAB HG002 / GRCh38 (chr20 and chr1)

Beyond the synthetic case, the same hap.py build was run against vcfclick
on **real GIAB data**, in two independent 5-Mb windows: the HG002 v4.2.1
benchmark truth over GRCh38 `chr20:1,000,000–6,000,000` (8,288 variants)
and `chr1:1,000,000–6,000,000` (7,984 variants), each inside the real
high-confidence BED, reference streamed from the GIAB
`no_alt_analysis_set`. The query is a controlled perturbation of that
truth (~4% of variants dropped → known FN, ~2% of genotypes flipped →
known errors) so both tools score the same real, messy call set.
Per-stratum, `ALL` filter:

| region | stratum | metric | hap.py (xcmp) | vcfclick | Δ |
|---|---|---|---|---|---|
| chr20 | SNP | recall / prec | 0.9402 / 0.9796 | 0.9414 / 0.9795 | +0.001 / −0.000 |
| chr20 | SNP | FN / FP | 405 / 133 | 410 / 138 | +5 / +5 |
| chr20 | INDEL | recall / prec | 0.9344 / 0.9747 | 0.9453 / 0.9791 | +0.011 / +0.005 |
| chr20 | INDEL | FN / FP | 75 / 29 | 76 / 28 | +1 / −1 |
| chr1 | SNP | recall / prec | 0.9414 / 0.9803 | 0.9423 / 0.9800 | +0.001 / −0.000 |
| chr1 | SNP | FN / FP | 403 / 130 | 410 / 137 | +7 / +7 |
| chr1 | INDEL | recall / prec | 0.9258 / 0.9715 | 0.9355 / 0.9761 | +0.010 / +0.005 |
| chr1 | INDEL | FN / FP | 56 / 21 | 59 / 21 | +3 / +0 |

Both chromosomes reproduce the same picture. **SNP concordance is
essentially exact** (recall and precision within 0.1%; FN/FP within ~7 of
~6,900 variants). **INDEL agrees within ~1%**, vcfclick reading slightly
higher on both — attributable to counting convention: vcfclick decomposes
multiallelics into more biallelic rows. (These numbers predate the
confident-region gating fix described below; the version-0.8.0 columns
here undercounted UNK — vcfclick reported UNK = 0 where hap.py bucketed a
few hundred out-of-region calls as not-assessable. The gating fix aligns
vcfclick's UNK with hap.py's and shifts these matched out-of-region calls
out of TP into UNK, moving the metrics still closer to hap.py.)

### Independent caller (GIAB T2T-Q100 assembly) vs v4.2.1 truth

The controlled runs above perturb the truth, so they don't stress
cross-representation matching at scale. This one does: the query is a
genuinely **independent HG002 call set** — the assembly-based
`GRCh38_HG2-T2TQ100-V1.1` dipcall (a different method entirely, and
*phased*) — scored against the v4.2.1 truth over the same `chr20:1–6M`
confident regions with both tools.

| stratum | metric | hap.py (xcmp) | vcfclick v0.8.0 | vcfclick (fixed) | Δ fixed |
|---|---|---|---|---|---|
| SNP | recall / precision | 1.0000 / 1.0000 | 0.9993 / 0.9990 | 0.9999 / 0.9996 | −0.000 / −0.000 |
| INDEL | recall / precision | 1.0000 / 0.9991 | 0.9604 / 0.9603 | 0.9975 / 0.9975 | −0.002 / −0.002 |

**Both strata now match hap.py to ~0.2%** on a fully independent, phased
assembly caller. Getting there took chasing the INDEL gap through three
successive causes — the first two real bugs, the third the actual driver:

1. **Not normalization.** vcfclick's `left_align` matches `bcftools norm
   -f` exactly on all 1,214 query indels (0 divergences); after identical
   normalization 1,211 / 1,217 indels share an exact key.
2. **Het-alt multiallelic indels** (`1/2` spelled as two `0/1` records)
   were scored `am`, and the `haplotype` engine *skipped* `am` rows.
   Fixed by feeding `am` rows through the diplotype-equivalence check
   (representation-equal het-alts resolve to TP; genuine zygosity errors
   stay `am`). INDEL recall 0.85 → 0.96.
3. **Confident-region gating (the real driver).** The residual after (2)
   was *not* an algorithm limit — the remaining ~50 indels were all
   **outside the confident BED**, where GIAB truth is unknown and hap.py
   marks them UNK. `classify`'s "both sides present" branch never checked
   `in_conf` (only the single-sided branches did), so out-of-region
   matches were scored as FN+FP instead of UNK. Gating that branch lifts
   INDEL recall 0.96 → **0.9975** and — the same bug — turns the earlier
   "vcfclick UNK = 0 vs hap.py ~370" discrepancy into a match (vcfclick
   now reports UNK ≈ hap.py). The residual is 3 INDEL FN of 1,217.

**Bug found and fixed by this run.** The independent phased caller also
surfaced a genotype-matching bug: the keyed verdict compared the raw GT
string, so a phased query `1|0` never matched an unphased truth `0/1` and
every shared heterozygote scored as a spurious mismatch (8,142
bcftools-shared variants → 0 TP before the fix). Genotype comparison is
now phase-insensitive (a real `0/1`-vs-`1/1` zygosity error still
mismatches). Every prior test used matching phase, so only real
independent data caught it.

**Bottom line.** On a fully independent, phased assembly caller, vcfclick's
`haplotype` engine tracks hap.py to **~0.2% on both SNPs and INDELs**,
after fixing the three issues this exercise surfaced (phase-insensitive
genotypes, het-alt rescue, and confident-region gating on matched pairs).
The residual is a handful of individual variants, not a systematic gap.
Reproduce with the `jmcdani20/hap.py:v0.3.12` image; hap.py/vcfeval are
dev-only oracles, never runtime dependencies.

The self-benchmark invariant (truth == query ⇒ Recall = Precision =
F1 = 1.0 for every stratum) is asserted in
`tests/benchmark/test_pipeline.py`; the module unit tests live under
`tests/benchmark/`.
