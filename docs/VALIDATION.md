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

The naive rule over-calls de novo by **~12×**. `tests/fixtures/giab/denovo_trio.vcf.gz`
encodes 7 of these sites (the 4 confident + 3 no-call) with real HG002
`GT:GQ:DP:AD`. On that fixture:

- `vcfclick db trio --category denovo` **(ingested `--keep-reference`)**
  returns **exactly the 4** confident sites and **excludes the 3**
  no-call ones.
- The **same data without `--keep-reference`** returns **0** — vcfclick
  refuses to call de novo when it cannot prove a parent is 0/0.

## 2. External cross-check (`bcftools +mendelian2`)

Run independently on the same fixture, `bcftools +mendelian2 -m e` flags
as Mendelian-erroneous **exactly the 4 sites** vcfclick calls de novo —
and, like vcfclick, does not flag the 3 no-call sites (a missing parent
genotype is undecidable). The two tools agree site-for-site. This is
asserted in CI (`test_giab_denovo_matches_bcftools_mendelian`, skipped
where bcftools is absent).

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
