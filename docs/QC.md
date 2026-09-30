# Sample QC (`vcfclick db qc`)

`vcfclick db qc <name>` reports per-sample quality-control metrics over a
cohort, computed in one pass against the genotypes table. It works on
both backends and needs no external data.

```bash
vcfclick db qc my-cohort
vcfclick db qc my-cohort --format json   # for pipelines
```

| Column | Meaning |
|---|---|
| `variants` | stored non-reference calls for the sample |
| `het` / `hom` | heterozygous (`gt=1`) / homozygous-alt (`gt=2`) counts |
| `het/hom` | het-to-hom-alt ratio (a standard genotyping sanity check) |
| `ti/tv` | transition/transversion ratio over SNVs (≈2.0–2.1 for WGS, ≈3 for exomes; low values suggest false positives) |
| `chrX-het` | heterozygous fraction on chromosome X |
| `sex` | sex inferred from chrX heterozygosity (`male` if low, `female` if high), `*` when it disagrees with the pedigree |

## Sex check

Males are hemizygous on the non-PAR X, so their chrX heterozygous
fraction is near zero; females sit near 0.3–0.6. When a pedigree is
loaded (`vcfclick db ped`), the inferred sex is compared to the declared
sex and a mismatch is flagged (`female*`) with a warning — the classic
signal of a sample swap or mislabel. Samples with too few chrX calls to
decide are reported as `unknown`.

## Relatedness

`vcfclick db relatedness <name>` estimates kinship for every pair of
samples in an ingestion with the KING-robust estimator (Manichaikul et
al. 2010), over biallelic SNVs where both samples are called:

```
kinship = (sites both het  -  2 x sites with opposite homozygotes) / (het sites of A + het sites of B)
```

| Kinship | Relationship |
|---|---|
| > 0.354 | duplicate sample or MZ twin |
| 0.177 – 0.354 | first degree: `parent-child` if opposite homozygotes (IBS0) are < 0.5% of sites, else `full-siblings` |
| 0.088 – 0.177 | second degree |
| 0.044 – 0.088 | third degree |
| below | unrelated |

By default only related pairs are listed; `--all` lists every pair. When a
pedigree is loaded, every declared parent is checked against the
genotypes, and a declared parent who isn't a first-degree `parent-child`
match is flagged `mismatch` (sample swap, mislabel or non-paternity).

Two limits come from the data:

- **It needs genome-wide markers.** Within one gene or a small region,
  markers are inherited together, so unrelated people who share a common
  haplotype look like duplicates. Below 1,000 SNVs or a 10 Mb span, pairs
  are left unclassified with a warning (`--force` classifies anyway).
- **No-calls.** For ingestions made with `--keep-reference`, an absent
  genotype is a real no-call and is skipped for that pair. Otherwise only
  non-reference calls were stored, so absent is read as `0/0`; missing
  calls then look hom-ref and pull kinship slightly down. Joint-called VCFs
  have few, so relationship calls remain reliable.

Kinship is computed with dense matrix products, one ingestion at a time,
and SNVs are thinned evenly to `--max-sites` (default 200,000) to bound
memory. Cohorts above 20,000 samples are refused.

To cross-check relatedness and the sex call against somalier, see
[Validation](VALIDATION.md#relatedness-and-sex-cross-check-with-somalier-optional).

## What it does not report

The genotypes table is **sparse** — only non-reference calls are stored,
so a sample's `0/0` and `./.` are indistinguishable by absence. Genotype
*missingness* / call rate therefore cannot be computed here; QC reports
only the metrics the stored non-reference calls support honestly.
