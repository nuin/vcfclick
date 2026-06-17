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

## What it does not report

The genotypes table is **sparse** — only non-reference calls are stored,
so a sample's `0/0` and `./.` are indistinguishable by absence. Genotype
*missingness* / call rate therefore cannot be computed here; QC reports
only the metrics the stored non-reference calls support honestly.
