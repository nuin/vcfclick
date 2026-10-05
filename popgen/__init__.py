"""Population-genetics statistics over a vcfclick database.

Layers, mirroring the rest of vcfclick:

  * `popgen.estimators` — pure numpy estimators (θ_W, π, Tajima's D,
    Fay & Wu's H, Hudson F_ST, SFS projection). No database access.
  * `popgen.ancestral`  — `info_AA` normalisation and polarisation.
  * `popgen.counts`     — per-site, per-group allele counts via SQL, in
    both backend dialects.
  * `popgen.analysis`   — site filters + the `summary` / `sfs` / `fst` /
    `windows` computations that `vcfclick db popgen` prints.

See docs/POPGEN.md for the definitions and the choices made.
"""
