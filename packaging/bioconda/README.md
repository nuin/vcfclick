# bioconda recipe for vcfclick

`meta.yaml` is the conda-build recipe ready for submission to
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes).

## Why bioconda

Genomics labs install via `conda install -c bioconda <tool>`. A tool that
isn't in bioconda is invisible to a large fraction of the field. The
recipe here uses the DuckDB backend introduced in vcfclick v0.3.0 —
chDB is not available on conda-forge and the upstream maintainer has
declared it impractical to package (see
[chdb-io/chdb#189](https://github.com/chdb-io/chdb/issues/189)).

## Submitting a new version

1. Wait for the new vcfclick release to be on PyPI.
2. Compute the sdist SHA:
   ```bash
   curl -L https://pypi.io/packages/source/v/vcfclick/vcfclick-${VERSION}.tar.gz \
     | shasum -a 256
   ```
3. Edit `meta.yaml`: bump `version` and replace the `sha256:` placeholder.
4. Open a PR against `bioconda/bioconda-recipes`:
   ```bash
   git clone https://github.com/bioconda/bioconda-recipes
   mkdir -p bioconda-recipes/recipes/vcfclick
   cp packaging/bioconda/meta.yaml bioconda-recipes/recipes/vcfclick/
   cd bioconda-recipes
   # Open PR; the bioconda bot validates and reviews.
   ```

## Local validation

```bash
# Build the recipe against the bioconda channel
conda install -n base -c conda-forge conda-build
conda-build -c conda-forge -c bioconda packaging/bioconda/

# Or via the bioconda bootstrap:
bioconda-utils build --packages vcfclick --mulled-test
```

## Backend selection on conda installs

The conda install pulls only DuckDB. Auto-detect (in `storage.backend()`)
picks DuckDB whenever chdb isn't importable, so users do not need to set
`VCFCLICK_BACKEND` explicitly.

Users who want the chDB backend on top of a conda install can layer it:
```bash
conda install -c bioconda vcfclick
pip install chdb       # adds the chDB backend; opt-in
export VCFCLICK_BACKEND=chdb
```
