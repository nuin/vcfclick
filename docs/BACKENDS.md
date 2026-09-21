# Backends

vcfclick has one logical schema and two local storage engines:

- **chDB**: embedded ClickHouse engine. This is the default when `chdb`
  imports successfully and is the best fit for cohort-scale local
  work.
- **DuckDB**: embedded single-file engine. This is useful for
  Bioconda/conda environments, constrained installs, and lightweight
  local workflows.

Both backends expose the same main tables: `variants`, `genotypes`,
`samples`, and `ingestions`. Some operations are still chDB-first,
especially `vcfclick db stats`, because it uses ClickHouse-specific
catalog and Map functions.

## Select A Backend

Set `VCFCLICK_BACKEND` before running a command:

```bash
VCFCLICK_BACKEND=chdb vcfclick db list
VCFCLICK_BACKEND=duckdb vcfclick db list
```

If the variable is unset, vcfclick auto-detects:

1. use chDB if the `chdb` Python package is importable;
2. otherwise fall back to DuckDB.

Use the same backend consistently for a database name. The backends use
different on-disk formats, even though they live under the same
`VCFCLICK_HOME/dbs/<name>/` naming scheme.

## Where Data Lives

By default:

```text
~/.vcfclick/dbs/<name>/
```

Change that with:

```bash
export VCFCLICK_HOME=/path/to/vcfclick-home
```

This is the easiest way to keep separate workspaces, for example one
home for production-like chDB cohorts and one for DuckDB test imports.

## chDB

chDB is the default install path from PyPI because `chdb` is a runtime
dependency of the package. Note that chDB ships only as a PyPI wheel —
it is **not packaged for conda**, so a `conda install` cannot pull it
in. conda / Bioconda installs run on DuckDB (see below).

Use it when:

- the cohort is large enough that ClickHouse-style scans matter;
- you want `vcfclick db stats`;
- you are using the default PyPI install on macOS or Linux.

Example:

```bash
export VCFCLICK_BACKEND=chdb
vcfclick db create study
vcfclick db ingest study cohort.vcf.gz --cohort cases
vcfclick db stats study
```

### macOS: why chDB is version-pinned

On macOS, `pyproject.toml` pins `chdb==4.2.0` and `chdb-core==26.5.0`
rather than tracking the latest release. This is a workaround for an
upstream packaging bug, not a preference.

chDB's prebuilt `_chdb.abi3.so` places its `__LINKEDIT` symbol string
table on a 4-byte boundary in most published builds. macOS 26+ (Darwin
27) tightened dyld's Mach-O checks to require **8-byte** alignment, so
those builds fail to load at import:

```
ImportError: dlopen(.../chdb/_chdb.abi3.so):
  mis-aligned LINKEDIT string pool, fileOffset=0x12C00DE4
```

The reported offset is exactly the binary's `LC_SYMTAB.stroff`. You can
check any build yourself:

```bash
otool -l .../chdb/_chdb.abi3.so | grep -A5 LC_SYMTAB | grep stroff
# stroff % 8 == 0  -> loads;  != 0  -> dyld refuses
```

Measured across published wheels (arm64):

| package | `stroff % 8` | loads on macOS 26+ |
|---|---|---|
| chdb 4.1.8 (bundled engine) | 4 | no |
| chdb-core 26.3.0 | 4 | no |
| **chdb-core 26.5.0** | **0** | **yes** |
| chdb-core 26.7.3 | 4 | no |

`chdb` 4.2+ ships no macOS wheel of its own — the engine comes from
`chdb-core` — and 4.2.0 is the matching Python package (4.2.1 hits an
unrelated circular import on Python 3.14; 4.3+ require a misaligned
core). Linux is unaffected and tracks `chdb>=4.1.8` as before.

**Lift the pin** once upstream links the string table aligned: drop the
two `sys_platform == 'darwin'` entries from `pyproject.toml`.

## DuckDB

DuckDB is selected explicitly:

```bash
export VCFCLICK_BACKEND=duckdb
vcfclick db create study
vcfclick db ingest study cohort.vcf.gz --cohort cases
```

Use it when:

- you are installing through Bioconda or a conda-first workflow;
- you want a simple single-file backend;
- you are testing backend-portable behavior.

Current caveat: `vcfclick db stats` is not yet implemented for DuckDB.
Use `vcfclick db info` for basic counts.

## Move Data Between Backends

Do not copy backend directories between engines. Use Parquet:

```bash
# Export from the source backend.
VCFCLICK_BACKEND=chdb vcfclick db dump study --out study-dump

# Import into the target backend.
VCFCLICK_BACKEND=duckdb vcfclick db create study-duck
VCFCLICK_BACKEND=duckdb vcfclick db ingest-parquet study-duck study-dump \
  --cohort cases \
  --ingest-id imported-from-chdb
```

The same path also works in the other direction.

## Bundles Are Backend-Neutral

`vcfclick db push` creates a tarball containing Parquet files:

```bash
vcfclick db push study study.tar.gz
```

`vcfclick db pull` restores those files into the active backend:

```bash
VCFCLICK_BACKEND=duckdb vcfclick db pull study study.tar.gz
```

That makes bundles suitable for sharing demo cohorts, moving data
between machines, and migrating between backends.

## Choosing For A Lab

For individual labs and core facilities, the practical default is:

- use chDB for normal local cohort work;
- use DuckDB when packaging constraints require it;
- share data as Parquet bundles, not backend directories;
- set `VCFCLICK_HOME` explicitly on shared workstations so users know
  where cohort data is stored.
