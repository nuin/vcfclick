# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo", "vcfclick[web]", "pandas", "matplotlib"]
# ///
"""vcfclick demo as a reactive marimo notebook.

Run it (vcfclick + marimo are pulled in by the inline script metadata):

    uvx marimo edit --sandbox examples/vcfclick_demo.py     # edit
    uvx marimo run  --sandbox examples/vcfclick_demo.py     # app view

Or, in an env where `vcfclick` is already installed:

    pip install marimo pandas matplotlib
    marimo edit examples/vcfclick_demo.py

Unlike the .ipynb sibling this is a plain-Python, git-friendly, reactive
notebook — but it needs a real Python kernel (vcfclick's chdb/cyvcf2 cannot run
in a browser-WASM build).
"""

import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # vcfclick — a reactive demo

    [vcfclick](https://github.com/nuin/vcfclick) turns VCF cohorts into
    local, queryable SQL databases — trio/family analysis, gnomAD rarity
    filtering, sample QC, and a GATK3-`CombineVariants` reimplementation.

    This notebook **downloads its own data** and runs end to end. Trio
    analysis is validated against the **GIAB Ashkenazi trio** benchmark
    ([docs/VALIDATION.md](https://github.com/nuin/vcfclick/blob/main/docs/VALIDATION.md)).
    """)
    return


@app.cell
def _():
    import json
    import os
    import pathlib
    import subprocess
    import urllib.request

    # Throwaway sandbox so we never touch a real ~/.vcfclick
    os.environ["VCFCLICK_HOME"] = "/tmp/vcfclick-marimo"
    os.environ["VCFCLICK_ANNOTATIONS_DB"] = "/tmp/vcfclick-marimo/ann.duckdb"
    DATA = pathlib.Path("/tmp/vcfclick-marimo/data")
    DATA.mkdir(parents=True, exist_ok=True)

    def vcf(*args: str) -> str:
        """Run a vcfclick command, returning its combined stdout+stderr."""
        r = subprocess.run(["vcfclick", *args], capture_output=True, text=True)
        return (r.stdout + r.stderr).strip()

    return DATA, json, urllib, vcf


@app.cell
def _(DATA, mo, urllib, vcf):
    _RAW = "https://raw.githubusercontent.com/nuin/vcfclick/main/tests/fixtures/"
    _FILES = [
        "tiny.vcf.gz", "qc_sex.vcf.gz", "qc_sex.ped", "gnomad_cftr.vcf.gz",
        "callset_a.vcf", "callset_b.vcf",
        "giab/cftr_trio.vcf.gz", "giab/cftr_trio.ped",
        "giab/denovo_trio.vcf.gz", "giab/denovo_trio.ped",
    ]
    for _f in _FILES:
        urllib.request.urlretrieve(_RAW + _f, DATA / _f.split("/")[-1])
    mo.md(f"`vcfclick {vcf('--version')}` · downloaded **{len(_FILES)}** example VCFs.")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 1. A cohort in, SQL out

    Ingest a VCF, query it. Variants/genotypes/samples land in an embedded
    ClickHouse (chDB) database; queries run in-process.
    """)
    return


@app.cell
def _(DATA, mo, vcf):
    vcf("db", "create", "demo")
    vcf("db", "ingest", "demo", str(DATA / "tiny.vcf.gz"),
        "--cohort", "study", "--ingest-id", "v1", "--serial")
    mo.md(
        "```\n"
        + vcf("db", "query", "demo",
              "SELECT chrom, pos, ref, alt FROM variants ORDER BY pos")
        + "\n```"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Trio de novo — and why it's *defensible*

    De novo means the child carries a variant **neither parent has** —
    which needs the parents to be *provably* hom-reference. The sparse
    genotype table can't tell `0/0` from `./.`, so de novo requires a
    `--keep-reference` ingest that stores confident parent `0/0`.

    The fixture is the real GIAB trio at 7 chr20 sites: 4 confident de
    novos (both parents in GIAB's high-confidence BED), 3 with a no-call
    parent.
    """)
    return


@app.cell
def _(DATA, mo, vcf):
    vcf("db", "create", "dn")
    vcf("db", "ingest", "dn", str(DATA / "denovo_trio.vcf.gz"),
        "--cohort", "fam", "--ingest-id", "g1", "--serial", "--keep-reference")
    vcf("db", "ped", "dn", str(DATA / "denovo_trio.ped"))
    mo.md(
        "**With `--keep-reference` → exactly 4 confident sites:**\n\n```\n"
        + vcf("db", "trio", "dn", "--proband", "HG002",
              "--category", "denovo", "--max-af", "1.0")
        + "\n```"
    )
    return


@app.cell
def _(DATA, mo, vcf):
    # Same data, naive sparse ingest (no --keep-reference) → de novo is
    # undecidable, so it correctly returns 0 rather than guessing.
    vcf("db", "create", "dn_naive")
    vcf("db", "ingest", "dn_naive", str(DATA / "denovo_trio.vcf.gz"),
        "--cohort", "fam", "--ingest-id", "g1", "--serial")
    vcf("db", "ped", "dn_naive", str(DATA / "denovo_trio.ped"))
    _out = vcf("db", "trio", "dn_naive", "--proband", "HG002",
               "--category", "denovo", "--max-af", "1.0")
    mo.md("**Without it → 0** (won't guess where a parent is a no-call):\n\n```\n"
          + [ln for ln in _out.splitlines() if "candidates" in ln][0] + "\n```")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Compound-het + gnomAD rarity filtering

    A real CFTR slice reproduces every model; CFTR is a genuine
    compound-het candidate gene. Loading real gnomAD frequencies and adding
    `--gnomad-max-af` then drops the variants that are actually common.
    """)
    return


@app.cell
def _(DATA, mo, vcf):
    vcf("db", "create", "fam")
    vcf("db", "ingest", "fam", str(DATA / "cftr_trio.vcf.gz"),
        "--cohort", "trio", "--ingest-id", "g1", "--serial", "--keep-reference")
    vcf("db", "ped", "fam", str(DATA / "cftr_trio.ped"))
    vcf("annotations", "load-gnomad", str(DATA / "gnomad_cftr.vcf.gz"))
    # CFTR coords (normally `vcfclick annotations load` pulls all of GENCODE):
    from annotations.db import get_connection as _gc

    _c = _gc()
    _c.execute(
        "INSERT INTO refseq_genes VALUES ('CFTR','chr7',117480025,117668665,'+','1080','x')"
    )
    _c.close()
    mo.md(
        "**all models**\n```\n"
        + vcf("db", "trio", "fam", "--proband", "HG002")
        + "\n```\n**+ gnomAD popmax < 0.01**\n```\n"
        + vcf("db", "trio", "fam", "--proband", "HG002", "--gnomad-max-af", "0.01")
        + "\n```"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Per-sample QC + a chrX sex check

    `db qc` infers sex from chrX heterozygosity and flags it against the
    pedigree — the sample-swap signal. Read straight into a DataFrame
    (marimo renders it interactively).
    """)
    return


@app.cell
def _(DATA, json, vcf):
    import pandas as pd

    vcf("db", "create", "qc")
    vcf("db", "ingest", "qc", str(DATA / "qc_sex.vcf.gz"),
        "--cohort", "c", "--ingest-id", "i1", "--serial")
    vcf("db", "ped", "qc", str(DATA / "qc_sex.ped"))
    qc_df = pd.DataFrame(json.loads(vcf("db", "qc", "qc", "--format", "json")))[
        ["sample_id", "het", "hom_alt", "het_hom_ratio", "ti_tv",
         "chrx_het_frac", "inferred_sex", "sex_mismatch"]
    ]
    qc_df
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 5. Combine call sets — the GATK3 `CombineVariants` GATK4 removed

    `combine` unions call sets that may *share* samples, resolving overlaps
    by input priority and annotating each record with `set=` provenance.
    """)
    return


@app.cell
def _(DATA, mo, vcf):
    vcf("combine", str(DATA / "callset_a.vcf"), str(DATA / "callset_b.vcf"),
        "-o", "/tmp/vcfclick-marimo/combined.vcf", "--name", "first", "--name", "second")
    _body = "\n".join(
        ln for ln in open("/tmp/vcfclick-marimo/combined.vcf") if not ln.startswith("##")
    )
    mo.md("```\n" + _body + "\n```")
    return


@app.cell(hide_code=True)
def _():
    import matplotlib.pyplot as plt

    # Across a 5 Mb chr20 GIAB region the naive rule calls 50 de novos;
    # GIAB's BED confirms only 4 (see docs/VALIDATION.md).
    _fig, _ax = plt.subplots(figsize=(5, 3))
    _ax.bar(["naive\nde novo", "confirmed\n(GIAB BED)"], [50, 4],
            color=["#f59e0b", "#16a34a"])
    _ax.set_ylabel("de-novo calls in 5 Mb")
    _ax.set_title("--keep-reference avoids ~46 false positives")
    for _i, _v in enumerate([50, 4]):
        _ax.text(_i, _v + 1, str(_v), ha="center", fontweight="bold")
    plt.tight_layout()
    _ax
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Next steps

    - **Web UI:** `vcfclick web demo` — a local SQL explorer + NL→SQL +
      trio/combine panels.
    - **Docs:** [README](https://github.com/nuin/vcfclick) ·
      [Trio](https://github.com/nuin/vcfclick/blob/main/docs/TRIO.md) ·
      [Validation](https://github.com/nuin/vcfclick/blob/main/docs/VALIDATION.md)
    """)
    return


if __name__ == "__main__":
    app.run()
