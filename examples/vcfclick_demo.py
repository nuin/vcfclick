# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo", "vcfclick[web]", "pandas", "matplotlib"]
# ///
"""vcfclick — an interactive demo as a reactive marimo notebook.

Run it (vcfclick + marimo are pulled in by the inline script metadata):

    uvx marimo edit --sandbox examples/vcfclick_demo.py     # edit
    uvx marimo run  --sandbox examples/vcfclick_demo.py     # app view

Or, in an env where `vcfclick` is already installed:

    pip install marimo pandas matplotlib
    marimo edit examples/vcfclick_demo.py

The controls (SQL editor, trio dropdown/slider, --keep-reference toggle,
min-callsets slider) drive the analysis live — change one and the
dependent cells re-run. Needs a real Python kernel (vcfclick's
chdb/cyvcf2 cannot run in a browser-WASM build).
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
    # vcfclick — an interactive demo

    [vcfclick](https://github.com/nuin/vcfclick) turns VCF cohorts into
    local, queryable SQL databases — trio/family analysis, gnomAD rarity
    filtering, sample QC, and a GATK3-`CombineVariants` reimplementation.

    This notebook **downloads its own data** and is **interactive**: edit
    the SQL, change the inheritance model, drag the gnomAD-frequency
    slider, flip `--keep-reference` — the results below re-compute live.
    Trio analysis is validated against the **GIAB Ashkenazi trio**
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

    os.environ["VCFCLICK_HOME"] = "/tmp/vcfclick-marimo"
    os.environ["VCFCLICK_ANNOTATIONS_DB"] = "/tmp/vcfclick-marimo/ann.duckdb"
    DATA = pathlib.Path("/tmp/vcfclick-marimo/data")
    DATA.mkdir(parents=True, exist_ok=True)

    def vcf(*args: str) -> str:
        """Run a vcfclick command, returning its combined stdout+stderr."""
        r = subprocess.run(["vcfclick", *args], capture_output=True, text=True)
        return (r.stdout + r.stderr).strip()

    return DATA, json, urllib, vcf


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Download the example data
    """)
    return


@app.cell
def _(DATA, mo, urllib, vcf):
    import pandas as pd

    RAW = "https://raw.githubusercontent.com/nuin/vcfclick/main/tests/fixtures/"
    FILES = [
        "tiny.vcf.gz",
        "qc_sex.vcf.gz",
        "qc_sex.ped",
        "gnomad_cftr.vcf.gz",
        "callset_a.vcf",
        "callset_b.vcf",
        "giab/cftr_trio.vcf.gz",
        "giab/cftr_trio.ped",
        "giab/denovo_trio.vcf.gz",
        "giab/denovo_trio.ped",
    ]
    _rows = []
    with mo.status.progress_bar(
        total=len(FILES), title="Downloading from GitHub"
    ) as _bar:
        for _f in FILES:
            _name = _f.split("/")[-1]
            _bar.update(subtitle=_name)
            _dest = DATA / _name
            urllib.request.urlretrieve(RAW + _f, _dest)
            _rows.append({"file": _name, "KB": round(_dest.stat().st_size / 1024, 1)})
    files_df = pd.DataFrame(_rows)
    data_ready = True  # downstream cells depend on this so they wait for the data

    mo.vstack(
        [
            mo.md(
                f"Downloaded **{len(FILES)}** example VCFs · `vcfclick {vcf('--version')}`"
            ),
            files_df,
        ]
    )
    return data_ready, pd


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 1. SQL, live

    Ingest a VCF, then **edit the query** below — the result re-runs as you
    type.
    """)
    return


@app.cell
def _(DATA, data_ready, vcf):
    if data_ready:
        vcf("db", "create", "demo")
        vcf(
            "db",
            "ingest",
            "demo",
            str(DATA / "tiny.vcf.gz"),
            "--cohort",
            "study",
            "--ingest-id",
            "v1",
            "--serial",
        )
    demo_db = "demo"
    return (demo_db,)


@app.cell
def _(mo):
    sql = mo.ui.text_area(
        value="SELECT chrom, pos, ref, alt FROM variants ORDER BY pos",
        label="SQL",
        full_width=True,
        rows=2,
    )
    sql
    return (sql,)


@app.cell
def _(demo_db, mo, pd, sql, vcf):
    import io as _io

    _out = vcf("db", "query", demo_db, sql.value, "--format", "CSVWithNames")
    try:
        _df = pd.read_csv(_io.StringIO(_out))
        _view = mo.ui.table(_df, selection=None) if len(_df) else mo.md("_0 rows_")
    except Exception:
        _view = mo.md("```\n" + _out + "\n```")  # surface SQL errors as text
    _view
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Trio analysis — pick a model, drag the rarity slider

    On the real GIAB CFTR trio. `comphet` groups two rare hets in the same
    gene; the slider applies a **gnomAD popmax** rarity cut (1.0 = no
    filter). Change either control and the candidate counts update.
    """)
    return


@app.cell
def _(DATA, data_ready, vcf):
    if data_ready:
        vcf("db", "create", "fam")
        vcf(
            "db",
            "ingest",
            "fam",
            str(DATA / "cftr_trio.vcf.gz"),
            "--cohort",
            "trio",
            "--ingest-id",
            "g1",
            "--serial",
            "--keep-reference",
        )
        vcf("db", "ped", "fam", str(DATA / "cftr_trio.ped"))
        vcf("annotations", "load-gnomad", str(DATA / "gnomad_cftr.vcf.gz"))
        from annotations.db import get_connection as _gc

        _c = _gc()
        _c.execute(
            "INSERT INTO refseq_genes VALUES ('CFTR','chr7',117480025,117668665,'+','1080','x')"
        )
        _c.close()
    fam_db = "fam"
    return (fam_db,)


@app.cell
def _(mo):
    category = mo.ui.dropdown(
        options=["all", "denovo", "recessive", "dominant", "comphet"],
        value="all",
        label="--category",
    )
    gnomad = mo.ui.slider(
        start=0.0,
        stop=1.0,
        step=0.01,
        value=1.0,
        label="--gnomad-max-af",
        show_value=True,
    )
    mo.hstack([category, gnomad], justify="start", gap=2)
    return category, gnomad


@app.cell
def _(category, fam_db, gnomad, mo, vcf):
    _args = [
        "db",
        "trio",
        fam_db,
        "--proband",
        "HG002",
        "--category",
        category.value,
        "--gnomad-max-af",
        str(gnomad.value),
    ]
    mo.md(f"`vcfclick {' '.join(_args[1:])}`\n\n```\n{vcf(*_args)}\n```")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Why de novo needs `--keep-reference`

    Toggle it. The genotype table is sparse, so a parent absent at a site
    is `0/0` **or** a no-call `./.`. With `--keep-reference` (parent `0/0`
    stored) vcfclick proves the **4** confident de novos; without it,
    absence is undecidable so it returns **0** rather than guess. Same real
    GIAB data, both ways.
    """)
    return


@app.cell
def _(DATA, data_ready, vcf):
    if data_ready:
        for _name, _kr in [("dn", True), ("dn_naive", False)]:
            vcf("db", "create", _name)
            _a = [
                "db",
                "ingest",
                _name,
                str(DATA / "denovo_trio.vcf.gz"),
                "--cohort",
                "fam",
                "--ingest-id",
                "g1",
                "--serial",
            ]
            if _kr:
                _a.append("--keep-reference")
            vcf(*_a)
            vcf("db", "ped", _name, str(DATA / "denovo_trio.ped"))
    dn_dbs = {"on": "dn", "off": "dn_naive"}
    return (dn_dbs,)


@app.cell
def _(mo):
    keepref = mo.ui.switch(value=True, label="ingest with `--keep-reference`")
    keepref
    return (keepref,)


@app.cell
def _(dn_dbs, keepref, mo, vcf):
    _db = dn_dbs["on" if keepref.value else "off"]
    _out = vcf(
        "db",
        "trio",
        _db,
        "--proband",
        "HG002",
        "--category",
        "denovo",
        "--max-af",
        "1.0",
    )
    _n = next((ln for ln in _out.splitlines() if "candidates" in ln), "")
    mo.md(f"**`--keep-reference` = {keepref.value}** → {_n or _out}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Per-sample QC + chrX sex check

    `db qc` infers sex from chrX heterozygosity and flags it against the
    pedigree — the sample-swap signal (`SWAP`).
    """)
    return


@app.cell
def _(DATA, data_ready, json, pd, vcf):
    if data_ready:
        vcf("db", "create", "qc")
        vcf(
            "db",
            "ingest",
            "qc",
            str(DATA / "qc_sex.vcf.gz"),
            "--cohort",
            "c",
            "--ingest-id",
            "i1",
            "--serial",
        )
        vcf("db", "ped", "qc", str(DATA / "qc_sex.ped"))
    pd.DataFrame(json.loads(vcf("db", "qc", "qc", "--format", "json")))[
        [
            "sample_id",
            "het",
            "hom_alt",
            "het_hom_ratio",
            "ti_tv",
            "chrx_het_frac",
            "inferred_sex",
            "sex_mismatch",
        ]
    ]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 5. Combine call sets (GATK3 `CombineVariants`)

    Union call sets that may *share* samples, resolving overlaps by input
    priority with `set=` provenance. Drag `--min-callsets` to keep only
    consensus sites.
    """)
    return


@app.cell
def _(mo):
    min_callsets = mo.ui.slider(
        start=1, stop=2, step=1, value=1, label="--min-callsets", show_value=True
    )
    min_callsets
    return (min_callsets,)


@app.cell
def _(DATA, data_ready, min_callsets, mo, vcf):
    if data_ready:
        vcf(
            "combine",
            str(DATA / "callset_a.vcf"),
            str(DATA / "callset_b.vcf"),
            "-o",
            "/tmp/vcfclick-marimo/combined.vcf",
            "--name",
            "first",
            "--name",
            "second",
            "--min-callsets",
            str(min_callsets.value),
        )
    _body = "\n".join(
        ln
        for ln in open("/tmp/vcfclick-marimo/combined.vcf")
        if not ln.startswith("##")
    )
    mo.md("```\n" + _body + "\n```")
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
