"""SQL-native stratification of concordance results against the annotation store.

This is the capability hap.py structurally lacks: join per-variant TP/FP/FN
against gnomAD / ClinVar / genes and report recall/precision per stratum.
"""

from __future__ import annotations

import duckdb
import pyarrow as pa

from benchmark.stratify_db import (
    stratify_by_clinvar,
    stratify_by_gene,
    stratify_by_gnomad,
)


def _ann_store(path):
    con = duckdb.connect(str(path))
    con.execute(
        "CREATE TABLE gnomad_af (chrom VARCHAR, pos UINTEGER, ref VARCHAR, "
        "alt VARCHAR, af DOUBLE, af_grpmax DOUBLE)"
    )
    con.execute(
        "INSERT INTO gnomad_af VALUES "
        "('chr1',100,'A','G',0.0001,0.0002),"  # rare
        "('chr1',200,'C','T',0.20,0.25)"  # common
    )
    con.execute(
        "CREATE TABLE clinvar_variants (chrom VARCHAR, pos UINTEGER, ref VARCHAR, "
        "alt VARCHAR, clin_sig VARCHAR, review_status VARCHAR, clinvar_id VARCHAR, "
        "condition VARCHAR)"
    )
    con.execute(
        "INSERT INTO clinvar_variants VALUES "
        "('chr1',100,'A','G','Pathogenic','','VCV1','x')"
    )
    con.execute(
        "CREATE TABLE refseq_genes (gene_symbol VARCHAR, chrom VARCHAR, "
        "start_pos UINTEGER, end_pos UINTEGER, strand VARCHAR, refseq_id VARCHAR, "
        "description VARCHAR)"
    )
    con.execute("INSERT INTO refseq_genes VALUES ('GENEA','chr1',50,150,'+','1','x')")
    con.close()
    return str(path)


def _conc():
    # rare variant chr1:100 A>G: a TP (truth+query both TP)
    # common variant chr1:200 C>T: truth FN (missed) + no query row
    # a novel (not in gnomad) chr1:300 G>A: query FP
    rows = [
        dict(
            side="truth",
            filter_view="ALL",
            chrom="chr1",
            pos=100,
            ref="A",
            alt="G",
            vtype="SNP",
            bd="TP",
        ),
        dict(
            side="query",
            filter_view="ALL",
            chrom="chr1",
            pos=100,
            ref="A",
            alt="G",
            vtype="SNP",
            bd="TP",
        ),
        dict(
            side="truth",
            filter_view="ALL",
            chrom="chr1",
            pos=200,
            ref="C",
            alt="T",
            vtype="SNP",
            bd="FN",
        ),
        dict(
            side="query",
            filter_view="ALL",
            chrom="chr1",
            pos=300,
            ref="G",
            alt="A",
            vtype="SNP",
            bd="FP",
        ),
    ]
    return pa.Table.from_pylist(rows)


def test_stratify_by_gnomad_af_bins(tmp_path):
    ann = _ann_store(tmp_path / "ann.duckdb")
    out = {r["stratum"]: r for r in stratify_by_gnomad(_conc(), ann)}
    # rare bin: the chr1:100 TP
    assert out["rare"]["truth_tp"] == 1 and out["rare"]["recall"] == 1.0
    # common bin: the chr1:200 FN
    assert out["common"]["truth_fn"] == 1 and out["common"]["recall"] == 0.0
    # novel (no gnomad row): the chr1:300 FP
    assert out["novel"]["query_fp"] == 1


def test_stratify_by_clinvar_significance(tmp_path):
    ann = _ann_store(tmp_path / "ann.duckdb")
    out = {r["stratum"]: r for r in stratify_by_clinvar(_conc(), ann)}
    # chr1:100 A>G is Pathogenic and a TP
    assert out["Pathogenic"]["truth_tp"] == 1
    assert out["Pathogenic"]["precision"] == 1.0


def test_stratify_by_gene(tmp_path):
    ann = _ann_store(tmp_path / "ann.duckdb")
    out = {r["stratum"]: r for r in stratify_by_gene(_conc(), ann)}
    # chr1:100 (TP) and chr1:200 (FN) fall within GENEA [50,150]? 200 is outside.
    assert out["GENEA"]["truth_tp"] == 1  # chr1:100
    assert out["GENEA"]["truth_fn"] == 0  # chr1:200 is outside 50..150


def test_stratify_by_regions(tmp_path):
    from benchmark.stratify_db import stratify_by_regions

    # BED is 0-based half-open. pos100 -> 0-based 99 in [50,150); pos300 -> 299 in [250,350).
    lc = tmp_path / "lowcomplex.bed"
    lc.write_text("chr1\t50\t150\n")
    sd = tmp_path / "segdup.bed"
    sd.write_text("chr1\t250\t350\n")

    out = {
        r["stratum"]: r
        for r in stratify_by_regions(
            _conc(), {"lowcomplex": str(lc), "segdup": str(sd)}
        )
    }
    assert out["lowcomplex"]["truth_tp"] == 1  # chr1:100 TP
    assert out["segdup"]["query_fp"] == 1  # chr1:300 FP
    assert out["none"]["truth_fn"] == 1  # chr1:200 in neither region
