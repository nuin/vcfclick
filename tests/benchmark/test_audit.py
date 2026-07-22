"""Auditability: per-error (FN/FP) variants with their annotation context —
the "why did this caller miss these, and do they matter?" view."""

from __future__ import annotations

import duckdb
import pyarrow as pa

from benchmark.audit import annotated_errors


def _ann(path):
    con = duckdb.connect(str(path))
    con.execute(
        "CREATE TABLE gnomad_af (chrom VARCHAR, pos UINTEGER, ref VARCHAR, "
        "alt VARCHAR, af DOUBLE, af_grpmax DOUBLE)"
    )
    con.execute("INSERT INTO gnomad_af VALUES ('chr1',200,'C','T',0.20,0.25)")
    con.execute(
        "CREATE TABLE clinvar_variants (chrom VARCHAR, pos UINTEGER, ref VARCHAR, "
        "alt VARCHAR, clin_sig VARCHAR, review_status VARCHAR, clinvar_id VARCHAR, "
        "condition VARCHAR)"
    )
    con.execute(
        "INSERT INTO clinvar_variants VALUES "
        "('chr1',200,'C','T','Pathogenic','','VCV9','cond')"
    )
    con.execute(
        "CREATE TABLE refseq_genes (gene_symbol VARCHAR, chrom VARCHAR, "
        "start_pos UINTEGER, end_pos UINTEGER, strand VARCHAR, refseq_id VARCHAR, "
        "description VARCHAR)"
    )
    con.execute("INSERT INTO refseq_genes VALUES ('GENEB','chr1',180,220,'+','2','x')")
    con.close()
    return str(path)


def _conc():
    return pa.Table.from_pylist(
        [
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
    )


def test_annotated_false_negatives_carry_context(tmp_path):
    ann = _ann(tmp_path / "ann.duckdb")
    fns = annotated_errors(_conc(), ann, kind="FN")
    assert len(fns) == 1
    r = fns[0]
    assert (r["chrom"], r["pos"], r["ref"], r["alt"]) == ("chr1", 200, "C", "T")
    assert r["gene"] == "GENEB"
    assert r["clin_sig"] == "Pathogenic"
    assert abs(r["af"] - 0.20) < 1e-9


def test_annotated_false_positives(tmp_path):
    ann = _ann(tmp_path / "ann.duckdb")
    fps = annotated_errors(_conc(), ann, kind="FP")
    assert len(fps) == 1
    assert fps[0]["pos"] == 300 and fps[0]["gene"] is None  # unannotated FP
