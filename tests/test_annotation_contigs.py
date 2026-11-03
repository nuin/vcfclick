"""ClinVar and gnomAD lookups match a contig however the cohort spells it.

The store holds `chr17`, `chrX`, `chrM` (the loaders write the canonical
`storage.contigs.canonical_chrom` name). A cohort from an Ensembl-style VCF
says `17`, `X`, `MT`. Every lookup tries each alias and answers under the
caller's spelling, so the `--gnomad-max-af` gate and the ClinVar column work
for both.
"""

from __future__ import annotations

import json


def _seed(rows_clinvar=(), rows_gnomad=()):
    from annotations.db import get_connection

    conn = get_connection()
    for r in rows_clinvar:
        conn.execute(
            "INSERT INTO clinvar_variants VALUES (?, ?, ?, ?, ?, 'reviewed', 'V1', 'c')",
            list(r),
        )
    for r in rows_gnomad:
        conn.execute("INSERT INTO gnomad_af VALUES (?, ?, ?, ?, ?, ?)", list(r))
    conn.close()


def test_point_lookups_take_any_contig_alias(isolated_annotation_db):
    from annotations.db import clinvar_lookup, gnomad_af

    _seed(
        [("chr17", 150, "A", "G", "Pathogenic"), ("chrM", 73, "A", "G", "Benign")],
        [("chr17", 150, "A", "G", 0.001, 0.002)],
    )
    c = clinvar_lookup("17", 150, "A", "G")
    assert c is not None and c.clin_sig == "Pathogenic" and c.chrom == "17"
    assert clinvar_lookup("MT", 73, "A", "G").clin_sig == "Benign"
    assert clinvar_lookup("chr17", 150, "A", "G").chrom == "chr17"
    g = gnomad_af("17", 150, "A", "G")
    assert g is not None and g.af_grpmax == 0.002 and g.chrom == "17"


def test_batch_lookups_key_by_the_callers_spelling(isolated_annotation_db):
    from annotations.db import clinvar_many, gnomad_af_many

    _seed(
        [("chr17", 150, "A", "G", "Pathogenic"), ("chrX", 9, "C", "T", "Benign")],
        [("chr17", 150, "A", "G", 0.001, 0.002), ("chrX", 9, "C", "T", 0.3, 0.4)],
    )
    keys = [("17", 150, "A", "G"), ("X", 9, "C", "T"), ("chr17", 150, "A", "G")]
    g = gnomad_af_many(keys)
    assert set(g) == set(keys)
    assert g[("X", 9, "C", "T")].af == 0.3
    assert g[("X", 9, "C", "T")].chrom == "X"
    c = clinvar_many(keys)
    assert set(c) == set(keys)
    assert c[("17", 150, "A", "G")].clin_sig == "Pathogenic"


def test_an_older_store_with_bare_contigs_still_answers(isolated_annotation_db):
    """gnomAD loaded before the loader canonicalised (a bare-contig file)."""
    from annotations.db import gnomad_af, gnomad_af_many

    _seed(rows_gnomad=[("7", 100, "A", "G", 0.01, 0.02)])
    assert gnomad_af("chr7", 100, "A", "G").af == 0.01
    assert ("chr7", 100, "A", "G") in gnomad_af_many([("chr7", 100, "A", "G")])


def test_gnomad_loader_writes_canonical_contigs(tmp_path, isolated_annotation_db):
    from annotations.db import read_connection
    from annotations.loaders.gnomad import load

    vcf = tmp_path / "bare.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=7>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "7\t100\t.\tA\tG\t.\tPASS\tAF=0.5\n"
    )
    load(vcf, replace=True)
    conn = read_connection()
    try:
        assert conn.execute("SELECT chrom FROM gnomad_af").fetchall() == [("chr7",)]
    finally:
        conn.close()


def test_clinvar_contigs_are_canonical():
    from annotations.loaders.clinvar import _normalise_contig

    assert _normalise_contig("17") == "chr17"
    assert _normalise_contig("MT") == "chrM"
    assert _normalise_contig("X") == "chrX"
    assert _normalise_contig("chr7") == "chr7"
    # an unplaced contig is not given a `chr` prefix it never had
    assert _normalise_contig("NW_009646201.1") == "NW_009646201.1"


def test_gnomad_coverage_reads_contig_aliases():
    from cli.effective_gates import VariantStats, _gnomad_coverage

    variants = VariantStats(n=10, with_af=0, chroms=(("17", 6), ("X", 4)))
    share, cohort = _gnomad_coverage(variants, {"chr17", "chrX"})
    assert share == 1.0 and cohort == ["17", "X"]


def _positions(doc: dict) -> set[int]:
    return {
        int(c["pos"])
        for m in doc["models"].values()
        for c in m.get("candidates", []) or []
    }


def test_trio_gnomad_gate_judges_bare_contigs(tmp_path, isolated_annotation_db):
    """A cohort with `chr1` and gnomAD loaded as `1` (an older store): the
    gate drops the common variant at chr1:300 instead of keeping every
    candidate as absent."""
    from tests.test_trio import _setup, _trio

    home = tmp_path / "home"
    _setup(home)
    _seed(rows_gnomad=[("1", 300, "G", "A", 0.5, 0.5)])
    args = ("--category", "dominant", "--max-af", "0.3", "--format", "json")
    base = json.loads(_trio(home, *args).stdout)
    gated = json.loads(_trio(home, *args, "--gnomad-max-af", "0.01").stdout)
    assert 300 in _positions(base)
    assert 300 not in _positions(gated)
