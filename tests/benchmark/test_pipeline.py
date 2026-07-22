"""Integration tests for benchmark/pipeline.py and cli/benchmark.py.

Hand-built truth/query VCFs over the tiny.fa coordinates (chr1=CAAAAT,
chr2=ACGTACGT, chrM=GGGGCCCC). Self-benchmark (truth==query) must score 1.0;
a degraded query recovers exact injected FP/FN counts.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from benchmark.constants import FILTER_ALL, VT_COMPLEX, VT_INDEL, VT_SNP
from benchmark.pipeline import exact_vs_normalized_delta, run_benchmark

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "benchmark"
REF = str(FIXTURES / "tiny.fa")
CONF = str(FIXTURES / "conf.bed")

_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=6>\n"
    "##contig=<ID=chr2,length=8>\n"
    "##contig=<ID=chrM,length=8>\n"
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)

# Two SNPs on chr1 (pos2 A>G transition, pos5 A>C transversion) + one deletion
# on chr2 (pos3 GT>G). All three fall inside conf.bed.
_TRUTH = [
    "chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    "chr1\t5\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n",
    "chr2\t3\t.\tGT\tG\t50\tPASS\t.\tGT\t0/1\n",
]

# Drop the pos5 SNP (→ FN) and inject a pos3 A>T SNP (→ FP).
_DEGRADED = [
    "chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    "chr1\t3\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\n",
    "chr2\t3\t.\tGT\tG\t50\tPASS\t.\tGT\t0/1\n",
]


def _write_vcf(path: Path, records: list[str]) -> str:
    path.write_text(_HEADER + "".join(records))
    return str(path)


def _row(summary: list[dict], vtype: str, filter_view: str) -> dict:
    for r in summary:
        if r["Type"] == vtype and r["Filter"] == filter_view:
            return r
    raise AssertionError(f"no summary row for {vtype}/{filter_view}")


def test_self_benchmark_is_perfect(tmp_path):
    truth = _write_vcf(tmp_path / "truth.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _TRUTH)
    outdir = tmp_path / "out"

    res = run_benchmark(truth, query, REF, str(outdir), regions=CONF)
    summary = res["summary"]

    for vt in (VT_SNP, VT_INDEL):
        row = _row(summary, vt, FILTER_ALL)
        assert row["recall"] == pytest.approx(1.0)
        assert row["precision"] == pytest.approx(1.0)
        assert row["f1"] == pytest.approx(1.0)


def test_degraded_query_exact_counts(tmp_path):
    truth = _write_vcf(tmp_path / "truth.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _DEGRADED)
    outdir = tmp_path / "out"

    res = run_benchmark(truth, query, REF, str(outdir), regions=CONF)
    summary = res["summary"]

    snp = _row(summary, VT_SNP, FILTER_ALL)
    assert snp["truth_tp"] == 1
    assert snp["truth_fn"] == 1
    assert snp["query_tp"] == 1
    assert snp["query_fp"] == 1
    assert snp["recall"] == pytest.approx(0.5)
    assert snp["precision"] == pytest.approx(0.5)

    indel = _row(summary, VT_INDEL, FILTER_ALL)
    assert indel["truth_tp"] == 1
    assert indel["query_tp"] == 1
    assert indel["query_fp"] == 0
    assert indel["recall"] == pytest.approx(1.0)
    assert indel["precision"] == pytest.approx(1.0)


def test_reports_written(tmp_path):
    truth = _write_vcf(tmp_path / "truth.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _TRUTH)
    outdir = tmp_path / "out"

    run_benchmark(truth, query, REF, str(outdir), regions=CONF)

    for name in (
        "summary.csv",
        "vcfclick.summary.csv",
        "run_meta.json",
        "metrics.json",
        "index.html",
    ):
        assert os.path.exists(outdir / name), f"missing {name}"


# chr2 = ACGTACGT. Truth MNP AC>GT at pos1; query the two atomized SNPs A>G@1,
# C>T@2. Keyed-different (normalized FP+FN) but sequence-equivalent (haplotype TP).
_MNP_TRUTH = ["chr2\t1\t.\tAC\tGT\t50\tPASS\t.\tGT\t1/1\n"]
_MNP_QUERY = [
    "chr2\t1\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\n",
    "chr2\t2\t.\tC\tT\t50\tPASS\t.\tGT\t1/1\n",
]


def test_haplotype_engine_rescues_mnp_split(tmp_path):
    truth = _write_vcf(tmp_path / "truth.vcf", _MNP_TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _MNP_QUERY)

    norm = run_benchmark(
        truth, query, REF, str(tmp_path / "n"), regions=None, engine="normalized"
    )
    assert _row(norm["summary"], VT_SNP, FILTER_ALL)["query_fp"] == 2

    hap = run_benchmark(
        truth, query, REF, str(tmp_path / "h"), regions=None, engine="haplotype"
    )
    snp = _row(hap["summary"], VT_SNP, FILTER_ALL)
    assert snp["query_fp"] == 0
    assert snp["query_tp"] == 2


def test_haplotype_over_budget_not_reclassified(tmp_path):
    # Drive classify_haplotype directly with a small budget so the 3-member
    # cluster (1 truth + 2 query) is over budget and is never reclassified TP.
    import dataclasses

    from benchmark import stratify
    from benchmark.pipeline import _parse_side
    from benchmark.reconcile import classify_haplotype
    from benchmark.reference import Reference

    truth = _write_vcf(tmp_path / "truth.vcf", _MNP_TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _MNP_QUERY)
    ref = Reference(REF)
    excluded = {"symbolic": 0, "ref_mismatch": 0}

    def _prep(path, side):
        recs = _parse_side(path, side, ref, "error", excluded, False)
        recs = [dataclasses.replace(r, in_conf=True) for r in recs]
        return stratify.tag(recs)

    tr = _prep(truth, "truth")
    qr = _prep(query, "query")

    rows = classify_haplotype(tr, qr, FILTER_ALL, ref.fetch, max_cluster=2)
    query_tp = sum(
        1 for r in rows if r.side == "query" and r.bd == "TP" and r.vtype == VT_SNP
    )
    assert query_tp != 2  # over-budget cluster is not reclassified TP
    assert not any(r.bd == "TP" for r in rows)  # counted too-complex, never TP


def test_exact_engine_misses_right_shifted_indel(tmp_path):
    # chr1 = CAAAAT. Truth is a right-shifted 1bp deletion (AT>T at pos5); query
    # is its left-aligned form (CA>C at pos1). Normalized matches; exact misses.
    truth = _write_vcf(
        tmp_path / "truth.vcf", ["chr1\t5\t.\tAT\tT\t50\tPASS\t.\tGT\t0/1\n"]
    )
    query = _write_vcf(
        tmp_path / "query.vcf", ["chr1\t1\t.\tCA\tC\t50\tPASS\t.\tGT\t0/1\n"]
    )

    norm = run_benchmark(
        truth, query, REF, str(tmp_path / "n"), regions=None, engine="normalized"
    )
    n_indel = _row(norm["summary"], VT_INDEL, FILTER_ALL)
    assert n_indel["truth_tp"] == 1 and n_indel["query_fp"] == 0

    exact = run_benchmark(
        truth, query, REF, str(tmp_path / "e"), regions=None, engine="exact"
    )
    e_indel = _row(exact["summary"], VT_INDEL, FILTER_ALL)
    assert e_indel["truth_fn"] == 1 and e_indel["query_fp"] == 1


def test_exact_vs_normalized_delta_positive_on_shifted_indel(tmp_path):
    truth = _write_vcf(
        tmp_path / "truth.vcf", ["chr1\t5\t.\tAT\tT\t50\tPASS\t.\tGT\t0/1\n"]
    )
    query = _write_vcf(
        tmp_path / "query.vcf", ["chr1\t1\t.\tCA\tC\t50\tPASS\t.\tGT\t0/1\n"]
    )
    delta = exact_vs_normalized_delta(truth, query, REF, regions=None)
    assert delta[VT_INDEL]["recall"] > 0
    assert delta[VT_INDEL]["precision"] > 0


def test_cli_verb_runs(tmp_path):
    from click.testing import CliRunner

    from cli.benchmark import benchmark

    truth = _write_vcf(tmp_path / "truth.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _TRUTH)
    outdir = tmp_path / "out"
    res = CliRunner().invoke(
        benchmark,
        [
            "--truth",
            truth,
            "--query",
            query,
            "--ref",
            REF,
            "--regions",
            CONF,
            "-o",
            str(outdir),
        ],
    )
    assert res.exit_code == 0, res.output
    assert (outdir / "summary.csv").exists()
    assert "recall=1.0000" in res.output


def test_cli_haplotype_engine_runs(tmp_path):
    from click.testing import CliRunner

    from cli.benchmark import benchmark

    truth = _write_vcf(tmp_path / "truth.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _TRUTH)
    res = CliRunner().invoke(
        benchmark,
        [
            "--truth",
            truth,
            "--query",
            query,
            "--ref",
            REF,
            "-o",
            str(tmp_path / "out"),
            "--engine",
            "haplotype",
        ],
    )
    assert res.exit_code == 0, res.output
    assert "recall=1.0000" in res.output


def test_no_regions_marks_all_confident(tmp_path):
    # Query-only variant with no BED → in_conf True everywhere → scored FP,
    # never UNK. Recall stays 1.0 (truth all matched); precision drops.
    truth = _write_vcf(tmp_path / "truth.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _DEGRADED)
    res = run_benchmark(truth, query, REF, str(tmp_path / "out"), regions=None)
    snp = _row(res["summary"], VT_SNP, FILTER_ALL)
    assert snp["query_unk"] == 0
    assert snp["query_fp"] == 1


def test_homref_query_record_is_not_an_fp(tmp_path):
    # A 0/0 (homref) query record is not a variant call; it must not score as FP.
    truth = _write_vcf(tmp_path / "truth.vcf", _TRUTH)
    query = _write_vcf(
        tmp_path / "query.vcf",
        _TRUTH + ["chr2\t5\t.\tA\tC\t50\tPASS\t.\tGT\t0/0\n"],
    )
    res = run_benchmark(truth, query, REF, str(tmp_path / "out"), regions=None)
    snp = _row(res["summary"], VT_SNP, FILTER_ALL)
    assert snp["query_fp"] == 0


def test_chrm_haploid_matches_diploid(tmp_path):
    # chrM is always haploid: a haploid truth `1` must match a diploid `1/1`
    # query as a TP (chrX/Y are deliberately NOT auto-equated).
    truth = _write_vcf(
        tmp_path / "truth.vcf", ["chrM\t5\t.\tC\tA\t50\tPASS\t.\tGT\t1\n"]
    )
    query = _write_vcf(
        tmp_path / "query.vcf", ["chrM\t5\t.\tC\tA\t50\tPASS\t.\tGT\t1/1\n"]
    )
    res = run_benchmark(truth, query, REF, str(tmp_path / "out"), regions=None)
    snp = _row(res["summary"], VT_SNP, FILTER_ALL)
    assert snp["truth_tp"] == 1 and snp["query_tp"] == 1
    assert snp["truth_fn"] == 0 and snp["query_fp"] == 0


def test_complex_variants_surfaced_in_run_meta(tmp_path):
    # An MNP (AC>GT) types as COMPLEX (BVT UNK); summary.csv stays strict
    # SNP/INDEL, but the count must be visible in run_meta (never silently dropped).
    rec = ["chr2\t1\t.\tAC\tGT\t50\tPASS\t.\tGT\t1/1\n"]
    truth = _write_vcf(tmp_path / "truth.vcf", rec)
    query = _write_vcf(tmp_path / "query.vcf", rec)
    res = run_benchmark(truth, query, REF, str(tmp_path / "out"), regions=None)
    assert res["run_meta"]["unsummarized_types"].get(VT_COMPLEX, 0) >= 1


def test_strict_requires_regions(tmp_path):
    from benchmark.reference import BenchmarkError

    t = _write_vcf(tmp_path / "t.vcf", _TRUTH)
    q = _write_vcf(tmp_path / "q.vcf", _TRUTH)
    with pytest.raises(BenchmarkError):
        run_benchmark(t, q, REF, str(tmp_path / "out"), regions=None, strict=True)


def test_parquet_format_written(tmp_path):
    t = _write_vcf(tmp_path / "t.vcf", _TRUTH)
    q = _write_vcf(tmp_path / "q.vcf", _TRUTH)
    outdir = tmp_path / "out"
    run_benchmark(t, q, REF, str(outdir), regions=CONF, report_formats=["parquet"])
    assert (outdir / "benchmark.parquet").exists()


def test_conf_containment_full_makes_straddling_indel_unk(tmp_path):
    # chr2 conf region is [0,4); a pos4 TA>T deletion spans 0-based [3,5) — its
    # start is inside but its span straddles the edge.
    q = _write_vcf(tmp_path / "q.vcf", ["chr2\t4\t.\tTA\tT\t50\tPASS\t.\tGT\t0/1\n"])
    t = _write_vcf(tmp_path / "t.vcf", [])
    start = run_benchmark(
        t, q, REF, str(tmp_path / "os"), regions=CONF, conf_containment="start"
    )
    full = run_benchmark(
        t, q, REF, str(tmp_path / "of"), regions=CONF, conf_containment="full"
    )
    assert _row(start["summary"], VT_INDEL, FILTER_ALL)["query_fp"] == 1
    assert _row(full["summary"], VT_INDEL, FILTER_ALL)["query_unk"] == 1


def test_decompose_mnp_flag_splits_into_snps(tmp_path):
    rec = ["chr2\t1\t.\tAC\tGT\t50\tPASS\t.\tGT\t1/1\n"]  # MNP AC>GT
    t = _write_vcf(tmp_path / "t.vcf", rec)
    q = _write_vcf(tmp_path / "q.vcf", rec)
    off = run_benchmark(
        t, q, REF, str(tmp_path / "o1"), regions=None, decompose_mnp=False
    )
    assert off["run_meta"]["unsummarized_types"].get(VT_COMPLEX, 0) >= 1
    on = run_benchmark(
        t, q, REF, str(tmp_path / "o2"), regions=None, decompose_mnp=True
    )
    assert _row(on["summary"], VT_SNP, FILTER_ALL)["truth_tp"] == 2


def _invoke(args):
    from click.testing import CliRunner

    from cli.benchmark import benchmark

    return CliRunner().invoke(benchmark, args)


def test_cli_pass_only_headline(tmp_path):
    t = _write_vcf(tmp_path / "t.vcf", _TRUTH)
    q = _write_vcf(tmp_path / "q.vcf", _TRUTH)
    res = _invoke(
        [
            "--truth",
            t,
            "--query",
            q,
            "--ref",
            REF,
            "--regions",
            CONF,
            "-o",
            str(tmp_path / "out"),
            "--pass-only",
        ]
    )
    assert res.exit_code == 0, res.output
    assert "PASS" in res.output and " ALL " not in res.output


def test_cli_parquet_format(tmp_path):
    t = _write_vcf(tmp_path / "t.vcf", _TRUTH)
    q = _write_vcf(tmp_path / "q.vcf", _TRUTH)
    outdir = tmp_path / "out"
    res = _invoke(
        [
            "--truth",
            t,
            "--query",
            q,
            "--ref",
            REF,
            "--regions",
            CONF,
            "-o",
            str(outdir),
            "--report-formats",
            "parquet",
        ]
    )
    assert res.exit_code == 0, res.output
    assert (outdir / "benchmark.parquet").exists()


def test_cli_strict_without_regions_errors(tmp_path):
    t = _write_vcf(tmp_path / "t.vcf", _TRUTH)
    q = _write_vcf(tmp_path / "q.vcf", _TRUTH)
    res = _invoke(
        [
            "--truth",
            t,
            "--query",
            q,
            "--ref",
            REF,
            "-o",
            str(tmp_path / "out"),
            "--strict",
        ]
    )
    assert res.exit_code == 1  # ClickException, not a usage error (2)
    assert "strict" in res.output.lower()


def test_stratify_writes_annotation_csvs(tmp_path, monkeypatch):
    # End-to-end: benchmark then stratify the concordance against a gnomAD store.
    import duckdb

    annp = tmp_path / "ann.duckdb"
    con = duckdb.connect(str(annp))
    con.execute(
        "CREATE TABLE gnomad_af (chrom VARCHAR, pos UINTEGER, ref VARCHAR, "
        "alt VARCHAR, af DOUBLE, af_grpmax DOUBLE)"
    )
    con.execute("INSERT INTO gnomad_af VALUES ('chr1',2,'A','G',0.0001,0.0002)")
    con.close()
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(annp))

    truth = _write_vcf(tmp_path / "truth.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _TRUTH)
    outdir = tmp_path / "out"
    run_benchmark(truth, query, REF, str(outdir), regions=CONF, stratify=["gnomad"])
    assert (outdir / "stratified_gnomad.csv").exists()
    text = (outdir / "stratified_gnomad.csv").read_text()
    assert "stratum" in text and "recall" in text


def test_roc_writes_tsv(tmp_path):
    truth = _write_vcf(tmp_path / "t.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "q.vcf", _TRUTH)
    outdir = tmp_path / "out"
    run_benchmark(truth, query, REF, str(outdir), regions=CONF, roc=True)
    roc = outdir / "roc.tsv"
    assert roc.exists()
    head = roc.read_text().splitlines()[0].split("\t")
    assert head == ["Type", "Threshold", "TP", "FP", "Recall", "Precision"]


def test_strat_region_writes_csv(tmp_path):
    bed = tmp_path / "lc.bed"
    bed.write_text("chr1\t0\t6\n")  # covers all of chr1 (len 6)
    truth = _write_vcf(tmp_path / "t.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "q.vcf", _TRUTH)
    outdir = tmp_path / "out"
    run_benchmark(
        truth,
        query,
        REF,
        str(outdir),
        regions=CONF,
        strat_regions={"lowcomplex": str(bed)},
    )
    out = outdir / "stratified_regions.csv"
    assert out.exists()
    text = out.read_text()
    assert "lowcomplex" in text and "stratum" in text
