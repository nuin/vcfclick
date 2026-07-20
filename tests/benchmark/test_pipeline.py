"""Integration tests for benchmark/pipeline.py and cli/benchmark.py.

Hand-built truth/query VCFs over the tiny.fa coordinates (chr1=CAAAAT,
chr2=ACGTACGT, chrM=GGGGCCCC). Self-benchmark (truth==query) must score 1.0;
a degraded query recovers exact injected FP/FN counts.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from benchmark.constants import FILTER_ALL, VT_INDEL, VT_SNP
from benchmark.pipeline import run_benchmark
from benchmark.reconcile import UnsupportedFeatureError

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


def test_haplotype_engine_unsupported(tmp_path):
    truth = _write_vcf(tmp_path / "truth.vcf", _TRUTH)
    query = _write_vcf(tmp_path / "query.vcf", _TRUTH)
    with pytest.raises(UnsupportedFeatureError):
        run_benchmark(truth, query, REF, str(tmp_path / "out"), engine="haplotype")


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


def test_cli_haplotype_errors_cleanly(tmp_path):
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
    assert res.exit_code != 0
    assert "haplotype" in res.output


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
    # An MNP (AC>GT) types as COMPLEX; summary.csv stays strict SNP/INDEL, but
    # the COMPLEX count must be visible in run_meta (never silently dropped).
    rec = ["chr2\t1\t.\tAC\tGT\t50\tPASS\t.\tGT\t1/1\n"]
    truth = _write_vcf(tmp_path / "truth.vcf", rec)
    query = _write_vcf(tmp_path / "query.vcf", rec)
    res = run_benchmark(truth, query, REF, str(tmp_path / "out"), regions=None)
    assert res["run_meta"]["unsummarized_types"].get("COMPLEX", 0) >= 1


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
    assert off["run_meta"]["unsummarized_types"].get("COMPLEX", 0) >= 1
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
