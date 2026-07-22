"""Multi-caller / cohort-scale benchmarking — the axis hap.py (one-pair) lacks.

Benchmark several callers against one truth in a single pass, then answer
cross-caller questions: per-caller recall/precision, and which variants a
caller misses that others catch.
"""

from __future__ import annotations

from pathlib import Path

from benchmark.cohort import (
    benchmark_callers,
    per_caller_metrics,
    variants_missed_by,
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "benchmark"
REF = str(FIXTURES / "tiny.fa")

_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=6>\n"
    "##contig=<ID=chr2,length=8>\n"
    '##FILTER=<ID=PASS,Description="x">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
)

_SNP = "chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
_INDEL = "chr2\t3\t.\tGT\tG\t50\tPASS\t.\tGT\t0/1\n"
_FP = "chr1\t4\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n"


def _vcf(path, records):
    path.write_text(_HEADER + "".join(records))
    return str(path)


def test_multi_caller_per_metrics_and_missed(tmp_path):
    truth = _vcf(tmp_path / "truth.vcf", [_SNP, _INDEL])
    callers = {
        "good": _vcf(tmp_path / "good.vcf", [_SNP, _INDEL]),  # perfect
        "no_indel": _vcf(tmp_path / "no_indel.vcf", [_SNP]),  # misses the indel
        "fp": _vcf(tmp_path / "fp.vcf", [_SNP, _INDEL, _FP]),  # one extra FP SNP
    }
    frame = benchmark_callers(truth, REF, callers, regions=None)
    m = {(r["caller"], r["vtype"]): r for r in per_caller_metrics(frame)}

    assert m[("good", "SNP")]["recall"] == 1.0
    assert m[("good", "INDEL")]["recall"] == 1.0
    assert m[("no_indel", "INDEL")]["recall"] == 0.0  # missed the indel
    assert m[("fp", "SNP")]["precision"] == 0.5  # one TP + one FP SNP

    # which variants does 'no_indel' miss that others catch?
    missed = variants_missed_by(frame, "no_indel", min_others=1)
    assert ("chr2", 3, "GT", "G") in missed
    assert ("chr1", 2, "A", "G") not in missed  # everyone got the SNP


def test_regression_history_append_and_load(tmp_path):
    from benchmark.cohort import append_run, load_history

    hist = tmp_path / "history.csv"
    append_run(
        str(hist),
        "v1",
        [{"caller": "gatk", "vtype": "SNP", "recall": 0.90, "precision": 0.98}],
    )
    append_run(
        str(hist),
        "v2",
        [{"caller": "gatk", "vtype": "SNP", "recall": 0.94, "precision": 0.99}],
    )
    rows = load_history(str(hist))
    v = {r["label"]: r for r in rows}
    assert v["v1"]["recall"] == "0.9"
    assert v["v2"]["recall"] == "0.94"
    assert len(rows) == 2  # one row per (label, caller, vtype)


def test_benchmark_cohort_cli(tmp_path):
    from click.testing import CliRunner

    from cli.benchmark import benchmark_cohort

    truth = _vcf(tmp_path / "truth.vcf", [_SNP, _INDEL])
    a = _vcf(tmp_path / "a.vcf", [_SNP, _INDEL])
    b = _vcf(tmp_path / "b.vcf", [_SNP])
    out = tmp_path / "out"
    res = CliRunner().invoke(
        benchmark_cohort,
        [
            "--truth",
            truth,
            "--ref",
            REF,
            "-o",
            str(out),
            "--caller",
            f"good={a}",
            "--caller",
            f"bad={b}",
            "--history",
            str(tmp_path / "h.csv"),
            "--label",
            "v1",
        ],
    )
    assert res.exit_code == 0, res.output
    assert (out / "per_caller_metrics.csv").exists()
    assert "bad" in res.output and "misses" in res.output
    assert (tmp_path / "h.csv").exists()
