from __future__ import annotations

import csv
import json
import os

from benchmark.constants import FILTER_ALL, FILTER_PASS, VT_INDEL, VT_SNP
from benchmark.report import write_reports

CANONICAL_COLUMNS = [
    "Type",
    "Filter",
    "TRUTH.TOTAL",
    "TRUTH.TP",
    "TRUTH.FN",
    "QUERY.TOTAL",
    "QUERY.TP",
    "QUERY.FP",
    "QUERY.UNK",
    "FP.gt",
    "FP.al",
    "METRIC.Recall",
    "METRIC.Precision",
    "METRIC.Frac_NA",
    "METRIC.F1_Score",
    "TRUTH.TOTAL.TiTv_ratio",
    "QUERY.TOTAL.TiTv_ratio",
    "TRUTH.TOTAL.het_hom_ratio",
    "QUERY.TOTAL.het_hom_ratio",
]


def _agg():
    """Known agg dict: SNP PASS 9 TP / 1 FN truth, 9 TP / 1 FP / 2 UNK query."""
    return {
        (FILTER_PASS, VT_SNP): {
            "truth_tp": 9,
            "truth_fn": 1,
            "query_tp": 9,
            "query_fp": 1,
            "query_unk": 2,
        },
        (FILTER_ALL, VT_SNP): {
            "truth_tp": 10,
            "truth_fn": 0,
            "query_tp": 10,
            "query_fp": 2,
            "query_unk": 2,
        },
        (FILTER_PASS, VT_INDEL): {
            "truth_tp": 3,
            "truth_fn": 1,
            "query_tp": 3,
            "query_fp": 0,
            "query_unk": 0,
        },
        (FILTER_ALL, VT_INDEL): {
            "truth_tp": 4,
            "truth_fn": 0,
            "query_tp": 4,
            "query_fp": 1,
            "query_unk": 1,
        },
    }


def _meta():
    return {"engine": "normalized", "ref": "tiny.fa", "regions": "conf.bed"}


def _read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.reader(fh))


def test_summary_csv_header_is_canonical_no_engine(tmp_path):
    write_reports(_agg(), _meta(), str(tmp_path), ["csv"])
    rows = _read_csv(os.path.join(str(tmp_path), "summary.csv"))
    assert rows[0] == CANONICAL_COLUMNS
    assert "Engine" not in rows[0]


def test_summary_csv_values_correct(tmp_path):
    write_reports(_agg(), _meta(), str(tmp_path), ["csv"])
    rows = _read_csv(os.path.join(str(tmp_path), "summary.csv"))
    header = rows[0]
    idx = {c: i for i, c in enumerate(header)}
    by_key = {(r[idx["Type"]], r[idx["Filter"]]): r for r in rows[1:]}

    def val(row, col):
        return row[idx[col]]

    sp = by_key[("SNP", "PASS")]
    assert val(sp, "TRUTH.TOTAL") == "10"  # 9 + 1
    assert val(sp, "TRUTH.TP") == "9"
    assert val(sp, "TRUTH.FN") == "1"
    assert val(sp, "QUERY.TOTAL") == "12"  # 9 + 1 + 2
    assert val(sp, "QUERY.FP") == "1"
    assert val(sp, "QUERY.UNK") == "2"
    assert float(val(sp, "METRIC.Recall")) == 0.9
    assert float(val(sp, "METRIC.Precision")) == 0.9
    assert abs(float(val(sp, "METRIC.Frac_NA")) - (2 / 12)) < 1e-9

    indel = by_key[("INDEL", "ALL")]
    assert val(indel, "TRUTH.TOTAL") == "4"
    assert float(val(indel, "METRIC.Recall")) == 1.0
    assert abs(float(val(indel, "METRIC.Precision")) - (4 / 5)) < 1e-9

    # hap.py columns present; agg has no subtype/blt so ratios are NaN and FP
    # is all-allele (FP.gt = 0, FP.al = QUERY.FP).
    assert val(sp, "FP.gt") == "0"
    assert val(sp, "FP.al") == "1"
    assert val(sp, "TRUTH.TOTAL.TiTv_ratio") == "NaN"


def test_summary_csv_has_four_rows(tmp_path):
    write_reports(_agg(), _meta(), str(tmp_path), ["csv"])
    rows = _read_csv(os.path.join(str(tmp_path), "summary.csv"))
    assert len(rows) == 5  # header + SNP/INDEL x PASS/ALL


def test_vcfclick_summary_first_column_is_engine(tmp_path):
    write_reports(_agg(), _meta(), str(tmp_path), ["csv"])
    rows = _read_csv(os.path.join(str(tmp_path), "vcfclick.summary.csv"))
    assert rows[0][0] == "Engine"
    assert rows[0][1:] == CANONICAL_COLUMNS
    assert rows[1][0] == "normalized"


def test_run_meta_json_round_trips(tmp_path):
    write_reports(_agg(), _meta(), str(tmp_path), ["json"])
    with open(os.path.join(str(tmp_path), "run_meta.json")) as fh:
        loaded = json.load(fh)
    assert loaded == _meta()


def test_metrics_json_written(tmp_path):
    write_reports(_agg(), _meta(), str(tmp_path), ["json"])
    with open(os.path.join(str(tmp_path), "metrics.json")) as fh:
        loaded = json.load(fh)
    assert isinstance(loaded, list)
    snp_pass = next(r for r in loaded if r["Type"] == "SNP" and r["Filter"] == "PASS")
    assert snp_pass["TRUTH.TOTAL"] == 10
    assert snp_pass["METRIC.Recall"] == 0.9


def test_index_html_is_self_contained(tmp_path):
    write_reports(_agg(), _meta(), str(tmp_path), ["html"])
    html = open(os.path.join(str(tmp_path), "index.html")).read()
    assert "http://" not in html
    assert "https://" not in html
    assert "src=" not in html
    assert "normalized" in html  # engine stamped
    assert "not hap.py-comparable" in html  # INDEL caveat


def test_formats_subset_honored(tmp_path):
    write_reports(_agg(), _meta(), str(tmp_path), ["json"])
    files = set(os.listdir(str(tmp_path)))
    assert "summary.csv" not in files
    assert "index.html" not in files
    assert "run_meta.json" in files
