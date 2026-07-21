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
    "METRIC.Recall",
    "METRIC.Precision",
    "METRIC.Frac_NA",
    "METRIC.F1_Score",
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
    by_key = {(r[0], r[1]): r for r in rows[1:]}

    snp_pass = by_key[("SNP", "PASS")]
    assert snp_pass[2] == "10"  # TRUTH.TOTAL = 9 + 1
    assert snp_pass[3] == "9"  # TRUTH.TP
    assert snp_pass[4] == "1"  # TRUTH.FN
    assert snp_pass[5] == "12"  # QUERY.TOTAL = 9 + 1 + 2
    assert snp_pass[6] == "9"  # QUERY.TP
    assert snp_pass[7] == "1"  # QUERY.FP
    assert snp_pass[8] == "2"  # QUERY.UNK
    assert float(snp_pass[9]) == 0.9  # Recall = 9/10
    assert float(snp_pass[10]) == 0.9  # Precision = 9/10
    assert abs(float(snp_pass[11]) - (2 / 12)) < 1e-9  # Frac_NA = 2/12
    assert abs(float(snp_pass[12]) - 0.9) < 1e-9  # F1 = 2*.9*.9/1.8

    indel_all = by_key[("INDEL", "ALL")]
    assert indel_all[2] == "4"  # TRUTH.TOTAL
    assert float(indel_all[9]) == 1.0  # Recall = 4/4
    assert abs(float(indel_all[10]) - (4 / 5)) < 1e-9  # Precision = 4/5


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
