"""Tests for scripts/compare_somalier.py — pure file comparison, no
database and no somalier binary needed."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

pytest.importorskip("numpy")

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "compare_somalier.py"
_spec = importlib.util.spec_from_file_location("compare_somalier", SCRIPT)
cs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cs)

PAIRS_HDR = "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\texpected_relatedness\n"
SAMPLES_HDR = "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tX_het\tX_hom_alt\n"


def _pair(a, b, kinship, rel="unrelated", ibs0=0.02, n=1000):
    return {
        "sample_a": a,
        "sample_b": b,
        "kinship": kinship,
        "ibs0": ibs0,
        "n_sites": n,
        "relationship": rel,
    }


@pytest.fixture
def files(tmp_path):
    rel = [
        {
            "ingest_id": "g1",
            "pairs": [
                _pair("HG002", "HG003", 0.25, "parent-child", 0.0),
                _pair("HG002", "HG004", 0.25, "parent-child", 0.0),
                _pair("HG003", "HG004", 0.0, "unrelated"),
            ],
        }
    ]
    qc = [
        {"sample_id": "HG002", "inferred_sex": "male"},
        {"sample_id": "HG003", "inferred_sex": "male"},
        {"sample_id": "HG004", "inferred_sex": "female"},
    ]
    # somalier: relatedness ≈ 2 x kinship, ibs0 a count over n; pair order
    # reversed on purpose to prove pairs are matched order-insensitively.
    pairs = PAIRS_HDR + (
        "HG003\tHG002\t0.49\t0\t500\t1000\t-1\n"
        "HG002\tHG004\t0.51\t1\t500\t1000\t-1\n"
        "HG003\tHG004\t0.01\t40\t100\t1000\t-1\n"
    )
    samples = SAMPLES_HDR + (
        "f\tHG002\t0\t0\t1\t2\t300\n"
        "f\tHG003\t0\t0\t1\t5\t280\n"
        "f\tHG004\t0\t0\t2\t150\t160\n"
    )
    paths = {
        "rel": tmp_path / "rel.json",
        "qc": tmp_path / "qc.json",
        "pairs": tmp_path / "somalier.pairs.tsv",
        "samples": tmp_path / "somalier.samples.tsv",
    }
    paths["rel"].write_text(json.dumps(rel))
    paths["qc"].write_text(json.dumps(qc))
    paths["pairs"].write_text(pairs)
    paths["samples"].write_text(samples)
    return paths


def _argv(p):
    return [
        "--vcfclick-relatedness", str(p["rel"]),
        "--somalier-pairs", str(p["pairs"]),
        "--vcfclick-qc", str(p["qc"]),
        "--somalier-samples", str(p["samples"]),
    ]  # fmt: skip


def test_somalier_relatedness_is_rescaled_to_kinship(files):
    so = cs.somalier_pairs(files["pairs"])
    assert so[("HG002", "HG003")]["kinship"] == pytest.approx(0.245)
    assert so[("HG002", "HG003")]["relationship"] == "parent-child"
    assert so[("HG003", "HG004")]["relationship"] == "unrelated"


def test_somalier_sex_from_x_het_fraction(files):
    assert cs.somalier_sex(files["samples"]) == {
        "HG002": "male",
        "HG003": "male",
        "HG004": "female",
    }


def test_agreement_exits_zero(files, capsys):
    assert cs.main(_argv(files)) == 0
    assert "no disagreements" in capsys.readouterr().out


def test_disagreement_exits_one(files, capsys):
    # somalier sees HG004 as male and HG003/HG004 as full siblings
    files["samples"].write_text(
        SAMPLES_HDR
        + "f\tHG002\t0\t0\t1\t2\t300\nf\tHG003\t0\t0\t1\t5\t280\nf\tHG004\t0\t0\t2\t3\t300\n"
    )
    files["pairs"].write_text(
        PAIRS_HDR
        + "HG003\tHG002\t0.49\t0\t500\t1000\t-1\n"
        + "HG002\tHG004\t0.51\t1\t500\t1000\t-1\n"
        + "HG003\tHG004\t0.50\t30\t400\t1000\t-1\n"
    )
    assert cs.main(_argv(files)) == 1
    err = capsys.readouterr().err
    assert "HG003/HG004" in err and "full-siblings" in err
    assert "HG004: sex vcfclick female, somalier male" in err


def test_insufficient_data_is_not_a_disagreement(files):
    rel = json.loads(files["rel"].read_text())
    rel[0]["pairs"][2] = _pair("HG003", "HG004", None, "insufficient-data", None, 0)
    files["rel"].write_text(json.dumps(rel))
    assert cs.main(_argv(files)) == 0


def test_relatives_vcfclick_did_not_report_are_a_disagreement(files, capsys):
    # vcfclick output without --all: only the related pairs are listed.
    # somalier also sees HG003/HG004 as first-degree relatives.
    rel = json.loads(files["rel"].read_text())
    del rel[0]["pairs"][2]
    files["rel"].write_text(json.dumps(rel))
    assert cs.main(_argv(files)) == 0  # somalier: unrelated → nothing to flag

    files["pairs"].write_text(
        PAIRS_HDR
        + "HG003\tHG002\t0.49\t0\t500\t1000\t-1\n"
        + "HG002\tHG004\t0.51\t1\t500\t1000\t-1\n"
        + "HG003\tHG004\t0.50\t30\t400\t1000\t-1\n"
    )
    assert cs.main(_argv(files)) == 1
    assert "HG003/HG004: somalier full-siblings, not reported by vcfclick" in (
        capsys.readouterr().err
    )
