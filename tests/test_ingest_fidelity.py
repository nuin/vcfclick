"""Row-building fidelity: what one VCF record becomes in the tables.

The unit tests here read small plain-text VCFs with cyvcf2 in-process
(no storage session) and call the row builders in ingest.vcf_rows
directly. The end-to-end tests go through the CLI on both backends.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from cyvcf2 import VCF

from ingest._arrow import (
    GENOTYPES_ARROW_SCHEMA,
    GENOTYPES_COLUMNS,
    VARIANTS_ARROW_SCHEMA,
    VARIANTS_COLUMNS,
    write_parquet,
)
from ingest.vcf_rows import build_record_rows
from tests.conftest import bgzip_vcf, run_cli

HEADER = """\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000000>
##INFO=<ID=SB,Number=4,Type=Integer,Description="Strand bias">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Min depth">
##FORMAT=<ID=SBX,Number=4,Type=Integer,Description="Strand counts">
##FORMAT=<ID=FX,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=FS,Number=1,Type=String,Description="A string">
"""


def _write_vcf(path: Path, samples: list[str], records: list[str]) -> Path:
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    body = "\n".join([cols + "\t" + "\t".join(samples), *records]) + "\n"
    path.write_text(HEADER + body)
    return path


def _records(path: Path):
    vcf = VCF(str(path))
    return list(vcf.samples), list(vcf)


def _rows(path: Path, keep_reference: bool = False):
    samples, variants = _records(path)
    from ingest.routing import classify_header

    extra = classify_header(VCF(str(path)))["extra_format"]
    out = []
    for v in variants:
        out.append(build_record_rows(v, samples, extra, "i1", keep_reference))
    return samples, out


def _gt_dict(genotype_rows: list[list]) -> dict[str, dict]:
    idx = {c: i for i, c in enumerate(GENOTYPES_COLUMNS)}
    return {r[idx["sample_id"]]: dict(zip(GENOTYPES_COLUMNS, r)) for r in genotype_rows}


def _variant_dict(variant_row: list) -> dict:
    return dict(zip(VARIANTS_COLUMNS, variant_row))


# ─── finding 6: numeric FORMAT overflow is stored as the numpy repr ────


def test_numeric_format_extra_is_comma_joined(tmp_path):
    """MIN_DP=5 is '5', SBX=1,2,3,4 is '1,2,3,4', a missing value is '.',
    and a shorter vector drops cyvcf2's end-of-vector padding."""
    vcf = _write_vcf(
        tmp_path / "x.vcf",
        ["A", "B", "C"],
        [
            "chr1\t100\t.\tA\tG\t50\tPASS\tDP=9\tGT:MIN_DP:SBX:FX:FS"
            "\t0/1:5:1,2,3,4:1.5:foo\t1/1:.:.:.:.\t0/1:7:1,2:2.5:bar"
        ],
    )
    _, rows = _rows(vcf)
    gts = _gt_dict(rows[0][1])
    assert gts["A"]["format_extra"] == {
        "MIN_DP": "5",
        "SBX": "1,2,3,4",
        "FX": "1.5",
        "FS": "foo",
    }
    assert gts["B"]["format_extra"] == {"MIN_DP": ".", "SBX": ".", "FX": ".", "FS": "."}
    assert gts["C"]["format_extra"]["SBX"] == "1,2"
    assert gts["C"]["format_extra"]["FX"] == "2.5"


@pytest.mark.parametrize("backend", ["duckdb", "chdb"])
def test_numeric_format_extra_queryable_after_ingest(vcfclick_home, tmp_path, backend):
    """CAST(format_extra['MIN_DP'] AS INT) works on a stored value."""
    src = _write_vcf(
        tmp_path / "fx.vcf",
        ["A", "B"],
        ["chr1\t100\t.\tA\tG\t50\tPASS\tDP=9\tGT:MIN_DP:SBX\t0/1:5:1,2,3,4\t1/1:.:."],
    )
    gz = bgzip_vcf(src, tmp_path / "fx.vcf.gz")
    run_cli(vcfclick_home, backend, "db", "create", "fx")
    run_cli(vcfclick_home, backend, "db", "ingest", "fx", str(gz), "--serial")
    r = run_cli(
        vcfclick_home,
        backend,
        "db",
        "query",
        "fx",
        "SELECT sample_id, format_extra['MIN_DP'], format_extra['SBX'] "
        "FROM genotypes ORDER BY sample_id",
        "--format",
        "JSONCompact",
    )
    assert json.loads(r.stdout)["data"] == [["A", "5", "1,2,3,4"], ["B", ".", "."]]


# ─── finding 36: a multi-value INFO field aimed at a scalar column ──────


def test_multi_value_info_scalar_goes_to_info_extra(tmp_path):
    """INFO SB declared Number=4 (GATK) cannot fill the float info_SB
    column: the tuple lands in info_extra and the row still writes."""
    vcf = _write_vcf(
        tmp_path / "sb.vcf",
        ["A"],
        ["chr1\t100\t.\tA\tG\t50\tPASS\tSB=1,2,3,4;AC=1,0;DP=9\tGT\t0/1"],
    )
    _, rows = _rows(vcf)
    variant = _variant_dict(rows[0][0])
    assert variant["info_SB"] is None
    assert variant["info_AC"] is None
    assert variant["info_DP"] == 9
    assert variant["info_extra"] == {"SB": "1,2,3,4", "AC": "1,0"}
    write_parquet([rows[0][0]], VARIANTS_ARROW_SCHEMA, tmp_path / "v.parquet")


# ─── finding 37: FORMAT integers above the uint16 column range ──────────


def test_wide_format_integers_clamp_to_uint16(tmp_path):
    """DP=70000 and AD=80000,5 are clamped to 65535, so the batch writes
    instead of aborting the ingest with an ArrowInvalid."""
    vcf = _write_vcf(
        tmp_path / "dp.vcf",
        ["A"],
        [
            "chr1\t100\t.\tA\tG\t50\tPASS\tDP=9\tGT:DP:AD:PL:PS"
            "\t0|1:70000:80000,5:0,70000,99:200000000"
        ],
    )
    _, rows = _rows(vcf)
    gt = _gt_dict(rows[0][1])["A"]
    assert (gt["dp"], gt["ad_ref"], gt["ad_alt"]) == (65535, 65535, 5)
    assert (gt["pl_ref_ref"], gt["pl_ref_alt"], gt["pl_alt_alt"]) == (0, 65535, 99)
    assert gt["ps"] == 200000000  # ps is UInt32: a phase set keeps its value
    write_parquet(rows[0][1], GENOTYPES_ARROW_SCHEMA, tmp_path / "g.parquet")


# ─── finding 40: samples but no GT ──────────────────────────────────────


def test_record_without_gt_builds_no_genotype_rows(tmp_path):
    """FORMAT=DP only: called counts are "not recorded" and there are no
    genotype rows, instead of a raw cyvcf2 exception."""
    vcf = _write_vcf(
        tmp_path / "nogt.vcf",
        ["A"],
        ["chr1\t200\t.\tA\tG\t50\tPASS\tDP=9\tDP\t10"],
    )
    _, rows = _rows(vcf)
    variant_row, genotype_rows, missing_rows = rows[0]
    variant = _variant_dict(variant_row)
    assert (variant["n_called"], variant["an_called"], variant["ac_called"]) == (
        None,
        None,
        None,
    )
    assert genotype_rows == []
    assert missing_rows == []


# ─── finding 39: a record with no ALT allele ────────────────────────────


def test_no_alt_record_names_the_remedy(tmp_path):
    """ALT '.' is not a multi-allelic site: the error says the record has
    no ALT allele and names the bcftools filter that drops such records.
    A real multi-allelic record keeps the norm remedy."""
    from ingest.vcf_rows import require_single_alt

    vcf = _write_vcf(
        tmp_path / "alt.vcf",
        ["A"],
        [
            "chr1\t100\t.\tA\t.\t50\tPASS\tDP=9\tGT\t0/0",
            "chr1\t200\t.\tA\tG,T\t50\tPASS\tDP=9\tGT\t1/2",
        ],
    )
    _, (ref_only, multi) = _records(vcf)
    with pytest.raises(ValueError) as e:
        require_single_alt(ref_only, "in.vcf.gz")
    msg = str(e.value)
    assert "chr1:100" in msg and "no ALT allele" in msg
    assert "bcftools view -e 'ALT=\".\"'" in msg
    assert "Multi-allelic" not in msg
    with pytest.raises(ValueError, match="Multi-allelic site at chr1:200"):
        require_single_alt(multi, "in.vcf.gz")
    assert "bcftools norm -m -" in str(
        pytest.raises(ValueError, require_single_alt, multi, "in.vcf.gz").value
    )
