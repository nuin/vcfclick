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

from ingest._arrow import GENOTYPES_COLUMNS, VARIANTS_COLUMNS
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
