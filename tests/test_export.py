"""Tests for `vcfclick db export` — a database (or a slice of it) back to VCF.

The core check is a round trip: ingest a VCF, export it, ingest the export
into a second database, and require identical variants and genotypes.
"""

from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
HAS_HTSLIB = shutil.which("bcftools") and shutil.which("tabix")


def _vc(home: Path, *args: str, ok: bool = True):
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    if ok:
        assert r.returncode == 0, f"{' '.join(args)} failed:\n{r.stderr}"
    return r


def _db(home: Path, vcf: Path, name: str, *extra: str) -> None:
    _vc(home, "db", "create", name)
    _vc(
        home,
        "db",
        "ingest",
        name,
        str(vcf),
        "--cohort",
        "c",
        "--ingest-id",
        "i1",
        "--serial",
        *extra,
    )


def _q(home: Path, name: str, sql: str) -> list[list]:
    return json.loads(
        _vc(home, "db", "query", name, sql, "--format", "JSONCompact").stdout
    )["data"]


def _records(path: Path) -> list[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return [line.rstrip("\n").split("\t") for line in f if not line.startswith("#")]


def _header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return [line.rstrip("\n") for line in f if line.startswith("#")]


def _export(home: Path, out: Path, *args: str, name: str = "src") -> Path:
    _vc(home, "db", "export", name, "-o", str(out), *args)
    return out


# --- round trip -----------------------------------------------------------------


def test_round_trip_is_lossless(vcfclick_home, tiny_vcf, tmp_path):
    _db(vcfclick_home, tiny_vcf, "src")
    out = _export(vcfclick_home, tmp_path / "out.vcf.gz")
    _db(vcfclick_home, out, "dst")
    vcols = "chrom, pos, ref, alt, qual, filter, info_AC, info_AF, info_AN, info_DP"
    order = "ORDER BY chrom, pos, ref, alt"
    assert _q(vcfclick_home, "src", f"SELECT {vcols} FROM variants {order}") == _q(
        vcfclick_home, "dst", f"SELECT {vcols} FROM variants {order}"
    )
    gcols = "chrom, pos, ref, alt, sample_id, gt, gq, dp, ad_ref, ad_alt"
    gorder = "ORDER BY chrom, pos, ref, alt, sample_id"
    src = _q(vcfclick_home, "src", f"SELECT {gcols} FROM genotypes {gorder}")
    assert src and src == _q(
        vcfclick_home, "dst", f"SELECT {gcols} FROM genotypes {gorder}"
    )


@pytest.mark.skipif(not HAS_HTSLIB, reason="bcftools/tabix not installed")
def test_output_is_valid_bgzf_vcf(vcfclick_home, tiny_vcf, tmp_path):
    _db(vcfclick_home, tiny_vcf, "src")
    out = _export(vcfclick_home, tmp_path / "out.vcf.gz")
    subprocess.run(["tabix", "-p", "vcf", str(out)], check=True)
    region = subprocess.run(
        ["bcftools", "view", "-H", "-r", "chr1:200-600", str(out)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    assert [line.split("\t")[1] for line in region] == ["250", "500"]
    lint = subprocess.run(
        ["bcftools", "view", str(out)], capture_output=True, text=True
    )
    assert (
        lint.returncode == 0 and "[W::" not in lint.stderr and "[E::" not in lint.stderr
    ), lint.stderr


def test_plain_vcf_and_stdout(vcfclick_home, tiny_vcf, tmp_path):
    _db(vcfclick_home, tiny_vcf, "src")
    plain = _export(vcfclick_home, tmp_path / "out.vcf")
    assert len(_records(plain)) == 5
    stdout = _vc(vcfclick_home, "db", "export", "src", "-o", "-").stdout
    assert stdout.startswith("##fileformat=VCFv4.2")
    assert len([ln for ln in stdout.splitlines() if not ln.startswith("#")]) == 5


# --- genotypes the sparse store can't distinguish ---------------------------------


def test_absent_genotypes_default_to_ref_with_a_header_note(
    vcfclick_home, tiny_vcf, tmp_path
):
    _db(vcfclick_home, tiny_vcf, "src")
    out = _export(vcfclick_home, tmp_path / "out.vcf")
    header = _header(out)
    assert any(h.startswith("##vcfclick_note=") and "0/0" in h for h in header)
    first = _records(out)[0]  # chr1:100, S1 was 0/0 in the source
    assert first[9].split(":")[0] == "0/0"
    nocall = _export(vcfclick_home, tmp_path / "nc.vcf", "--absent-as", "nocall")
    assert _records(nocall)[0][9].split(":")[0] == "./."


def test_keep_reference_ingest_keeps_no_calls(vcfclick_home, tmp_path):
    vcf = tmp_path / "nc.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n"
        "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t./.\t0/1\n"
    )
    _db(vcfclick_home, vcf, "src", "--keep-reference")
    rec = _records(_export(vcfclick_home, tmp_path / "out.vcf"))[0]
    assert [c.split(":")[0] for c in rec[9:]] == ["0/0", "./.", "0/1"]


# --- slicing ------------------------------------------------------------------------


def test_region_samples_sites_only_and_where(vcfclick_home, tiny_vcf, tmp_path):
    _db(vcfclick_home, tiny_vcf, "src")
    reg = _records(
        _export(vcfclick_home, tmp_path / "r.vcf", "--region", "chr1:200-600")
    )
    assert [r[1] for r in reg] == ["250", "500"]
    two = _records(
        _export(
            vcfclick_home,
            tmp_path / "r2.vcf",
            "--region",
            "chr1:1-150",
            "--region",
            "chr1:800-1000",
        )
    )
    assert [r[1] for r in two] == ["100", "900"]

    sub = _export(vcfclick_home, tmp_path / "s.vcf", "--samples", "S3,S2")
    assert _header(sub)[-1].split("\t")[9:] == ["S2", "S3"]
    assert all(len(r) == 11 for r in _records(sub))

    sites = _records(_export(vcfclick_home, tmp_path / "so.vcf", "--sites-only"))
    assert all(len(r) == 8 for r in sites) and len(sites) == 5
    assert _header(tmp_path / "so.vcf")[-1].split("\t")[-1] == "INFO"

    common = _records(
        _export(vcfclick_home, tmp_path / "w.vcf", "--where", "info_AF > 0.4")
    )
    assert [r[1] for r in common] == ["500", "750"]


def test_regions_bed(vcfclick_home, tiny_vcf, tmp_path):
    _db(vcfclick_home, tiny_vcf, "src")
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t199\t500\n")  # BED is 0-based half-open: covers 200..500
    assert [
        r[1]
        for r in _records(
            _export(vcfclick_home, tmp_path / "b.vcf", "--regions-bed", str(bed))
        )
    ] == ["250", "500"]


def test_gene(vcfclick_home, tiny_vcf, tmp_path, monkeypatch):
    ann = tmp_path / "ann.duckdb"
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(ann))
    import annotations.db as adb

    conn = adb.get_connection()
    conn.execute(
        "INSERT INTO refseq_genes VALUES ('GENEX', 'chr1', 200, 800, '+', '1', 'test')"
    )
    conn.close()
    _db(vcfclick_home, tiny_vcf, "src")
    assert [
        r[1]
        for r in _records(_export(vcfclick_home, tmp_path / "g.vcf", "--gene", "GENEX"))
    ] == ["250", "500", "750"]


def test_pass_only(vcfclick_home, tmp_path):
    vcf = tmp_path / "f.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
        '##FILTER=<ID=LowQual,Description="low">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"
        "chr1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "chr1\t200\t.\tC\tT\t5\tLowQual\t.\tGT\t0/1\n"
        "chr1\t300\t.\tG\tA\t40\t.\t.\tGT\t1/1\n"
    )
    _db(vcfclick_home, vcf, "src")
    allv = _export(vcfclick_home, tmp_path / "all.vcf")
    recs = _records(allv)
    # Ingest stores PASS and unset (.) alike (NULL), so both come back as "."
    # and the header says so; a real filter name survives.
    assert [(r[1], r[2], r[6]) for r in recs] == [
        ("100", "rs1", "."),
        ("200", ".", "LowQual"),
        ("300", ".", "."),
    ]
    header = _header(allv)
    assert any(h.startswith("##FILTER=<ID=LowQual") for h in header)
    assert any(h.startswith("##vcfclick_note_filter=") for h in header)
    assert [
        r[1]
        for r in _records(_export(vcfclick_home, tmp_path / "p.vcf", "--pass-only"))
    ] == ["100", "300"]


# --- errors ---------------------------------------------------------------------------


def test_multiple_ingestions_need_an_ingest_id(vcfclick_home, tiny_vcf, tmp_path):
    _db(vcfclick_home, tiny_vcf, "src")
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "src",
        str(tiny_vcf),
        "--cohort",
        "c2",
        "--ingest-id",
        "i2",
        "--serial",
    )
    r = _vc(
        vcfclick_home, "db", "export", "src", "-o", str(tmp_path / "x.vcf"), ok=False
    )
    assert (
        r.returncode != 0
        and "--ingest-id" in r.stderr
        and "i1" in r.stderr
        and "i2" in r.stderr
    )
    assert (
        len(_records(_export(vcfclick_home, tmp_path / "i2.vcf", "--ingest-id", "i2")))
        == 5
    )


def test_bad_inputs_are_reported(vcfclick_home, tiny_vcf, tmp_path):
    _db(vcfclick_home, tiny_vcf, "src")
    for args, msg in [
        (["--samples", "S1,NOPE"], "NOPE"),
        (["--region", "chr1:abc"], "region"),
        (["--where", "1=1; DROP TABLE variants"], "single"),
    ]:
        r = _vc(
            vcfclick_home,
            "db",
            "export",
            "src",
            "-o",
            str(tmp_path / "x.vcf"),
            *args,
            ok=False,
        )
        assert r.returncode != 0 and msg in r.stderr, (args, r.stderr)
