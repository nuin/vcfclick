"""Tests for `vcfclick db ingest-batch`.

Two fixtures stand in for per-sample VCFs:
  tests/fixtures/per_sample_a.vcf.gz  — sample SAMPLE_A, 3 variants
  tests/fixtures/per_sample_b.vcf.gz  — sample SAMPLE_B, 3 variants

For continue-on-error the existing multi-allelic fixture is reused:
  tests/fixtures/multiallelic.vcf.gz  — fails on the in-stream check
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
FIXTURES = Path(__file__).parent / "fixtures"
PER_SAMPLE_A = FIXTURES / "per_sample_a.vcf.gz"
PER_SAMPLE_B = FIXTURES / "per_sample_b.vcf.gz"
MULTI = FIXTURES / "multiallelic.vcf.gz"


def _vc(
    home: Path, *args: str, expect_failure: bool = False
) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    if expect_failure:
        assert r.returncode != 0, f"expected failure but got rc=0:\n{r.stdout}"
    else:
        assert r.returncode == 0, (
            f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
            f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
        )
    return r


def _query(home: Path, db: str, sql: str) -> str:
    """TabSeparated output so assertions check the value directly, not
    the engine's box-drawing render."""
    if "FORMAT " not in sql.upper():
        sql = sql + " FORMAT TabSeparated"
    return _vc(home, "db", "query", db, sql).stdout


# ─── derive_ingest_id helper ─────────────────────────────────────────


def test_derive_ingest_id_strips_vcf_extensions():
    """Filename-derivation: HG00096.vcf.gz → HG00096."""
    from cli.db import _derive_ingest_id

    assert _derive_ingest_id(Path("/x/HG00096.vcf.gz")) == "HG00096"
    assert _derive_ingest_id(Path("/x/HG00096.vcf.bgz")) == "HG00096"
    assert _derive_ingest_id(Path("/x/HG00096.vcf")) == "HG00096"
    # Files without a known VCF extension fall back to .stem
    assert _derive_ingest_id(Path("/x/raw_file.txt")) == "raw_file"


# ─── --from-dir happy path ───────────────────────────────────────────


def test_from_dir_ingests_all_vcfs_in_directory(vcfclick_home, tmp_path):
    """Two VCFs in a directory → two ingestions under one cohort,
    each with ingest_id derived from the filename stem."""
    staging = tmp_path / "per_sample"
    staging.mkdir()
    shutil.copy(PER_SAMPLE_A, staging)
    shutil.copy(PER_SAMPLE_B, staging)
    shutil.copy(PER_SAMPLE_A.with_suffix(".gz.tbi"), staging)
    shutil.copy(PER_SAMPLE_B.with_suffix(".gz.tbi"), staging)

    _vc(vcfclick_home, "db", "create", "demo")
    _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--from-dir",
        str(staging),
        "--cohort",
        "study1",
    )

    # 2 ingestions in the catalog
    out = _query(vcfclick_home, "demo", "SELECT count() FROM ingestions")
    assert out.strip() == "2"

    # ingest_ids derived from filenames
    out = _query(
        vcfclick_home,
        "demo",
        "SELECT ingest_id FROM ingestions ORDER BY ingest_id FORMAT TSV",
    )
    assert out.strip().splitlines() == ["per_sample_a", "per_sample_b"]

    # Both cohort labels are "study1"
    out = _query(
        vcfclick_home,
        "demo",
        "SELECT DISTINCT cohort FROM samples FORMAT TSV",
    )
    assert out.strip() == "study1"


# ─── --manifest happy path ───────────────────────────────────────────


def test_manifest_with_custom_ingest_ids_and_per_row_cohort(vcfclick_home, tmp_path):
    """Manifest carries its own sample_id and cohort columns; CLI
    --cohort is only a fallback (not used when row supplies its own)."""
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample_id\tvcf_path\tcohort\n"
        f"case_001\t{PER_SAMPLE_A}\tcase\n"
        f"ctrl_001\t{PER_SAMPLE_B}\tcontrol\n"
    )

    _vc(vcfclick_home, "db", "create", "demo")
    _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--manifest",
        str(manifest),
        "--cohort",
        "fallback",
    )

    out = _query(
        vcfclick_home,
        "demo",
        "SELECT ingest_id, cohort FROM samples ORDER BY ingest_id FORMAT TSV",
    )
    rows = [line.split("\t") for line in out.strip().splitlines()]
    assert rows == [["case_001", "case"], ["ctrl_001", "control"]]


def test_manifest_falls_back_to_cli_cohort_when_column_absent(vcfclick_home, tmp_path):
    """A manifest with no `cohort` column uses the --cohort flag."""
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(f"vcf_path\n{PER_SAMPLE_A}\n{PER_SAMPLE_B}\n")

    _vc(vcfclick_home, "db", "create", "demo")
    _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--manifest",
        str(manifest),
        "--cohort",
        "shared",
    )

    out = _query(
        vcfclick_home,
        "demo",
        "SELECT DISTINCT cohort FROM samples FORMAT TSV",
    )
    assert out.strip() == "shared"


def test_manifest_paths_resolved_relative_to_manifest_dir(vcfclick_home, tmp_path):
    """nf-core convention: relative paths in the manifest resolve
    against the manifest file's directory."""
    staging = tmp_path / "vcfs"
    staging.mkdir()
    shutil.copy(PER_SAMPLE_A, staging)
    shutil.copy(PER_SAMPLE_A.with_suffix(".gz.tbi"), staging)

    manifest = tmp_path / "samples.tsv"
    manifest.write_text("vcf_path\nvcfs/per_sample_a.vcf.gz\n")

    _vc(vcfclick_home, "db", "create", "demo")
    _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--manifest",
        str(manifest),
        "--cohort",
        "x",
    )

    out = _query(vcfclick_home, "demo", "SELECT count() FROM ingestions")
    assert out.strip() == "1"


# ─── continue-on-error ───────────────────────────────────────────────


def test_continue_on_error_skips_bad_files_keeps_good(vcfclick_home, tmp_path):
    """A failing file (multi-allelic) leaves its writes rolled back
    and the batch continues. Summary lists the failure, exit code
    is non-zero."""
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample_id\tvcf_path\n"
        f"good_a\t{PER_SAMPLE_A}\n"
        f"bad_one\t{MULTI}\n"
        f"good_b\t{PER_SAMPLE_B}\n"
    )

    _vc(vcfclick_home, "db", "create", "demo")
    r = _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--manifest",
        str(manifest),
        "--cohort",
        "test",
        expect_failure=True,  # any failure → rc != 0
    )

    # Summary numbers in stdout
    assert "total:    3" in r.stdout
    assert "ingested: 2" in r.stdout
    assert "failed:   1" in r.stdout
    assert "bad_one" in r.stdout
    assert "Multi-allelic" in r.stdout  # the underlying error

    # The two good files made it in; the bad one left no rows
    out = _query(
        vcfclick_home,
        "demo",
        "SELECT ingest_id FROM ingestions ORDER BY ingest_id FORMAT TSV",
    )
    assert out.strip().splitlines() == ["good_a", "good_b"]

    out = _query(
        vcfclick_home,
        "demo",
        "SELECT count() FROM variants WHERE ingest_id = 'bad_one'",
    )
    assert out.strip() == "0"


# ─── argv validation ─────────────────────────────────────────────────


def test_from_dir_and_manifest_are_mutually_exclusive(vcfclick_home, tmp_path):
    _vc(vcfclick_home, "db", "create", "demo")
    manifest = tmp_path / "m.tsv"
    manifest.write_text("vcf_path\n")
    r = _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--from-dir",
        str(tmp_path),
        "--manifest",
        str(manifest),
        "--cohort",
        "x",
        expect_failure=True,
    )
    assert "mutually exclusive" in r.stderr


def test_neither_from_dir_nor_manifest_errors(vcfclick_home):
    _vc(vcfclick_home, "db", "create", "demo")
    r = _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--cohort",
        "x",
        expect_failure=True,
    )
    assert "--from-dir or --manifest" in r.stderr


def test_from_dir_requires_cohort(vcfclick_home, tmp_path):
    _vc(vcfclick_home, "db", "create", "demo")
    staging = tmp_path / "v"
    staging.mkdir()
    shutil.copy(PER_SAMPLE_A, staging)
    r = _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--from-dir",
        str(staging),
        expect_failure=True,
    )
    assert "--cohort is required" in r.stderr


def test_manifest_missing_vcf_path_column_errors(vcfclick_home, tmp_path):
    _vc(vcfclick_home, "db", "create", "demo")
    manifest = tmp_path / "m.tsv"
    manifest.write_text("sample_id\tcohort\nA\tcase\n")  # no vcf_path column
    r = _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--manifest",
        str(manifest),
        "--cohort",
        "x",
        expect_failure=True,
    )
    assert "vcf_path" in r.stderr


def test_manifest_rejects_duplicate_ingest_ids(vcfclick_home, tmp_path):
    _vc(vcfclick_home, "db", "create", "demo")
    manifest = tmp_path / "m.tsv"
    manifest.write_text(
        f"sample_id\tvcf_path\ndup\t{PER_SAMPLE_A}\ndup\t{PER_SAMPLE_B}\n"
    )
    r = _vc(
        vcfclick_home,
        "db",
        "ingest-batch",
        "demo",
        "--manifest",
        str(manifest),
        "--cohort",
        "x",
        expect_failure=True,
    )
    assert "duplicate ingest_id" in r.stderr
