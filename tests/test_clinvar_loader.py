"""Tests for the ClinVar VCF loader.

Exercises the library-level loader API against a small ClinVar-shaped
fixture (5 valid alleles spread across 4 records: one bi-allelic, one
multi-allelic with two ALTs, one MT-contig, and one record with ALT
'.' that must be skipped).

The DuckDB store is monkey-patched per-test so we don't touch the
user's real annotations/annotations.duckdb.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
CLINVAR_FIXTURE = Path(__file__).parent / "fixtures" / "clinvar_mini.vcf.gz"


@pytest.fixture
def isolated_annotation_db(tmp_path, monkeypatch):
    """Redirect annotations.db.DUCKDB_PATH to a temp file per test."""
    import annotations.db as adb

    monkeypatch.setattr(adb, "DUCKDB_PATH", tmp_path / "test_annotations.duckdb")
    yield tmp_path / "test_annotations.duckdb"


def test_load_counts_decomposed_alleles(isolated_annotation_db):
    """Fixture has 4 records, one of which is multi-allelic (2 ALTs) and
    one of which has ALT='.' (skip). Expected rows = 3 + 2 - 1 = 5."""
    from annotations.loaders.clinvar import load

    n = load(CLINVAR_FIXTURE)
    assert n == 5


def test_contigs_normalised_to_chr_prefix(isolated_annotation_db):
    """ClinVar's bare numeric '17' → 'chr17', 'MT' → 'chrM'."""
    from annotations.db import get_connection
    from annotations.loaders.clinvar import load

    load(CLINVAR_FIXTURE)
    rows = get_connection().execute(
        "SELECT DISTINCT chrom FROM clinvar_variants ORDER BY chrom"
    ).fetchall()
    contigs = {r[0] for r in rows}
    assert contigs == {"chr17", "chrM"}


def test_lookup_returns_pathogenic(isolated_annotation_db):
    from annotations.db import clinvar_lookup
    from annotations.loaders.clinvar import load

    load(CLINVAR_FIXTURE)
    r = clinvar_lookup("chr17", 43044295, "C", "T")
    assert r is not None
    assert r.clin_sig == "Pathogenic"
    assert "Hereditary" in r.condition
    assert r.clinvar_id == "1001"


def test_multiallelic_record_decomposes_to_separate_rows(isolated_annotation_db):
    """Pos 43044300 G→A,C must produce two distinct lookups, each
    returning the same shared CLN values."""
    from annotations.db import clinvar_lookup
    from annotations.loaders.clinvar import load

    load(CLINVAR_FIXTURE)
    ga = clinvar_lookup("chr17", 43044300, "G", "A")
    gc = clinvar_lookup("chr17", 43044300, "G", "C")
    assert ga is not None and gc is not None
    assert ga.alt == "A"
    assert gc.alt == "C"
    # Pipe-joined multi-value CLNSIG preserved as-is.
    assert "Likely_pathogenic" in ga.clin_sig
    assert "Uncertain_significance" in ga.clin_sig
    assert ga.clinvar_id == gc.clinvar_id  # same VCV accession


def test_alt_dot_record_is_skipped(isolated_annotation_db):
    """A record with ALT='.' cannot compose into the (chrom,pos,ref,alt)
    primary key — must be silently dropped, not fabricated."""
    from annotations.db import clinvar_lookup, get_connection
    from annotations.loaders.clinvar import load

    load(CLINVAR_FIXTURE)
    # The pos 43044500 record had ALT='.' — no row in either form should exist.
    assert clinvar_lookup("chr17", 43044500, "A", ".") is None
    assert clinvar_lookup("chr17", 43044500, "A", "") is None
    rows = get_connection().execute(
        "SELECT count() FROM clinvar_variants WHERE pos = 43044500"
    ).fetchone()
    assert rows[0] == 0


def test_mt_contig_lookup(isolated_annotation_db):
    """ClinVar's 'MT' contig must be looked up via 'chrM' (not 'chrMT')."""
    from annotations.db import clinvar_lookup
    from annotations.loaders.clinvar import load

    load(CLINVAR_FIXTURE)
    assert clinvar_lookup("chrM", 8993, "T", "G") is not None
    assert clinvar_lookup("chrMT", 8993, "T", "G") is None


def test_replace_mode_clears_prior_rows(isolated_annotation_db):
    """Default replace=True wipes the table before re-loading."""
    from annotations.db import get_connection
    from annotations.loaders.clinvar import load

    load(CLINVAR_FIXTURE)
    load(CLINVAR_FIXTURE)  # second load with replace=True (default)
    rows = get_connection().execute("SELECT count() FROM clinvar_variants").fetchone()
    assert rows[0] == 5  # still 5, not 10


def test_cli_command_registered():
    """`vcfclick annotations load-clinvar --help` must list the new command."""
    r = subprocess.run(
        [VCFCLICK_BIN, "annotations", "load-clinvar", "--help"],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, f"--help failed:\n{r.stdout}\n{r.stderr}"
    assert "--vcf" in r.stdout
    assert "--keep-existing" in r.stdout
    assert "ClinVar" in r.stdout
