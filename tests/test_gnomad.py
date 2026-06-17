"""Tests for the gnomAD allele-frequency annotation.

tests/fixtures/gnomad_cftr.vcf.gz holds REAL gnomAD v4.1 genomes AF /
AF_grpmax for the CFTR golden alleles (source URL in the header). The
loader/lookup are validated against those published values.
"""

from __future__ import annotations

from pathlib import Path

from annotations.db import GnomadAF

FIX = Path(__file__).parent / "fixtures"
GNOMAD_VCF = FIX / "gnomad_cftr.vcf.gz"


def test_gnomad_loader_and_lookup(tmp_path, monkeypatch):
    """Loading the real gnomAD slice and looking a variant up returns its
    published AF and popmax."""
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    from annotations import gnomad_af
    from annotations.loaders.gnomad import load

    assert load(GNOMAD_VCF, replace=True) == 5

    g = gnomad_af("chr7", 117488888, "A", "C")
    assert g is not None
    assert abs(g.af - 0.548299) < 1e-4
    assert abs(g.af_grpmax - 0.827770) < 1e-4
    assert g.popmax == g.af_grpmax  # popmax prefers the group-max AF

    rare = gnomad_af("chr7", 117480471, "A", "G")
    assert abs(rare.popmax - 0.009072) < 1e-4  # genuinely rare


def test_gnomad_lookup_absent_is_none(tmp_path, monkeypatch):
    """A locus outside the loaded slice returns None — the caller treats
    that as rare, not as AF 0."""
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    from annotations import gnomad_af
    from annotations.loaders.gnomad import load

    load(GNOMAD_VCF, replace=True)
    assert gnomad_af("chr7", 1, "A", "C") is None
    assert gnomad_af("chr7", 117488888, "A", "T") is None  # wrong alt


def test_popmax_falls_back_to_overall_af():
    """popmax uses AF_grpmax when present, else the overall AF."""
    assert GnomadAF("chr1", 1, "A", "G", 0.02, 0.05).popmax == 0.05
    assert GnomadAF("chr1", 1, "A", "G", 0.02, None).popmax == 0.02
    assert GnomadAF("chr1", 1, "A", "G", None, None).popmax is None
