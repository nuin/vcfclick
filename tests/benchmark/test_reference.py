from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pyfaidx")

from benchmark.reference import (
    ContigError,
    Reference,
    build_alias_map,
    canonical_contig,
)

FASTA = Path(__file__).resolve().parent.parent / "fixtures" / "benchmark" / "tiny.fa"


def test_fetch_is_zero_based_half_open_uppercase():
    ref = Reference(FASTA)
    assert ref.fetch("chr1", 0, 1) == "C"
    assert ref.fetch("chr1", 1, 5) == "AAAA"
    assert ref.fetch("chr1", 0, 6) == "CAAAAT"


def test_contig_len():
    ref = Reference(FASTA)
    assert ref.contig_len("chr1") == 6
    assert ref.contig_len("chr2") == 8


def test_build_alias_map_chr_and_mt():
    amap = build_alias_map(["chr1", "chr2", "chrM"])
    assert amap["1"] == "chr1"
    assert amap["chr1"] == "chr1"
    assert amap["MT"] == "chrM"
    assert amap["M"] == "chrM"
    assert amap["chrM"] == "chrM"


def test_canonical_contig_resolves_alias():
    known = {"chr1", "chrM"}
    assert canonical_contig("1", known) == "chr1"
    assert canonical_contig("MT", known) == "chrM"
    assert canonical_contig("chr1", known) == "chr1"
    assert canonical_contig("chr9", known) is None


def test_fetch_via_alias():
    ref = Reference(FASTA)
    # "1" resolves to "chr1"; "MT" resolves to "chrM"
    assert ref.fetch("1", 0, 1) == "C"
    assert ref.fetch("MT", 0, 4) == "GGGG"


def test_validate_length_mismatch_raises():
    ref = Reference(FASTA)
    ref.validate_length("chr1", 6)  # ok, no raise
    with pytest.raises(ContigError):
        ref.validate_length("chr1", 999)
