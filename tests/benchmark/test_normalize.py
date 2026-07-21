from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pyfaidx")

from benchmark.normalize import (
    canonical_gt,
    decompose_mnp,
    left_align,
    split_multiallelic,
    trim,
)
from benchmark.reference import Reference

FASTA = Path(__file__).resolve().parent.parent / "fixtures" / "benchmark" / "tiny.fa"


# ---- trim: VCF minimal representation (no reference movement) ----


def test_trim_snp_passthrough():
    assert trim(5, "A", "G") == (5, "A", "G")


def test_trim_right_shared_suffix():
    assert trim(100, "CTT", "CT") == (100, "CT", "C")


def test_trim_both_sides_adjusts_pos():
    assert trim(100, "GCAT", "GCT") == (101, "CA", "C")


def test_trim_keeps_one_base_each():
    # deletion already minimal; must not empty an allele
    assert trim(4, "AA", "A") == (4, "AA", "A")


# ---- left_align: shift indels leftward through a repeat, vs the reference ----
# tiny.fa chr1 = CAAAAT  (1-based: C1 A2 A3 A4 A5 T6)


def _fetch():
    return Reference(FASTA).fetch


def test_left_align_deletion_rolls_to_leftmost_anchor():
    # delete one A of the run, written at the right → leftmost representation
    assert left_align(_fetch(), "chr1", 4, "AA", "A") == (1, "CA", "C")


def test_left_align_insertion_rolls_to_leftmost_anchor():
    assert left_align(_fetch(), "chr1", 5, "A", "AA") == (1, "C", "CA")


def test_left_align_snp_is_unchanged():
    assert left_align(_fetch(), "chr1", 2, "A", "G") == (2, "A", "G")


def test_left_align_is_idempotent():
    once = left_align(_fetch(), "chr1", 4, "AA", "A")
    assert left_align(_fetch(), "chr1", *once) == once


# ---- split_multiallelic: -m - style split with GT remap ----


def test_split_biallelic_passthrough():
    rows = split_multiallelic(100, "A", ["C"], (0, 1))
    assert len(rows) == 1
    assert rows[0].alt == "C"
    assert rows[0].gt == (0, 1)
    assert rows[0].other_alt_present is False


def test_split_het_alt_remaps_gt_and_flags_other_alt():
    rows = split_multiallelic(100, "A", ["C", "G"], (1, 2))
    assert [r.alt for r in rows] == ["C", "G"]
    assert rows[0].gt == (1, 0) and rows[0].other_alt_present is True
    assert rows[1].gt == (0, 1) and rows[1].other_alt_present is True
    # sibling alleles share a locus id
    assert rows[0].locus_id == rows[1].locus_id


# ---- canonical_gt: order-independent unphased, missing, phased ----


def test_canonical_gt_sorts_unphased():
    assert canonical_gt((1, 0)) == "0/1"
    assert canonical_gt((2, 1)) == "1/2"
    assert canonical_gt((1, 1)) == "1/1"


def test_canonical_gt_missing_allele():
    assert canonical_gt((1, -1)) == "./1"


def test_canonical_gt_phased_keeps_order():
    assert canonical_gt((1, 0), phased=True) == "1|0"


# ---- decompose_mnp: split equal-length multibase substitutions into SNPs ----


def test_decompose_mnp_splits_into_snps():
    assert decompose_mnp(10, "AT", "GC") == [(10, "A", "G"), (11, "T", "C")]


def test_decompose_mnp_skips_matching_bases():
    assert decompose_mnp(10, "AC", "GC") == [(10, "A", "G")]


def test_decompose_mnp_passthrough_non_mnp():
    assert decompose_mnp(10, "A", "G") == [(10, "A", "G")]  # SNP
    assert decompose_mnp(10, "AT", "A") == [(10, "AT", "A")]  # indel, not MNP


# ---- left_align never returns an empty allele ----


def test_left_align_rejects_empty_allele():
    with pytest.raises(ValueError):
        left_align(_fetch(), "chr1", 1, "", "A")
