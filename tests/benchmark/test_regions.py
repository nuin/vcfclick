from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("numpy")

from benchmark.model import NormRecord
from benchmark.regions import ConfRegions

BED = Path(__file__).resolve().parent.parent / "fixtures" / "benchmark" / "conf.bed"


def _oracle(intervals, chrom, pos1):
    """Brute-force start-containment: start0 <= pos1-1 < end0."""
    pos0 = pos1 - 1
    return any(c == chrom and s <= pos0 < e for c, s, e in intervals)


# ---- contains: membership vs brute-force oracle ----


def test_contains_matches_oracle():
    intervals = [("chr1", 10, 20), ("chr1", 25, 30), ("chr2", 0, 5)]
    cr = ConfRegions(intervals)
    for chrom in ("chr1", "chr2", "chrX"):
        for pos1 in range(0, 35):
            assert cr.contains(chrom, pos1) == _oracle(intervals, chrom, pos1), (
                chrom,
                pos1,
            )


# ---- boundary: start inclusive, end exclusive (0-based interval, 1-based POS) ----


def test_boundaries_start_inclusive_end_exclusive():
    # interval [10, 20) in 0-based → covers 0-based 10..19 → 1-based POS 11..20
    cr = ConfRegions([("chr1", 10, 20)])
    assert cr.contains("chr1", 10) is False  # pos0 = 9  < 10
    assert cr.contains("chr1", 11) is True  # pos0 = 10 == start
    assert cr.contains("chr1", 20) is True  # pos0 = 19 last inside
    assert cr.contains("chr1", 21) is False  # pos0 = 20 == end (exclusive)


# ---- merge overlapping AND adjacent intervals ----


def test_merge_overlapping():
    # [0,5) and [3,8) overlap → [0,8): 1-based POS 1..8 all inside, 9 outside
    cr = ConfRegions([("chr1", 0, 5), ("chr1", 3, 8)])
    assert all(cr.contains("chr1", p) for p in range(1, 9))
    assert cr.contains("chr1", 9) is False


def test_merge_adjacent():
    # [0,5) and [5,10) touch → [0,10); no gap at the seam (POS 5 and 6 both in)
    cr = ConfRegions([("chr1", 0, 5), ("chr1", 5, 10)])
    assert cr.contains("chr1", 5) is True
    assert cr.contains("chr1", 6) is True
    assert cr.contains("chr1", 10) is True
    assert cr.contains("chr1", 11) is False


def test_non_adjacent_keeps_gap():
    # [0,5) and [6,10) leave a 1-base gap at 0-based 5 → 1-based POS 6 outside
    cr = ConfRegions([("chr1", 0, 5), ("chr1", 6, 10)])
    assert cr.contains("chr1", 5) is True  # pos0 4 inside first
    assert cr.contains("chr1", 6) is False  # pos0 5 in the gap
    assert cr.contains("chr1", 7) is True  # pos0 6 inside second


def test_unordered_input_is_sorted():
    cr = ConfRegions([("chr1", 25, 30), ("chr1", 0, 5), ("chr1", 10, 15)])
    assert cr.contains("chr1", 3) is True
    assert cr.contains("chr1", 12) is True
    assert cr.contains("chr1", 27) is True
    assert cr.contains("chr1", 8) is False


# ---- empty / unknown contig ----


def test_empty_regions_contains_nothing():
    cr = ConfRegions([])
    assert cr.contains("chr1", 1) is False


def test_unknown_contig_not_contained():
    cr = ConfRegions([("chr1", 0, 10)])
    assert cr.contains("chrZ", 5) is False


# ---- from_bed: parse the fixture BED ----


def test_from_bed_parses_and_merges():
    cr = ConfRegions.from_bed(BED)
    # chr1 rows [1,3)+[3,5) adjacent → merged [1,5): 1-based POS 2..5 inside
    assert cr.contains("chr1", 1) is False  # pos0 0 < 1
    assert cr.contains("chr1", 2) is True
    assert cr.contains("chr1", 5) is True
    assert cr.contains("chr1", 6) is False
    # chr2 [0,4): 1-based POS 1..4 inside
    assert cr.contains("chr2", 1) is True
    assert cr.contains("chr2", 4) is True
    assert cr.contains("chr2", 5) is False


# ---- tag: set in_conf on the normalized locus start ----


def _rec(chrom, pos, in_conf=False):
    return NormRecord(
        side="query",
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="C",
        gt="0/1",
        is_pass=True,
        other_alt_present=False,
        locus_id=(pos, "A", ("C",)),
        in_conf=in_conf,
    )


def test_tag_sets_in_conf_and_copies():
    cr = ConfRegions([("chr1", 10, 20)])
    recs = [_rec("chr1", 15), _rec("chr1", 5), _rec("chrZ", 15)]
    out = cr.tag(recs)
    assert [r.in_conf for r in out] == [True, False, False]
    # returns copies; originals untouched
    assert [r.in_conf for r in recs] == [False, False, False]
    assert out[0] is not recs[0]
    # only in_conf changed
    assert out[0].pos == 15 and out[0].chrom == "chr1"


def test_full_containment_rejects_boundary_straddling_indel():
    from benchmark.model import NormRecord

    # region covers 0-based [2,5) => 1-based POS 3,4,5
    rec = NormRecord(
        side="query",
        chrom="chr1",
        pos=5,
        ref="AA",
        alt="A",
        gt="0/1",
        is_pass=True,
        other_alt_present=False,
        locus_id=(5, "AA", ("A",)),
    )
    start = ConfRegions([("chr1", 2, 5)]).tag([rec], containment="start")[0]
    full = ConfRegions([("chr1", 2, 5)]).tag([rec], containment="full")[0]
    assert start.in_conf is True  # POS 5 start is inside
    assert full.in_conf is False  # span 0-based [4,6) exceeds region end 5
