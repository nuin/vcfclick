"""Unit tests for benchmark/haplotype.py — local-haplotype reconciliation core.

Uses the tiny.fa fixture (chr2 = ACGTACGT) so the canonical n:m fixture from the
P2 spec (truth MNP AC>GT vs query A>G + C>T) is exercised against a real
reference slice.
"""

from __future__ import annotations

from pathlib import Path

from benchmark import haplotype
from benchmark.model import NormRecord
from benchmark.reference import Reference

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "benchmark"
FETCH = Reference(str(FIXTURES / "tiny.fa")).fetch


def _rec(side, chrom, pos, ref, alt, gt="1/1"):
    return NormRecord(
        side=side,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        gt=gt,
        is_pass=True,
        other_alt_present=False,
        locus_id=(chrom, pos),
    )


# --- apply_haplotype ---------------------------------------------------------


def test_apply_haplotype_snp():
    # chr2 = ACGTACGT; window 1-based [1,2] = "AC"; apply A>G@1 -> "GC".
    r = _rec("query", "chr2", 1, "A", "G")
    assert haplotype.apply_haplotype(FETCH, "chr2", 1, 2, [r]) == "GC"


def test_apply_haplotype_mnp():
    r = _rec("truth", "chr2", 1, "AC", "GT")
    assert haplotype.apply_haplotype(FETCH, "chr2", 1, 2, [r]) == "GT"


def test_apply_haplotype_skips_homref():
    # A 0/0 record carries no alt allele and must not be applied.
    r = _rec("query", "chr2", 1, "A", "G", gt="0/0")
    assert haplotype.apply_haplotype(FETCH, "chr2", 1, 2, [r]) == "AC"


def test_apply_haplotype_right_to_left_indel():
    # A length-changing deletion left of a SNP: applying left-first would corrupt
    # the SNP offset. window chr2 [1,4] = "ACGT".
    # del GT>G @3 (removes T4), snp A>C ... use del at 2: CG>C removes G3.
    dele = _rec("truth", "chr2", 2, "CG", "C")  # ACGT -> A C T
    snp = _rec("truth", "chr2", 4, "T", "A")  # T4 -> A
    # Expected: apply right-to-left. T4->A gives ACGA; CG>C @2 gives A C A = "ACA".
    assert haplotype.apply_haplotype(FETCH, "chr2", 1, 4, [dele, snp]) == "ACA"


# --- cluster -----------------------------------------------------------------


def test_cluster_merges_within_flank():
    a = _rec("truth", "chr2", 1, "A", "G")
    b = _rec("truth", "chr2", 6, "C", "A")
    assert len(haplotype.cluster([a, b], flank=30)) == 1


def test_cluster_splits_beyond_flank():
    a = _rec("truth", "chr2", 1, "A", "G")
    b = _rec("truth", "chr2", 6, "C", "A")
    clusters = haplotype.cluster([a, b], flank=2)
    assert len(clusters) == 2


def test_cluster_isolates_chroms():
    a = _rec("truth", "chr1", 2, "A", "G")
    b = _rec("truth", "chr2", 2, "C", "T")
    clusters = haplotype.cluster([a, b], flank=1000)
    assert len(clusters) == 2


# --- haplotype_equivalent ----------------------------------------------------


def test_haplotype_equivalent_canonical_mnp_vs_two_snps():
    # The canonical n:m fixture: truth AC>GT (1/1) == query A>G + C>T (both 1/1).
    truth = [_rec("truth", "chr2", 1, "AC", "GT")]
    query = [
        _rec("query", "chr2", 1, "A", "G"),
        _rec("query", "chr2", 2, "C", "T"),
    ]
    assert haplotype.haplotype_equivalent(FETCH, truth, query, flank=30) is True


def test_haplotype_equivalent_non_equivalent():
    # Query only spells the first SNP; sequences differ (GT vs GC).
    truth = [_rec("truth", "chr2", 1, "AC", "GT")]
    query = [_rec("query", "chr2", 1, "A", "G")]
    assert haplotype.haplotype_equivalent(FETCH, truth, query, flank=30) is False


def test_haplotype_equivalent_het_phase_mismatch():
    # Truth hom AC>GT vs two het SNPs: no phasing reproduces GT/GT -> no match.
    truth = [_rec("truth", "chr2", 1, "AC", "GT", gt="1/1")]
    query = [
        _rec("query", "chr2", 1, "A", "G", gt="0/1"),
        _rec("query", "chr2", 2, "C", "T", gt="0/1"),
    ]
    assert haplotype.haplotype_equivalent(FETCH, truth, query, flank=30) is False


def test_haplotype_equivalent_het_cis_match():
    # Both sides het, both alts in cis (0/1 sharing one haplotype) reproduce the
    # same diplotype as a single het MNP.
    truth = [_rec("truth", "chr2", 1, "AC", "GT", gt="0/1")]
    query = [
        _rec("query", "chr2", 1, "A", "G", gt="0/1"),
        _rec("query", "chr2", 2, "C", "T", gt="0/1"),
    ]
    assert haplotype.haplotype_equivalent(FETCH, truth, query, flank=30) is True


def test_haplotype_equivalent_representation_shifted_deletion():
    # chr1 = CAAAAT: deleting one A from the run is representation-ambiguous.
    # truth left-aligns it as CA>C @1; query spells the same net deletion as
    # AA>A @4. Both yield CAAAT over the shared window, so they are equivalent.
    truth = [_rec("truth", "chr1", 1, "CA", "C", gt="1/1")]
    query = [_rec("query", "chr1", 4, "AA", "A", gt="1/1")]
    assert haplotype.haplotype_equivalent(FETCH, truth, query, flank=30) is True


def test_apply_haplotype_none_on_overlapping_records():
    # Overlapping same-side records cannot be co-applied -> infeasible (None), so
    # a pathological pair can never coincidentally match the other side (a false TP).
    r1 = _rec("truth", "chr1", 2, "AA", "G")
    r2 = _rec("truth", "chr1", 2, "A", "T")
    assert haplotype.apply_haplotype(FETCH, "chr1", 2, 3, [r1, r2]) is None


def test_overlapping_records_rejected_even_when_splice_looks_valid():
    # ref chr1[0:2] = "CA". pos1 CA>T and pos2 A>AA have OVERLAPPING footprints
    # ([0,2) and [1,2)); co-applying them is infeasible and must yield None, so a
    # pathological overlap can never coincidentally match another side (false TP).
    r1 = _rec("truth", "chr1", 1, "CA", "T")
    r2 = _rec("truth", "chr1", 2, "A", "AA")
    assert haplotype.apply_haplotype(FETCH, "chr1", 1, 2, [r1, r2]) is None


def test_phased_trans_hets_not_matched_to_cis_mnp():
    # Truth is two phased hets in TRANS: A>G 0|1 and (next base) C>T 1|0 — the two
    # alts are on OPPOSITE haplotypes, so neither haplotype carries both. A cis het
    # MNP that puts both alts on one haplotype must NOT be judged equivalent.
    # chr2 = ACGTACGT: pos1 A>G, pos2 C>T.
    truth = [
        _rec("truth", "chr2", 1, "A", "G", gt="0|1"),
        _rec("truth", "chr2", 2, "C", "T", gt="1|0"),
    ]
    query = [_rec("query", "chr2", 1, "AC", "GT", gt="0/1")]  # cis: both on one hap
    assert haplotype.haplotype_equivalent(FETCH, truth, query) is False
