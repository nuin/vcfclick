from __future__ import annotations

from benchmark.constants import (
    BLT_HAPLOID,
    BLT_HET,
    BLT_HETALT,
    BLT_HOMALT,
    BLT_HOMREF,
    BLT_NOCALL,
    VT_COMPLEX,
    VT_INDEL,
    VT_NOCALL,
    VT_SNP,
)
from benchmark.model import NormRecord
from benchmark.stratify import blt, subtype, tag, vtype

# ---- vtype ----


def test_vtype_snp():
    assert vtype("A", "G") == VT_SNP


def test_vtype_indel_insertion():
    assert vtype("A", "AT") == VT_INDEL


def test_vtype_indel_deletion():
    assert vtype("AT", "A") == VT_INDEL


def test_vtype_complex_equal_len_multibase():
    assert vtype("AT", "GC") == VT_COMPLEX


def test_vtype_nocall_variants():
    assert vtype("A", ".") == VT_NOCALL
    assert vtype("A", "") == VT_NOCALL
    assert vtype("A", "*") == VT_NOCALL
    assert vtype("A", None) == VT_NOCALL


# ---- subtype: SNP transition vs transversion ----


def test_subtype_ti_all_transition_pairs():
    assert subtype("A", "G") == "ti"
    assert subtype("G", "A") == "ti"
    assert subtype("C", "T") == "ti"
    assert subtype("T", "C") == "ti"


def test_subtype_tv_transversion():
    assert subtype("A", "C") == "tv"


# ---- subtype: indel size bins, insertion and deletion ----


def test_subtype_insertion_bins():
    assert subtype("A", "AT") == "I1_5"  # +1
    assert subtype("A", "A" + "T" * 7) == "I6_15"  # +7
    assert subtype("A", "A" + "T" * 20) == "I16_PLUS"  # +20


def test_subtype_deletion_bins():
    assert subtype("AT", "A") == "D1_5"  # -1
    assert subtype("A" + "T" * 7, "A") == "D6_15"  # -7
    assert subtype("A" + "T" * 20, "A") == "D16_PLUS"  # -20


def test_subtype_bin_boundaries():
    assert subtype("A", "A" + "T" * 5) == "I1_5"  # +5 upper edge
    assert subtype("A", "A" + "T" * 6) == "I6_15"  # +6 lower edge
    assert subtype("A", "A" + "T" * 15) == "I6_15"  # +15 upper edge
    assert subtype("A", "A" + "T" * 16) == "I16_PLUS"  # +16 lower edge


def test_subtype_complex_and_nocall_empty():
    assert subtype("AT", "GC") == ""
    assert subtype("A", ".") == ""


# ---- blt: biallelic locus (zygosity) type ----


def test_blt_homalt():
    assert blt("1/1") == BLT_HOMALT


def test_blt_het():
    assert blt("0/1") == BLT_HET


def test_blt_hetalt():
    assert blt("1/2") == BLT_HETALT


def test_blt_homref():
    assert blt("0/0") == BLT_HOMREF


def test_blt_haploid():
    assert blt("1") == BLT_HAPLOID
    assert blt("0") == BLT_HAPLOID


def test_blt_nocall():
    assert blt(".") == BLT_NOCALL
    assert blt("./.") == BLT_NOCALL


def test_blt_phased_separator():
    assert blt("0|1") == BLT_HET
    assert blt("1|1") == BLT_HOMALT


# ---- tag: fill vtype/subtype/blt on records ----


def _rec(ref, alt, gt):
    return NormRecord(
        side="query",
        chrom="chr1",
        pos=1,
        ref=ref,
        alt=alt,
        gt=gt,
        is_pass=True,
        other_alt_present=False,
        locus_id=(1, ref, (alt,)),
    )


def test_tag_fills_all_fields():
    recs = [_rec("A", "G", "0/1"), _rec("A", "AT", "1/1")]
    out = tag(recs)
    assert out[0].vtype == VT_SNP
    assert out[0].subtype == "ti"
    assert out[0].blt == BLT_HET
    assert out[1].vtype == VT_INDEL
    assert out[1].subtype == "I1_5"
    assert out[1].blt == BLT_HOMALT


def test_tag_is_pure():
    recs = [_rec("A", "G", "0/1")]
    tag(recs)
    assert recs[0].vtype == ""  # original unchanged
