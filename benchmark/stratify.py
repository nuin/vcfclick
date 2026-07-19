"""Variant typing for `vcfclick benchmark`.

Assigns each normalized biallelic call a variant type (SNP/INDEL/COMPLEX/NOCALL),
a subtype (ti/tv for SNPs, signed size bin for indels), and a biallelic locus
(zygosity) type from its canonical genotype. Pure; operates on trimmed alleles.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Sequence

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

_NOCALL_ALTS = frozenset({".", "", "*"})
_TRANSITIONS = frozenset({frozenset("AG"), frozenset("CT")})


def vtype(ref: str, alt: str) -> str:
    """Classify a biallelic record by allele lengths."""
    if alt is None or alt in _NOCALL_ALTS:
        return VT_NOCALL
    if len(ref) == 1 and len(alt) == 1 and ref != alt:
        return VT_SNP
    if len(ref) != len(alt):
        return VT_INDEL
    return VT_COMPLEX


def _size_bin(n: int) -> str:
    """Absolute length-change bin label suffix."""
    if n <= 5:
        return "1_5"
    if n <= 15:
        return "6_15"
    return "16_PLUS"


def subtype(ref: str, alt: str) -> str:
    """Subtype within a vtype: ti/tv for SNPs, I/D size bin for indels, else ''."""
    vt = vtype(ref, alt)
    if vt == VT_SNP:
        return "ti" if frozenset((ref, alt)) in _TRANSITIONS else "tv"
    if vt == VT_INDEL:
        prefix = "I" if len(alt) > len(ref) else "D"
        return prefix + _size_bin(abs(len(alt) - len(ref)))
    return ""


def blt(gt: str) -> str:
    """Map a canonical genotype string to a biallelic locus (zygosity) type."""
    parts = gt.replace("|", "/").split("/")
    if any(p == "." for p in parts):
        return BLT_NOCALL
    alleles = [int(p) for p in parts]
    if len(alleles) == 1:
        return BLT_HAPLOID
    uniq = set(alleles)
    if uniq == {0}:
        return BLT_HOMREF
    if 0 in uniq:
        return BLT_HET
    if len(uniq) == 1:
        return BLT_HOMALT
    return BLT_HETALT


def tag(records: Sequence[NormRecord]) -> list[NormRecord]:
    """Return records with vtype/subtype/blt filled from their alleles and GT."""
    return [
        replace(
            r,
            vtype=vtype(r.ref, r.alt),
            subtype=subtype(r.ref, r.alt),
            blt=blt(r.gt),
        )
        for r in records
    ]
