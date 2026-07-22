"""Shared record schema for `vcfclick benchmark` — the contract every module
builds against. Pure data; no behavior.

Flow:  cyvcf2 → NormRecord (per normalized biallelic call, one per side) →
regions/stratify tag `in_conf`/`vtype`/`subtype`/`blt` → reconcile emits
ClassifiedRow (once per filter view) → aggregate counts → metrics → report.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class NormRecord:
    """A normalized biallelic call on one side, before classification.

    `gt` is the canonical genotype string (see normalize.canonical_gt).
    `locus_id` links sibling alleles from one original multiallelic site so a
    genotype verdict can be taken on the reconstructed per-locus allele multiset.
    Region/type tags default empty and are filled by regions.py / stratify.py.
    """

    side: str  # "truth" | "query"
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: str
    is_pass: bool
    other_alt_present: bool
    locus_id: tuple
    in_conf: bool = False
    vtype: str = ""  # constants.VT_*
    subtype: str = ""  # ti/tv | indel size bin | ""
    blt: str = ""  # constants.BLT_*
    qual: float = 0.0  # query QUAL (for ROC threshold sweeps)


@dataclass(frozen=True)
class ClassifiedRow:
    """One classified variant for one filter view. Aggregation groups over
    (filter_view, vtype, subtype, side, bd)."""

    side: str  # "truth" | "query"
    filter_view: str  # constants.FILTER_ALL | FILTER_PASS
    chrom: str
    pos: int
    ref: str
    alt: str
    vtype: str
    subtype: str
    in_conf: bool
    bd: str  # constants.BD_*
    bk: str  # constants.BK_*
    blt: str
    qual: float = 0.0
