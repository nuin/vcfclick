"""Annotation service. Single tier, single license.

    from annotations import position_for_gene, gene_at, clinvar_lookup

The transcript / exon / CDS / UTR hierarchy lives in
annotations/transcripts.py (Phase 2; not yet implemented).
"""

from annotations.db import (
    GeneRange,
    ClinVarRecord,
    position_for_gene,
    gene_at,
    clinvar_lookup,
)

__all__ = [
    "GeneRange",
    "ClinVarRecord",
    "position_for_gene",
    "gene_at",
    "clinvar_lookup",
]
