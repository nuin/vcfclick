"""Ancestral-allele (`info_AA`) normalisation and site polarisation.

`variants.info_AA` is stored verbatim from the VCF. Conventions differ
between sources — the 1000 Genomes releases write `AA=G|||` (the base,
then indel context fields), Ensembl EPO files use upper case for a
high-confidence call and lower case for a low-confidence one, and `.`,
`-`, `N` or an empty string mean unknown. `normalise_aa` reduces all of
those to (allele, high_confidence) or None.
"""

from __future__ import annotations

from dataclasses import dataclass

ANCESTRAL_MODES = ("aa", "aa-high", "ref", "none")

_UNKNOWN = {"", ".", "-", "?"}


@dataclass(frozen=True)
class AncestralAllele:
    """A normalised ancestral allele: upper-cased sequence + confidence."""

    allele: str
    high_confidence: bool


def normalise_aa(value: str | None) -> AncestralAllele | None:
    """Normalise a raw INFO/AA value.

    Takes the text before the first `|`, strips whitespace; returns None
    for unknown (`.`, `-`, `?`, empty, or all-`N`). Upper case means high
    confidence, lower case low confidence (Ensembl EPO convention); a
    mixed-case value counts as low confidence.
    """
    if value is None:
        return None
    token = str(value).split("|", 1)[0].strip()
    if token in _UNKNOWN or set(token.upper()) == {"N"}:
        return None
    if not token.isalpha():
        return None
    return AncestralAllele(allele=token.upper(), high_confidence=token.isupper())


# Orientation of a polarised site: is the ancestral allele REF or ALT?
ANCESTRAL_IS_REF = 0
ANCESTRAL_IS_ALT = 1


def polarise(ref: str, alt: str, aa: str | None, mode: str) -> int | None:
    """Which allele is ancestral at a site, or None if it can't be said.

    mode `aa`      — use INFO/AA at high or low confidence;
    mode `aa-high` — use INFO/AA only when upper case (high confidence);
    mode `ref`     — assume REF is ancestral everywhere;
    mode `none`    — never polarise (folded statistics only).

    With `aa`/`aa-high`, a site is polarised only when the ancestral
    allele equals REF (ancestral = REF) or ALT (derived = REF, so counts
    are flipped); any other value leaves it unpolarised.
    """
    if mode not in ANCESTRAL_MODES:
        raise ValueError(f"unknown ancestral mode {mode!r}")
    if mode == "none":
        return None
    if mode == "ref":
        return ANCESTRAL_IS_REF
    norm = normalise_aa(aa)
    if norm is None or (mode == "aa-high" and not norm.high_confidence):
        return None
    if norm.allele == ref.upper():
        return ANCESTRAL_IS_REF
    if norm.allele == alt.upper():
        return ANCESTRAL_IS_ALT
    return None
