"""INFO/AA normalisation and site polarisation (popgen.ancestral)."""

from __future__ import annotations

import pytest

from popgen.ancestral import (
    ANCESTRAL_IS_ALT,
    ANCESTRAL_IS_REF,
    AncestralAllele,
    normalise_aa,
    polarise,
)


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("A", AncestralAllele("A", True)),
        ("a", AncestralAllele("A", False)),
        ("G|||", AncestralAllele("G", True)),  # 1000 Genomes style
        ("t|-|-|", AncestralAllele("T", False)),
        (" C ", AncestralAllele("C", True)),
        ("ACT", AncestralAllele("ACT", True)),  # indel ancestral sequence
        ("aCt", AncestralAllele("ACT", False)),  # mixed case: low confidence
        (".", None),
        ("-", None),
        ("?", None),
        ("", None),
        ("N", None),
        ("n", None),
        ("NN", None),
        ("|A", None),  # nothing before the first |
        ("1", None),
        (None, None),
    ],
)
def test_normalise_aa(raw, expected):
    assert normalise_aa(raw) == expected


@pytest.mark.parametrize(
    "ref,alt,aa,mode,expected",
    [
        ("A", "G", "A", "aa", ANCESTRAL_IS_REF),
        ("A", "G", "g", "aa", ANCESTRAL_IS_ALT),
        ("A", "G", "g", "aa-high", None),  # low confidence rejected
        ("A", "G", "G", "aa-high", ANCESTRAL_IS_ALT),
        ("A", "G", "C", "aa", None),  # matches neither allele
        ("A", "G", ".", "aa", None),
        ("A", "G", None, "aa", None),
        ("A", "G", None, "ref", ANCESTRAL_IS_REF),
        ("A", "G", "G", "ref", ANCESTRAL_IS_REF),  # ref mode ignores AA
        ("A", "G", "A", "none", None),
        ("AT", "A", "at", "aa", ANCESTRAL_IS_REF),
    ],
)
def test_polarise(ref, alt, aa, mode, expected):
    assert polarise(ref, alt, aa, mode) == expected


def test_polarise_rejects_unknown_mode():
    with pytest.raises(ValueError):
        polarise("A", "G", "A", "outgroup")
