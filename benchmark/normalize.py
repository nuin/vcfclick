"""Reference-based VCF normalization for `vcfclick benchmark`.

Pure, backend-neutral. `trim` is minimal-representation parsimony; `left_align`
adds reference-based left-shifting (Tan et al. 2015) with an extend-on-demand
buffer (no fixed window). Positions are 1-based; `fetch` is 0-based half-open.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

Fetch = Callable[[str, int, int], str]


def trim(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """VCF minimal representation: drop shared suffix then shared prefix.

    Never empties an allele (keeps ≥1 base each); prefix trimming advances POS.
    """
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def left_align(
    fetch: Fetch, chrom: str, pos: int, ref: str, alt: str
) -> tuple[int, str, str]:
    """Left-align and trim an indel against the reference.

    Rolls the variant leftward through a repeat to its leftmost representation,
    fetching further-left reference bases on demand until it stabilises or the
    contig start is reached. Idempotent.
    """
    while True:
        if len(ref) == 0 or len(alt) == 0:
            if pos <= 1:  # at contig start — cannot extend further
                break
            base = fetch(chrom, pos - 2, pos - 1)  # base immediately left of POS
            ref, alt, pos = base + ref, base + alt, pos - 1
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        else:
            break
    # Restore parsimony on the left (keep ≥1 base each).
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if not ref or not alt:
        raise ValueError(
            f"normalization produced an empty allele at {chrom}:{pos} "
            "(malformed input: no left anchor at the contig start)"
        )
    return pos, ref, alt


def decompose_mnp(pos: int, ref: str, alt: str) -> list[tuple[int, str, str]]:
    """Split an equal-length multibase substitution (MNP) into per-position SNPs.

    Only differing positions become SNPs. Non-MNPs (SNPs, indels) pass through
    unchanged. OFF by default in the pipeline: atomizing an unphased het MNP
    loses cis/trans phase, which a keyed comparer cannot recover (a P2 concern).
    """
    if len(ref) != len(alt) or len(ref) <= 1:
        return [(pos, ref, alt)]
    return [(pos + i, ref[i], alt[i]) for i in range(len(ref)) if ref[i] != alt[i]]


def atomize(pos: int, ref: str, alt: str) -> list[tuple[int, str, str]]:
    """Split an allele into primitives: SNPs plus at most one indel.

    Trims to the minimal representation first, then:
      * equal-length multibase (MNP) -> one SNP per differing position;
      * a pure indel (one side is the single anchor base) -> unchanged, it is
        already primitive;
      * a genuine complex allele (both sides multibase, lengths differ) -> the
        leading substitutions as SNPs, then the length change as one anchored
        indel.

    This converges callers that pack substitutions into one record with callers
    that emit them separately. It does NOT converge differently *anchored*
    spellings of the same indel — that is reference left-alignment's job
    (`left_align`), which `--reference` applies first.
    """
    pos, ref, alt = trim(pos, ref, alt)
    if len(ref) == len(alt):
        if len(ref) == 1:
            return [(pos, ref, alt)]
        return [(pos + i, ref[i], alt[i]) for i in range(len(ref)) if ref[i] != alt[i]]
    if len(ref) == 1 or len(alt) == 1:
        return [(pos, ref, alt)]  # pure indel: already a primitive
    # Complex: substitutions over the shared span, then the indel. The last
    # shared position stays as the indel's anchor base.
    n = min(len(ref), len(alt))
    out = [(pos + i, ref[i], alt[i]) for i in range(n - 1) if ref[i] != alt[i]]
    out.append((pos + n - 1, ref[n - 1 :], alt[n - 1 :]))
    return out


@dataclass(frozen=True)
class AltRow:
    """One biallelic record after splitting a multiallelic site."""

    pos: int
    ref: str
    alt: str
    gt: tuple[int, ...]
    other_alt_present: bool
    locus_id: tuple


def split_multiallelic(
    pos: int, ref: str, alts: Sequence[str], gt: tuple[int, ...]
) -> list[AltRow]:
    """Split into biallelic rows, remapping GT `-m -` style.

    For alt index `i`, GT alleles equal to `i` become 1 and all other non-missing
    alleles become 0; `other_alt_present` records that a *different* alt allele was
    carried at the site (so a het-alt `1/2` is never mistaken for a plain `0/1`).
    """
    locus_id = (pos, ref, tuple(alts))
    rows: list[AltRow] = []
    for i, alt in enumerate(alts, start=1):
        remapped = tuple(-1 if a < 0 else (1 if a == i else 0) for a in gt)
        other_alt = any(a > 0 and a != i for a in gt)
        rows.append(AltRow(pos, ref, alt, remapped, other_alt, locus_id))
    return rows


def canonical_gt(alleles: tuple[int, ...], phased: bool = False) -> str:
    """Render a genotype; unphased genotypes are order-independent (sorted)."""
    ordered = alleles if phased else tuple(sorted(alleles))
    sep = "|" if phased else "/"
    return sep.join("." if a < 0 else str(a) for a in ordered)
