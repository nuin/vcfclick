"""Local-haplotype reconciliation core for `vcfclick benchmark` (P2).

The vcfeval/xcmp-style residual pass: replay each side's variants onto a shared
reference window to synthesize alternate haplotype sequences, then declare two
clusters equivalent when they yield a common (multiset of) haplotype sequence(s).
This resolves representation-different but sequence-equivalent calls that P1's
keyed match miscounts as FP+FN (an MNP `AC>GT` vs two SNPs `A>G`+`C>T`, a complex
indel vs its atomized spelling, a shifted homopolymer deletion).

SQL-free and backend-neutral (pure Python). Coordinates: `apply_haplotype` takes
a 1-based inclusive window; internal footprints are 0-based half-open.
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from typing import Callable, Sequence

from benchmark.model import NormRecord

Fetch = Callable[[str, int, int], str]  # (chrom, start0, end0) -> uppercased bases

_MAX_HAPLOTYPE_PATHS = 4096


def _alleles(r: NormRecord) -> list[str]:
    """Called genotype indices, dropping missing/empty slots."""
    return [p for p in r.gt.replace("|", "/").split("/") if p not in (".", "")]


def _carried(r: NormRecord) -> bool:
    """True when the record's genotype places its alt on at least one haplotype."""
    return "1" in _alleles(r)


def _is_het(r: NormRecord) -> bool:
    """Alt on exactly one haplotype (e.g. 0/1)."""
    a = _alleles(r)
    return "0" in a and "1" in a


def _is_hom(r: NormRecord) -> bool:
    """Alt on both haplotypes (1/1, or a haploid `1` pre-expansion)."""
    a = _alleles(r)
    return "1" in a and "0" not in a


def _is_phased(r: NormRecord) -> bool:
    return "|" in r.gt


def _phased_alt_slot(r: NormRecord) -> int:
    """Haplotype index (0 or 1) carrying the alt in a phased het `a|b`."""
    return r.gt.split("|").index("1")


def _footprint(r: NormRecord) -> tuple[int, int]:
    """0-based half-open reference span the record occupies."""
    start0 = r.pos - 1
    return (start0, start0 + len(r.ref))


def _window1(records: Sequence[NormRecord]) -> tuple[int, int]:
    """Tight 1-based inclusive window spanning every record's footprint."""
    start1 = min(r.pos for r in records)
    end1 = max(r.pos + len(r.ref) - 1 for r in records)
    return start1, end1


def apply_haplotype(
    fetch: Fetch,
    chrom: str,
    win_start1: int,
    win_end1: int,
    records: Sequence[NormRecord],
) -> str | None:
    """Alt sequence over the 1-based inclusive window `[win_start1, win_end1]`.

    Applies every carried (genotype contains a `1`) record onto the reference
    slice, splicing right-to-left (highest POS first) so a length-changing indel
    never invalidates a not-yet-applied leftward variant's offset. Returns None if
    a record's REF disagrees with the window or two carried records overlap
    (infeasible), so pathological input can never coincidentally match.
    """
    carried = [r for r in records if _carried(r)]
    return _apply_variants(fetch, chrom, win_start1, win_end1, carried)


def _apply_variants(
    fetch: Fetch,
    chrom: str,
    win_start1: int,
    win_end1: int,
    variants: Sequence[NormRecord],
) -> str | None:
    """Apply an explicit variant set to one haplotype; None if infeasible.

    Unlike `apply_haplotype` this trusts the caller's phasing (no genotype
    filtering) and returns None when a variant's REF does not match the window
    (defensive) or when two co-applied variants overlap.
    """
    # Reject overlapping footprints up front — after length-changing splices the
    # REF check alone can't tell an overlap from a legitimate edit (codex).
    prev_start = None
    for r in sorted(variants, key=_footprint, reverse=True):
        start0, end0 = _footprint(r)
        if prev_start is not None and end0 > prev_start:
            return None
        prev_start = start0

    seq = fetch(chrom, win_start1 - 1, win_end1)
    for r in sorted(variants, key=lambda r: r.pos, reverse=True):
        off = r.pos - win_start1
        if off < 0 or seq[off : off + len(r.ref)] != r.ref:
            return None
        seq = seq[:off] + r.alt + seq[off + len(r.ref) :]
    return seq


def cluster(records: Sequence[NormRecord], flank: int) -> list[list[NormRecord]]:
    """Group records whose flank-padded footprints overlap, per chromosome.

    Sweep-line merge: sorted by footprint start, a record joins the open cluster
    when its start lies within `flank` bp of the running cluster end. `flank`
    governs only clustering distance, never the replay window.
    """
    by_chrom: dict[str, list[NormRecord]] = defaultdict(list)
    for r in records:
        by_chrom[r.chrom].append(r)

    clusters: list[list[NormRecord]] = []
    for chrom in sorted(by_chrom):
        recs = sorted(by_chrom[chrom], key=_footprint)
        cur: list[NormRecord] = []
        cur_end = 0
        for r in recs:
            s0, e0 = _footprint(r)
            if not cur:
                cur, cur_end = [r], e0
            elif s0 - flank <= cur_end:
                cur.append(r)
                cur_end = max(cur_end, e0)
            else:
                clusters.append(cur)
                cur, cur_end = [r], e0
        if cur:
            clusters.append(cur)
    return clusters


def _enumerate_diplotypes(
    fetch: Fetch,
    chrom: str,
    win_start1: int,
    win_end1: int,
    members: Sequence[NormRecord],
    max_paths: int,
) -> set[tuple[str, str]] | None:
    """Feasible unordered diplotype (hapA, hapB) sequence pairs for one side.

    Hom-alt records sit on both haplotypes. A *phased* het is fixed to its
    declared haplotype (its alt slot in `a|b`), preserving cis/trans relations so
    a trans `0|1`+`1|0` pair is never collapsed onto one haplotype and matched to
    a cis MNP. Only *unphased* hets are enumerated across haplotypes. Returns None
    when the enumerated count exceeds `max_paths` (too complex).
    """
    homs = [r for r in members if _carried(r) and _is_hom(r)]
    phased = [r for r in members if _carried(r) and _is_het(r) and _is_phased(r)]
    unphased = [r for r in members if _carried(r) and _is_het(r) and not _is_phased(r)]
    if 2 ** len(unphased) > max_paths:
        return None
    phased_a = [r for r in phased if _phased_alt_slot(r) == 0]
    phased_b = [r for r in phased if _phased_alt_slot(r) == 1]
    pin = not phased  # only free to pin the frame when no phased het fixes it
    dips: set[tuple[str, str]] = set()
    for assignment in itertools.product("AB", repeat=len(unphased)):
        if pin and unphased and assignment[0] != "A":  # dedupe the A/B swap
            continue
        hap_a = homs + phased_a + [h for h, x in zip(unphased, assignment) if x == "A"]
        hap_b = homs + phased_b + [h for h, x in zip(unphased, assignment) if x == "B"]
        a = _apply_variants(fetch, chrom, win_start1, win_end1, hap_a)
        b = _apply_variants(fetch, chrom, win_start1, win_end1, hap_b)
        if a is None or b is None:  # overlapping same-haplotype variants
            continue
        dips.add(tuple(sorted((a, b))))
    return dips


def haplotype_equivalent(
    fetch: Fetch,
    cluster_truth: Sequence[NormRecord],
    cluster_query: Sequence[NormRecord],
    flank: int = 30,
    *,
    max_paths: int = _MAX_HAPLOTYPE_PATHS,
) -> bool:
    """True iff truth and query describe a common diplotype over the shared window.

    The window spans both sides' footprints, so identical flanking sequence
    cancels on each side. The all-hom case (canonical `1/1`) reduces to a single
    haplotype-sequence comparison. When hets are present, feasible diplotypes are
    enumerated (bounded by `max_paths`) and the sides match only when they share
    a diplotype; an over-budget enumeration yields False (never a misclassify).

    `flank` is accepted for a uniform call signature but does not affect the
    equivalence verdict (the window is tight).
    """
    _ = flank
    recs = list(cluster_truth) + list(cluster_query)
    if not recs:
        return True
    chrom = recs[0].chrom
    win_start1, win_end1 = _window1(recs)

    has_het = any(_is_het(r) for r in recs)
    if not has_het:
        t_seq = apply_haplotype(fetch, chrom, win_start1, win_end1, cluster_truth)
        q_seq = apply_haplotype(fetch, chrom, win_start1, win_end1, cluster_query)
        return t_seq is not None and t_seq == q_seq

    td = _enumerate_diplotypes(
        fetch, chrom, win_start1, win_end1, cluster_truth, max_paths
    )
    qd = _enumerate_diplotypes(
        fetch, chrom, win_start1, win_end1, cluster_query, max_paths
    )
    if td is None or qd is None:
        return False
    return bool(td & qd)
