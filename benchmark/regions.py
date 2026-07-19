"""Confident-region (BED) membership for `vcfclick benchmark`.

Pure, backend-neutral. Intervals are 0-based half-open, matching BED. Membership
uses "start containment": a 1-based VCF POS is in-region iff its 0-based start
lies within some interval. Overlapping/adjacent intervals are merged into
disjoint, sorted per-contig arrays for O(log n) `numpy.searchsorted` lookup.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from benchmark.model import NormRecord


def _merge(intervals: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent [start, end) pairs into disjoint, sorted ones."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:  # overlap or touch → extend
            prev_start, prev_end = merged[-1]
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


class ConfRegions:
    """Disjoint per-contig confident intervals with start-containment lookup."""

    def __init__(self, intervals: Iterable[tuple[str, int, int]]) -> None:
        by_contig: dict[str, list[tuple[int, int]]] = {}
        for chrom, start0, end0 in intervals:
            by_contig.setdefault(chrom, []).append((int(start0), int(end0)))
        self._starts: dict[str, np.ndarray] = {}
        self._ends: dict[str, np.ndarray] = {}
        for chrom, ivals in by_contig.items():
            merged = _merge(ivals)
            self._starts[chrom] = np.array([s for s, _ in merged], dtype=np.int32)
            self._ends[chrom] = np.array([e for _, e in merged], dtype=np.int32)

    @classmethod
    def from_bed(cls, path: str | Path) -> "ConfRegions":
        """Parse a BED file (chrom, start0, end0); skip blanks/comments/headers."""
        intervals: list[tuple[str, int, int]] = []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split("\t")
                intervals.append((fields[0], int(fields[1]), int(fields[2])))
        return cls(intervals)

    def contains(self, chrom: str, pos1: int) -> bool:
        """True iff 1-based POS is in a confident interval (start containment)."""
        starts = self._starts.get(chrom)
        if starts is None:
            return False
        pos0 = pos1 - 1
        idx = int(np.searchsorted(starts, pos0, side="right")) - 1
        return idx >= 0 and pos0 < int(self._ends[chrom][idx])

    def tag(self, records: Iterable[NormRecord]) -> list[NormRecord]:
        """Return copies with `in_conf` set from the normalized locus start POS."""
        return [
            dataclasses.replace(r, in_conf=self.contains(r.chrom, r.pos))
            for r in records
        ]
