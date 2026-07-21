"""Reference-FASTA access and contig validation for `vcfclick benchmark`.

Backend-neutral. Owns 0-based/1-based discipline and chrom-name aliasing so
truth, query, and BED contigs all resolve to the FASTA's naming.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

_MITO_ALIASES = ("M", "MT", "chrM", "chrMT")
_MITO_KEYS = {"m", "mt", "chrm", "chrmt"}


class BenchmarkError(Exception):
    """Recoverable, user-facing benchmark error."""


class ContigError(BenchmarkError):
    """Contig naming or length mismatch against the reference."""


def build_alias_map(names: Iterable[str]) -> dict[str, str]:
    """Map every alias spelling of a contig to its canonical FASTA name.

    Handles `chr` prefix add/strip and the mitochondrial `M`/`MT`/`chrM`/`chrMT`
    family. Canonical names always map to themselves and win over aliases.
    """
    amap: dict[str, str] = {}
    for canon in names:
        low = canon.lower()
        if low.startswith("chr"):
            amap.setdefault(canon[3:], canon)
        else:
            amap.setdefault("chr" + canon, canon)
        if low in _MITO_KEYS:
            for alias in _MITO_ALIASES:
                amap.setdefault(alias, canon)
    # Canonical names win over any alias collision.
    for canon in names:
        amap[canon] = canon
    return amap


def canonical_contig(name: str, known: Iterable[str]) -> str | None:
    """Resolve `name` to a known FASTA contig name, or None if unresolvable."""
    return build_alias_map(known).get(name)


class Reference:
    """Thin wrapper over a pyfaidx FASTA with alias-aware access."""

    def __init__(self, path: str | Path) -> None:
        from pyfaidx import Fasta

        self._fa = Fasta(str(path))
        self._alias = build_alias_map(list(self._fa.keys()))

    def _canon(self, chrom: str) -> str:
        canon = self._alias.get(chrom)
        if canon is None:
            raise ContigError(f"contig {chrom!r} is not in the reference")
        return canon

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        """Reference bases for 0-based half-open [start0, end0), uppercased."""
        return str(self._fa[self._canon(chrom)][start0:end0]).upper()

    def contig_len(self, chrom: str) -> int:
        return len(self._fa[self._canon(chrom)])

    @property
    def contigs(self) -> set[str]:
        return set(self._fa.keys())

    def validate_length(self, chrom: str, length: int) -> None:
        """Raise if `length` disagrees with the reference contig length."""
        actual = self.contig_len(chrom)
        if actual != length:
            raise ContigError(
                f"contig {chrom!r} length {length} != reference length {actual}"
            )
