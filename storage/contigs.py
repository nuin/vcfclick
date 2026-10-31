"""Contig names: the one place every alias set is defined.

Pure (no SQL, no storage session), so it can be imported anywhere. VCFs in
the wild name the same chromosome `chrX`, `X`, `chr23` or `23`; each
feature that needs "the X chromosome" or "the autosomes" asks this module
instead of keeping its own list.

The alias sets are the human conventions vcfclick supports. They are
compared case-sensitively and in full (a name such as `chrX_KI270880v1_alt`
is none of X, Y, PAR or mitochondrial, and is not an autosome either).

`popgen.analysis.is_autosome` is a different, deliberately species-agnostic
classifier ("anything unrecognised is an autosome") that its own tests pin;
this module is strict: only the human autosomes 1-22 count.
"""

from __future__ import annotations

import re

# Ordered tuples: SQL built from them is stable text.
X_NAMES: tuple[str, ...] = ("chrX", "X", "chr23", "23")
Y_NAMES: tuple[str, ...] = ("chrY", "Y", "chr24", "24")
# Pseudo-autosomal labels some pipelines give the PAR as a contig of its own.
PAR_NAMES: tuple[str, ...] = ("chrXY", "XY", "chr25", "25", "PAR", "chrPAR")
MITO_NAMES: tuple[str, ...] = ("chrM", "M", "MT", "chrMT", "chr26", "26")

NON_AUTOSOMAL: frozenset[str] = frozenset(X_NAMES + Y_NAMES + PAR_NAMES + MITO_NAMES)

_AUTOSOME_RE = re.compile(r"^(?:chr)?([1-9]|1[0-9]|2[0-2])$")


def is_x(chrom: str) -> bool:
    return chrom in X_NAMES


def is_y(chrom: str) -> bool:
    return chrom in Y_NAMES


def is_par_label(chrom: str) -> bool:
    """A contig that IS the pseudo-autosomal region (chrXY, PAR, ...). The
    PAR coordinates inside chrX/chrY are a separate, build-aware question."""
    return chrom in PAR_NAMES


def is_mito(chrom: str) -> bool:
    return chrom in MITO_NAMES


def is_sex_or_mt(chrom: str) -> bool:
    """X, Y, a PAR label or the mitochondrion."""
    return chrom in NON_AUTOSOMAL


def is_autosome(chrom: str) -> bool:
    """chr1..chr22 or 1..22 (no leading zeros). Everything else is not: the
    sex and mitochondrial aliases, unplaced, alt, decoy and viral contigs."""
    return _AUTOSOME_RE.match(chrom) is not None


def autosome_names() -> tuple[str, ...]:
    """Every accepted spelling of the autosomes, for `chrom IN (...)`:
    `1`..`22` then `chr1`..`chr22`."""
    bare = tuple(str(n) for n in range(1, 23))
    return bare + tuple(f"chr{n}" for n in bare)


def canonical_chrom(chrom: str) -> str:
    """The `chr`-prefixed GRCh38 name the annotation store uses for a contig,
    whatever the VCF calls it: `17` and `chr17` are `chr17`; every X, Y and
    mitochondrial alias above is `chrX`, `chrY` and `chrM`. A name that is none
    of these (an alt or unplaced contig, a PAR label) is returned unchanged."""
    if chrom in X_NAMES:
        return "chrX"
    if chrom in Y_NAMES:
        return "chrY"
    if chrom in MITO_NAMES:
        return "chrM"
    m = _AUTOSOME_RE.match(chrom)
    return f"chr{m.group(1)}" if m else chrom


def chrom_aliases(chrom: str) -> tuple[str, ...]:
    """Every spelling of the contig `chrom` names, for `chrom IN (...)`: the
    whole alias group for X, Y and the mitochondrion, `17` and `chr17` for an
    autosome, the name alone otherwise. Ordered, so SQL built from it is
    stable text."""
    canonical = canonical_chrom(chrom)
    if canonical == "chrX":
        return X_NAMES
    if canonical == "chrY":
        return Y_NAMES
    if canonical == "chrM":
        return MITO_NAMES
    if _AUTOSOME_RE.match(canonical):
        return (canonical[3:], canonical)
    return (chrom,)
