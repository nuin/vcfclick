"""chrX regions: non-PAR, PAR and the build that decides between them.

Pure (no click, no session): the intervals, the build rule and the SQL
predicate built from them. `npx_intervals(build)` is the only definition of
hemizygous territory in vcfclick; everything else on chrX is PAR and is
treated like an autosome.

    GRCh38   PAR1 10,001-2,781,479      PAR2 155,701,383-156,030,895
    GRCh37   PAR1 60,001-2,699,520      PAR2 154,931,044-155,260,560

The build is never assumed. A chrX POS above 155,270,560 (the GRCh37 chrX
length) proves GRCh38; otherwise it is `unknown`, and the two bands where
the builds disagree (`AMBIGUOUS_BAND`) are left unclassified: excluded from
the X-linked model and from the X guard, and counted. CHM13/T2T is not
supported (its chrX never exceeds the GRCh37 length, so it resolves to
`unknown`).

Classification is by POS: an indel that straddles a PAR edge is classified
by its start.
"""

from __future__ import annotations

from storage.contigs import PAR_NAMES
from storage.contigs import X_NAMES as X_NAMES

GRCH37, GRCH38, UNKNOWN = "GRCh37", "GRCh38", "unknown"
BUILDS: tuple[str, ...] = (GRCH37, GRCH38)
BUILD_CHOICES: tuple[str, ...] = ("auto", GRCH37, GRCH38)

# The GRCh37 chrX length. A position above it exists only on GRCh38.
GRCH37_CHRX_LENGTH = 155_270_560

Intervals = tuple[tuple[int, int], ...]

# Closed intervals of non-PAR X (NPX), where a male is hemizygous.
NPX_INTERVALS: dict[str, Intervals] = {
    GRCH38: ((2_781_480, 155_701_382),),
    GRCH37: ((2_699_521, 154_931_043), (155_260_561, 155_270_560)),
    # the intersection of the two builds' NPX
    UNKNOWN: ((2_781_480, 154_931_043), (155_260_561, 155_701_382)),
}

# Where the builds disagree: NPX on one, PAR on the other. Only an unknown
# build has them.
AMBIGUOUS_BAND: dict[str, Intervals] = {
    GRCH38: (),
    GRCH37: (),
    UNKNOWN: ((2_699_521, 2_781_479), (154_931_044, 155_260_560)),
}

# Sex contigs vcfclick does not interpret: the PAR given as a contig of its
# own (storage.contigs.PAR_NAMES). Data on them is not X data.
UNSUPPORTED_SEX_CONTIGS: tuple[str, ...] = PAR_NAMES


def npx_intervals(build: str) -> Intervals:
    return NPX_INTERVALS[build]


def ambiguous_band(build: str) -> Intervals:
    return AMBIGUOUS_BAND[build]


def _inside(pos: int, intervals: Intervals) -> bool:
    return any(lo <= pos <= hi for lo, hi in intervals)


def region(chrom: str, pos: int, build: str) -> str:
    """`npx`, `par`, `ambiguous` (an unknown build's bands) or `other` (not
    chrX)."""
    if chrom not in X_NAMES:
        return "other"
    if _inside(pos, NPX_INTERVALS[build]):
        return "npx"
    if _inside(pos, AMBIGUOUS_BAND[build]):
        return "ambiguous"
    return "par"


def is_npx(chrom: str, pos: int, build: str) -> bool:
    return region(chrom, pos, build) == "npx"


def is_hemizygous(sex: str | None, chrom: str, pos: int, build: str) -> bool:
    """Is a call of a person of this sex at (chrom, pos) haploid? True only
    for a male on sure non-PAR X."""
    return sex == "male" and is_npx(chrom, pos, build)


def resolve_build(option: str, max_x_pos: int | None) -> tuple[str, str, bool]:
    """(build, build_source, conflict) from the `--build` option and the
    largest chrX position in the data.

    An explicit build is taken as given (source `option`); GRCh37 with a
    position above the GRCh37 chrX length is a conflict. `auto` resolves to
    GRCh38 (source `data`) when such a position exists, else `unknown`."""
    beyond_37 = max_x_pos is not None and int(max_x_pos) > GRCH37_CHRX_LENGTH
    if option in BUILDS:
        return option, "option", option == GRCH37 and beyond_37
    if option != "auto":
        raise ValueError(f"unknown build {option!r}")
    if beyond_37:
        return GRCH38, "data", False
    return UNKNOWN, "unknown", False


def npx_sql(build: str, *, chrom: str = "chrom", pos: str = "pos") -> str:
    """SQL boolean: (chrom, pos) is sure non-PAR X under `build`."""
    from storage.sql import chrom_in_sql, in_intervals_sql

    return (
        f"({chrom_in_sql(chrom, X_NAMES)} AND "
        f"{in_intervals_sql(pos, NPX_INTERVALS[build])})"
    )


def band_sql(build: str, *, chrom: str = "chrom", pos: str = "pos") -> str:
    """SQL boolean: (chrom, pos) lies in an ambiguous band of `build`."""
    from storage.sql import chrom_in_sql, in_intervals_sql

    return (
        f"({chrom_in_sql(chrom, X_NAMES)} AND "
        f"{in_intervals_sql(pos, AMBIGUOUS_BAND[build])})"
    )
