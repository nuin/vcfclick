"""Pure service layer for the vcfclick Textual UI.

This module intentionally imports no Textual symbols. The TUI and tests both
use these functions so parsing, SQL generation, and backend behavior stay
testable without a terminal event loop.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import re
from typing import Any, Literal


class TuiServiceError(Exception):
    """Recoverable user-facing service error."""

    code = "service_error"


class LocusInputError(TuiServiceError):
    """Raised when a gene/range input cannot be interpreted."""

    code = "invalid_locus"


class DatabaseError(TuiServiceError):
    """Raised for missing or invalid active database state."""

    code = "database_error"


class AnnotationUnavailableError(TuiServiceError):
    """Raised when annotation lookup cannot answer a request."""

    code = "annotation_unavailable"


class UnsupportedFeatureError(TuiServiceError):
    """Raised when a backend does not support a TUI operation yet."""

    code = "unsupported_feature"


@dataclass(frozen=True)
class ParsedLocus:
    kind: Literal["gene", "range"]
    label: str
    chrom: str | None
    start_pos: int | None
    end_pos: int | None
    gene_symbol: str | None


@dataclass(frozen=True)
class ResolvedLocus:
    label: str
    chrom: str
    start_pos: int
    end_pos: int
    gene_symbol: str | None = None
    source: Literal["gene", "range"] = "range"


@dataclass(frozen=True)
class QueryResult:
    sql: str
    columns: list[str]
    rows: list[list[Any]]
    row_count: int


@dataclass(frozen=True)
class DatabaseSummary:
    name: str
    path: str
    size_bytes: int
    variants: int | None
    genotypes: int | None
    samples: int | None
    ingestions: int | None


@dataclass(frozen=True)
class LocusSummary:
    locus: ResolvedLocus
    counts: QueryResult
    cohorts: QueryResult
    quality: QueryResult
    preview: QueryResult


_RANGE_RE = re.compile(
    r"^(?P<chrom>[A-Za-z0-9_.-]+):(?P<start>[0-9,]+)-(?P<end>[0-9,]+)$"
)
_GENE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")
_BARE_CHROM_RE = re.compile(r"^chr(?:[0-9]+|[XYM]|MT)$", re.IGNORECASE)


def parse_locus_input(raw: str) -> ParsedLocus:
    """Parse a user-entered gene symbol or `chrom:start-end` range."""
    text = raw.strip()
    if not text:
        raise LocusInputError("Enter a gene symbol or chrom:start-end range.")

    range_match = _RANGE_RE.fullmatch(text)
    if range_match:
        chrom = range_match.group("chrom")
        start_pos = int(range_match.group("start").replace(",", ""))
        end_pos = int(range_match.group("end").replace(",", ""))
        if start_pos < 1 or end_pos < start_pos:
            raise LocusInputError("Range must use positive coordinates with start <= end.")
        label = f"{chrom}:{start_pos}-{end_pos}"
        return ParsedLocus(
            kind="range",
            label=label,
            chrom=chrom,
            start_pos=start_pos,
            end_pos=end_pos,
            gene_symbol=None,
        )

    if _GENE_RE.fullmatch(text) and not _BARE_CHROM_RE.fullmatch(text):
        symbol = text.upper()
        return ParsedLocus(
            kind="gene",
            label=symbol,
            chrom=None,
            start_pos=None,
            end_pos=None,
            gene_symbol=symbol,
        )

    raise LocusInputError("Use a gene symbol or range like chr17:43044295-43125483.")
