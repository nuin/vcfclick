"""Result values and recoverable errors shared by TUI services."""

from dataclasses import dataclass
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
