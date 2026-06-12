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


def list_database_names() -> list[str]:
    """Return local database names sorted by the storage layer."""
    from storage import list_dbs

    return list_dbs()


def validate_database(name: str) -> str:
    """Return `name` when it exists, otherwise raise a recoverable error."""
    from storage import db_path

    path = db_path(name)
    if not path.exists():
        raise DatabaseError(f"Database {name!r} does not exist.")
    return name


def _query_json(name: str, sql: str) -> QueryResult:
    """Execute SQL and normalize chDB/DuckDB JSONCompact output."""
    from storage import get_session

    validate_database(name)
    sess = get_session(name)
    raw = sess.query(sql, "JSONCompact").bytes().decode()
    parsed = json.loads(raw)
    columns = [m["name"] for m in parsed.get("meta", [])]
    rows = parsed.get("data", [])
    return QueryResult(
        sql=sql,
        columns=columns,
        rows=rows,
        row_count=len(rows),
    )


def execute_sql(name: str, sql: str) -> QueryResult:
    """Run user-entered SQL against a named database."""
    cleaned = sql.strip()
    if not cleaned:
        raise TuiServiceError("Enter a SQL query.")
    return _query_json(name, cleaned)


def _scalar_count(name: str, table: str) -> int | None:
    from storage import count_expr

    result = _query_json(name, f"SELECT {count_expr()} AS n FROM {table}")
    if not result.rows:
        return None
    return int(result.rows[0][0])


def database_summary(name: str) -> DatabaseSummary:
    """Return basic DB metadata for the Operations mode."""
    from storage import db_disk_size, db_path

    validate_database(name)
    path = db_path(name)
    counts: dict[str, int | None] = {}
    for table in ("variants", "genotypes", "samples", "ingestions"):
        try:
            counts[table] = _scalar_count(name, table)
        except Exception:
            counts[table] = None

    return DatabaseSummary(
        name=name,
        path=str(path),
        size_bytes=db_disk_size(name),
        variants=counts["variants"],
        genotypes=counts["genotypes"],
        samples=counts["samples"],
        ingestions=counts["ingestions"],
    )


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
            raise LocusInputError(
                "Range must use positive coordinates with start <= end."
            )
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
