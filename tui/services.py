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

    try:
        path = db_path(name)
    except ValueError as exc:
        raise DatabaseError(str(exc)) from exc
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
    try:
        return _query_json(name, cleaned)
    except TuiServiceError:
        raise
    except Exception as exc:
        raise TuiServiceError(f"Unable to execute SQL: {exc}") from exc


_TABLE_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _table_exists(name: str, table: str) -> bool:
    """Return whether `table` exists in the explicitly named database."""
    from storage import backend, get_session

    if not _TABLE_RE.fullmatch(table):
        raise TuiServiceError(f"Unsafe table name: {table!r}")

    sess = get_session(name)
    if backend() == "duckdb":
        sql = (
            "SELECT count(*) FROM information_schema.tables "
            f"WHERE table_name = '{table}'"
        )
    else:
        sql = (
            "SELECT count() FROM system.tables "
            f"WHERE database = currentDatabase() AND name = '{table}'"
        )
    result = sess.query(sql, "CSV").bytes().decode().strip()
    last = result.splitlines()[-1].strip() if result else "0"
    return int(last) > 0


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
            if not _table_exists(name, table):
                counts[table] = None
                continue
            counts[table] = _scalar_count(name, table)
        except TuiServiceError:
            raise
        except Exception as exc:
            raise TuiServiceError(
                f"Unable to summarize database {name!r}: {exc}"
            ) from exc

    return DatabaseSummary(
        name=name,
        path=str(path),
        size_bytes=db_disk_size(name),
        variants=counts["variants"],
        genotypes=counts["genotypes"],
        samples=counts["samples"],
        ingestions=counts["ingestions"],
    )


def stats_summary(name: str, top: int = 20) -> dict[str, Any]:
    """Return stats payload where supported by the active backend."""
    from storage import backend, get_session

    validate_database(name)
    if backend() == "duckdb":
        raise UnsupportedFeatureError("Stats are not implemented on DuckDB yet.")

    from cli.db_stats import _stats_payload

    try:
        return _stats_payload(get_session(name), top)
    except TuiServiceError:
        raise
    except Exception as exc:
        raise TuiServiceError(f"Unable to summarize stats for {name!r}: {exc}") from exc


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


def resolve_locus(parsed: ParsedLocus) -> ResolvedLocus:
    """Resolve parsed user input into concrete coordinates."""
    if parsed.kind == "range":
        assert parsed.chrom is not None
        assert parsed.start_pos is not None
        assert parsed.end_pos is not None
        return ResolvedLocus(
            label=parsed.label,
            chrom=parsed.chrom,
            start_pos=parsed.start_pos,
            end_pos=parsed.end_pos,
            gene_symbol=None,
            source="range",
        )

    assert parsed.gene_symbol is not None
    try:
        import annotations

        gene = annotations.position_for_gene(parsed.gene_symbol)
    except Exception as exc:
        raise AnnotationUnavailableError(
            f"Gene annotations are unavailable: {exc}"
        ) from exc

    if gene is None:
        raise AnnotationUnavailableError(f"Gene {parsed.gene_symbol!r} was not found.")

    return ResolvedLocus(
        label=parsed.gene_symbol,
        chrom=gene.chrom,
        start_pos=int(gene.start_pos),
        end_pos=int(gene.end_pos),
        gene_symbol=gene.gene_symbol,
        source="gene",
    )


def _locus_where(locus: ResolvedLocus, alias: str | None = None) -> str:
    from storage import sql_quote_str

    prefix = f"{alias}." if alias else ""
    return (
        f"{prefix}chrom = {sql_quote_str(locus.chrom)} "
        f"AND {prefix}pos BETWEEN {int(locus.start_pos)} AND {int(locus.end_pos)}"
    )


def _quality_sql(locus: ResolvedLocus) -> str:
    from storage import count_expr

    where = _locus_where(locus)
    return (
        "SELECT "
        f"{count_expr()} AS genotype_rows, "
        "count(gq) AS rows_with_gq, "
        "count(dp) AS rows_with_dp "
        f"FROM genotypes WHERE {where}"
    )


def build_locus_summary(name: str, locus: ResolvedLocus) -> LocusSummary:
    """Run the v1 summary query set for a resolved locus."""
    from storage import count_expr

    where_v = _locus_where(locus, "v")
    where_g = _locus_where(locus, "g")

    counts_sql = (
        "SELECT "
        f"(SELECT {count_expr()} FROM ("
        "SELECT DISTINCT v.ingest_id, v.chrom, v.pos, v.ref, v.alt "
        f"FROM variants v WHERE {where_v}"
        ") variant_rows) AS variants, "
        f"(SELECT {count_expr()} FROM ("
        "SELECT DISTINCT g.ingest_id, g.sample_id "
        f"FROM genotypes g WHERE {where_g}"
        ") carrier_rows) AS carrier_samples, "
        f"(SELECT {count_expr()} FROM genotypes g WHERE {where_g}) "
        "AS non_ref_genotype_rows"
    )
    cohorts_sql = (
        "SELECT cohort, "
        f"{count_expr()} AS samples "
        "FROM (SELECT DISTINCT ingest_id, sample_id, cohort FROM samples) sample_rows "
        "GROUP BY cohort ORDER BY samples DESC, cohort"
    )
    preview_sql = (
        "SELECT chrom, pos, ref, alt, vcf_id, qual, filter, info_AF, info_AC "
        "FROM variants "
        f"WHERE {_locus_where(locus)} "
        "ORDER BY chrom, pos, ref, alt LIMIT 50"
    )

    try:
        return LocusSummary(
            locus=locus,
            counts=_query_json(name, counts_sql),
            cohorts=_query_json(name, cohorts_sql),
            quality=_query_json(name, _quality_sql(locus)),
            preview=_query_json(name, preview_sql),
        )
    except TuiServiceError:
        raise
    except Exception as exc:
        raise TuiServiceError(f"Unable to summarize locus: {exc}") from exc
