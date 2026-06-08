"""Ingest data into an active chDB store from a Parquet dump.

This is the inverse of `vcfclick db dump`. The same Parquet files that
fall out of a dump are the wire format that comes back in here — so a
dump on machine A followed by an ingest-parquet on machine B reconstructs
the cohort byte-for-byte (modulo the new ingest_id + cohort labels the
caller supplies).

Beyond round-tripping vcfclick → vcfclick, this is the "be a citizen of
the columnar stack" path: any external tool that can write Parquet
matching the schemas in ingest/_arrow.py (DuckDB, polars, Spark, ...)
can land data here without ever touching cyvcf2. The schemas are the
public interchange contract.

The ingest_id and cohort columns of the source files are NOT honoured.
They get rewritten to the values the caller supplies via --ingest-id
and --cohort. That preserves the model used everywhere else in
vcfclick: an ingest is one atomic upload identified by one ingest_id,
and replays of the same ingest_id replace the prior data.

Atomicity model mirrors the VCF path. Phase 1 reads each source file's
Parquet footer and rejects schema-mismatches before touching chDB.
Phase 2 holds the per-(DB, ingest_id) file lock, deletes the prior
rows under this ingest_id, and runs the INSERT...SELECT FROM file()
calls in one go. A failure in Phase 1 leaves the prior data intact;
a failure in Phase 2 rolls back to empty under this ingest_id.
"""

from __future__ import annotations

import json
import logging
import uuid
from pathlib import Path

import pyarrow.parquet as pq

from ingest._arrow import (
    GENOTYPES_ARROW_SCHEMA,
    INGESTIONS_ARROW_SCHEMA,
    SAMPLES_ARROW_SCHEMA,
    TABLE_COLUMNS,
    VARIANTS_ARROW_SCHEMA,
    column_list_sql,
)
from storage import (
    count_expr,
    get_session,
    ingest_id_lock,
    insert_via_parquet,
    parquet_file_expr,
    rollback_ingest,
    sql_quote_str,
    validate_ingest_id,
)

log = logging.getLogger(__name__)


_EXPECTED_SCHEMAS = {
    "variants": VARIANTS_ARROW_SCHEMA,
    "genotypes": GENOTYPES_ARROW_SCHEMA,
    "samples": SAMPLES_ARROW_SCHEMA,
}

# Server-side default column present on every ReplacingMergeTree table
# (the version column, `DEFAULT now()` in the SQL DDL). `db dump`
# includes it in the output Parquet; ingest tolerates it on input but
# does NOT carry it through — chDB re-defaults it on the new INSERT.
# Allow it as a known extra; reject anything else extra as a typo.
_SERVER_DEFAULT_COLUMNS = frozenset({"ingested_at"})


def _validate_parquet_schema(path: Path, table: str) -> None:
    """Reject the file early if its column set doesn't agree with the
    locked Arrow schema for `table`. Types are left to chDB to enforce
    on INSERT — it raises a precise error for incompatible types.

    Allows but ignores `ingested_at` so that vcfclick's own dumps round
    trip cleanly (the dump emits it, the ingest re-defaults it on
    re-insert).
    """
    expected = {f.name for f in _EXPECTED_SCHEMAS[table]}
    found = set(pq.read_schema(path).names)
    missing = expected - found
    extra = found - expected - _SERVER_DEFAULT_COLUMNS
    problems = []
    if missing:
        problems.append(f"missing columns: {sorted(missing)}")
    if extra:
        problems.append(f"unexpected columns: {sorted(extra)}")
    if problems:
        raise ValueError(
            f"{path} doesn't conform to the {table} schema:\n  " + "\n  ".join(problems)
        )


def _import_with_override(
    table: str,
    src: Path,
    ingest_id: str,
    cohort: str | None = None,
) -> None:
    """INSERT into chDB from a Parquet file, rewriting ingest_id (and
    cohort, for the samples table) to the caller's values via the
    SELECT list. chDB reads the file in place — no restaging needed.
    """
    cols = TABLE_COLUMNS[table]
    select_exprs = []
    for c in cols:
        if c == "ingest_id":
            select_exprs.append(f"{sql_quote_str(ingest_id)} AS ingest_id")
        elif c == "cohort" and cohort is not None:
            select_exprs.append(f"{sql_quote_str(cohort)} AS cohort")
        else:
            select_exprs.append(f'"{c}"')
    sess = get_session()
    sess.query(
        f"INSERT INTO {table} ({column_list_sql(cols)}) "
        f"SELECT {', '.join(select_exprs)} "
        f"FROM {parquet_file_expr(str(src))}"
    )


def _count_under_ingest(table: str, ingest_id: str) -> int:
    sess = get_session()
    raw = (
        sess.query(
            f"SELECT {count_expr()} FROM {table} WHERE ingest_id = "
            f"{sql_quote_str(ingest_id)} FORMAT JSONCompact"
        )
        .bytes()
        .decode()
    )
    return int(json.loads(raw)["data"][0][0])


def _distinct_samples_from_genotypes_pq(path: Path) -> list[str]:
    sess = get_session()
    raw = (
        sess.query(
            f"SELECT DISTINCT sample_id FROM {parquet_file_expr(str(path))} "
            f"ORDER BY sample_id FORMAT JSONCompact"
        )
        .bytes()
        .decode()
    )
    return [r[0] for r in json.loads(raw)["data"]]


def ingest_from_parquet(
    dump_dir: str,
    cohort: str,
    ingest_id: str | None = None,
) -> str:
    """Ingest a Parquet dump into the active database under a new label.

    `dump_dir` is a directory containing at least `variants.parquet`,
    and optionally `genotypes.parquet` and `samples.parquet`. The
    schemas must match the ones in ingest/_arrow.py — call them via
    `from ingest._arrow import VARIANTS_ARROW_SCHEMA, ...` to read the
    canonical types and column orders.

    Sample handling, in order of precedence:
      * samples.parquet present → imported as-is, with ingest_id and
        cohort columns rewritten to the caller's values.
      * samples.parquet missing but genotypes.parquet present → the
        sample list is derived via `SELECT DISTINCT sample_id` against
        the genotypes file. `sex` is left NULL.
      * Neither present → no samples row is written. Valid for a
        variants-only cohort summary (e.g. an external AF table).
    """
    if ingest_id is None:
        ingest_id = str(uuid.uuid4())
    validate_ingest_id(ingest_id)

    dump = Path(dump_dir)
    if not dump.is_dir():
        raise ValueError(f"{dump} is not a directory")

    variants_pq = dump / "variants.parquet"
    genotypes_pq = dump / "genotypes.parquet"
    samples_pq = dump / "samples.parquet"

    if not variants_pq.exists():
        raise FileNotFoundError(
            f"{dump} is missing variants.parquet (required for any ingest)"
        )

    log.info("[parquet-ingest] dir:       %s", dump)
    log.info("[parquet-ingest] ingest_id: %s", ingest_id)
    log.info("[parquet-ingest] cohort:    %s", cohort)

    _validate_parquet_schema(variants_pq, "variants")
    if genotypes_pq.exists():
        _validate_parquet_schema(genotypes_pq, "genotypes")
    if samples_pq.exists():
        _validate_parquet_schema(samples_pq, "samples")

    with ingest_id_lock(ingest_id):
        return _ingest_parquet_locked(
            variants_pq,
            genotypes_pq if genotypes_pq.exists() else None,
            samples_pq if samples_pq.exists() else None,
            cohort,
            ingest_id,
            dump,
        )


def _ingest_parquet_locked(
    variants_pq: Path,
    genotypes_pq: Path | None,
    samples_pq: Path | None,
    cohort: str,
    ingest_id: str,
    dump: Path,
) -> str:
    # Local import to avoid pulling the cyvcf2-heavy vcf_load module
    # at import time of this lighter Parquet path.
    from ingest.vcf_load import _ensure_schema

    _ensure_schema()

    commit_started = False
    try:
        # Once we call rollback_ingest, prior rows under this ingest_id
        # are gone. Mark the boundary so the except arm only rolls back
        # if we crossed it.
        commit_started = True
        rollback_ingest(ingest_id)

        if samples_pq is not None:
            _import_with_override("samples", samples_pq, ingest_id, cohort=cohort)
        elif genotypes_pq is not None:
            sample_ids = _distinct_samples_from_genotypes_pq(genotypes_pq)
            insert_via_parquet(
                "samples",
                SAMPLES_ARROW_SCHEMA,
                [
                    {
                        "ingest_id": ingest_id,
                        "sample_id": sid,
                        "cohort": cohort,
                        "sex": None,
                    }
                    for sid in sample_ids
                ],
            )

        _import_with_override("variants", variants_pq, ingest_id)
        if genotypes_pq is not None:
            _import_with_override("genotypes", genotypes_pq, ingest_id)

        n_variants = _count_under_ingest("variants", ingest_id)
        n_samples = _count_under_ingest("samples", ingest_id)
        insert_via_parquet(
            "ingestions",
            INGESTIONS_ARROW_SCHEMA,
            [
                {
                    "ingest_id": ingest_id,
                    "cohort": cohort,
                    "vcf_path": f"parquet://{dump.resolve()}",
                    "n_variants": n_variants,
                    "n_samples": n_samples,
                }
            ],
        )
    except BaseException:
        if commit_started:
            log.warning(
                "[parquet-ingest] failed mid-commit — rolling back %s", ingest_id
            )
            try:
                rollback_ingest(ingest_id)
            except Exception as rb_err:
                log.error("[parquet-ingest] rollback FAILED: %s", rb_err)
        raise

    log.info("[parquet-ingest] done. %d variants, %d samples.", n_variants, n_samples)
    return ingest_id
