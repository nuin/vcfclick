"""Prepare a VCF export selection and fetch bounded variant windows."""

from __future__ import annotations

import json
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from .vcf_context import ExportContext, ExportError, ExportRequest, Region, SiteCalls

_WINDOW_VARIANTS = 50_000
_REGION_RE = re.compile(r"^([^:\s]+)(?::([\d,]+)(?:-([\d,]+))?)?$")


def parse_region(text: str) -> tuple[str, int | None, int | None]:
    m = _REGION_RE.match(text.strip())
    if not m:
        raise ExportError(
            f"invalid region {text!r}; use chr1, chr1:1000 or chr1:1000-2000"
        )
    chrom, a, b = m.group(1), m.group(2), m.group(3)
    start = int(a.replace(",", "")) if a else None
    end = int(b.replace(",", "")) if b else start
    if start is not None and end is not None and end < start:
        raise ExportError(f"invalid region {text!r}: end before start")
    return chrom, start, end


def read_bed(path: str) -> list[tuple[str, int, int]]:
    """BED intervals (0-based, half-open) as 1-based inclusive (chrom, start, end)."""
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) < 3:
            raise ExportError(f"{path}: BED line needs chrom, start, end: {line!r}")
        out.append((f[0], int(f[1]) + 1, int(f[2])))
    return out


def _q(s: str) -> str:
    return "'" + str(s).replace("\\", "\\\\").replace("'", "\\'") + "'"


def _region_sql(regions: list[tuple[str, int | None, int | None]]) -> str:
    from storage.gene_query import chrom_aliases

    parts = []
    for chrom, start, end in regions:
        names = ", ".join(_q(c) for c in chrom_aliases(chrom))
        clause = f"chrom IN ({names})"
        if start is not None:
            clause += f" AND pos BETWEEN {int(start)} AND {int(end)}"
        parts.append(f"({clause})")
    return "(" + " OR ".join(parts) + ")"


def _chrom_key(c: str) -> tuple:
    bare = c[3:] if c.lower().startswith("chr") else c
    if bare.isdigit():
        return (0, int(bare), "")
    return (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(bare.upper(), 3), bare)


def _rows(sess, sql: str) -> tuple[list[str], list[list]]:
    d = json.loads(sess.query(sql, "JSONCompact").bytes().decode())
    return [m["name"] for m in d.get("meta", [])], d.get("data", [])


def _select_ingestion(sess, request: ExportRequest) -> str:
    ingests = [
        row[0]
        for row in _rows(
            sess, "SELECT ingest_id FROM ingestions ORDER BY ingested_at, ingest_id"
        )[1]
    ]
    if not ingests:
        raise ExportError(f"{request.db_name!r} has no ingestions to export")
    if request.ingest_id is None:
        if len(ingests) > 1:
            raise ExportError(
                f"{request.db_name!r} has {len(ingests)} ingestions ({', '.join(ingests)}); "
                "export one at a time with --ingest-id"
            )
        return ingests[0]
    if request.ingest_id not in ingests:
        raise ExportError(
            f"no ingestion {request.ingest_id!r} in {request.db_name!r}; have: {', '.join(ingests)}"
        )
    return request.ingest_id


def _select_samples(sess, ingest_id: str, request: ExportRequest) -> list[str]:
    all_samples = [
        row[0]
        for row in _rows(
            sess,
            f"SELECT DISTINCT sample_id FROM samples WHERE ingest_id = {_q(ingest_id)} ORDER BY sample_id",
        )[1]
    ]
    if request.samples:
        unknown = [sample for sample in request.samples if sample not in all_samples]
        if unknown:
            raise ExportError(
                f"unknown sample(s) in {ingest_id!r}: {', '.join(unknown)}"
            )
        chosen = set(request.samples)
        all_samples = [sample for sample in all_samples if sample in chosen]
    return [] if request.sites_only else all_samples


def _variant_where(ingest_id: str, request: ExportRequest) -> str:
    conds = [f"ingest_id = {_q(ingest_id)}"]
    if request.regions:
        conds.append(_region_sql(request.regions))
    if request.where:
        if ";" in request.where:
            raise ExportError("--where must be a single SQL condition (no ';')")
        conds.append(f"({request.where})")
    if request.pass_only:
        conds.append("(filter IS NULL OR filter IN ('PASS', '.'))")
    return " AND ".join(conds)


def prepare_export(sess, request: ExportRequest) -> ExportContext:
    from storage.db import typed_columns_sql

    ingest_id = _select_ingestion(sess, request)
    samples = _select_samples(sess, ingest_id, request)
    keep_ref = (
        int(
            _rows(
                sess,
                f"SELECT count(*) FROM genotypes WHERE ingest_id = {_q(ingest_id)} AND gt = 0",
            )[1][0][0]
        )
        > 0
    )
    absent_gt = "./." if (keep_ref or request.absent_as == "nocall") else "0/0"
    variant_where = _variant_where(ingest_id, request)
    typed = _rows(sess, typed_columns_sql("variants"))[1]
    info_cols = [
        (column, typ, int(flag))
        for column, typ, flag in typed
        if column.startswith("info_") and column != "info_extra"
    ]
    chrom_stats = _rows(
        sess,
        f"SELECT chrom, min(pos), max(pos), count(*) FROM variants WHERE {variant_where} GROUP BY chrom",
    )[1]
    chrom_stats.sort(key=lambda row: _chrom_key(row[0]))
    return ExportContext(
        request=request,
        ingest_id=ingest_id,
        samples=samples,
        keep_ref=keep_ref,
        absent_gt=absent_gt,
        variant_where=variant_where,
        info_cols=info_cols,
        has_extra=any(column == "info_extra" for column, _, _ in typed),
        chrom_stats=chrom_stats,
    )


def _windows(context: ExportContext) -> Iterator[tuple[str, int, int]]:
    # Window by chromosome and position to bound memory for large databases.
    for chrom, lo, hi, count in context.chrom_stats:
        lo, hi, count = int(lo), int(hi), int(count)
        n_windows = max(1, -(-count // _WINDOW_VARIANTS))
        width = (hi - lo) // n_windows + 1
        for w in range(n_windows):
            yield chrom, lo + w * width, min(hi, lo + (w + 1) * width - 1)


def _window_calls(
    sess, context: ExportContext, region: Region
) -> dict[tuple, SiteCalls]:
    calls: dict[tuple, SiteCalls] = {}
    if not context.samples:
        return calls
    chrom, start, end = region
    sample_filter = ""
    if context.request.samples:
        sample_filter = (
            " AND sample_id IN (" + ", ".join(_q(s) for s in context.samples) + ")"
        )
    _, rows = _rows(
        sess,
        "SELECT pos, ref, alt, sample_id, gt, gq, dp, ad_ref, ad_alt, ft FROM genotypes "
        f"WHERE ingest_id = {_q(context.ingest_id)} AND chrom = {_q(chrom)} "
        f"AND pos BETWEEN {start} AND {end}{sample_filter}",
    )
    for pos, ref, alt, sid, *values in rows:
        calls.setdefault((int(pos), ref, alt), {})[sid] = values
    return calls


def iter_records(
    sess, context: ExportContext
) -> Iterator[tuple[dict[str, Any], SiteCalls]]:
    for chrom, start, end in _windows(context):
        window = (
            f"{context.variant_where} AND chrom = {_q(chrom)} "
            f"AND pos BETWEEN {start} AND {end}"
        )
        cols, rows = _rows(
            sess,
            f"SELECT {', '.join(context.select_columns)} FROM variants WHERE {window} ORDER BY pos, ref, alt",
        )
        if not rows:
            continue
        calls = _window_calls(sess, context, (chrom, start, end))
        for row in rows:
            rec = dict(zip(cols, row))
            yield rec, calls.get((int(rec["pos"]), rec["ref"], rec["alt"]), {})
