"""Prepare a VCF export selection and fetch bounded variant windows."""

from __future__ import annotations

import json
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from storage.sql import count_expr

from .vcf_context import ExportContext, ExportError, ExportRequest, Region, SiteCalls

# A window is one variants query plus one genotypes query held in memory.
# Its size follows the genotype rows it may return, n_variants x n_samples
# (the bound under --keep-reference; sparse ingests return fewer).
_WINDOW_ROWS = 1_000_000
_WINDOW_VARIANTS = 50_000  # without genotype rows (--sites-only), and the cap
_WINDOW_MIN_VARIANTS = 200
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
    missing_recorded = bool(samples) and _missing_recorded(sess, ingest_id)
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
        qual_type=next((typ for column, typ, _ in typed if column == "qual"), "Float32"),
        missing_recorded=missing_recorded,
        selective=bool(request.regions or request.where or request.pass_only),
    )


def _missing_recorded(sess, ingest_id: str) -> bool:
    """True when the ingestion has at least one `missing_genotypes` row.

    A database from before that table, or an ingest run with
    --no-record-missing, has none: a no-call then cannot be told from an
    absent hom-ref call and the absent rule applies to both.
    """
    from storage.db import table_exists

    if not table_exists("missing_genotypes"):
        return False
    n = _rows(
        sess,
        f"SELECT {count_expr()} FROM missing_genotypes WHERE ingest_id = {_q(ingest_id)}",
    )[1][0][0]
    return int(n) > 0


def _variants_per_window(n_samples: int) -> int:
    if n_samples <= 0:
        return _WINDOW_VARIANTS
    return max(_WINDOW_MIN_VARIANTS, min(_WINDOW_VARIANTS, _WINDOW_ROWS // n_samples))


def _windows(sess, context: ExportContext) -> Iterator[tuple[str, int, int]]:
    """Position windows per chromosome, each holding about `per` selected
    variants. The cuts are the positions at every `per`-th variant in
    (pos, ref, alt) order, so a cluster (an exome's genes) is split by
    count; a window ends at the last position of its last variant, so the
    other alleles at that position stay with it."""
    per = _variants_per_window(len(context.samples))
    for chrom, lo, hi, count in context.chrom_stats:
        lo, hi, count = int(lo), int(hi), int(count)
        if count <= per:
            yield chrom, lo, hi
            continue
        cuts = _rows(
            sess,
            "SELECT DISTINCT pos FROM (SELECT pos, row_number() OVER (ORDER BY pos, ref, alt) AS rn "
            f"FROM variants WHERE {context.variant_where} AND chrom = {_q(chrom)}) AS ranked "
            f"WHERE rn % {per} = 0 ORDER BY pos",
        )[1]
        start = lo
        for (cut,) in cuts:
            cut = int(cut)
            if cut < start:
                continue
            yield chrom, start, cut
            start = cut + 1
        if start <= hi:
            yield chrom, start, hi


def _window_calls(
    sess, context: ExportContext, region: Region
) -> dict[tuple, SiteCalls]:
    calls: dict[tuple, SiteCalls] = {}
    if not context.samples:
        return calls
    chrom, start, end = region
    span = f"chrom = {_q(chrom)} AND pos BETWEEN {start} AND {end}"
    where = f"g.ingest_id = {_q(context.ingest_id)} AND g.{span}"
    if context.request.samples:
        where += " AND g.sample_id IN (" + ", ".join(_q(s) for s in context.samples) + ")"
    # --where, --pass-only and regions select variants; fetch only their
    # calls, so a 1% selection does not read the other 99% of genotypes.
    join = ""
    if context.selective:
        join = (
            " INNER JOIN (SELECT pos, ref, alt FROM variants "
            f"WHERE {context.variant_where} AND {span}) AS v "
            "ON g.pos = v.pos AND g.ref = v.ref AND g.alt = v.alt"
        )
    _, rows = _rows(
        sess,
        "SELECT g.pos, g.ref, g.alt, g.sample_id, g.gt, g.gq, g.dp, g.ad_ref, g.ad_alt, g.ft "
        f"FROM genotypes AS g{join} WHERE {where}",
    )
    for pos, ref, alt, sid, *values in rows:
        calls.setdefault((int(pos), ref, alt), {})[sid] = values
    if context.missing_recorded:
        _, rows = _rows(
            sess,
            "SELECT g.pos, g.ref, g.alt, g.sample_id "
            f"FROM missing_genotypes AS g{join} WHERE {where}",
        )
        for pos, ref, alt, sid in rows:
            calls.setdefault((int(pos), ref, alt), {}).setdefault(sid, None)
    return calls


def iter_records(
    sess, context: ExportContext
) -> Iterator[tuple[dict[str, Any], SiteCalls]]:
    for chrom, start, end in _windows(sess, context):
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
