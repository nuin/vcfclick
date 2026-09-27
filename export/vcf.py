"""Export a database, or a slice of it, back to VCF (`vcfclick db export`).

One ingestion per file: different ingestions carry different samples and
site lists, so merging them would invent genotypes for samples that were
never genotyped at a site.

Genotypes live in a sparse table, so reconstructing a VCF needs one rule
for samples with no row at a site:
- ingestion made with --keep-reference: 0/0 calls are stored, so absent
  means a no-call and is written ./.;
- otherwise only non-reference calls were stored and a hom-ref call cannot
  be told apart from a no-call. Absent is written 0/0 by default
  (--absent-as nocall writes ./. instead), and the header says so.

Exported per genotype: GT, GQ, DP, AD, FT. Phase is not preserved: the
store keeps a phased flag but not which haplotype carries the ALT, so all
calls are written unphased. PL/GL are not exported. INFO is rebuilt from
the typed info_* columns plus the info_extra map. FILTER PASS and unset (.)
are both stored as NULL at ingest, so both come back as ".".

.vcf.gz output is BGZF (tabix-indexable); .vcf is plain text; "-" is stdout.
"""

from __future__ import annotations

import json
import re
import struct
import sys
import tempfile
import zlib
from datetime import datetime, timezone
from pathlib import Path
from typing import IO

# Known INFO fields: (Number, Type). Anything else falls back to its SQL type.
_INFO_SPEC = {
    "AC": ("A", "Integer"),
    "AF": ("A", "Float"),
    "AN": ("1", "Integer"),
    "DP": ("1", "Integer"),
    "MQ": ("1", "Float"),
    "MQ0": ("1", "Integer"),
    "NS": ("1", "Integer"),
    "BQ": ("1", "Float"),
    "SB": ("1", "Float"),
    "END": ("1", "Integer"),
    "CIGAR": ("A", "String"),
    "AA": ("1", "String"),
    "QD": ("1", "Float"),
    "FS": ("1", "Float"),
    "SOR": ("1", "Float"),
    "MQRankSum": ("1", "Float"),
    "ReadPosRankSum": ("1", "Float"),
    "ExcessHet": ("1", "Float"),
    "InbreedingCoeff": ("1", "Float"),
    "MLEAC": ("A", "Integer"),
    "MLEAF": ("A", "Float"),
    "BaseQRankSum": ("1", "Float"),
    "ClippingRankSum": ("1", "Float"),
    "AD": ("R", "Integer"),
}
_FORMAT_SPEC = {
    "GT": ("1", "String", "Genotype"),
    "GQ": ("1", "Integer", "Genotype quality"),
    "DP": ("1", "Integer", "Read depth"),
    "AD": ("R", "Integer", "Allelic depths for the ref and alt alleles"),
    "FT": ("1", "String", "Per-sample genotype filter"),
}
_GT = {0: "0/0", 1: "0/1", 2: "1/1", -1: "./1"}
_WINDOW_VARIANTS = 50_000
_REGION_RE = re.compile(r"^([^:\s]+)(?::([\d,]+)(?:-([\d,]+))?)?$")


class ExportError(ValueError):
    """A user-facing problem with the export request."""


# --- BGZF ---------------------------------------------------------------------------

_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


class BgzfWriter:
    """Minimal BGZF writer (the blocked gzip htslib and tabix read)."""

    _MAX_INPUT = 65280  # leaves room for the block header in the 64 KiB limit

    def __init__(self, fh: IO[bytes]):
        self._fh = fh
        self._buf = bytearray()

    def write(self, text: str) -> None:
        self._buf += text.encode()
        while len(self._buf) >= self._MAX_INPUT:
            self._block(bytes(self._buf[: self._MAX_INPUT]))
            del self._buf[: self._MAX_INPUT]

    def _block(self, data: bytes) -> None:
        comp = zlib.compressobj(6, zlib.DEFLATED, -15)
        cdata = comp.compress(data) + comp.flush()
        if len(cdata) > 65536 - 26:  # incompressible: split and retry
            half = len(data) // 2
            self._block(data[:half])
            self._block(data[half:])
            return
        bsize = 18 + len(cdata) + 8 - 1
        header = (
            b"\x1f\x8b\x08\x04"
            + b"\x00\x00\x00\x00"
            + b"\x00\xff"
            + struct.pack("<H", 6)
            + b"BC"
            + struct.pack("<HH", 2, bsize)
        )
        self._fh.write(
            header
            + cdata
            + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
        )

    def close(self) -> None:
        if self._buf:
            self._block(bytes(self._buf))
            self._buf.clear()
        self._fh.write(_BGZF_EOF)
        self._fh.close()


# --- request parsing ---------------------------------------------------------------


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


# --- value formatting --------------------------------------------------------------


def _num(v) -> str:
    if isinstance(v, float) or (isinstance(v, str) and "." in v):
        f = float(v)
        return str(int(f)) if f.is_integer() else f"{f:.6g}"
    return str(v)


def _rows(sess, sql: str) -> tuple[list[str], list[list]]:
    d = json.loads(sess.query(sql, "JSONCompact").bytes().decode())
    return [m["name"] for m in d.get("meta", [])], d.get("data", [])


def _map_items(v) -> list[tuple[str, str]]:
    if not v:
        return []
    if isinstance(v, dict):
        return list(v.items())
    if isinstance(v, list):  # DuckDB can hand maps back as key/value pairs
        return [
            (p["key"], p["value"]) if isinstance(p, dict) else (p[0], p[1]) for p in v
        ]
    return []


# --- the export -------------------------------------------------------------------


def export_vcf(
    sess,
    db_name: str,
    out: str,
    *,
    ingest_id: str | None = None,
    samples: list[str] | None = None,
    regions: list[tuple[str, int | None, int | None]] | None = None,
    where: str | None = None,
    pass_only: bool = False,
    sites_only: bool = False,
    absent_as: str = "ref",
    version: str = "",
) -> dict:
    """Write the VCF. Returns a summary dict (variants, samples, path)."""
    from storage.db import typed_columns_sql

    ingests = [
        r[0]
        for r in _rows(
            sess, "SELECT ingest_id FROM ingestions ORDER BY ingested_at, ingest_id"
        )[1]
    ]
    if not ingests:
        raise ExportError(f"{db_name!r} has no ingestions to export")
    if ingest_id is None:
        if len(ingests) > 1:
            raise ExportError(
                f"{db_name!r} has {len(ingests)} ingestions ({', '.join(ingests)}); "
                "export one at a time with --ingest-id"
            )
        ingest_id = ingests[0]
    elif ingest_id not in ingests:
        raise ExportError(
            f"no ingestion {ingest_id!r} in {db_name!r}; have: {', '.join(ingests)}"
        )
    iq = _q(ingest_id)

    all_samples = [
        r[0]
        for r in _rows(
            sess,
            f"SELECT DISTINCT sample_id FROM samples WHERE ingest_id = {iq} ORDER BY sample_id",
        )[1]
    ]
    if samples:
        unknown = [s for s in samples if s not in all_samples]
        if unknown:
            raise ExportError(
                f"unknown sample(s) in {ingest_id!r}: {', '.join(unknown)}"
            )
        chosen = set(samples)
        out_samples = [s for s in all_samples if s in chosen]
    else:
        out_samples = all_samples
    if sites_only:
        out_samples = []

    keep_ref = (
        int(
            _rows(
                sess,
                f"SELECT count(*) FROM genotypes WHERE ingest_id = {iq} AND gt = 0",
            )[1][0][0]
        )
        > 0
    )
    absent_gt = "./." if (keep_ref or absent_as == "nocall") else "0/0"

    conds = [f"ingest_id = {iq}"]
    if regions:
        conds.append(_region_sql(regions))
    if where:
        if ";" in where:
            raise ExportError("--where must be a single SQL condition (no ';')")
        conds.append(f"({where})")
    if pass_only:
        conds.append("(filter IS NULL OR filter IN ('PASS', '.'))")
    vwhere = " AND ".join(conds)

    typed = _rows(sess, typed_columns_sql("variants"))[1]
    info_cols = [
        (c, t, int(flag))
        for c, t, flag in typed
        if c.startswith("info_") and c != "info_extra"
    ]
    has_extra = any(c == "info_extra" for c, _, _ in typed)
    select = (
        ["chrom", "pos", "vcf_id", "ref", "alt", "qual", "filter"]
        + [c for c, _, _ in info_cols]
        + (["info_extra"] if has_extra else [])
    )

    # Window by chromosome, and by position within big chromosomes, so memory
    # stays bounded no matter how large the database is.
    chrom_stats = _rows(
        sess,
        f"SELECT chrom, min(pos), max(pos), count(*) FROM variants WHERE {vwhere} GROUP BY chrom",
    )[1]
    chrom_stats.sort(key=lambda r: _chrom_key(r[0]))

    used_info: dict[str, str] = {}  # key -> SQL type (for fallback header types)
    used_extra: set[str] = set()
    used_filters: set[str] = set()
    used_fmt = ["GT"]
    fmt_seen = {"GQ": False, "DP": False, "AD": False, "FT": False}
    sample_index = {s: i for i, s in enumerate(out_samples)}
    n_variants = 0

    tmp = tempfile.TemporaryFile("w+", encoding="utf-8")
    body_fmt_placeholder = (
        "\x00FMT\x00"  # FORMAT column filled in once all fields are known
    )
    try:
        for chrom, lo, hi, count in chrom_stats:
            lo, hi, count = int(lo), int(hi), int(count)
            n_windows = max(1, -(-count // _WINDOW_VARIANTS))
            width = (hi - lo) // n_windows + 1
            for w in range(n_windows):
                a, b = lo + w * width, min(hi, lo + (w + 1) * width - 1)
                win = f"{vwhere} AND chrom = {_q(chrom)} AND pos BETWEEN {a} AND {b}"
                cols, vrows = _rows(
                    sess,
                    f"SELECT {', '.join(select)} FROM variants WHERE {win} ORDER BY pos, ref, alt",
                )
                if not vrows:
                    continue
                calls: dict[tuple, dict[str, list]] = {}
                if out_samples:
                    sfilter = ""
                    if samples:
                        sfilter = (
                            " AND sample_id IN ("
                            + ", ".join(_q(s) for s in out_samples)
                            + ")"
                        )
                    _, grows = _rows(
                        sess,
                        "SELECT pos, ref, alt, sample_id, gt, gq, dp, ad_ref, ad_alt, ft FROM genotypes "
                        f"WHERE ingest_id = {iq} AND chrom = {_q(chrom)} AND pos BETWEEN {a} AND {b}{sfilter}",
                    )
                    for pos, ref, alt, sid, gt, gq, dp, adr, ada, ft in grows:
                        calls.setdefault((int(pos), ref, alt), {})[sid] = [
                            gt,
                            gq,
                            dp,
                            adr,
                            ada,
                            ft,
                        ]
                for row in vrows:
                    rec = dict(zip(cols, row))
                    n_variants += 1
                    tmp.write(
                        _record(
                            rec,
                            info_cols,
                            has_extra,
                            calls.get((int(rec["pos"]), rec["ref"], rec["alt"]), {}),
                            out_samples,
                            sample_index,
                            absent_gt,
                            used_info,
                            used_extra,
                            used_filters,
                            fmt_seen,
                            body_fmt_placeholder,
                        )
                    )

        used_fmt += [f for f in ("GQ", "DP", "AD", "FT") if fmt_seen[f]]
        header = _header(
            db_name=db_name,
            ingest_id=ingest_id,
            version=version,
            chroms=[r[0] for r in chrom_stats],
            used_info=used_info,
            used_extra=used_extra,
            used_filters=used_filters,
            used_fmt=used_fmt if out_samples else [],
            samples=out_samples,
            keep_ref=keep_ref,
            absent_gt=absent_gt,
            request={
                "regions": regions,
                "where": where,
                "pass_only": pass_only,
                "samples": samples,
                "sites_only": sites_only,
            },
        )
        _write(out, header, tmp, body_fmt_placeholder, used_fmt)
    finally:
        tmp.close()
    return {
        "path": out,
        "ingest_id": ingest_id,
        "variants": n_variants,
        "samples": len(out_samples),
        "absent_written_as": absent_gt,
    }


def _record(
    rec,
    info_cols,
    has_extra,
    site_calls,
    out_samples,
    sample_index,
    absent_gt,
    used_info,
    used_extra,
    used_filters,
    fmt_seen,
    fmt_placeholder,
) -> str:
    info = []
    ad_ref = ad_alt = None
    for col, typ, is_flag in info_cols:
        v = rec.get(col)
        key = col[5:]
        if key == "AD_ref":
            ad_ref = v
            continue
        if key == "AD_alt":
            ad_alt = v
            continue
        if is_flag:
            if v not in (None, 0, "0", False):
                info.append(key)
                used_info[key] = "Flag"
            continue
        if v is None:
            continue
        info.append(f"{key}={_num(v)}")
        used_info[key] = typ
    if ad_ref is not None and ad_alt is not None:
        info.append(f"AD={_num(ad_ref)},{_num(ad_alt)}")
        used_info["AD"] = "Integer"
    if has_extra:
        for k, v in _map_items(rec.get("info_extra")):
            info.append(k if v in ("", None) else f"{k}={v}")
            used_extra.add(k)

    filt = rec.get("filter") or "."
    if filt not in ("PASS", "."):
        for f in filt.split(";"):
            used_filters.add(f)
    qual = "." if rec.get("qual") is None else _num(rec["qual"])
    fields = [
        rec["chrom"],
        str(int(rec["pos"])),
        rec.get("vcf_id") or ".",
        rec["ref"],
        rec["alt"],
        qual,
        filt,
        ";".join(info) if info else ".",
    ]
    if not out_samples:
        return "\t".join(fields) + "\n"

    cols = [absent_gt] * len(out_samples)
    for sid, (gt, gq, dp, adr, ada, ft) in site_calls.items():
        j = sample_index.get(sid)
        if j is None:
            continue
        parts = {
            "GT": _GT.get(int(gt), "./."),
            "GQ": ".",
            "DP": ".",
            "AD": ".",
            "FT": ".",
        }
        if gq is not None:
            parts["GQ"] = str(int(gq))
            fmt_seen["GQ"] = True
        if dp is not None:
            parts["DP"] = str(int(dp))
            fmt_seen["DP"] = True
        if adr is not None and ada is not None:
            parts["AD"] = f"{int(adr)},{int(ada)}"
            fmt_seen["AD"] = True
        if ft is not None:
            parts["FT"] = ft
            fmt_seen["FT"] = True
        # Stash all fields; the final FORMAT decides which are written.
        cols[j] = "\x01".join(parts[k] for k in ("GT", "GQ", "DP", "AD", "FT"))
    return "\t".join(fields + [fmt_placeholder] + cols) + "\n"


def _header(
    *,
    db_name,
    ingest_id,
    version,
    chroms,
    used_info,
    used_extra,
    used_filters,
    used_fmt,
    samples,
    keep_ref,
    absent_gt,
    request,
) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "##fileformat=VCFv4.2",
        f"##source=vcfclick {version}".rstrip(),
        f"##vcfclick_export=<Database={db_name},IngestID={ingest_id},Date={now}>",
    ]
    shown = {k: v for k, v in request.items() if v}
    if shown:
        lines.append("##vcfclick_export_filters=" + json.dumps(shown, default=str))
    if samples:
        if keep_ref:
            lines.append(
                '##vcfclick_note="Ingested with --keep-reference: stored 0/0 calls are written 0/0; '
                'absent genotypes are no-calls, written ./."'
            )
        else:
            lines.append(
                f'##vcfclick_note="Only non-reference genotypes were stored at ingest, so a hom-ref '
                f"call cannot be told from a no-call: absent genotypes are written {absent_gt}. "
                'Phase is not preserved (all calls unphased)."'
            )
    lines.append(
        '##vcfclick_note_filter="FILTER PASS and unset (.) are stored alike at ingest, so both '
        'are written as ."'
    )
    lines += [f"##contig=<ID={c}>" for c in chroms]
    for f in sorted(used_filters):
        lines.append(f'##FILTER=<ID={f},Description="From the source VCF">')
    for key in sorted(used_info):
        typ = used_info[key]
        if typ == "Flag":
            lines.append(
                f'##INFO=<ID={key},Number=0,Type=Flag,Description="From the source VCF">'
            )
            continue
        number, vtype = _INFO_SPEC.get(key, ("1", _vcf_type(typ)))
        lines.append(
            f'##INFO=<ID={key},Number={number},Type={vtype},Description="From the source VCF">'
        )
    for key in sorted(used_extra - set(used_info)):
        lines.append(
            f'##INFO=<ID={key},Number=.,Type=String,Description="From the source VCF (untyped)">'
        )
    for f in used_fmt:
        n, t, d = _FORMAT_SPEC[f]
        lines.append(f'##FORMAT=<ID={f},Number={n},Type={t},Description="{d}">')
    cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if samples:
        cols += ["FORMAT", *samples]
    return "\n".join(lines + ["\t".join(cols)]) + "\n"


def _vcf_type(sql_type: str) -> str:
    t = sql_type.lower()
    if "int" in t:
        return "Integer"
    if "float" in t or "double" in t or "real" in t or "decimal" in t:
        return "Float"
    return "String"


def _write(out: str, header: str, tmp, placeholder: str, used_fmt: list[str]) -> None:
    keep = [i for i, f in enumerate(("GT", "GQ", "DP", "AD", "FT")) if f in used_fmt]
    fmt = ":".join(used_fmt)

    def finish(line: str) -> str:
        if placeholder not in line:
            return line
        fixed, _, rest = line.partition(placeholder)
        out_cols = []
        for cell in rest.rstrip("\n").lstrip("\t").split("\t"):
            if "\x01" in cell:
                parts = cell.split("\x01")
                cell = ":".join(parts[i] for i in keep)
            out_cols.append(cell)
        return fixed + fmt + "\t" + "\t".join(out_cols) + "\n"

    tmp.seek(0)
    if out == "-":
        sys.stdout.write(header)
        for line in tmp:
            sys.stdout.write(finish(line))
        sys.stdout.flush()
        return
    path = Path(out)
    if path.suffix == ".gz":
        w = BgzfWriter(path.open("wb"))
        w.write(header)
        for line in tmp:
            w.write(finish(line))
        w.close()
    else:
        with path.open("w", encoding="utf-8") as f:
            f.write(header)
            for line in tmp:
                f.write(finish(line))
