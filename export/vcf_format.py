"""Render VCF records and the header derived from their observed fields."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any

from .vcf_context import FORMAT_PLACEHOLDER, ExportContext, SiteCalls

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


def _num(v) -> str:
    if isinstance(v, float) or (isinstance(v, str) and "." in v):
        f = float(v)
        return str(int(f)) if f.is_integer() else f"{f:.6g}"
    return str(v)


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


def record(rec: dict[str, Any], site_calls: SiteCalls, context: ExportContext) -> str:
    info = []
    ad_ref = ad_alt = None
    for col, typ, is_flag in context.info_cols:
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
                context.used_info[key] = "Flag"
            continue
        if v is None:
            continue
        info.append(f"{key}={_num(v)}")
        context.used_info[key] = typ
    if ad_ref is not None and ad_alt is not None:
        info.append(f"AD={_num(ad_ref)},{_num(ad_alt)}")
        context.used_info["AD"] = "Integer"
    if context.has_extra:
        for k, v in _map_items(rec.get("info_extra")):
            info.append(k if v in ("", None) else f"{k}={v}")
            context.used_extra.add(k)

    filt = rec.get("filter") or "."
    if filt not in ("PASS", "."):
        for f in filt.split(";"):
            context.used_filters.add(f)
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
    if not context.samples:
        return "\t".join(fields) + "\n"

    cols = _genotypes(site_calls, context)
    return "\t".join(fields + [FORMAT_PLACEHOLDER] + cols) + "\n"


def _genotypes(site_calls: SiteCalls, context: ExportContext) -> list[str]:
    cols = [context.absent_gt] * len(context.samples)
    for sid, (gt, gq, dp, adr, ada, ft) in site_calls.items():
        j = context.sample_index.get(sid)
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
            context.fmt_seen.add("GQ")
        if dp is not None:
            parts["DP"] = str(int(dp))
            context.fmt_seen.add("DP")
        if adr is not None and ada is not None:
            parts["AD"] = f"{int(adr)},{int(ada)}"
            context.fmt_seen.add("AD")
        if ft is not None:
            parts["FT"] = ft
            context.fmt_seen.add("FT")
        # Stash all fields; the final FORMAT decides which are written.
        cols[j] = "\x01".join(parts[k] for k in ("GT", "GQ", "DP", "AD", "FT"))
    return cols


def header(context: ExportContext, used_fmt: list[str]) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "##fileformat=VCFv4.2",
        f"##source=vcfclick {context.request.version}".rstrip(),
        f"##vcfclick_export=<Database={context.request.db_name},IngestID={context.ingest_id},Date={now}>",
    ]
    shown = {k: v for k, v in context.request.filters().items() if v}
    if shown:
        lines.append("##vcfclick_export_filters=" + json.dumps(shown, default=str))
    if context.samples:
        if context.keep_ref:
            lines.append(
                '##vcfclick_note="Ingested with --keep-reference: stored 0/0 calls are written 0/0; '
                'absent genotypes are no-calls, written ./."'
            )
        else:
            lines.append(
                f'##vcfclick_note="Only non-reference genotypes were stored at ingest, so a hom-ref '
                f"call cannot be told from a no-call: absent genotypes are written {context.absent_gt}. "
                'Phase is not preserved (all calls unphased)."'
            )
    lines.append(
        '##vcfclick_note_filter="FILTER PASS and unset (.) are stored alike at ingest, so both '
        'are written as ."'
    )
    lines += [f"##contig=<ID={c}>" for c, *_ in context.chrom_stats]
    for f in sorted(context.used_filters):
        lines.append(f'##FILTER=<ID={f},Description="From the source VCF">')
    for key in sorted(context.used_info):
        typ = context.used_info[key]
        if typ == "Flag":
            lines.append(
                f'##INFO=<ID={key},Number=0,Type=Flag,Description="From the source VCF">'
            )
            continue
        number, vtype = _INFO_SPEC.get(key, ("1", _vcf_type(typ)))
        lines.append(
            f'##INFO=<ID={key},Number={number},Type={vtype},Description="From the source VCF">'
        )
    for key in sorted(context.used_extra - set(context.used_info)):
        lines.append(
            f'##INFO=<ID={key},Number=.,Type=String,Description="From the source VCF (untyped)">'
        )
    for f in used_fmt:
        n, t, d = _FORMAT_SPEC[f]
        lines.append(f'##FORMAT=<ID={f},Number={n},Type={t},Description="{d}">')
    cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if context.samples:
        cols += ["FORMAT", *context.samples]
    return "\n".join(lines + ["\t".join(cols)]) + "\n"


def _vcf_type(sql_type: str) -> str:
    t = sql_type.lower()
    if "int" in t:
        return "Integer"
    if "float" in t or "double" in t or "real" in t or "decimal" in t:
        return "Float"
    return "String"
