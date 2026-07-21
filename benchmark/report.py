"""Report writers for `vcfclick benchmark`: strict GA4GH-ish summary.csv, a
vcfclick.summary.csv with an Engine column, JSON provenance/metrics, and a
self-contained HTML. All ratios computed in Python from integer counts."""

from __future__ import annotations

import csv
import html
import json
import os

from benchmark.constants import FILTER_ALL, FILTER_PASS, VT_INDEL, VT_SNP
from benchmark.metrics import metrics_from_counts

# Strict canonical column order (no Engine column — stays hap.py-parser-safe).
SUMMARY_COLUMNS = [
    "Type",
    "Filter",
    "TRUTH.TOTAL",
    "TRUTH.TP",
    "TRUTH.FN",
    "QUERY.TOTAL",
    "QUERY.TP",
    "QUERY.FP",
    "QUERY.UNK",
    "METRIC.Recall",
    "METRIC.Precision",
    "METRIC.Frac_NA",
    "METRIC.F1_Score",
]

# One row per (Type, Filter) in this deterministic order.
_ROW_ORDER = [
    (VT_SNP, FILTER_PASS),
    (VT_SNP, FILTER_ALL),
    (VT_INDEL, FILTER_PASS),
    (VT_INDEL, FILTER_ALL),
]


def _summary_rows(agg: dict) -> list[dict]:
    """Build one summary record per (Type, Filter), values from int counts."""
    rows = []
    for vtype, filter_view in _ROW_ORDER:
        c = agg.get((filter_view, vtype), {})
        truth_tp = c.get("truth_tp", 0)
        truth_fn = c.get("truth_fn", 0)
        query_tp = c.get("query_tp", 0)
        query_fp = c.get("query_fp", 0)
        query_unk = c.get("query_unk", 0)
        m = metrics_from_counts(truth_tp, truth_fn, query_tp, query_fp, query_unk)
        rows.append(
            {
                "Type": vtype,
                "Filter": filter_view,
                "TRUTH.TOTAL": truth_tp + truth_fn,
                "TRUTH.TP": truth_tp,
                "TRUTH.FN": truth_fn,
                "QUERY.TOTAL": query_tp + query_fp + query_unk,
                "QUERY.TP": query_tp,
                "QUERY.FP": query_fp,
                "QUERY.UNK": query_unk,
                "METRIC.Recall": m.recall,
                "METRIC.Precision": m.precision,
                "METRIC.Frac_NA": m.frac_na,
                "METRIC.F1_Score": m.f1,
            }
        )
    return rows


def _write_summary_csv(rows: list[dict], path: str) -> None:
    """Strict summary.csv: canonical columns, no Engine."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(SUMMARY_COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in SUMMARY_COLUMNS])


def _write_vcfclick_summary_csv(rows: list[dict], engine: str, path: str) -> None:
    """vcfclick.summary.csv: leading Engine column, then canonical columns."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Engine", *SUMMARY_COLUMNS])
        for r in rows:
            w.writerow([engine, *[r[c] for c in SUMMARY_COLUMNS]])


def _write_html(rows: list[dict], engine: str, path: str) -> None:
    """Self-contained HTML: inline CSS, no external refs, headline cards + table."""
    cards = []
    for r in rows:
        caveat = (
            " <span class='caveat'>INDEL: not hap.py-comparable (lower bound)</span>"
            if r["Type"] == VT_INDEL
            else ""
        )
        cards.append(
            "<div class='card'>"
            f"<h3>{html.escape(r['Type'])} / {html.escape(r['Filter'])}</h3>"
            f"<p>Recall {r['METRIC.Recall']:.4f} &middot; "
            f"Precision {r['METRIC.Precision']:.4f} &middot; "
            f"F1 {r['METRIC.F1_Score']:.4f}</p>"
            f"{caveat}</div>"
        )

    header_cells = "".join(f"<th>{html.escape(c)}</th>" for c in SUMMARY_COLUMNS)
    body_rows = []
    for r in rows:
        cells = []
        for c in SUMMARY_COLUMNS:
            v = r[c]
            cells.append(
                f"<td>{v:.4f}</td>" if isinstance(v, float) else f"<td>{v}</td>"
            )
        body_rows.append("<tr>" + "".join(cells) + "</tr>")

    doc = (
        "<main>"
        f"<h1>vcfclick benchmark <small>engine={html.escape(engine)}</small></h1>"
        "<p class='note'>INDEL: not hap.py-comparable; numbers are lower bounds.</p>"
        f"<section class='cards'>{''.join(cards)}</section>"
        f"<table><thead><tr>{header_cells}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody></table>"
        "</main>"
    )
    style = (
        "<style>"
        "body{font-family:system-ui,sans-serif;margin:2rem;color:#111;background:#fff}"
        "h1 small{font-weight:400;color:#666;font-size:.6em}"
        ".cards{display:flex;flex-wrap:wrap;gap:1rem}"
        ".card{border:1px solid #ccc;border-radius:8px;padding:1rem;min-width:12rem}"
        ".caveat,.note{color:#a15c00;font-size:.85em}"
        "table{border-collapse:collapse;margin-top:1.5rem;font-size:.85em}"
        "th,td{border:1px solid #ccc;padding:.3rem .5rem;text-align:right}"
        "th:first-child,td:first-child{text-align:left}"
        "@media(prefers-color-scheme:dark){body{background:#111;color:#eee}"
        "th,td,.card{border-color:#444}}"
        "</style>"
    )
    with open(path, "w") as fh:
        fh.write(style + doc)


def _write_parquet(classified, path: str) -> None:
    """Write the classified frame to parquet; caller only invokes when provided."""
    import pyarrow.parquet as pq

    pq.write_table(classified, path)


def write_reports(
    agg: dict,
    run_meta: dict,
    outdir: str,
    formats: list[str],
    classified=None,
) -> None:
    """Write the requested report formats into outdir from an aggregate dict."""
    os.makedirs(outdir, exist_ok=True)
    fmts = set(formats)
    rows = _summary_rows(agg)
    engine = run_meta.get("engine", "")

    if "csv" in fmts:
        _write_summary_csv(rows, os.path.join(outdir, "summary.csv"))
        _write_vcfclick_summary_csv(
            rows, engine, os.path.join(outdir, "vcfclick.summary.csv")
        )

    if "json" in fmts:
        with open(os.path.join(outdir, "run_meta.json"), "w") as fh:
            json.dump(run_meta, fh, indent=2)
        with open(os.path.join(outdir, "metrics.json"), "w") as fh:
            json.dump(rows, fh, indent=2)

    if "html" in fmts:
        _write_html(rows, engine, os.path.join(outdir, "index.html"))

    if "parquet" in fmts and classified is not None:
        _write_parquet(classified, os.path.join(outdir, "benchmark.parquet"))
