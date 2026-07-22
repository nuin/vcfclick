"""Orchestration for `vcfclick benchmark` — the services.py-style pure API the
CLI verb calls. Wires reference/normalize/regions/stratify/reconcile/aggregate/
report together for the P1 `normalized` engine. SQL-free; no click here.

Flow: open reference → parse+decompose+normalize both sides into NormRecords →
tag confident regions and variant types → classify per filter view → aggregate →
write reports. Returns run_meta + a per-(Type,Filter) summary.
"""

from __future__ import annotations

import dataclasses
import logging

from benchmark import regions as regions_mod
from benchmark import stratify as _strat
from benchmark.aggregate import aggregate_counts
from benchmark.constants import (
    FILTER_ALL,
    FILTER_PASS,
    VT_INDEL,
    VT_SNP,
)
from benchmark.metrics import metrics_from_counts
from benchmark.model import NormRecord
from benchmark.normalize import (
    canonical_gt,
    decompose_mnp,
    left_align,
    split_multiallelic,
    trim,
)
from benchmark.reconcile import classify, classify_haplotype
from benchmark.reference import BenchmarkError, Reference, canonical_contig
from benchmark.report import write_reports

log = logging.getLogger(__name__)

ALL_FORMATS = ["csv", "json", "parquet", "html"]
_MITO_KEYS = {"m", "mt", "chrm", "chrmt"}


def _is_mito(chrom: str) -> bool:
    return chrom.lower() in _MITO_KEYS


def _excluded_alt(alt: str | None) -> bool:
    """Symbolic (`<DEL>`), breakend (`N[chr:pos[`), spanning-`*`, or empty alt —
    unrepresentable by the normalized engine; the whole locus is routed out."""
    if alt is None or alt in ("", ".", "*"):
        return True
    return alt[0] == "<" or "[" in alt or "]" in alt


def _parse_side(
    path: str,
    side: str,
    ref: Reference,
    on_ref_mismatch: str,
    excluded: dict,
    decompose: bool,
    align: bool = True,
) -> list[NormRecord]:
    """Parse one VCF into normalized biallelic NormRecords.

    `align=True` left-aligns indels against the reference (the `normalized` and
    `haplotype` engines); `align=False` keeps the minimal-trim representation
    without left-shifting (the internal `exact` diagnostic engine).
    """
    from cyvcf2 import VCF

    known = ref.contigs
    out: list[NormRecord] = []
    for var in VCF(path):
        chrom = canonical_contig(var.CHROM, known)
        if chrom is None:
            raise BenchmarkError(f"contig {var.CHROM!r} is not in the reference")

        alts = list(var.ALT)
        if not alts or any(_excluded_alt(a) for a in alts):
            excluded["symbolic"] += 1
            continue

        ref_allele = var.REF.upper()
        actual = ref.fetch(chrom, var.POS - 1, var.POS - 1 + len(ref_allele))
        if actual != ref_allele:
            if on_ref_mismatch == "skip":
                excluded["ref_mismatch"] += 1
                continue
            raise BenchmarkError(
                f"{side} {chrom}:{var.POS} REF {ref_allele!r} != reference "
                f"{actual!r} (use --on-ref-mismatch skip to drop)"
            )

        gts = var.genotypes
        gt = tuple(gts[0][:-1]) if gts else (-1,)
        phased = bool(gts[0][-1]) if gts else False
        is_pass = var.FILTER is None or var.FILTER in ("PASS", ".")
        qual = float(var.QUAL) if var.QUAL is not None else 0.0

        for altrow in split_multiallelic(var.POS, ref_allele, alts, gt):
            if 1 not in altrow.gt:
                continue  # sample carries no copy of this alt (homref / no-call)
            pieces = (
                decompose_mnp(altrow.pos, altrow.ref, altrow.alt)
                if decompose
                else [(altrow.pos, altrow.ref, altrow.alt)]
            )
            alleles = altrow.gt
            if _is_mito(chrom) and len(alleles) == 1:
                alleles = alleles * 2  # haploid 1 ≡ diploid 1/1 on chrM
            for ppos, pref, palt in pieces:
                if align:
                    pos, r, a = left_align(ref.fetch, chrom, ppos, pref, palt)
                else:
                    pos, r, a = trim(ppos, pref, palt)
                out.append(
                    NormRecord(
                        side=side,
                        chrom=chrom,
                        pos=pos,
                        ref=r,
                        alt=a,
                        gt=canonical_gt(alleles, phased),
                        is_pass=is_pass,
                        other_alt_present=altrow.other_alt_present,
                        locus_id=altrow.locus_id,
                        qual=qual,
                    )
                )
    return out


def _summary(agg: dict) -> list[dict]:
    """One record per (Type, Filter) with integer counts + Python-computed metrics."""
    rows: list[dict] = []
    for vtype in (VT_SNP, VT_INDEL):
        for fv in (FILTER_PASS, FILTER_ALL):
            c = agg.get((fv, vtype), {})
            tp_t = c.get("truth_tp", 0)
            fn_t = c.get("truth_fn", 0)
            tp_q = c.get("query_tp", 0)
            fp_q = c.get("query_fp", 0)
            unk_q = c.get("query_unk", 0)
            m = metrics_from_counts(tp_t, fn_t, tp_q, fp_q, unk_q)
            rows.append(
                {
                    "Type": vtype,
                    "Filter": fv,
                    "truth_tp": tp_t,
                    "truth_fn": fn_t,
                    "truth_total": tp_t + fn_t,
                    "query_tp": tp_q,
                    "query_fp": fp_q,
                    "query_unk": unk_q,
                    "query_total": tp_q + fp_q + unk_q,
                    "recall": m.recall,
                    "precision": m.precision,
                    "f1": m.f1,
                    "frac_na": m.frac_na,
                }
            )
    return rows


def classified_rows(
    truth: str,
    query: str,
    ref: str,
    regions: str | None = None,
    engine: str = "normalized",
    on_ref_mismatch: str = "error",
    strict: bool = False,
    conf_containment: str = "start",
    decompose_mnp: bool = False,
    excluded: dict | None = None,
) -> list:
    """Parse, normalize, region-tag, and classify both sides into ClassifiedRows.

    The reusable core of a benchmark run (no aggregation or reporting) — the
    multi-caller cohort layer builds its combined frame from this.
    """
    if engine not in ("normalized", "haplotype", "exact"):
        raise BenchmarkError(f"unknown engine {engine!r}")
    reference = Reference(ref)
    if excluded is None:
        excluded = {"symbolic": 0, "ref_mismatch": 0}
    align = engine != "exact"  # exact keys on trimmed reps without left-shifting

    truth_recs = _parse_side(
        truth, "truth", reference, on_ref_mismatch, excluded, decompose_mnp, align
    )
    query_recs = _parse_side(
        query, "query", reference, on_ref_mismatch, excluded, decompose_mnp, align
    )

    if regions is None:
        if strict:
            raise BenchmarkError("no --regions BED given and --strict is set")
        log.warning("no --regions BED given; treating every call as confident")
        truth_recs = [dataclasses.replace(r, in_conf=True) for r in truth_recs]
        query_recs = [dataclasses.replace(r, in_conf=True) for r in query_recs]
    else:
        conf = regions_mod.ConfRegions.from_bed(regions)
        truth_recs = conf.tag(truth_recs, containment=conf_containment)
        query_recs = conf.tag(query_recs, containment=conf_containment)

    truth_recs = _strat.tag(truth_recs)
    query_recs = _strat.tag(query_recs)

    if engine == "haplotype":
        rows = classify_haplotype(truth_recs, query_recs, FILTER_ALL, reference.fetch)
        rows += classify_haplotype(truth_recs, query_recs, FILTER_PASS, reference.fetch)
    else:
        rows = classify(truth_recs, query_recs, FILTER_ALL)
        rows += classify(truth_recs, query_recs, FILTER_PASS)
    return rows


def run_benchmark(
    truth: str,
    query: str,
    ref: str,
    outdir: str,
    regions: str | None = None,
    engine: str = "normalized",
    report_formats: list[str] | None = None,
    on_ref_mismatch: str = "error",
    strict: bool = False,
    conf_containment: str = "start",
    decompose_mnp: bool = False,
    stratify: list[str] | None = None,
    roc: bool = False,
) -> dict:
    """Benchmark `query` against `truth` over reference `ref`, writing reports.

    `engine="haplotype"` runs the P2 local-haplotype reconciliation; `engine=
    "exact"` is an internal diagnostic that keys on the trimmed (not left-aligned)
    representation. With `regions=None` every call is treated as confident (a
    warning, or a hard error under `strict`). `conf_containment` is "start" or
    "full"; `decompose_mnp` atomizes MNPs into SNPs (off by default). Returns
    `{"run_meta", "summary"}`.
    """
    formats = list(report_formats) if report_formats else list(ALL_FORMATS)
    excluded = {"symbolic": 0, "ref_mismatch": 0}
    rows = classified_rows(
        truth,
        query,
        ref,
        regions=regions,
        engine=engine,
        on_ref_mismatch=on_ref_mismatch,
        strict=strict,
        conf_containment=conf_containment,
        decompose_mnp=decompose_mnp,
        excluded=excluded,
    )

    agg = aggregate_counts(rows)
    summary = _summary(agg)

    # summary.csv is strict SNP/INDEL (hap.py shape); surface any other type
    # (COMPLEX/MNP, NOCALL) here so it is never silently dropped.
    unsummarized: dict[str, int] = {}
    for r in rows:
        if r.filter_view == FILTER_ALL and r.vtype not in (VT_SNP, VT_INDEL):
            unsummarized[r.vtype] = unsummarized.get(r.vtype, 0) + 1

    run_meta = {
        "engine": engine,
        "truth": truth,
        "query": query,
        "ref": ref,
        "regions": regions,
        "outdir": outdir,
        "on_ref_mismatch": on_ref_mismatch,
        "report_formats": formats,
        "excluded": excluded,
        "unsummarized_types": unsummarized,
    }

    classified = None
    if "parquet" in formats or stratify:
        import pyarrow as pa

        classified = pa.Table.from_pylist([dataclasses.asdict(r) for r in rows])
    write_reports(agg, run_meta, outdir, formats, classified=classified)

    if roc:
        import os

        from benchmark.roc import write_roc_tsv

        write_roc_tsv(rows, [VT_SNP, VT_INDEL], os.path.join(outdir, "roc.tsv"))
        run_meta["roc"] = "roc.tsv"

    if stratify:
        from annotations.db import _store_path
        from benchmark.stratify_db import write_stratified

        ann_path = str(_store_path())
        try:
            written = write_stratified(classified, ann_path, list(stratify), outdir)
            run_meta["stratified"] = [p.rsplit("/", 1)[-1] for p in written]
        except (
            Exception
        ) as exc:  # annotation store missing/unreadable — surface, don't crash
            log.warning("stratification skipped: %s", exc)
            run_meta["stratified"] = []

    return {"run_meta": run_meta, "summary": summary}


def exact_vs_normalized_delta(
    truth: str,
    query: str,
    ref: str,
    regions: str | None = None,
) -> dict:
    """Per-Type `normalized`-minus-`exact` metric deltas (the offline gap signal).

    Runs both engines on the same inputs and returns, for each Type in
    {SNP, INDEL}, the recall/precision/f1 improvement that reference-based
    left-alignment alone buys (ALL filter view). All-positive on a fixture whose
    only defect is representation shift. Self-contained (no external deps); writes
    reports to throwaway temp dirs.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as d_ex, tempfile.TemporaryDirectory() as d_nm:
        ex = run_benchmark(
            truth,
            query,
            ref,
            d_ex,
            regions=regions,
            engine="exact",
            report_formats=["csv"],
        )
        nm = run_benchmark(
            truth,
            query,
            ref,
            d_nm,
            regions=regions,
            engine="normalized",
            report_formats=["csv"],
        )

    def _by_type(summary: list[dict]) -> dict:
        return {r["Type"]: r for r in summary if r["Filter"] == FILTER_ALL}

    ex_t = _by_type(ex["summary"])
    nm_t = _by_type(nm["summary"])
    delta: dict = {}
    for vt in (VT_SNP, VT_INDEL):
        e, n = ex_t[vt], nm_t[vt]
        delta[vt] = {
            "recall": n["recall"] - e["recall"],
            "precision": n["precision"] - e["precision"],
            "f1": n["f1"] - e["f1"],
        }
    return delta
