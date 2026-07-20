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
from benchmark import stratify
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
)
from benchmark.reconcile import UnsupportedFeatureError, classify
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
) -> list[NormRecord]:
    """Parse one VCF into normalized biallelic NormRecords."""
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
                pos, r, a = left_align(ref.fetch, chrom, ppos, pref, palt)
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
) -> dict:
    """Benchmark `query` against `truth` over reference `ref`, writing reports.

    `engine="haplotype"` is a P2 feature and raises UnsupportedFeatureError.
    With `regions=None` every call is treated as confident (a warning, or a hard
    error under `strict`). `conf_containment` is "start" or "full"; `decompose_mnp`
    atomizes MNPs into SNPs (off by default). Returns `{"run_meta", "summary"}`.
    """
    if engine == "haplotype":
        raise UnsupportedFeatureError(
            "haplotype reconciliation is a P2 feature; use --engine normalized"
        )
    if engine != "normalized":
        raise BenchmarkError(f"unknown engine {engine!r}")

    formats = list(report_formats) if report_formats else list(ALL_FORMATS)
    reference = Reference(ref)
    excluded = {"symbolic": 0, "ref_mismatch": 0}

    truth_recs = _parse_side(
        truth, "truth", reference, on_ref_mismatch, excluded, decompose_mnp
    )
    query_recs = _parse_side(
        query, "query", reference, on_ref_mismatch, excluded, decompose_mnp
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

    truth_recs = stratify.tag(truth_recs)
    query_recs = stratify.tag(query_recs)

    rows = classify(truth_recs, query_recs, FILTER_ALL)
    rows += classify(truth_recs, query_recs, FILTER_PASS)

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
    if "parquet" in formats:
        import pyarrow as pa

        classified = pa.Table.from_pylist([dataclasses.asdict(r) for r in rows])
    write_reports(agg, run_meta, outdir, formats, classified=classified)
    return {"run_meta": run_meta, "summary": summary}
