"""P1 classification for `vcfclick benchmark` (the `normalized` engine crux).

Full-outer join of truth and query on the normalized key `(chrom,pos,ref,alt)`,
run independently *per filter view*. In the PASS view only the *query* is
filtered by FILTER (truth is always scored), so a query call failing FILTER is
absent there and its matched truth call becomes an FN — never a post-hoc tag.

The genotype verdict compares the reconstructed per-locus allele multiset via
`(gt, other_alt_present)`, so a het-alt `1/2` split row is never matched to a
plain `0/1`. Duplicate normalized keys within one side are routed to a counted
UNK bucket (`bd=N`) instead of being silently merged.

P1 is SQL-free and backend-neutral (pure Python).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import replace as _replace
from typing import Callable

from benchmark import haplotype as _hap
from benchmark.constants import (
    BD_FN,
    BD_FP,
    BD_N,
    BD_TP,
    BK_AM,
    BK_GM,
    BK_LM,
    BK_NONE,
    FILTER_ALL,
)
from benchmark.model import ClassifiedRow, NormRecord

Fetch = Callable[[str, int, int], str]  # (chrom, start0, end0) -> uppercased bases


class UnsupportedFeatureError(NotImplementedError):
    """A P2-only reconciliation feature was requested from the P1 engine."""


Key = tuple  # (chrom, pos, ref, alt)


def _key(r: NormRecord) -> Key:
    return (r.chrom, r.pos, r.ref, r.alt)


def _gt_key(r: NormRecord) -> tuple:
    """Locus-level genotype identity for the keyed match: the phase-insensitive
    allele multiset plus whether a sibling alt was carried (so `1/2`≠`0/1`).

    Genotype concordance ignores phase — a phased `1|0` and an unphased `0/1`
    are the same heterozygous genotype (assembly callers emit phased GTs; the
    haplotype engine still reads phase from `r.gt`, which is left intact).
    """
    alleles = tuple(sorted(r.gt.replace("|", "/").split("/")))
    return (alleles, r.other_alt_present)


def _present(records: list[NormRecord], filter_view: str) -> list[NormRecord]:
    if filter_view == FILTER_ALL:
        return list(records)
    return [r for r in records if r.is_pass]


def _partition(
    records: list[NormRecord],
) -> tuple[dict[Key, NormRecord], list[NormRecord]]:
    """Split into a unique-key map and the records whose key collides on a side."""
    buckets: dict[Key, list[NormRecord]] = defaultdict(list)
    for r in records:
        buckets[_key(r)].append(r)
    unique: dict[Key, NormRecord] = {}
    dups: list[NormRecord] = []
    for k, rs in buckets.items():
        if len(rs) == 1:
            unique[k] = rs[0]
        else:
            dups.extend(rs)
    return unique, dups


def _row(r: NormRecord, filter_view: str, bd: str, bk: str) -> ClassifiedRow:
    return ClassifiedRow(
        side=r.side,
        filter_view=filter_view,
        chrom=r.chrom,
        pos=r.pos,
        ref=r.ref,
        alt=r.alt,
        vtype=r.vtype,
        subtype=r.subtype,
        in_conf=r.in_conf,
        bd=bd,
        bk=bk,
        blt=r.blt,
    )


def classify(
    truth: list[NormRecord], query: list[NormRecord], filter_view: str
) -> list[ClassifiedRow]:
    """Classify one filter view into per-side `ClassifiedRow`s.

    Verdicts (see design doc truth table): both present + equal locus GT → TP/TP
    (`gm`); both present + GT differ → FN/FP (`am`, the het↔hom double penalty);
    truth-only in_conf → FN; query-only in_conf → FP; query-only out-of-conf → N
    (UNK); truth-only out-of-conf → dropped. Duplicate-key loci → N.

    `FP.al` (a query FP position-colliding with a differently-spelled truth FN)
    carries `bk=.` here and is derived downstream from these rows.
    """
    # hap.py filters only the query by FILTER; truth is the gold standard and is
    # always scored (a non-PASS truth call matched by a PASS query is still a TP).
    truth_p = list(truth)
    query_p = _present(query, filter_view)

    truth_map, truth_dups = _partition(truth_p)
    query_map, query_dups = _partition(query_p)

    rows: list[ClassifiedRow] = []

    # Key-uniqueness guard: collided keys never enter the join.
    for r in truth_dups + query_dups:
        rows.append(_row(r, filter_view, BD_N, BK_NONE))
    assert not ({_key(r) for r in truth_dups} & truth_map.keys())
    assert not ({_key(r) for r in query_dups} & query_map.keys())

    for key in sorted(set(truth_map) | set(query_map)):
        t = truth_map.get(key)
        q = query_map.get(key)
        if t is not None and q is not None:
            if _gt_key(t) == _gt_key(q):
                rows.append(_row(t, filter_view, BD_TP, BK_GM))
                rows.append(_row(q, filter_view, BD_TP, BK_GM))
            else:
                rows.append(_row(t, filter_view, BD_FN, BK_AM))
                rows.append(_row(q, filter_view, BD_FP, BK_AM))
        elif t is not None:
            if t.in_conf:
                rows.append(_row(t, filter_view, BD_FN, BK_NONE))
            # truth-only out-of-conf: dropped (never an FN)
        else:  # q is not None
            if q.in_conf:
                rows.append(_row(q, filter_view, BD_FP, BK_NONE))
            else:
                rows.append(_row(q, filter_view, BD_N, BK_NONE))
    return rows


# Design-doc alias (pipeline calls it `classify_keyed`).
classify_keyed = classify


def _ident(r) -> tuple:
    """Row/record identity used to splice reclassified verdicts back in."""
    return (r.side, r.chrom, r.pos, r.ref, r.alt)


def _full_diploid(r: NormRecord) -> bool:
    """True only for a fully-called diploid genotype (two present alleles).

    The haplotype rescue collapses a record's alt into a single applied sequence;
    that is faithful only for a genuine diploid genotype. A half-call (`./1`) or
    haploid (`1`) would otherwise be treated as hom-alt and falsely rescued
    against a real `1/1`, so such records are kept at their P1 verdict.
    """
    alleles = r.gt.replace("|", "/").split("/")
    return len(alleles) == 2 and "." not in alleles


def classify_haplotype(
    truth: list[NormRecord],
    query: list[NormRecord],
    filter_view: str,
    fetch: Fetch,
    flank: int = 20,
    max_cluster: int = 16,
) -> list[ClassifiedRow]:
    """P2 local-haplotype (`lm`) reconciliation over the P1 keyed residual.

    Runs the P1 keyed match, then rescues representation-different but
    sequence-equivalent calls that P1 miscounts as FP+FN: the in-conf unmatched
    truth (FN) and query (FP) records are clustered by `haplotype.cluster` and,
    per candidate cluster (both sides present), replayed onto the reference. A
    cluster whose two sides yield a common (di)plotype is reclassified TP/`lm`.
    Over-budget clusters (> `max_cluster` members) are left at their P1 FP/FN
    verdict (never rescued, never routed to `BD_N` — that would drop the truth FN
    from the recall denominator). All other rows pass through unchanged.
    """
    rows = classify(truth, query, filter_view)

    # Index the norm records feeding the residual (residual keys are unique per
    # side — duplicate keys were already routed to BD_N by classify).
    t_index: dict[Key, NormRecord] = {}
    for r in truth:
        t_index.setdefault(_key(r), r)
    q_index: dict[Key, NormRecord] = {}
    for r in _present(query, filter_view):
        q_index.setdefault(_key(r), r)

    fn_recs: list[NormRecord] = []
    fp_recs: list[NormRecord] = []
    for row in rows:
        if row.bk not in (BK_NONE, BK_AM) or not row.in_conf:
            continue  # unmatched (.) OR allele-match/GT-differ (am) residual is
            # eligible; the diplotype check gates it, so het-alt representation
            # differences rescue while genuine zygosity errors stay am.
        k = (row.chrom, row.pos, row.ref, row.alt)
        if row.side == "truth" and row.bd == BD_FN:
            rec = t_index.get(k)
            if rec is not None and _full_diploid(rec):
                fn_recs.append(rec)
        elif row.side == "query" and row.bd == BD_FP:
            rec = q_index.get(k)
            if rec is not None and _full_diploid(rec):
                fp_recs.append(rec)

    reclassified: set[tuple] = set()
    for clust in _hap.cluster(fn_recs + fp_recs, flank):
        c_truth = [r for r in clust if r.side == "truth"]
        c_query = [r for r in clust if r.side == "query"]
        if not c_truth or not c_query:
            continue  # single-sided residual: no counterpart, stays FP/FN
        if len(clust) > max_cluster:
            # Over budget: not rescued. Keep the P1 FP/FN verdict — never route to
            # BD_N, which would drop the truth FN from the recall denominator and
            # let a messy over-complex query outscore P1 (codex #8).
            continue
        if _hap.haplotype_equivalent(fetch, c_truth, c_query, flank):
            reclassified.update(_ident(r) for r in clust)

    if not reclassified:
        return rows

    out: list[ClassifiedRow] = []
    for row in rows:
        if _ident(row) in reclassified:
            out.append(_replace(row, bd=BD_TP, bk=BK_LM))
        else:
            out.append(row)
    return out
