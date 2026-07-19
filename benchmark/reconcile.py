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

from benchmark.constants import (
    BD_FN,
    BD_FP,
    BD_N,
    BD_TP,
    BK_AM,
    BK_GM,
    BK_NONE,
    FILTER_ALL,
)
from benchmark.model import ClassifiedRow, NormRecord


class UnsupportedFeatureError(NotImplementedError):
    """A P2-only reconciliation feature was requested from the P1 engine."""


Key = tuple  # (chrom, pos, ref, alt)


def _key(r: NormRecord) -> Key:
    return (r.chrom, r.pos, r.ref, r.alt)


def _gt_key(r: NormRecord) -> tuple:
    """Locus-level genotype identity: the remapped GT plus whether a sibling
    alt was carried at the locus (so `1/2`≠`0/1`)."""
    return (r.gt, r.other_alt_present)


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


def classify_haplotype(
    truth: list[NormRecord], query: list[NormRecord], filter_view: str
) -> list[ClassifiedRow]:
    """P2 haplotype (`lm`) reconciliation — not implemented in P1."""
    raise UnsupportedFeatureError(
        "haplotype reconciliation (BK=lm) is a P2 feature; use the normalized engine"
    )
