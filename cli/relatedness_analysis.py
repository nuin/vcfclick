"""KING-robust kinship calculation and pedigree checks."""

import json
from functools import partial

import click
import numpy as np

# KING kinship thresholds (powers of 2 between the expected values).
_DUPLICATE = 0.354
_FIRST = 0.177
_SECOND = 0.0884
_THIRD = 0.0442
# IBS0 proportion below which a first-degree pair reads as parent/child.
_PO_MAX_IBS0 = 0.005
_MAX_SAMPLES = 20000
# Kinship needs many independent markers. Within one gene or a small region,
# markers are inherited together (linkage), so unrelated people who share a
# common haplotype look identical. Below these, pairs are not classified.
_MIN_SPAN_BP = 10_000_000
_MIN_SITES = 1000
# A pair needs this many heterozygous sites between them to be estimable.
_MIN_PAIR_HETS = 20


def _pair(stats, samples, classify_pairs: bool, i: int, j: int) -> dict:
    n = int(stats["n_sites"][i, j])
    k = float(stats["kinship"][i, j])
    frac = stats["ibs0"][i, j] / n if n else float("nan")
    hets = int(stats["het_pair"][i, j])
    if hets < _MIN_PAIR_HETS:
        k = float("nan")
    return {
        "sample_a": samples[i],
        "sample_b": samples[j],
        "kinship": None if np.isnan(k) else round(k, 4),
        "ibs0": None if np.isnan(frac) else round(float(frac), 4),
        "n_sites": n,
        "relationship": classify(k, frac) if classify_pairs else "insufficient-data",
    }


def _related_pairs(stats, samples, pair, include_all, classify_pairs, min_kinship):
    pairs = []
    if stats is not None:
        n = len(samples)
        if include_all:
            candidates = [(i, j) for i in range(n) for j in range(i + 1, n)]
        elif classify_pairs:
            # Select in numpy first: building a record for every one of n^2/2
            # pairs just to drop almost all of them is the slow part.
            keep = (stats["kinship"] >= min_kinship) & (
                stats["het_pair"] >= _MIN_PAIR_HETS
            )
            candidates = [tuple(ij) for ij in np.argwhere(np.triu(keep, 1))]
        else:
            candidates = []
        for i, j in candidates:
            p = pair(int(i), int(j))
            related = p["relationship"] not in ("unrelated", "insufficient-data")
            if include_all or related:
                pairs.append(p)
        pairs.sort(key=lambda p: -(p["kinship"] if p["kinship"] is not None else -9))

    return pairs


def _pedigree_checks(pedigree, samples, stats, pair):
    idx = {s: i for i, s in enumerate(samples)}
    checks = []
    for child, parent, role in pedigree:
        if stats is None or child not in idx or parent not in idx:
            checks.append(
                {
                    "child": child,
                    "parent": parent,
                    "role": role,
                    "kinship": None,
                    "ibs0": None,
                    "observed": None,
                    "verdict": "not genotyped",
                }
            )
            continue
        p = pair(idx[child], idx[parent])
        verdict = (
            "insufficient data"
            if p["relationship"] == "insufficient-data"
            else "ok"
            if p["relationship"] == "parent-child"
            else "mismatch"
        )
        checks.append(
            {
                "child": child,
                "parent": parent,
                "role": role,
                "kinship": p["kinship"],
                "ibs0": p["ibs0"],
                "observed": p["relationship"],
                "verdict": verdict,
            }
        )

    return checks


def king_robust(g: np.ndarray) -> dict[str, np.ndarray]:
    """Pairwise KING-robust statistics.

    `g` is a (sites x samples) int8 matrix of alt-allele counts: 0, 1, 2,
    or -1 for a no-call. Returns samples x samples matrices: kinship, ibs0
    (count), hethet (count) and n_sites (sites where both are called).
    """
    called = (g >= 0).astype(np.float32)
    het = (g == 1).astype(np.float32)
    homref = (g == 0).astype(np.float32)
    homalt = (g == 2).astype(np.float32)

    hethet = het.T @ het
    ibs0 = homref.T @ homalt
    ibs0 = ibs0 + ibs0.T
    # het count of sample i over the sites where sample j is also called
    het_given_called = het.T @ called
    denom = het_given_called + het_given_called.T
    with np.errstate(divide="ignore", invalid="ignore"):
        kinship = np.where(denom > 0, (hethet - 2 * ibs0) / denom, np.nan)
    return {
        "het_pair": denom.astype(np.int64),
        "kinship": kinship,
        "ibs0": ibs0.astype(np.int64),
        "hethet": hethet.astype(np.int64),
        "n_sites": (called.T @ called).astype(np.int64),
    }


def classify(kinship: float, ibs0_frac: float) -> str:
    if np.isnan(kinship):
        return "insufficient-data"
    if kinship > _DUPLICATE:
        return "duplicate"
    if kinship > _FIRST:
        return "parent-child" if ibs0_frac < _PO_MAX_IBS0 else "full-siblings"
    if kinship > _SECOND:
        return "second-degree"
    if kinship > _THIRD:
        return "third-degree"
    return "unrelated"


def _rows(sess, sql: str) -> list[list]:
    return json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]


def _span_bp(sites: list[list]) -> int:
    """Total genomic span covered by the sites: sum over chromosomes of max - min."""
    by_chrom: dict[str, list[int]] = {}
    for chrom, pos, *_ in sites:
        lo_hi = by_chrom.setdefault(chrom, [int(pos), int(pos)])
        lo_hi[0] = min(lo_hi[0], int(pos))
        lo_hi[1] = max(lo_hi[1], int(pos))
    return sum(hi - lo for lo, hi in by_chrom.values())


def _load(
    sess, ingest_id: str, max_sites: int
) -> tuple[list[str], np.ndarray, str, int, int]:
    """Samples, (sites x samples) matrix, mode, total SNV sites, genomic span."""
    from cli.db_diff import _quote_str

    q = _quote_str(ingest_id)
    samples = [
        r[0]
        for r in _rows(
            sess,
            f"SELECT DISTINCT sample_id FROM samples WHERE ingest_id = {q} ORDER BY sample_id",
        )
    ]
    keep_ref = (
        int(
            _rows(
                sess, f"SELECT count(*) FROM genotypes WHERE ingest_id = {q} AND gt = 0"
            )[0][0]
        )
        > 0
    )
    snv = "length(ref) = 1 AND length(alt) = 1"
    sites = _rows(
        sess,
        f"SELECT chrom, pos, ref, alt FROM variants WHERE ingest_id = {q} AND {snv} ORDER BY chrom, pos, ref, alt",
    )
    total = len(sites)
    span = _span_bp(sites)
    if total > max_sites:
        # Deterministic, evenly spread thinning keeps memory bounded.
        step = total / max_sites
        sites = [sites[int(i * step)] for i in range(max_sites)]
    site_index = {tuple(s): i for i, s in enumerate(sites)}
    sample_index = {s: j for j, s in enumerate(samples)}

    g = np.full((len(sites), len(samples)), -1 if keep_ref else 0, dtype=np.int8)
    min_gt = 0 if keep_ref else 1
    for chrom, pos, ref, alt, sample_id, gt in _rows(
        sess,
        f"SELECT chrom, pos, ref, alt, sample_id, gt FROM genotypes "
        f"WHERE ingest_id = {q} AND {snv} AND gt >= {min_gt}",
    ):
        i = site_index.get((chrom, pos, ref, alt))
        j = sample_index.get(sample_id)
        if i is not None and j is not None and int(gt) >= 0:
            g[i, j] = int(gt)
    return samples, g, "keep-reference" if keep_ref else "sparse", total, span


def _pedigree(sess, ingest_id: str) -> list[tuple[str, str, str]]:
    """(child, parent, role) for every declared parent in the pedigree."""
    from cli.db_diff import _quote_str
    from storage import table_exists

    if not table_exists("pedigree"):
        return []
    rows = _rows(
        sess,
        f"SELECT sample_id, father_id, mother_id FROM pedigree WHERE ingest_id = {_quote_str(ingest_id)}",
    )
    out = []
    for child, father, mother in rows:
        for parent, role in ((father, "father"), (mother, "mother")):
            if parent and parent != "0":
                out.append((child, parent, role))
    return out


def relatedness(
    sess,
    ingest_id: str,
    *,
    max_sites: int,
    include_all: bool,
    min_kinship: float,
    force: bool = False,
) -> dict:
    samples, g, mode, total, span = _load(sess, ingest_id, max_sites)
    warning = None
    if span < _MIN_SPAN_BP or total < _MIN_SITES:
        warning = (
            f"only {total:,} SNVs spanning {span / 1e6:.2f} Mb; kinship needs genome-wide markers "
            f"(at least {_MIN_SITES:,} SNVs over {_MIN_SPAN_BP // 1_000_000} Mb). Within a small region, "
            "unrelated people sharing a common haplotype look related"
            + (
                ". Classified anyway (--force); treat as unreliable."
                if force
                else ", so pairs are not classified."
            )
        )
    classify_pairs = warning is None or force
    if len(samples) > _MAX_SAMPLES:
        raise click.ClickException(
            f"{len(samples):,} samples in {ingest_id!r}; pairwise kinship above {_MAX_SAMPLES:,} "
            "samples needs more memory than this command uses."
        )
    stats = king_robust(g) if len(samples) and g.shape[0] else None

    pair = partial(_pair, stats, samples, classify_pairs)
    pairs = _related_pairs(
        stats, samples, pair, include_all, classify_pairs, min_kinship
    )

    checks = _pedigree_checks(_pedigree(sess, ingest_id), samples, stats, pair)

    return {
        "ingest_id": ingest_id,
        "mode": mode,
        "n_samples": len(samples),
        "n_sites": total,
        "sites_used": int(g.shape[0]),
        "span_bp": span,
        "warning": warning,
        "pairs": pairs,
        "pedigree_checks": checks,
    }
