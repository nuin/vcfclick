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
# Sites per chunk in king_robust: bounds the float32 indicator matrices.
_CHUNK_SITES = 10_000


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


def king_robust(
    g: np.ndarray, chunk_sites: int = _CHUNK_SITES
) -> dict[str, np.ndarray]:
    """Pairwise KING-robust statistics.

    `g` is a (sites x samples) int8 matrix of alt-allele counts: 0, 1, 2,
    or -1 for a no-call. Returns samples x samples matrices: kinship, ibs0
    (count), hethet (count) and n_sites (sites where both are called).

    The 0/1 indicator matrices are built `chunk_sites` rows at a time and
    their products summed into samples x samples accumulators, so peak
    memory is about chunk_sites x samples x 16 bytes plus the outputs.
    Each product is an exact integer count in float32 (below 2^24 sites,
    the --max-sites cap), so the sum equals the single-pass result.
    """
    n = g.shape[1]
    hethet = np.zeros((n, n), dtype=np.float32)
    ibs0 = np.zeros((n, n), dtype=np.float32)
    het_given_called = np.zeros((n, n), dtype=np.float32)
    both_called = np.zeros((n, n), dtype=np.float32)
    step = max(1, int(chunk_sites))
    for lo in range(0, g.shape[0], step):
        c = g[lo : lo + step]
        called = (c >= 0).astype(np.float32)
        het = (c == 1).astype(np.float32)
        hethet += het.T @ het
        # het count of sample i over the sites where sample j is also called
        het_given_called += het.T @ called
        both_called += called.T @ called
        del called
        homref = (c == 0).astype(np.float32)
        homalt = (c == 2).astype(np.float32)
        ibs0 += homref.T @ homalt
        del het, homref, homalt
    ibs0 = ibs0 + ibs0.T
    denom = het_given_called + het_given_called.T
    with np.errstate(divide="ignore", invalid="ignore"):
        kinship = np.where(denom > 0, (hethet - 2 * ibs0) / denom, np.nan)
    return {
        "het_pair": denom.astype(np.int64),
        "kinship": kinship,
        "ibs0": ibs0.astype(np.int64),
        "hethet": hethet.astype(np.int64),
        "n_sites": both_called.astype(np.int64),
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


def _arrow(sess, sql: str):
    from storage.db import query_arrow

    return query_arrow(sql, sess)


def _thin_keeps(k: int, total: int, max_sites: int) -> bool:
    """True when the site of 0-based rank k is one of the max_sites kept,
    the ranks floor(i * total / max_sites). Mirrors `_thin_sql`."""
    if total <= max_sites:
        return True
    r = (k * max_sites) % total
    return r == 0 or r > total - max_sites


def _thin_sql(col: str, total: int, max_sites: int) -> str:
    """SQL form of `_thin_keeps` on the rank column `col`. Integer
    arithmetic only, so both backends keep the same sites."""
    if total <= max_sites:
        return "1 = 1"
    m, t = int(max_sites), int(total)
    return f"(({col} * {m}) % {t} = 0 OR ({col} * {m}) % {t} > {t - m})"


def _kept_sites_cte(where: str, total: int, max_sites: int) -> str:
    """`WITH ... kept AS (chrom, pos, ref, alt, i)`: the thinned SNV sites,
    with i their 0-based row in the genotype matrix."""
    return (
        "WITH ranked AS (SELECT chrom, pos, ref, alt, "
        "row_number() OVER (ORDER BY chrom, pos, ref, alt) - 1 AS k "
        f"FROM variants WHERE {where}), "
        "kept AS (SELECT chrom, pos, ref, alt, "
        "row_number() OVER (ORDER BY k) - 1 AS i "
        f"FROM ranked WHERE {_thin_sql('k', total, max_sites)}) "
    )


def _kept_join(alias: str) -> str:
    return " AND ".join(
        f"{alias}.{c} = kept.{c}" for c in ("chrom", "pos", "ref", "alt")
    )


def _fill(g: np.ndarray, table, samples: list[str], value: int | None) -> None:
    """Write `table`'s rows (i, sample_id[, gt]) into g: gt, or `value`."""
    import pyarrow as pa
    import pyarrow.compute as pc

    if table.num_rows == 0:
        return
    rows = np.asarray(table.column("i").to_numpy(), dtype=np.int64)
    ids = pc.cast(table.column("sample_id"), pa.string())
    cols = pc.fill_null(pc.index_in(ids, value_set=pa.array(samples, pa.string())), -1)
    cols = np.asarray(cols.to_numpy(), dtype=np.int64)
    if value is None:
        vals = np.asarray(table.column("gt").to_numpy(), dtype=np.int64)
        ok = (cols >= 0) & (vals >= 0)
    else:
        vals = np.full(len(rows), value, dtype=np.int64)
        ok = cols >= 0
    g[rows[ok], cols[ok]] = vals[ok].astype(np.int8)


def _genotype_matrix(
    sess, q: str, cte: str, n_sites: int, samples: list[str], keep_ref: bool
) -> np.ndarray:
    """(sites x samples) int8 alt-allele counts over the kept sites, -1 for a
    no-call. Only the kept sites' rows leave the database."""
    g = np.full((n_sites, len(samples)), -1 if keep_ref else 0, dtype=np.int8)
    min_gt = 0 if keep_ref else 1
    _fill(
        g,
        _arrow(
            sess,
            f"{cte}SELECT kept.i AS i, g.sample_id AS sample_id, g.gt AS gt "
            f"FROM genotypes g INNER JOIN kept ON {_kept_join('g')} "
            f"WHERE g.ingest_id = {q} AND g.gt >= {min_gt}",
        ),
        samples,
        None,
    )
    if not keep_ref:
        _mark_no_calls(sess, q, cte, g, samples)
    return g


def _mark_no_calls(sess, q: str, cte: str, g: np.ndarray, samples: list[str]) -> None:
    """Sparse mode reads an absent genotype as 0/0; the full no-calls an
    ingest recorded in missing_genotypes are set to -1 so they are skipped
    for every pair. An ingest made with --no-record-missing has none."""
    from storage import table_exists

    if not table_exists("missing_genotypes"):
        return
    _fill(
        g,
        _arrow(
            sess,
            f"{cte}SELECT kept.i AS i, m.sample_id AS sample_id "
            f"FROM missing_genotypes m INNER JOIN kept ON {_kept_join('m')} "
            f"WHERE m.ingest_id = {q}",
        ),
        samples,
        -1,
    )


def _load(
    sess, ingest_id: str, max_sites: int
) -> tuple[list[str], np.ndarray, str, int, int]:
    """Samples, (sites x samples) matrix, mode, total SNV sites, genomic span.

    Refuses a cohort above _MAX_SAMPLES before anything is allocated.
    """
    from cli.db_diff import _quote_str
    from storage.sql import to_int64

    q = _quote_str(ingest_id)
    samples = [
        r[0]
        for r in _rows(
            sess,
            f"SELECT DISTINCT sample_id FROM samples WHERE ingest_id = {q} ORDER BY sample_id",
        )
    ]
    if len(samples) > _MAX_SAMPLES:
        raise click.ClickException(
            f"{len(samples):,} samples in {ingest_id!r}; pairwise kinship above {_MAX_SAMPLES:,} "
            "samples needs more memory than this command uses."
        )
    keep_ref = bool(
        _rows(
            sess,
            f"SELECT 1 FROM genotypes WHERE ingest_id = {q} AND gt = 0 LIMIT 1",
        )
    )
    where = f"ingest_id = {q} AND length(ref) = 1 AND length(alt) = 1"
    total = int(_rows(sess, f"SELECT count(*) FROM variants WHERE {where}")[0][0])
    span = int(
        _rows(
            sess,
            f"SELECT {to_int64('coalesce(sum(s), 0)')} FROM ("
            f"SELECT max(pos) - min(pos) AS s FROM variants WHERE {where} "
            "GROUP BY chrom) t",
        )[0][0]
    )
    n_sites = min(total, max_sites)
    cte = _kept_sites_cte(where, total, max_sites)
    g = _genotype_matrix(sess, q, cte, n_sites, samples, keep_ref)
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
