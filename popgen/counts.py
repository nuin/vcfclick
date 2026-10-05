"""Per-site, per-group allele counts from a vcfclick database.

The SQL here only aggregates; every statistic is computed in Python from
the counts it returns. All statements are plain SQL that both backends
accept (no `countIf`, no `FINAL`, `count(*)` rather than `count()`,
aggregates cast to BIGINT), and results come back as Arrow, so chDB and
DuckDB return identical numbers without any text parsing.

What is counted, per site (chrom, pos, ref, alt) of one ingestion:

  * from `variants`: FILTER, INFO/AA and the per-site called counts
    n_called / an_called / ac_called (NULL in older databases);
  * from `genotypes` (sparse: called non-reference genotypes only):
    the ALT dosage sum(gt) and the heterozygote count per group;
  * from `missing_genotypes`: the number of fully missing calls per group.

Each of the two sample-level tables is scanned ONCE per chromosome: the
scan groups by (site, label) with a LEFT JOIN to the panel, so samples
without a label fall into one extra "unlabelled" bucket and the
cohort-wide totals (needed for the exactness check and the MAF filter)
are the sum over all buckets. Work and memory are bounded by one
chromosome at a time: the caller (`analysis.prepare`) filters each
chromosome and keeps only the retained sites, in compact integer types.

A group's called haplotypes at a diploid autosomal site are
2 × (group size − missing in group). That derivation is exact for
diploid, fully called or fully missing genotypes; `analysis` checks it
against the per-site an_called / ac_called and drops the sites where it
is not (partial calls such as ./1, haploid or polyploid calls).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pyarrow as pa

from storage import sql_quote_str

GROUP_COLUMNS = ("population", "super_population")
ALL_GROUP = "all"
UNLABELLED = ""  # bucket label of samples without a panel label


@dataclass(frozen=True)
class Region:
    """A 1-based inclusive interval; start/end None = whole chromosome."""

    chrom: str
    start: int | None = None
    end: int | None = None


@dataclass(frozen=True)
class SchemaFeatures:
    """Which popgen-relevant parts of the schema this database has."""

    called_columns: bool
    missing_table: bool
    populations_table: bool


def schema_features() -> SchemaFeatures:
    from storage import table_columns, table_exists

    return SchemaFeatures(
        called_columns="n_called" in table_columns("variants"),
        missing_table=table_exists("missing_genotypes"),
        populations_table=table_exists("populations"),
    )


@dataclass
class ChromCounts:
    """Columnar per-site data of one chromosome (sites ordered by pos,
    ref, alt). Bucket columns of the 2-D arrays follow the `buckets`
    passed to `fetch_chrom`; the last bucket is UNLABELLED."""

    chrom: str
    pos: np.ndarray
    ref: list[str]
    alt: list[str]
    filter: list[str | None]
    aa: list[str | None]
    # Per-site called counts from `variants`; -1 where NULL.
    n_called: np.ndarray
    an_called: np.ndarray
    ac_called: np.ndarray
    # (n_sites, n_buckets) sums from genotypes / missing_genotypes.
    alt_by_bucket: np.ndarray
    het_by_bucket: np.ndarray
    missing_by_bucket: np.ndarray

    def __len__(self) -> int:
        return len(self.pos)

    @property
    def alt_total(self) -> np.ndarray:
        return self.alt_by_bucket.sum(axis=1)

    @property
    def het_total(self) -> np.ndarray:
        return self.het_by_bucket.sum(axis=1)

    @property
    def missing_total(self) -> np.ndarray:
        return self.missing_by_bucket.sum(axis=1)


def _arrow(sql: str) -> pa.Table:
    from storage import query_arrow

    return query_arrow(sql)


def _strings(table: pa.Table, col: str) -> list:
    return table.column(col).to_pylist()


def _ints(table: pa.Table, col: str, null: int = -1) -> np.ndarray:
    arr = table.column(col)
    if arr.null_count == len(arr):  # an all-NULL column may be typed `null`
        return np.full(len(arr), null, dtype=np.int64)
    arr = arr.cast(pa.int64())  # UInt32 columns cannot hold the -1 marker
    return arr.fill_null(null).to_numpy(zero_copy_only=False)


def region_sql(regions: list[Region], prefix: str = "") -> str:
    """WHERE fragment for `regions` (empty = no restriction). Matches both
    chr-prefixed and bare chromosome names."""
    from storage.gene_query import chrom_aliases

    if not regions:
        return ""
    parts = []
    for r in regions:
        names = ", ".join(sql_quote_str(c) for c in chrom_aliases(r.chrom))
        clause = f"{prefix}chrom IN ({names})"
        if r.start is not None:
            end = r.end if r.end is not None else r.start
            clause += f" AND {prefix}pos BETWEEN {int(r.start)} AND {int(end)}"
        parts.append(f"({clause})")
    return " AND (" + " OR ".join(parts) + ")"


def ingestions(sess) -> list[str]:
    t = _arrow("SELECT DISTINCT ingest_id FROM ingestions ORDER BY ingest_id")
    return _strings(t, "ingest_id")


def samples(sess, ingest_id: str) -> list[str]:
    t = _arrow(
        "SELECT DISTINCT sample_id FROM samples WHERE ingest_id = "
        f"{sql_quote_str(ingest_id)} ORDER BY sample_id"
    )
    return _strings(t, "sample_id")


def _labels_sql(ingest_id: str, by: str) -> str:
    """sample_id -> label for the samples CURRENTLY in the ingestion. The
    join to `samples` keeps a stale panel row (a sample re-ingested away)
    from inflating a group."""
    if by not in GROUP_COLUMNS:
        raise ValueError(f"unknown grouping {by!r}")
    iid = sql_quote_str(ingest_id)
    return (
        f"SELECT DISTINCT p.sample_id AS sample_id, p.{by} AS grp "
        "FROM populations AS p INNER JOIN samples AS s "
        "ON s.ingest_id = p.ingest_id AND s.sample_id = p.sample_id "
        f"WHERE p.ingest_id = {iid} AND p.{by} IS NOT NULL"
    )


def labels(sess, ingest_id: str, by: str) -> dict[str, str]:
    """sample_id -> group label (`population` or `super_population`) for
    the samples of `ingest_id` that have one and are still ingested."""
    t = _arrow(_labels_sql(ingest_id, by) + " ORDER BY sample_id")
    return dict(zip(_strings(t, "sample_id"), _strings(t, "grp"), strict=True))


def chromosomes(sess, ingest_id: str, regions: list[Region]) -> list[tuple[str, int]]:
    """(stored chromosome name, number of in-scope sites)."""
    t = _arrow(
        "SELECT chrom, CAST(count(*) AS BIGINT) AS n FROM variants WHERE "
        f"ingest_id = {sql_quote_str(ingest_id)}{region_sql(regions)} GROUP BY chrom"
    )
    return list(zip(_strings(t, "chrom"), t.column("n").to_pylist(), strict=True))


_SITE_ORDER = "ORDER BY pos, ref, alt"


def fetch_chrom(
    sess,
    ingest_id: str,
    chrom: str,
    regions: list[Region],
    by: str | None,
    buckets: list[str],
    features: SchemaFeatures,
) -> ChromCounts:
    """Load the in-scope sites of one chromosome with per-bucket counts.

    `by` None means no grouping: every sample is in the UNLABELLED bucket
    and `buckets` is [UNLABELLED]. Otherwise `buckets` lists the group
    labels followed by UNLABELLED.
    """
    iid = sql_quote_str(ingest_id)
    cq = sql_quote_str(chrom)
    where = f"ingest_id = {iid} AND chrom = {cq}{region_sql(regions)}"
    called = (
        "n_called, an_called, ac_called"
        if features.called_columns
        else "NULL AS n_called, NULL AS an_called, NULL AS ac_called"
    )
    # Site index = rank in (pos, ref, alt) order. The aggregates below
    # rebuild the same index in SQL, so their rows land on the right site
    # with an array store instead of a Python key lookup per row.
    index = f"CAST(row_number() OVER ({_SITE_ORDER}) - 1 AS BIGINT) AS idx"
    sites = _arrow(
        f"SELECT {index}, pos, ref, alt, filter, info_AA, {called} "
        f"FROM variants WHERE {where} ORDER BY idx"
    )
    n = sites.num_rows
    nb = len(buckets)
    out = ChromCounts(
        chrom=chrom,
        pos=_ints(sites, "pos"),
        ref=_strings(sites, "ref"),
        alt=_strings(sites, "alt"),
        filter=_strings(sites, "filter"),
        aa=_strings(sites, "info_AA"),
        n_called=_ints(sites, "n_called"),
        an_called=_ints(sites, "an_called"),
        ac_called=_ints(sites, "ac_called"),
        alt_by_bucket=np.zeros((n, nb), dtype=np.int64),
        het_by_bucket=np.zeros((n, nb), dtype=np.int64),
        missing_by_bucket=np.zeros((n, nb), dtype=np.int64),
    )
    if n == 0:
        return out

    column = {label: i for i, label in enumerate(buckets)}
    site_cte = f"s AS (SELECT {index}, pos, ref, alt FROM variants WHERE {where})"
    if by is None:
        ctes = f"WITH {site_cte}"
        label_join, grp = "", f"'{UNLABELLED}'"
    else:
        ctes = f"WITH {site_cte}, lab AS ({_labels_sql(ingest_id, by)})"
        label_join = "LEFT JOIN lab ON lab.sample_id = t.sample_id"
        # chDB fills an unmatched LEFT JOIN column with '' and DuckDB with
        # NULL; coalesce makes both the UNLABELLED bucket.
        grp = f"coalesce(lab.grp, '{UNLABELLED}')"
    twhere = f"t.ingest_id = {iid} AND t.chrom = {cq}{region_sql(regions, 't.')}"
    join = (
        "INNER JOIN s ON s.pos = t.pos AND s.ref = t.ref AND s.alt = t.alt "
        + label_join
    )

    def scatter(sql: str, targets: list[np.ndarray]) -> None:
        t = _arrow(sql)
        if t.num_rows == 0:
            return
        idx = _ints(t, "idx")
        col = np.array(
            [column.get(g, nb - 1) for g in _strings(t, "grp")], dtype=np.int64
        )
        for target, name in zip(targets, ("v1", "v2"), strict=False):
            target[idx, col] = _ints(t, name, 0)

    scatter(
        f"{ctes} SELECT s.idx AS idx, {grp} AS grp, "
        "CAST(sum(t.gt) AS BIGINT) AS v1, "
        "CAST(sum(CASE WHEN t.gt = 1 THEN 1 ELSE 0 END) AS BIGINT) AS v2 "
        f"FROM genotypes AS t {join} WHERE {twhere} AND t.gt > 0 "
        f"GROUP BY s.idx, {grp}",
        [out.alt_by_bucket, out.het_by_bucket],
    )
    if features.missing_table:
        scatter(
            f"{ctes} SELECT s.idx AS idx, {grp} AS grp, "
            "CAST(count(*) AS BIGINT) AS v1 "
            f"FROM missing_genotypes AS t {join} WHERE {twhere} "
            f"GROUP BY s.idx, {grp}",
            [out.missing_by_bucket],
        )
    return out
