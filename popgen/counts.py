"""Per-site, per-group allele counts from a vcfclick database.

The SQL here only aggregates; every statistic is computed in Python from
the counts it returns. All statements are plain SQL that both backends
accept (no `countIf`, no `FINAL`, `count(*)` rather than `count()`), so
chDB and DuckDB return identical numbers by construction.

What is counted, per site (chrom, pos, ref, alt) of one ingestion:

  * from `variants`: FILTER, INFO/AA and the per-site called counts
    n_called / an_called / ac_called (NULL in older databases);
  * from `genotypes` (sparse: called non-reference genotypes only):
    the ALT dosage sum(gt) and the heterozygote count, cohort-wide and
    per group (joined to `populations`);
  * from `missing_genotypes`: the number of fully missing calls,
    cohort-wide and per group.

A group's called haplotypes at a diploid autosomal site are then
2 × (group size − missing in group). That derivation is exact for
diploid, fully called or fully missing genotypes; `analysis` checks it
against the per-site an_called / ac_called and drops the sites where it
is not (partial calls such as ./1, haploid or polyploid calls).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from storage import sql_quote_str

GROUP_COLUMNS = ("population", "super_population")
ALL_GROUP = "all"


@dataclass(frozen=True)
class Region:
    """A 1-based inclusive interval; start/end None = whole chromosome."""

    chrom: str
    start: int | None = None
    end: int | None = None


@dataclass
class SiteCounts:
    """Columnar per-site data for one ingestion and its groups."""

    chrom: list[str]
    pos: np.ndarray
    ref: list[str]
    alt: list[str]
    filter: list[str | None]
    aa: list[str | None]
    # Per-site called counts from `variants`; -1 where NULL.
    n_called: np.ndarray
    an_called: np.ndarray
    ac_called: np.ndarray
    # Cohort-wide (every sample of the ingestion), derived from the
    # sparse genotypes + missing_genotypes tables.
    alt_total: np.ndarray
    missing_total: np.ndarray
    n_samples: int
    # Per group: name -> array over sites.
    group_size: dict[str, int] = field(default_factory=dict)
    group_alt: dict[str, np.ndarray] = field(default_factory=dict)
    group_het: dict[str, np.ndarray] = field(default_factory=dict)
    group_missing: dict[str, np.ndarray] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.chrom)


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


def _tsv(sess, sql: str) -> list[list[str | None]]:
    raw = sess.query(sql, "TSV").bytes().decode()
    rows = []
    for line in raw.splitlines():
        if not line:
            continue
        rows.append([None if f == "\\N" else f for f in line.split("\t")])
    return rows


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
    rows = _tsv(sess, "SELECT DISTINCT ingest_id FROM ingestions ORDER BY ingest_id")
    return [r[0] for r in rows]


def samples(sess, ingest_id: str) -> list[str]:
    rows = _tsv(
        sess,
        "SELECT DISTINCT sample_id FROM samples WHERE ingest_id = "
        f"{sql_quote_str(ingest_id)} ORDER BY sample_id",
    )
    return [r[0] for r in rows]


def labels(sess, ingest_id: str, by: str) -> dict[str, str]:
    """sample_id -> group label (`population` or `super_population`) for
    the samples of `ingest_id` that have one."""
    if by not in GROUP_COLUMNS:
        raise ValueError(f"unknown grouping {by!r}")
    rows = _tsv(
        sess,
        f"SELECT sample_id, {by} FROM populations WHERE ingest_id = "
        f"{sql_quote_str(ingest_id)} AND {by} IS NOT NULL ORDER BY sample_id",
    )
    return {r[0]: r[1] for r in rows}


def any_ancestral(sess, ingest_id: str, regions: list[Region]) -> bool:
    """Does any in-scope site carry INFO/AA?"""
    rows = _tsv(
        sess,
        "SELECT count(*) FROM variants WHERE ingest_id = "
        f"{sql_quote_str(ingest_id)} AND info_AA IS NOT NULL{region_sql(regions)}",
    )
    return bool(rows) and int(rows[0][0]) > 0


_KEY = "chrom, pos, ref, alt"


def _site_aggregate(
    sess, sql: str, index: dict[tuple, int], n_sites: int, n_values: int
) -> dict[str, list[np.ndarray]]:
    """Run an aggregate returning (chrom, pos, ref, alt, group, v1..vk)
    and scatter it into per-group arrays aligned with `index`."""
    out: dict[str, list[np.ndarray]] = {}
    for row in _tsv(sess, sql):
        key = (row[0], int(row[1]), row[2], row[3])
        i = index.get(key)
        if i is None:  # genotype row without a variants row: ignore
            continue
        arrays = out.setdefault(
            row[4], [np.zeros(n_sites, dtype=np.int64) for _ in range(n_values)]
        )
        for a, v in zip(arrays, row[5:], strict=True):
            a[i] = int(v)
    return out


def fetch_counts(
    sess,
    ingest_id: str,
    regions: list[Region],
    by: str | None,
    features: SchemaFeatures,
    group_labels: dict[str, str],
) -> SiteCounts:
    """Load every in-scope site of `ingest_id` with its counts.

    `by` None means a single cohort-wide group named "all"; otherwise the
    samples are grouped by that `populations` column (`group_labels` maps
    sample_id -> label, used for the group sizes).
    """
    iid = sql_quote_str(ingest_id)
    where = f"ingest_id = {iid}{region_sql(regions)}"
    called = (
        "n_called, an_called, ac_called"
        if features.called_columns
        else "NULL, NULL, NULL"
    )
    site_rows = _tsv(
        sess,
        f"SELECT {_KEY}, filter, info_AA, {called} FROM variants "
        f"WHERE {where} ORDER BY {_KEY}",
    )
    index = {(r[0], int(r[1]), r[2], r[3]): i for i, r in enumerate(site_rows)}
    n = len(site_rows)

    def ints(col: int) -> np.ndarray:
        return np.array(
            [-1 if r[col] is None else int(r[col]) for r in site_rows], dtype=np.int64
        )

    sample_ids = samples(sess, ingest_id)
    counts = SiteCounts(
        chrom=[r[0] for r in site_rows],
        pos=np.array([int(r[1]) for r in site_rows], dtype=np.int64),
        ref=[r[2] for r in site_rows],
        alt=[r[3] for r in site_rows],
        filter=[r[4] for r in site_rows],
        aa=[r[5] for r in site_rows],
        n_called=ints(6),
        an_called=ints(7),
        ac_called=ints(8),
        alt_total=np.zeros(n, dtype=np.int64),
        missing_total=np.zeros(n, dtype=np.int64),
        n_samples=len(sample_ids),
    )
    if n == 0:
        return counts

    geno_values = "sum(gt), sum(CASE WHEN gt = 1 THEN 1 ELSE 0 END)"
    total = _site_aggregate(
        sess,
        f"SELECT {_KEY}, '{ALL_GROUP}', {geno_values} FROM genotypes "
        f"WHERE {where} AND gt > 0 GROUP BY {_KEY}",
        index,
        n,
        2,
    ).get(ALL_GROUP)
    het_total = np.zeros(n, dtype=np.int64)
    if total is not None:
        counts.alt_total, het_total = total
    if features.missing_table:
        miss = _site_aggregate(
            sess,
            f"SELECT {_KEY}, '{ALL_GROUP}', count(*) FROM missing_genotypes "
            f"WHERE {where} GROUP BY {_KEY}",
            index,
            n,
            1,
        ).get(ALL_GROUP)
        if miss is not None:
            counts.missing_total = miss[0]

    if by is None:
        counts.group_size = {ALL_GROUP: len(sample_ids)}
        counts.group_alt[ALL_GROUP] = counts.alt_total
        counts.group_het[ALL_GROUP] = het_total
        counts.group_missing[ALL_GROUP] = counts.missing_total
        return counts

    for label in group_labels.values():
        counts.group_size[label] = counts.group_size.get(label, 0) + 1
    gwhere = f"g.ingest_id = {iid}{region_sql(regions, 'g.')}"
    join = (
        "INNER JOIN populations AS p "
        "ON p.ingest_id = g.ingest_id AND p.sample_id = g.sample_id"
    )
    gkey = "g.chrom, g.pos, g.ref, g.alt"
    per_group = _site_aggregate(
        sess,
        f"SELECT {gkey}, p.{by}, sum(g.gt), "
        "sum(CASE WHEN g.gt = 1 THEN 1 ELSE 0 END) "
        f"FROM genotypes AS g {join} WHERE {gwhere} AND g.gt > 0 "
        f"AND p.{by} IS NOT NULL GROUP BY {gkey}, p.{by}",
        index,
        n,
        2,
    )
    per_group_missing = (
        _site_aggregate(
            sess,
            f"SELECT {gkey}, p.{by}, count(*) FROM missing_genotypes AS g {join} "
            f"WHERE {gwhere} AND p.{by} IS NOT NULL GROUP BY {gkey}, p.{by}",
            index,
            n,
            1,
        )
        if features.missing_table
        else {}
    )
    zeros = np.zeros(n, dtype=np.int64)
    for label in counts.group_size:
        alt, het = per_group.get(label, (zeros, zeros))
        counts.group_alt[label] = alt
        counts.group_het[label] = het
        counts.group_missing[label] = per_group_missing.get(label, [zeros])[0]
    return counts
