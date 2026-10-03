"""Inheritance queries and family resolution for trio analysis."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, NamedTuple

import click


class Trio(NamedTuple):
    """A resolved family: the cohort and the three sample ids."""

    ingest_id: str
    proband: str
    father: str
    mother: str


class Gates(NamedTuple):
    """The tunable quality / rarity thresholds, passed together."""

    min_gq: int
    min_dp: int
    max_af: float
    min_ab: float
    max_ab: float


def _sole_ingest_id(name: str) -> str | None:
    """Return the db's ingest_id if it has exactly one ingestion, else
    None. Lets `db ped` infer the target when unambiguous."""
    from storage import get_session

    sess = get_session(name)
    raw = (
        sess.query("SELECT DISTINCT ingest_id FROM ingestions FORMAT TabSeparated")
        .bytes()
        .decode()
    )
    ids = [s for s in raw.splitlines() if s.strip()]
    return ids[0] if len(ids) == 1 else None


def _quality_gate(alias: str, min_gq: int, min_dp: int) -> str:
    return (
        f"({alias}.gq IS NULL OR {alias}.gq >= {min_gq}) "
        f"AND ({alias}.dp IS NULL OR {alias}.dp >= {min_dp})"
    )


def _het_ab_gate(alias: str, min_ab: float, max_ab: float) -> str:
    # Allele balance for a het call: ad_alt / (ad_ref + ad_alt) should sit
    # near 0.5. Applies only to het rows (gt=1) with AD present.
    denom = f"({alias}.ad_ref + {alias}.ad_alt)"
    frac = f"({alias}.ad_alt * 1.0 / {denom})"
    return (
        f"({alias}.gt != 1 OR {alias}.ad_ref IS NULL OR {alias}.ad_alt IS NULL "
        f"OR {denom} = 0 OR {frac} BETWEEN {min_ab} AND {max_ab})"
    )


def _parent_joins(father: str, mother: str) -> str:
    from storage import sql_quote_str

    key = "f.ingest_id = g.ingest_id AND f.chrom = g.chrom AND f.pos = g.pos AND f.ref = g.ref AND f.alt = g.alt"
    mkey = key.replace("f.", "m.")
    return (
        f"INNER JOIN genotypes f ON {key} AND f.sample_id = {sql_quote_str(father)} "
        f"INNER JOIN genotypes m ON {mkey} AND m.sample_id = {sql_quote_str(mother)} "
        "LEFT JOIN variants v ON v.ingest_id = g.ingest_id AND v.chrom = g.chrom "
        "AND v.pos = g.pos AND v.ref = g.ref AND v.alt = g.alt"
    )


def _where(trio: Trio, gates: Gates, model: str) -> str:
    """Shared WHERE clause: the genotype model plus the quality / rarity
    gates, scoped to this proband and ingest."""
    from storage import sql_quote_str

    predicate = " AND ".join(
        [
            model,
            _quality_gate("g", gates.min_gq, gates.min_dp),
            _quality_gate("f", gates.min_gq, gates.min_dp),
            _quality_gate("m", gates.min_gq, gates.min_dp),
            _het_ab_gate("g", gates.min_ab, gates.max_ab),
            f"(v.info_AF IS NULL OR v.info_AF <= {gates.max_af})",
        ]
    )
    return (
        f"g.ingest_id = {sql_quote_str(trio.ingest_id)} "
        f"AND g.sample_id = {sql_quote_str(trio.proband)} AND {predicate}"
    )


def trio_sql(category: str, trio: Trio, gates: Gates, *, count_only: bool) -> str:
    """Build the SQL for one inheritance model. The proband row is `g`;
    parents are joined as `f`/`m`; `v` brings population AF."""
    from storage import count_expr

    if category == "denovo":
        model = "g.gt > 0 AND f.gt = 0 AND m.gt = 0"
    elif category == "recessive":
        model = "g.gt = 2 AND f.gt = 1 AND m.gt = 1"
    elif category == "dominant":
        model = "g.gt = 1 AND ((f.gt > 0 AND m.gt = 0) OR (f.gt = 0 AND m.gt > 0))"
    else:
        raise ValueError(f"unknown category {category!r}")

    joins = _parent_joins(trio.father, trio.mother)
    where = _where(trio, gates, model)
    if count_only:
        return f"SELECT {count_expr()} FROM genotypes g {joins} WHERE {where}"
    return (
        "SELECT g.chrom, g.pos, g.ref, g.alt, g.gt AS proband_gt, "
        "f.gt AS father_gt, m.gt AS mother_gt, v.info_AF AS af "
        f"FROM genotypes g {joins} WHERE {where} ORDER BY g.chrom, g.pos"
    )


def _comphet_sql(trio: Trio, gates: Gates) -> str:
    """Candidate variants for compound-het: each is a rare proband het
    inherited from exactly one parent (the dominant pattern), tagged with
    parent-of-origin. The gene grouping happens in Python — genes live in
    the annotation store, which can't be SQL-joined to the cohort."""
    model = "g.gt = 1 AND ((f.gt > 0 AND m.gt = 0) OR (f.gt = 0 AND m.gt > 0))"
    where = _where(trio, gates, model)
    origin = (
        "CASE WHEN f.gt > 0 AND m.gt = 0 THEN 'paternal' "
        "WHEN m.gt > 0 AND f.gt = 0 THEN 'maternal' END AS origin"
    )
    return (
        f"SELECT g.chrom, g.pos, g.ref, g.alt, v.info_AF AS af, {origin} "
        f"FROM genotypes g {_parent_joins(trio.father, trio.mother)} "
        f"WHERE {where} ORDER BY g.chrom, g.pos"
    )


def _comphet_genes(candidate_rows: list) -> dict:
    """Group origin-tagged candidate variants by gene, keeping only genes
    that carry BOTH a paternal and a maternal het (trans configuration →
    both gene copies hit). A variant overlapping several genes counts for
    each. Returns {gene_symbol: {"paternal": [...], "maternal": [...]}}."""
    from annotations import gene_at

    genes: dict = {}
    for chrom, pos, ref, alt, af, origin in candidate_rows:
        if origin not in ("paternal", "maternal"):
            continue
        for gr in gene_at(chrom, int(pos)):
            entry = genes.setdefault(gr.gene_symbol, {"paternal": [], "maternal": []})
            entry[origin].append((chrom, pos, ref, alt, af))
    return {sym: e for sym, e in genes.items() if e["paternal"] and e["maternal"]}


def _resolve_parents(sess, ingest_id: str, proband: str) -> tuple[str, str]:
    from storage import sql_quote_str

    raw = (
        sess.query(
            f"SELECT father_id, mother_id FROM pedigree WHERE ingest_id = "
            f"{sql_quote_str(ingest_id)} AND sample_id = {sql_quote_str(proband)} "
            "FORMAT JSONCompact"
        )
        .bytes()
        .decode()
    )
    data = json.loads(raw)["data"]
    if not data:
        raise click.ClickException(
            f"no pedigree entry for proband {proband!r}. Load a PED first: "
            f"vcfclick db ped <name> <file.ped>"
        )
    father, mother = data[0]
    if father in ("0", "", None) or mother in ("0", "", None):
        raise click.ClickException(
            f"proband {proband!r} is missing a parent in the pedigree "
            f"(father={father!r}, mother={mother!r}); trio analysis needs both."
        )
    return father, mother


def _has_reference_rows(sess) -> bool:
    from storage import count_expr

    raw = (
        sess.query(
            f"SELECT {count_expr()} FROM genotypes WHERE gt = 0 FORMAT JSONCompact"
        )
        .bytes()
        .decode()
    )
    return int(json.loads(raw)["data"][0][0]) > 0


def _gnomad_keep(chrom: str, pos, ref: str, alt: str, max_af: float) -> bool:
    """Keep a candidate whose gnomAD popmax AF is <= max_af, or that is
    absent from the loaded gnomAD slice — absence is treated as rare
    (the slice may simply not cover the locus), never as AF 0."""
    from annotations import gnomad_af

    g = gnomad_af(chrom, int(pos), ref, alt)
    return g is None or g.popmax is None or g.popmax <= max_af


@dataclass(frozen=True)
class TrioAnalysis:
    """Query a resolved trio using the same gates for every output format."""

    session: Any
    trio: Trio
    gates: Gates
    has_reference: bool
    gnomad_max_af: float | None

    def _gnomad(self, rows: list) -> list:
        if self.gnomad_max_af is None:
            return rows
        return [row for row in rows if _gnomad_keep(*row[:4], self.gnomad_max_af)]

    def detail_rows(self, category: str) -> list:
        sql = trio_sql(category, self.trio, self.gates, count_only=False)
        rows = json.loads(self.session.query(sql, "JSONCompact").bytes().decode())[
            "data"
        ]
        return self._gnomad(rows)

    def count(self, category: str) -> int:
        if self.gnomad_max_af is not None:
            return len(self.detail_rows(category))
        sql = trio_sql(category, self.trio, self.gates, count_only=True)
        return json.loads(self.session.query(sql, "JSONCompact").bytes().decode())[
            "data"
        ][0][0]

    def comphet_genes(self) -> dict:
        sql = _comphet_sql(self.trio, self.gates)
        rows = json.loads(self.session.query(sql, "JSONCompact").bytes().decode())[
            "data"
        ]
        return _comphet_genes(self._gnomad(rows))
