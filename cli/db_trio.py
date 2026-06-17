"""`vcfclick db ped` + `db trio` — family-based analysis.

A PED file declares who-is-whose-parent among already-ingested samples.
Load it, then `db trio` resolves a proband's parents and reports
candidate variants under each Mendelian inheritance model.

    vcfclick db ped  fam1 fam1.ped
    vcfclick db trio fam1 --proband CHILD

Inheritance models, with slivar-style genotype quality gates (GQ, depth,
allele balance) and population-AF rarity:

  * de novo     proband carries; BOTH parents provably hom-ref (gt=0).
                Requires a `--keep-reference` ingest — without stored
                parent 0/0 rows, "parent absent" is a no-call, not
                reference, so de novo is undecidable.
  * recessive   proband hom-alt; both parents heterozygous carriers.
                Works on a normal sparse ingest.
  * dominant    proband het; exactly one parent carries, the other
                provably hom-ref. Also needs --keep-reference (to prove
                the non-carrier parent is reference, not no-call).
  * comphet     two rare proband hets in the SAME gene, one inherited
                from each parent (trans → both gene copies hit). Needs
                gene annotations loaded (`vcfclick annotations load`) and
                --keep-reference (to prove each non-carrier parent is
                hom-ref). Reported per gene, not per variant.

Honest scope: this is candidate FILTERING, not variant calling. Quality
gates are only as strong as the FORMAT fields present (gq/dp/ad are
often absent on public joint-call releases, in which case the gates
pass-through). De novo confidence is bounded by genotype-level
evidence; PL-based Bayesian refinement is a future enhancement.
"""

from __future__ import annotations

import json
from typing import NamedTuple

import click

from cli.main import _set_db, db


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
    from storage import get_session, sql_quote_str  # noqa: F401

    sess = get_session(name)
    raw = (
        sess.query("SELECT DISTINCT ingest_id FROM ingestions FORMAT TabSeparated")
        .bytes()
        .decode()
    )
    ids = [s for s in raw.splitlines() if s.strip()]
    return ids[0] if len(ids) == 1 else None


@db.command(name="ped")
@click.argument("name")
@click.argument("ped_path", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--ingest-id",
    default=None,
    help="Ingest_id the pedigree's samples belong to. Inferred when the "
    "database has exactly one ingestion.",
)
def db_ped(name: str, ped_path: str, ingest_id: str | None) -> None:
    """Load family relationships from a PED/FAM file into NAME.

    The PED's individual ids must match sample ids already ingested
    under the target ingest_id (v1 assumes a joint-called trio, so all
    members share one ingest_id). Re-loading replaces the prior
    pedigree for that ingest_id.
    """
    from storage import db_path

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} does not exist.")

    _set_db(name)

    if ingest_id is None:
        ingest_id = _sole_ingest_id(name)
        if ingest_id is None:
            raise click.ClickException(
                "database has multiple ingestions; pass --ingest-id to say "
                "which cohort the pedigree applies to."
            )

    from ingest.pedigree import load_pedigree

    try:
        n = load_pedigree(ingest_id, ped_path)
    except ValueError as e:
        raise click.ClickException(str(e)) from e
    click.echo(f"loaded pedigree: {n} individuals under ingest_id={ingest_id}")


# ─────────────────────── trio analysis ───────────────────────

# Shared genotype quality gates (slivar-style), applied per sample alias.
# Lenient on NULL: gq/dp/ad are frequently absent on public joint-call
# VCFs, and a strict gate there would silently drop every row (see
# docs/SCHEMA.md). So the gate filters only where the data exists.


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


@db.command(name="trio")
@click.argument("name")
@click.option("--proband", required=True, help="Sample id of the affected child.")
@click.option(
    "--category",
    type=click.Choice(["denovo", "recessive", "dominant", "comphet", "all"]),
    default="all",
    show_default=True,
    help="Inheritance model. 'all' prints per-model candidate counts.",
)
@click.option("--min-gq", default=20, show_default=True, type=int)
@click.option("--min-dp", default=10, show_default=True, type=int)
@click.option(
    "--max-af",
    default=0.01,
    show_default=True,
    type=float,
    help="Keep variants with population info_AF <= this (rarity filter).",
)
@click.option("--min-ab", default=0.25, show_default=True, type=float)
@click.option("--max-ab", default=0.75, show_default=True, type=float)
@click.option(
    "--gnomad-max-af",
    default=None,
    type=float,
    help="Additionally drop candidates whose gnomAD popmax AF exceeds this "
    "(needs `vcfclick annotations load-gnomad`). Variants absent from the "
    "loaded gnomAD slice are kept as rare.",
)
def db_trio(
    name: str,
    proband: str,
    category: str,
    min_gq: int,
    min_dp: int,
    max_af: float,
    min_ab: float,
    max_ab: float,
    gnomad_max_af: float | None,
) -> None:
    """Report candidate variants under Mendelian inheritance models for
    a trio, with genotype quality gates and an AF rarity filter."""
    from storage import db_path, get_session

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} does not exist.")
    _set_db(name)
    sess = get_session(name)

    ingest_id = _sole_ingest_id(name)
    if ingest_id is None:
        raise click.ClickException(
            "database has multiple ingestions; trio analysis assumes one "
            "joint-called cohort. Re-ingest the trio as a single VCF "
            "(see `vcfclick merge`)."
        )
    father, mother = _resolve_parents(sess, ingest_id, proband)
    trio = Trio(ingest_id, proband, father, mother)
    gates = Gates(min_gq, min_dp, max_af, min_ab, max_ab)

    has_ref = _has_reference_rows(sess)
    needs_ref = {"denovo", "dominant", "comphet"}

    def _gnomad(rows: list) -> list:
        if gnomad_max_af is None:
            return rows
        return [r for r in rows if _gnomad_keep(r[0], r[1], r[2], r[3], gnomad_max_af)]

    def detail_rows(cat: str) -> list:
        """Candidate rows for a model, after the optional gnomAD filter."""
        sql = trio_sql(cat, trio, gates, count_only=False)
        rows = json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]
        return _gnomad(rows)

    def count(cat: str) -> int:
        if gnomad_max_af is not None:
            return len(detail_rows(cat))
        sql = trio_sql(cat, trio, gates, count_only=True)
        return json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"][0][0]

    def comphet_genes() -> dict:
        sql = _comphet_sql(trio, gates)
        rows = json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]
        return _comphet_genes(_gnomad(rows))

    click.echo(f"trio: proband={proband} father={father} mother={mother}")
    if not has_ref and (category in needs_ref or category == "all"):
        click.echo(
            "  note: this database has no stored hom-reference calls, so "
            "de-novo/dominant/comphet cannot prove a parent is 0/0 (vs "
            "no-call). Re-ingest with `--keep-reference` for those models.",
            err=True,
        )

    if category == "all":
        for cat in ["denovo", "recessive", "dominant"]:
            n = count(cat)
            blocked = (
                "" if has_ref or cat not in needs_ref else "  (needs --keep-reference)"
            )
            click.echo(f"  {cat:10s} {n:>6}{blocked}")
        genes = comphet_genes()
        blocked = "" if has_ref else "  (needs --keep-reference)"
        click.echo(f"  {'comphet':10s} {len(genes):>6}{blocked}  genes")
        return

    if category == "comphet":
        genes = comphet_genes()
        click.echo(f"\ncomphet candidate genes: {len(genes)}")
        for sym in sorted(genes):
            entry = genes[sym]
            click.echo(f"  {sym}")
            for origin in ("paternal", "maternal"):
                for chrom, pos, ref, alt, af in entry[origin]:
                    af_s = "NA" if af is None else f"{af}"
                    click.echo(f"    {origin:8s} {chrom}:{pos} {ref}>{alt}  AF={af_s}")
        return

    data = detail_rows(category)
    click.echo(f"\n{category} candidates: {len(data)}")
    for row in data:
        chrom, pos, ref, alt, pgt, fgt, mgt, af = row
        af_s = "NA" if af is None else f"{af}"
        click.echo(
            f"  {chrom}:{pos} {ref}>{alt}  proband_gt={pgt} "
            f"father_gt={fgt} mother_gt={mgt}  AF={af_s}"
        )
