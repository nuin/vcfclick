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

Honest scope: this is candidate FILTERING, not variant calling. Quality
gates are only as strong as the FORMAT fields present (gq/dp/ad are
often absent on public joint-call releases, in which case the gates
pass-through). De novo confidence is bounded by genotype-level
evidence; PL-based Bayesian refinement is a future enhancement.
"""

from __future__ import annotations

import json

import click

from cli.main import _set_db, db


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


def trio_sql(
    category: str,
    ingest_id: str,
    proband: str,
    father: str,
    mother: str,
    *,
    min_gq: int,
    min_dp: int,
    max_af: float,
    min_ab: float,
    max_ab: float,
    count_only: bool,
) -> str:
    """Build the SQL for one inheritance model. The proband row is `g`;
    parents are joined as `f`/`m`; `v` brings population AF."""
    from storage import count_expr, sql_quote_str

    iid = sql_quote_str(ingest_id)
    pro = sql_quote_str(proband)

    # Per-model genotype predicate.
    if category == "denovo":
        model = "g.gt > 0 AND f.gt = 0 AND m.gt = 0"
    elif category == "recessive":
        model = "g.gt = 2 AND f.gt = 1 AND m.gt = 1"
    elif category == "dominant":
        model = "g.gt = 1 AND ((f.gt > 0 AND m.gt = 0) OR (f.gt = 0 AND m.gt > 0))"
    else:
        raise ValueError(f"unknown category {category!r}")

    gates = " AND ".join(
        [
            model,
            _quality_gate("g", min_gq, min_dp),
            _quality_gate("f", min_gq, min_dp),
            _quality_gate("m", min_gq, min_dp),
            _het_ab_gate("g", min_ab, max_ab),
            f"(v.info_AF IS NULL OR v.info_AF <= {max_af})",
        ]
    )
    where = f"g.ingest_id = {iid} AND g.sample_id = {pro} AND {gates}"

    if count_only:
        return (
            f"SELECT {count_expr()} FROM genotypes g {_parent_joins(father, mother)} "
            f"WHERE {where}"
        )
    return (
        "SELECT g.chrom, g.pos, g.ref, g.alt, g.gt AS proband_gt, "
        "f.gt AS father_gt, m.gt AS mother_gt, v.info_AF AS af "
        f"FROM genotypes g {_parent_joins(father, mother)} "
        f"WHERE {where} ORDER BY g.chrom, g.pos"
    )


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


@db.command(name="trio")
@click.argument("name")
@click.option("--proband", required=True, help="Sample id of the affected child.")
@click.option(
    "--category",
    type=click.Choice(["denovo", "recessive", "dominant", "all"]),
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
def db_trio(
    name: str,
    proband: str,
    category: str,
    min_gq: int,
    min_dp: int,
    max_af: float,
    min_ab: float,
    max_ab: float,
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

    has_ref = _has_reference_rows(sess)
    needs_ref = {"denovo", "dominant"}

    def run(cat: str, count_only: bool):
        sql = trio_sql(
            cat,
            ingest_id,
            proband,
            father,
            mother,
            min_gq=min_gq,
            min_dp=min_dp,
            max_af=max_af,
            min_ab=min_ab,
            max_ab=max_ab,
            count_only=count_only,
        )
        return sess.query(sql, "JSONCompact").bytes().decode()

    click.echo(f"trio: proband={proband} father={father} mother={mother}")
    if not has_ref and (category in needs_ref or category == "all"):
        click.echo(
            "  note: this database has no stored hom-reference calls, so "
            "de-novo/dominant cannot prove a parent is 0/0 (vs no-call). "
            "Re-ingest with `--keep-reference` for those models.",
            err=True,
        )

    cats = ["denovo", "recessive", "dominant"] if category == "all" else [category]

    if category == "all":
        for cat in cats:
            n = json.loads(run(cat, count_only=True))["data"][0][0]
            blocked = (
                "" if has_ref or cat not in needs_ref else "  (needs --keep-reference)"
            )
            click.echo(f"  {cat:10s} {n:>6}{blocked}")
        return

    data = json.loads(run(category, count_only=False))["data"]
    click.echo(f"\n{category} candidates: {len(data)}")
    for row in data:
        chrom, pos, ref, alt, pgt, fgt, mgt, af = row
        af_s = "NA" if af is None else f"{af}"
        click.echo(
            f"  {chrom}:{pos} {ref}>{alt}  proband_gt={pgt} "
            f"father_gt={fgt} mother_gt={mgt}  AF={af_s}"
        )
