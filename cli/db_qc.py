"""`vcfclick db qc` — per-sample quality-control metrics.

Self-contained over the cohort's sparse genotypes table: heterozygous /
homozygous-alt counts and their ratio, the transition/transversion ratio
(SNVs only), and a chromosome-X heterozygosity sex check that is flagged
against the pedigree's declared sex when one is loaded.

The genotypes table is sparse (only non-reference calls), so genotype
*missingness* / call rate is not computable here — those need the no-call
information a sparse store discards. The metrics below are exactly the
ones the stored non-reference calls support honestly.
"""

from __future__ import annotations

import json

import click

from cli.main import _set_db, db
from cli.sex import FEMALE_MIN_HETFRAC, MALE_MAX_HETFRAC, MIN_CHRX_CALLS, infer_sex
from cli.xchrom import BUILD_CHOICES, npx_sql, resolve_build
from storage.contigs import X_NAMES
from storage.sql import sql_quote_str, str_list

# chrX names come from storage.contigs. By default the het fraction is taken
# over the whole X (PAR variants slightly inflate a male's het count and are
# rare relative to the whole X); --x-nonpar restricts it to sure non-PAR X
# (cli/xchrom.py), the region `db trio` infers sex on.

# The thresholds and the call live in cli/sex.py, shared with `db trio`.
# `_infer_sex` keeps its name here: scripts/compare_somalier.py imports it.
_MALE_MAX_HETFRAC = MALE_MAX_HETFRAC
_FEMALE_MIN_HETFRAC = FEMALE_MIN_HETFRAC
_MIN_CHRX_CALLS = MIN_CHRX_CALLS
_infer_sex = infer_sex


def _transition() -> str:
    return (
        "((ref = 'A' AND alt = 'G') OR (ref = 'G' AND alt = 'A') "
        "OR (ref = 'C' AND alt = 'T') OR (ref = 'T' AND alt = 'C'))"
    )


def _count(expr: str, alias: str) -> str:
    # Backend-agnostic conditional count (chDB + DuckDB both support this).
    return f"sum(CASE WHEN {expr} THEN 1 ELSE 0 END) AS {alias}"


def _qc_sql(x_builds: dict[str, str] | None = None) -> str:
    """The per-sample counts. `x_builds` ({ingest_id: build}, `--x-nonpar`)
    restricts each ingestion's chrX counts to sure non-PAR X of its own build;
    None counts the whole X."""
    snv = "length(ref) = 1 AND length(alt) = 1"
    ti = f"{snv} AND {_transition()}"
    tv = f"{snv} AND NOT {_transition()}"
    if x_builds is None:
        chrx = f"chrom IN ({str_list(X_NAMES)})"
    elif x_builds:
        chrx = (
            "("
            + " OR ".join(
                f"(ingest_id = {sql_quote_str(i)} AND {npx_sql(b)})"
                for i, b in sorted(x_builds.items())
            )
            + ")"
        )
    else:
        chrx = "(1 = 0)"
    cols = ", ".join(
        [
            "ingest_id",
            "sample_id",
            "count(*) AS variants",
            _count("gt = 1", "het"),
            _count("gt = 2", "hom_alt"),
            _count(ti, "ti"),
            _count(tv, "tv"),
            _count(f"{chrx} AND gt = 1", "chrx_het"),
            _count(f"{chrx} AND gt > 0", "chrx_calls"),
        ]
    )
    return (
        f"SELECT {cols} FROM genotypes WHERE gt > 0 "
        "GROUP BY ingest_id, sample_id ORDER BY ingest_id, sample_id"
    )


def _ped_sex(sess) -> dict:
    """{(ingest_id, sample_id): 'male'|'female'} from the pedigree, if any."""
    from storage import table_exists

    if not table_exists("pedigree"):
        return {}
    raw = (
        sess.query("SELECT ingest_id, sample_id, sex FROM pedigree FORMAT JSONCompact")
        .bytes()
        .decode()
    )
    return {
        (r[0], r[1]): r[2]
        for r in json.loads(raw)["data"]
        if r[2] in ("male", "female")
    }


def _qc_builds(sess, option: str) -> dict[str, tuple[str, str]]:
    """{ingest_id: (build, source)} for --x-nonpar: each ingestion's build
    from its own chrX positions (an ingestion with none is `unknown`). An
    explicit GRCh37 that an ingestion's data contradicts is an error."""
    raw = sess.query(
        "SELECT i.ingest_id, max(v.pos) FROM ingestions i LEFT JOIN variants v "
        f"ON v.ingest_id = i.ingest_id AND v.chrom IN ({str_list(X_NAMES)}) "
        "GROUP BY i.ingest_id ORDER BY i.ingest_id",
        "JSONCompact",
    )
    out = {}
    for ingest_id, max_pos in json.loads(raw.bytes().decode())["data"]:
        build, source, conflict = resolve_build(option, max_pos)
        if conflict:
            raise click.ClickException(
                f"--build GRCh37, but ingestion {ingest_id!r} has chrX positions "
                "beyond the GRCh37 chrX length (155,270,560): the data is not "
                "GRCh37."
            )
        out[ingest_id] = (build, source)
    return out


def _ratio(num: int, denom: int) -> str:
    return "inf" if denom == 0 else f"{num / denom:.2f}"


@db.command(name="qc")
@click.argument("name")
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["table", "json"]),
    default="table",
    show_default=True,
)
@click.option(
    "--x-nonpar",
    is_flag=True,
    help="Take the chrX het fraction on sure non-PAR X only (the region "
    "`db trio` infers sex on), so both commands report the same sex.",
)
@click.option(
    "--build",
    type=click.Choice(BUILD_CHOICES),
    default="auto",
    show_default=True,
    help="With --x-nonpar: the assembly whose PAR is left out. auto: GRCh38 "
    "when a chrX position exceeds the GRCh37 chrX length, else only the "
    "region that is non-PAR on both builds. CHM13/T2T is not supported.",
)
def db_qc(name: str, out_format: str, x_nonpar: bool, build: str) -> None:
    """Per-sample QC: het/hom ratio, Ti/Tv, and a chrX-het sex check."""
    from storage import db_path, get_session

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} does not exist.")
    _set_db(name)
    sess = get_session(name)

    x_builds = None
    if x_nonpar:
        resolved = _qc_builds(sess, build)
        x_builds = {i: b for i, (b, _) in resolved.items()}
        for ingest_id, (b, source) in resolved.items():
            click.echo(
                f"note: chrX het fraction on non-PAR X ({b} [{source}]) for "
                f"ingestion {ingest_id}",
                err=True,
            )
    raw = sess.query(_qc_sql(x_builds), "JSONCompact").bytes().decode()
    rows = json.loads(raw)["data"]
    ped = _ped_sex(sess)

    samples = []
    for (
        ingest_id,
        sample_id,
        variants,
        het,
        hom_alt,
        ti,
        tv,
        chrx_het,
        chrx_calls,
    ) in rows:
        het, hom_alt, ti, tv = int(het), int(hom_alt), int(ti), int(tv)
        chrx_het, chrx_calls = int(chrx_het), int(chrx_calls)
        inferred, frac = _infer_sex(chrx_het, chrx_calls)
        declared = ped.get((ingest_id, sample_id))
        mismatch = bool(
            declared and inferred in ("male", "female") and inferred != declared
        )
        samples.append(
            {
                "ingest_id": ingest_id,
                "sample_id": sample_id,
                "variants": int(variants),
                "het": het,
                "hom_alt": hom_alt,
                "het_hom_ratio": _ratio(het, hom_alt),
                "ti_tv": _ratio(ti, tv),
                "chrx_het_frac": None if frac is None else round(frac, 3),
                "inferred_sex": inferred,
                "pedigree_sex": declared,
                "sex_mismatch": mismatch,
            }
        )

    if out_format == "json":
        click.echo(json.dumps(samples, indent=2))
        return

    if not samples:
        click.echo("no genotype rows to QC.")
        return
    click.echo(
        f"{'sample':<16}{'variants':>10}{'het':>8}{'hom':>8}"
        f"{'het/hom':>9}{'ti/tv':>7}{'chrX-het':>10}{'sex':>10}"
    )
    for s in samples:
        xf = "n/a" if s["chrx_het_frac"] is None else f"{s['chrx_het_frac']:.2f}"
        sex = s["inferred_sex"] + ("*" if s["sex_mismatch"] else "")
        click.echo(
            f"{s['sample_id']:<16}{s['variants']:>10,}{s['het']:>8,}"
            f"{s['hom_alt']:>8,}{s['het_hom_ratio']:>9}{s['ti_tv']:>7}"
            f"{xf:>10}{sex:>11}"
        )
    if any(s["sex_mismatch"] for s in samples):
        click.echo(
            "\n* inferred chrX sex disagrees with the pedigree's declared sex "
            "(possible sample swap or mislabel).",
            err=True,
        )
