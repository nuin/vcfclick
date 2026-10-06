"""`vcfclick db popgen` — population-genetics statistics.

    vcfclick db popgen summary kg
    vcfclick db popgen sfs     kg --project 20 --format json
    vcfclick db popgen fst     kg --by super_population
    vcfclick db popgen windows kg --region 22:16e6-17e6 --window 100000 --fst

Per-site, per-group allele counts come from SQL (either backend); the
statistics are computed in Python (popgen/). Groups are the panel's
populations (`vcfclick db panel`), or the whole cohort as one group "all"
when no panel is loaded. See docs/POPGEN.md for definitions and choices.
"""

from __future__ import annotations

import functools
import json
from dataclasses import dataclass

import click

from cli.main import _set_db, db

_FORMATS = ("table", "tsv", "json")


@dataclass(frozen=True)
class CommonOptions:
    name: str
    ingest_id: str | None
    regions: tuple[str, ...]
    gene: str | None
    by: str
    include_indels: bool
    pass_only: bool
    min_call_rate: float
    maf: float
    ancestral: str | None
    include_sex_chroms: bool
    allow_untracked: bool
    out_format: str


def _common_options(default_format: str = "table"):
    """Click options shared by every popgen subcommand."""

    def decorate(f):
        options = [
            click.argument("name"),
            click.option(
                "--ingest-id",
                default=None,
                help="Ingestion to analyse (required when NAME has several).",
            ),
            click.option(
                "--region",
                "regions",
                multiple=True,
                help="chr, chr:pos or chr:start-end (repeatable).",
            ),
            click.option(
                "--gene",
                default=None,
                help="Restrict to a gene (needs `vcfclick annotations load`).",
            ),
            click.option(
                "--by",
                type=click.Choice(["population", "super_population", "all"]),
                default="population",
                show_default=True,
                help="Grouping. Falls back to one group 'all' when no panel is "
                "loaded for the ingestion.",
            ),
            click.option(
                "--include-indels",
                is_flag=True,
                help="Also use biallelic indels (default: biallelic SNVs only).",
            ),
            click.option(
                "--pass-only/--all-filters",
                default=True,
                show_default=True,
                help="Keep only sites whose FILTER is PASS or missing.",
            ),
            click.option(
                "--min-call-rate",
                type=click.FloatRange(0, 1),
                default=0.9,
                show_default=True,
                help="Minimum call rate a site needs in EVERY group.",
            ),
            click.option(
                "--maf",
                type=click.FloatRange(0, 0.5),
                default=0.0,
                show_default=True,
                help="Minimum cohort-wide minor-allele frequency (0 keeps "
                "monomorphic sites, which θ per site needs).",
            ),
            click.option(
                "--ancestral",
                type=click.Choice(["aa", "aa-high", "ref", "none"]),
                default=None,
                help="Polarisation: INFO/AA at any confidence (aa), "
                "high-confidence only (aa-high), REF as ancestral (ref), or "
                "folded statistics only (none). Default: aa when any site "
                "carries INFO/AA, else none.",
            ),
            click.option(
                "--include-sex-chroms",
                is_flag=True,
                help="Not supported yet: X/Y need per-sample ploidy.",
            ),
            click.option(
                "--allow-untracked",
                is_flag=True,
                help="Compute even when the ingestion has no missing-call "
                "tracking (ingested before it existed, or with "
                "--no-record-missing). Output then says "
                "missing_data_tracked: false.",
            ),
            click.option(
                "--format",
                "out_format",
                type=click.Choice(_FORMATS),
                default=default_format,
                show_default=True,
            ),
        ]
        for opt in reversed(options):
            f = opt(f)
        return f

    return decorate


def _project_option(f):
    return click.option(
        "--project",
        type=click.IntRange(2, None),
        default=None,
        help="SFS projection size in haplotypes (default: the smallest "
        "number of called haplotypes at any retained site, per group).",
    )(f)


def _split_common(kwargs: dict) -> tuple[CommonOptions, dict]:
    common = {
        k: kwargs.pop(k)
        for k in list(kwargs)
        if k in CommonOptions.__dataclass_fields__
    }
    return CommonOptions(**common), kwargs


def _resolve_regions(options: CommonOptions):
    from export.vcf import parse_region
    from popgen.analysis import is_autosome
    from popgen.counts import Region

    regions = []
    for text in options.regions:
        chrom, start, end = parse_region(text)
        regions.append(Region(chrom, start, end))
    if options.gene:
        from annotations.db import position_for_gene

        g = position_for_gene(options.gene)
        if g is None:
            raise click.ClickException(
                f"gene {options.gene.upper()!r} not found; load gene coordinates "
                "with `vcfclick annotations load`"
            )
        regions.append(Region(g.chrom, int(g.start_pos), int(g.end_pos)))
    bad = [r.chrom for r in regions if not is_autosome(r.chrom)]
    if bad:
        raise click.ClickException(
            f"{', '.join(sorted(set(bad)))}: sex chromosomes and MT are not "
            "supported by `db popgen` yet (autosomes only)"
        )
    return regions


def _resolve_ingest(sess, options: CommonOptions) -> str:
    from popgen.counts import ingestions

    have = ingestions(sess)
    if not have:
        raise click.ClickException(f"{options.name!r} has no ingestions")
    if options.ingest_id is None:
        if len(have) > 1:
            raise click.ClickException(
                f"{options.name!r} has {len(have)} ingestions ({', '.join(have)}); "
                "pick one with --ingest-id"
            )
        return have[0]
    if options.ingest_id not in have:
        raise click.ClickException(
            f"no ingestion {options.ingest_id!r} in {options.name!r}; "
            f"have: {', '.join(have)}"
        )
    return options.ingest_id


def _prepare(options: CommonOptions):
    from popgen.analysis import PopgenError, SiteFilters, prepare
    from storage import db_path, get_session

    if options.include_sex_chroms:
        raise click.ClickException(
            "--include-sex-chroms is not supported yet: X/Y statistics need "
            "per-sample ploidy (hemizygous males), which v1 does not model. "
            "Analyse autosomes, or subset females with a panel."
        )
    if not db_path(options.name).exists():
        raise click.ClickException(f"db {options.name!r} does not exist.")
    _set_db(options.name)
    sess = get_session(options.name)
    try:
        regions = _resolve_regions(options)
        ingest_id = _resolve_ingest(sess, options)
        prep = prepare(
            sess,
            ingest_id,
            regions,
            options.by,
            SiteFilters(
                include_indels=options.include_indels,
                pass_only=options.pass_only,
                min_call_rate=options.min_call_rate,
                maf=options.maf,
            ),
            options.ancestral,
            allow_untracked=options.allow_untracked,
        )
    except (PopgenError, ValueError) as e:
        raise click.ClickException(str(e)) from e
    for w in prep.warnings:
        click.echo(f"warning: {w}", err=True)
    return prep


def _header(prep, options: CommonOptions, command: str) -> dict:
    return {
        "command": command,
        "db": options.name,
        "ingest_id": prep.ingest_id,
        "grouping": prep.grouping,
        "groups": list(prep.groups),
        "ancestral": prep.ancestral,
        "missing_data_tracked": prep.missing_data_tracked,
        "filters": {
            "regions": list(options.regions),
            "gene": options.gene,
            "autosomes_only": True,
            "include_indels": options.include_indels,
            "pass_only": options.pass_only,
            "min_call_rate": options.min_call_rate,
            "maf": options.maf,
        },
        "sites": prep.report,
        "warnings": prep.warnings,
    }


# ─────────────────────────────── output ─────────────────────────────────


def _cell(v) -> str:
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def _emit_rows(rows: list[dict], columns: list[str], out_format: str) -> None:
    if out_format == "tsv":
        click.echo("\t".join(columns))
        for r in rows:
            click.echo("\t".join(_cell(r.get(c)) for c in columns))
        return
    cells = [[_cell(r.get(c)) for c in columns] for r in rows]
    widths = [
        max([len(c)] + [len(row[i]) for row in cells]) for i, c in enumerate(columns)
    ]
    click.echo("  ".join(c.rjust(w) for c, w in zip(columns, widths, strict=True)))
    for row in cells:
        click.echo("  ".join(v.rjust(w) for v, w in zip(row, widths, strict=True)))


def _emit_preamble(header: dict) -> None:
    s = header["sites"]
    dropped = ", ".join(f"{k}={v}" for k, v in s["dropped"].items() if v)
    click.echo(
        f"ingest_id: {header['ingest_id']}   grouping: {header['grouping']}   "
        f"ancestral: {header['ancestral']}   "
        f"missing_data_tracked: {str(header['missing_data_tracked']).lower()}"
    )
    click.echo(
        f"sites: {s['retained']} retained of {s['in_scope']} in scope"
        + (f" (dropped: {dropped})" if dropped else "")
        + (
            f"; polarised {s['polarised']}, unpolarised {s['unpolarised']}"
            if header["ancestral"] != "none"
            else ""
        )
    )
    if s["unlabelled_samples"]:
        click.echo(
            f"samples without a panel label (excluded): {s['unlabelled_samples']}"
        )
    click.echo()


def _round(v):
    """12 significant digits: hides last-bit float noise (log-space SFS
    projection) without touching any meaningful digit."""
    if isinstance(v, float):
        return float(f"{v:.12g}")
    if isinstance(v, dict):
        return {k: _round(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_round(x) for x in v]
    return v


def _emit(header: dict, payload: dict, out_format: str, table_fn) -> None:
    if out_format == "json":
        click.echo(json.dumps(_round({**header, **payload}), indent=2))
        return
    if out_format == "table":
        _emit_preamble(header)
    table_fn()


def _command(name: str, default_format: str = "table", project: bool = False):
    """Register a popgen subcommand with the shared options."""

    def decorate(f):
        @functools.wraps(f)
        def run(**kwargs):
            options, rest = _split_common(kwargs)
            prep = _prepare(options)
            if "project" in rest:
                from popgen.analysis import projection_warnings

                for w in projection_warnings(prep, rest["project"]):
                    prep.warnings.append(w)
                    click.echo(f"warning: {w}", err=True)
            return f(options, prep, **rest)

        cmd = _common_options(default_format)(run)
        if project:
            cmd = _project_option(cmd)
        return popgen.command(name=name)(cmd)

    return decorate


@db.group(name="popgen")
def popgen() -> None:
    """Population-genetics statistics (θ, π, Tajima's D, SFS, F_ST)."""


_SUMMARY_COLUMNS = [
    "group",
    "n_samples",
    "sites",
    "segregating_sites",
    "theta_w",
    "theta_w_per_site",
    "pi",
    "pi_per_site",
    "tajima_d",
    "fay_wu_h",
    "projection_n",
    "ho",
    "he",
    "f",
]


@_command("summary", project=True)
def popgen_summary(options: CommonOptions, prep, project: int | None) -> None:
    """Per group: S, θ_W, π, Tajima's D, Fay & Wu's H, Ho/He/F."""
    from popgen.analysis import summary

    rows = summary(prep, project)
    _emit(
        _header(prep, options, "summary"),
        {"results": rows},
        options.out_format,
        lambda: _emit_rows(rows, _SUMMARY_COLUMNS, options.out_format),
    )


@_command("sfs", project=True)
def popgen_sfs(options: CommonOptions, prep, project: int | None) -> None:
    """Per group site-frequency spectrum, projected for missing data."""
    from popgen.analysis import sfs

    rows = sfs(prep, project)

    def table() -> None:
        long = []
        for r in rows:
            for kind in ("unfolded", "folded"):
                for j, v in enumerate(r[kind] or []):
                    long.append(
                        {
                            "group": r["group"],
                            "spectrum": kind,
                            "projection_n": r["projection_n"],
                            "count": j,
                            "sites": v,
                        }
                    )
        _emit_rows(
            long,
            ["group", "spectrum", "projection_n", "count", "sites"],
            options.out_format,
        )
        if options.out_format == "table":
            for r in rows:
                click.echo(
                    f"\n{r['group']}: {r['sites_used']} sites projected to "
                    f"n={r['projection_n']} ({r['sites_dropped']} with fewer "
                    f"called haplotypes dropped); "
                    f"{r['polarised_sites_used']} polarised"
                )

    _emit(_header(prep, options, "sfs"), {"results": rows}, options.out_format, table)


@_command("fst")
@click.option(
    "--window", type=click.IntRange(1, None), default=None, help="Window size (bp)."
)
@click.option(
    "--step",
    type=click.IntRange(1, None),
    default=None,
    help="Window step (bp; default = window).",
)
def popgen_fst(
    options: CommonOptions, prep, window: int | None, step: int | None
) -> None:
    """Pairwise Hudson F_ST between groups (Bhatia et al. 2013)."""
    from popgen.analysis import fst

    if len(prep.groups) < 2:
        raise click.ClickException(
            "F_ST needs at least two groups; load a panel with `vcfclick db panel`"
        )
    result = fst(prep, window, step)

    def table() -> None:
        if options.out_format == "table" or not window:
            _emit_rows(
                result["pairs"],
                ["group1", "group2", "fst", "sites"],
                options.out_format,
            )
        if window:
            if options.out_format == "table":
                click.echo()
            cols = list(result["windows"][0]) if result["windows"] else []
            _emit_rows(result["windows"], cols, options.out_format)

    _emit(_header(prep, options, "fst"), result, options.out_format, table)


@_command("windows", default_format="tsv", project=True)
@click.option(
    "--window", type=click.IntRange(1, None), required=True, help="Window size (bp)."
)
@click.option(
    "--step",
    type=click.IntRange(1, None),
    default=None,
    help="Window step (bp; default = window).",
)
@click.option(
    "--fst", "with_fst", is_flag=True, help="Add pairwise Hudson F_ST columns."
)
def popgen_windows(
    options: CommonOptions,
    prep,
    project: int | None,
    window: int,
    step: int | None,
    with_fst: bool,
) -> None:
    """Sliding-window S, θ_W, π, Tajima's D per group (+ F_ST)."""
    from popgen.analysis import windows

    if with_fst and len(prep.groups) < 2:
        raise click.ClickException("--fst needs at least two groups (load a panel)")
    rows = windows(prep, window, step, with_fst, project)
    cols = list(rows[0]) if rows else ["chrom", "start", "end", "sites"]
    _emit(
        _header(prep, options, "windows"),
        {"window": window, "step": step or window, "results": rows},
        options.out_format,
        lambda: _emit_rows(rows, cols, options.out_format),
    )
