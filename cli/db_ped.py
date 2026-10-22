"""`vcfclick db ped` — load a pedigree: from a PED file, or from flags.

    vcfclick db ped fam1 fam1.ped
    vcfclick db ped fam1 --proband C1 --proband C2 --father F --mother M

A pedigree may hold several families and several generations. A PED file
describes the whole pedigree, so loading one replaces what is stored. The
flags describe a few people of ONE family; to add or change them without
touching anybody else, pass --merge: it rewrites only the people the flags
name (and only what the flags say about them) and leaves every other
individual as stored. With a PED file, --merge replaces the families the
file lists as wholes and says who it removed. A flag call without --merge
that would delete any individual it does not name fails instead, and says
so, with no exception: the programs that call `db ped` with flags (the
desktop apps) show nothing from a call that succeeds, and an unnamed person
may be an unaffected sibling whose phenotype is evidence the flags cannot
restate.

The flag grammar and its rules are in cli/pedigree_flags.py; the write and
its checks (all done before anything is deleted) are in ingest/pedigree.py.
"""

from __future__ import annotations

from dataclasses import dataclass

import click

from cli.main import _set_db, db
from cli.options import command_options
from cli.pedigree_flags import (
    DEFAULT_FAMILY,
    PedFlags,
    build_rows,
    keep_stored,
)
from cli.pedigree_load import pedigree_rows
from cli.trio_analysis import _sole_ingest_id

# Flags that only make sense without a PED file, for the "not both" check.
_FLAG_ONLY = (
    "--proband",
    "--father",
    "--mother",
    "--family",
    "--parents",
    "--affected",
    "--unaffected",
    "--unknown",
    "--sex",
    "--proband-sex",
)


@dataclass(frozen=True)
class PedigreeOptions:
    name: str
    ped_path: str | None
    ingest_id: str | None
    proband: tuple[str, ...]
    father: str | None
    mother: str | None
    family: str | None
    parents: tuple[str, ...]
    merge: bool
    affected: tuple[str, ...]
    unaffected: tuple[str, ...]
    unknown: tuple[str, ...]
    sex: tuple[str, ...]
    proband_sex: str

    def flags_given(self) -> list[str]:
        """The flag-only options present on the command line."""
        given = {
            "--proband": self.proband,
            "--father": self.father,
            "--mother": self.mother,
            "--family": self.family,
            "--parents": self.parents,
            "--affected": self.affected,
            "--unaffected": self.unaffected,
            "--unknown": self.unknown,
            "--sex": self.sex,
            "--proband-sex": self.proband_sex != "unknown",
        }
        return [flag for flag in _FLAG_ONLY if given[flag]]


def _by_family(rows: list[dict]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for r in rows:
        out.setdefault(r["family_id"], []).append(r["sample_id"])
    return {fam: sorted(members) for fam, members in sorted(out.items())}


def _families_line(rows: list[dict]) -> str:
    return "families: " + ", ".join(
        f"{fam}({len(members)})" for fam, members in _by_family(rows).items()
    )


def _removed_lines(removed: list[dict]) -> list[str]:
    """`removed from famB: K2, XB` for each family that lost individuals."""
    return [
        f"removed from {fam}: {', '.join(members)}"
        for fam, members in _by_family(removed).items()
    ]


def _refuse_to_drop_others(existing: list[dict], rows: list[dict]) -> None:
    """A flag call without --merge replaces the whole pedigree of the
    ingestion. If that would delete any individual the flags do not name,
    stop and say how to do what was probably meant. There is no exception
    for a small or single-family pedigree: the callers that send flags show
    nothing from a call that succeeds, and an unnamed sibling may carry an
    asserted phenotype that the flags cannot say again."""
    named = {r["sample_id"] for r in rows}
    others = [e for e in existing if e["sample_id"] not in named]
    if not others:
        return
    shown = []
    for fam, members in _by_family(others).items():
        head = ", ".join(members[:4]) + (", ..." if len(members) > 4 else "")
        shown.append(f"{fam}: {head}")
    declared = sorted({r["family_id"] for r in rows})
    # The samples this call names may already live in one family: then that
    # family is almost certainly the one meant (the desktop apps never send
    # --family, so the default fam1 would name a family that does not exist).
    homes = sorted({e["family_id"] for e in existing if e["sample_id"] in named})
    target = homes[0] if len(homes) == 1 else None
    if target is not None and declared != [target]:
        merge = f"`--family {target} --merge`"
    else:
        merge = "`--merge`"
    raise click.ClickException(
        f"the stored pedigree has {len(others)} individual(s) that this call "
        f"does not name ({'; '.join(shown)}), and replacing it would delete "
        f"them. Add {merge} to change only the people this call names and "
        "keep everyone else as stored, or load a PED file to replace the "
        "whole pedigree."
    )


def _warn_phenotype_changes(
    existing: list[dict], rows: list[dict], options: PedigreeOptions
) -> None:
    """Say so when a flag call overwrites a stored phenotype that someone
    asserted (affected or unaffected) with a different value the flags did
    not ask for: probands default to affected and everyone else to unknown."""
    asked = {*options.affected, *options.unaffected, *options.unknown}
    old = {e["sample_id"]: e["affected"] for e in existing}
    changes = [
        f"{r['sample_id']} ({old[r['sample_id']]} -> {r['affected'] or 'unknown'})"
        for r in rows
        if r["sample_id"] not in asked
        and old.get(r["sample_id"]) is not None
        and old[r["sample_id"]] != r["affected"]
    ]
    if changes:
        click.echo(
            "warning: replacing the stored phenotype of "
            + ", ".join(changes)
            + "; name them with --affected, --unaffected or --unknown to say "
            "what they are",
            err=True,
        )


@db.command(name="ped")
@click.argument("name")
@click.argument(
    "ped_path", required=False, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--ingest-id",
    default=None,
    help="Ingest_id the pedigree's samples belong to. Inferred when the "
    "database has exactly one ingestion.",
)
@click.option(
    "--proband",
    multiple=True,
    help="Without a PED file: an affected child's sample id. Repeat it for "
    "siblings who share the same parents.",
)
@click.option(
    "--father",
    default=None,
    help="Without a PED file: the proband's father (sample id; 0 when only "
    "the mother is known).",
)
@click.option(
    "--mother",
    default=None,
    help="Without a PED file: the proband's mother (sample id; 0 when only "
    "the father is known).",
)
@click.option(
    "--family",
    default=None,
    help=f"Without a PED file: the family id of every listed sample. "
    f"[default: {DEFAULT_FAMILY}]",
)
@click.option(
    "--parents",
    multiple=True,
    metavar="CHILD=FATHER,MOTHER",
    help="Without a PED file: the parents of another member, to declare "
    "grandparents or half siblings (e.g. F=PGF,PGM). Repeatable; either "
    "side may be 0.",
)
@click.option(
    "--merge",
    is_flag=True,
    help="Keep what this call does not name. With flags: rewrite only the "
    "people the flags name, and only what the flags say about them, and "
    "leave everybody else as stored. With a PED file: replace the families "
    "it lists as wholes (members it omits are removed, and listed) and keep "
    "the other families. Without it a PED file replaces the whole pedigree, "
    "and a flag call that would delete individuals it does not name fails.",
)
@click.option(
    "--affected",
    multiple=True,
    help="Without a PED file: mark a listed sample affected. Repeatable.",
)
@click.option(
    "--unaffected",
    multiple=True,
    help="Without a PED file: mark a listed sample unaffected. Repeatable.",
)
@click.option(
    "--unknown",
    multiple=True,
    help="Without a PED file: mark a listed sample's phenotype unknown. "
    "Repeatable. Probands default to affected, everyone else to unknown.",
)
@click.option(
    "--sex",
    multiple=True,
    metavar="ID=male|female|unknown",
    help="Without a PED file: a listed sample's sex. Repeatable. A father "
    "is male and a mother female unless this says otherwise.",
)
@click.option(
    "--proband-sex",
    type=click.Choice(["male", "female", "unknown"]),
    default="unknown",
    show_default=True,
    help="Without a PED file: the proband's sex (a single --proband only; "
    "use --sex for each sibling, or a PED file).",
)
@command_options(PedigreeOptions)
def db_ped(options: PedigreeOptions) -> None:
    """Load family relationships into NAME, from a PED/FAM file or as a
    family given by sample ids (--proband, repeatable for siblings, with
    --father and --mother; --parents adds grandparents and half siblings).

    The individual ids must match sample ids already ingested under the
    target ingest_id, and all members share that one ingest_id. A PED file
    replaces the stored pedigree for that ingest_id. Flags describe some
    people of one family; add --merge to change just them and keep everybody
    else as stored. Everything is checked before anything is deleted.
    """
    name, ped_path, ingest_id = options.name, options.ped_path, options.ingest_id
    from storage import db_path

    if not db_path(name).exists():
        raise click.ClickException(f"db {name!r} does not exist.")

    _set_db(name)

    if ingest_id is None:
        ingest_id = _sole_ingest_id(name)
        if ingest_id is None:
            raise click.ClickException(
                "database has multiple ingestion; pass --ingest-id to say "
                "which cohort the pedigree applies to."
            )

    given = options.flags_given()
    if ped_path is not None and given:
        raise click.ClickException(
            "give either a PED file or --proband/--father/--mother (with "
            f"--family, --parents and the like), not both (got {', '.join(given)})"
        )

    from ingest.pedigree import load_rows, parse_ped, removed_rows, resulting_rows
    from storage import get_session, table_exists

    if not table_exists("pedigree"):
        raise click.ClickException(
            "this database has no pedigree table (an older vcfclick created "
            "it); create a new database and re-ingest to load a pedigree."
        )
    sess = get_session(name)

    from_flags = ped_path is None
    try:
        # A PED file is the whole pedigree: only a merge or a flag call needs
        # to know what is stored.
        existing = (
            pedigree_rows(sess, ingest_id) if options.merge or from_flags else None
        )
        if from_flags:
            built = build_rows(
                PedFlags(
                    probands=options.proband,
                    father=options.father,
                    mother=options.mother,
                    family=options.family or DEFAULT_FAMILY,
                    parents=options.parents,
                    affected=options.affected,
                    unaffected=options.unaffected,
                    unknown=options.unknown,
                    sex=options.sex,
                    proband_sex=options.proband_sex,
                ),
                stored=existing if options.merge else None,
            )
            rows, warnings = built.rows, built.warnings
            if options.merge:
                rows = keep_stored(built, existing)
            else:
                _refuse_to_drop_others(existing, rows)
            _warn_phenotype_changes(existing, rows, options)
        else:
            rows, warnings = parse_ped(ped_path), []
        for warning in warnings:
            click.echo(f"warning: {warning}", err=True)
        n = load_rows(
            ingest_id,
            rows,
            merge=options.merge,
            existing=existing,
            by_sample=from_flags,
        )
        stored = resulting_rows(
            rows, existing, merge=options.merge, by_sample=from_flags
        )
    except ValueError as e:
        raise click.ClickException(str(e)) from e
    click.echo(f"loaded pedigree: {n} individuals under ingest_id={ingest_id}")
    click.echo(_families_line(stored))
    for line in _removed_lines(removed_rows(existing, stored)):
        click.echo(line)
