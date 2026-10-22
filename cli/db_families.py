"""`vcfclick db families` — the families in a database's pedigree (read-only).

    vcfclick db families fam1
    vcfclick db families fam1 --format json

Lists each family's members by generation, the probands (the one definition
`db trio` also uses), every trio and duo with its siblings, and anything
wrong or worth knowing about the pedigree. It reads tolerantly: a cycle or a
parent in another family that loaded earlier is reported as a problem, never
an error, and the command exits 0. `--format json` is what other programs
read; its shape is in cli/pedigree_report.py.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import click

from cli.main import _set_db, db
from cli.options import command_options
from cli.pedigree_load import NoPedigree, load_pedigree_graph
from cli.pedigree_report import families_document, families_text
from cli.trio_analysis import _sole_ingest_id
from storage.pedigree_graph import Pedigree


@dataclass(frozen=True)
class FamiliesOptions:
    name: str
    ingest_id: str | None
    out_format: str


@db.command(name="families")
@click.argument("name")
@click.option(
    "--ingest-id",
    default=None,
    help="Ingest_id whose pedigree to list. Inferred when the database has "
    "exactly one ingestion.",
)
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
    help="json: {ingest_id, families: [{family_id, members, probands, "
    "trios}], problems}.",
)
@command_options(FamiliesOptions)
def db_families(options: FamiliesOptions) -> None:
    """List the families in NAME's pedigree: members by generation, the
    probands `db trio` analyses, each trio or duo with its siblings, and
    problems in the pedigree. Read-only; exits 0 even when it finds
    problems."""
    from storage import db_path, get_session

    if not db_path(options.name).exists():
        raise click.ClickException(f"db {options.name!r} does not exist.")
    _set_db(options.name)
    ingest_id = options.ingest_id or _sole_ingest_id(options.name)
    if ingest_id is None:
        raise click.ClickException(
            "database has no single ingestion; pass --ingest-id to say "
            "which cohort's pedigree to list."
        )
    graph = Pedigree.from_rows([])
    try:
        graph = load_pedigree_graph(get_session(options.name), ingest_id, require=True)
        message = None
    except NoPedigree as e:
        message = e.message
    if options.out_format == "json":
        if message:
            click.echo(f"note: {message}", err=True)
        click.echo(json.dumps(families_document(graph, ingest_id), indent=2))
    elif message:
        click.echo(f"ingest_id {ingest_id}: {message}")
    else:
        click.echo(families_text(graph, ingest_id))
