"""`vcfclick db ped` (and, later, `db trio`) — family-based analysis.

A PED file declares who-is-whose-parent among already-ingested samples.
Load it, then trio analysis can resolve a proband's parents and compare
their genotypes at each site.

    vcfclick db ped fam1 fam1.ped
"""

from __future__ import annotations

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
