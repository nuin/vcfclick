"""Read the `pedigree` table into a `Pedigree` (storage/pedigree_graph.py).

This is the one place that SELECTs from `pedigree` for family structure:
`db families`, `db trio` and the web UI get theirs from here, and so will
every later feature. (`cli/trio_analysis._resolve_parents` keeps its own,
narrower query on purpose: it must keep resolving a proband the way it
always has, whatever else is wrong with the pedigree.)

Two things the table does not promise, handled here once:

  * The table may not exist (a database that predates it): `NoPedigree`
    says so in one sentence, instead of a SQL error.
  * `pedigree` is a ReplacingMergeTree on chDB, so until a background merge
    runs, a re-inserted sample can be read twice. Rows are collapsed to one
    per (ingest_id, sample_id), the most recently inserted winning.
"""

from __future__ import annotations

import json

import click

from storage.pedigree_graph import Pedigree

NO_PEDIGREE = (
    "no pedigree loaded. Load one with `vcfclick db ped <name> <file.ped>` "
    "(or `--proband/--father/--mother`)."
)

COLUMNS = ("sample_id", "family_id", "father_id", "mother_id", "sex", "affected")


class NoPedigree(click.ClickException):
    """The database has no pedigree to read."""


def _query_rows(sess, sql: str) -> list:
    return json.loads(sess.query(sql, "JSONCompact").bytes().decode())["data"]


def pedigree_rows(sess, ingest_id: str) -> list[dict]:
    """The pedigree rows of one ingestion, one per sample, in sample-id
    order, as dicts with the COLUMNS keys. Raises NoPedigree when the
    database has no `pedigree` table; an existing but empty table gives
    an empty list."""
    from storage import sql_quote_str, table_exists

    if not table_exists("pedigree"):
        raise NoPedigree(NO_PEDIGREE)
    data = _query_rows(
        sess,
        f"SELECT {', '.join(COLUMNS)} FROM pedigree "
        f"WHERE ingest_id = {sql_quote_str(ingest_id)} "
        "ORDER BY ingested_at, sample_id",
    )
    latest = {row[0]: dict(zip(COLUMNS, row, strict=True)) for row in data}
    return [latest[s] for s in sorted(latest)]


def known_samples(sess, ingest_id: str) -> set[str]:
    """The sample ids ingested under `ingest_id`."""
    from storage import sql_quote_str

    data = _query_rows(
        sess,
        "SELECT DISTINCT sample_id FROM samples "
        f"WHERE ingest_id = {sql_quote_str(ingest_id)}",
    )
    return {str(r[0]) for r in data}


def load_pedigree_graph(sess, ingest_id: str, *, require: bool = False) -> Pedigree:
    """The `Pedigree` of one ingestion, with `genotyped` set to the samples
    that were ingested. Tolerant: whatever loaded is returned, cycles and
    cross-family parents included (use `Pedigree.validate()` to see them).
    With `require`, a pedigree with no rows raises NoPedigree too."""
    rows = pedigree_rows(sess, ingest_id)
    if require and not rows:
        raise NoPedigree(NO_PEDIGREE)
    return Pedigree.from_rows(rows, genotyped=known_samples(sess, ingest_id))
