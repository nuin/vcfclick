"""`vcfclick` — small VCF databases, one per cohort.

Usage:
    vcfclick db create <name>
    vcfclick db list
    vcfclick db ingest <name> <vcf_path>
    vcfclick db query <name> "<sql>"
    vcfclick db info <name>
    vcfclick db dump <name> [--out <dir>]
    vcfclick db rm <name>

Each named DB is a self-contained chDB session at
$VCFCLICK_HOME/dbs/<name>/ (default ~/.vcfclick/dbs/<name>/). Setting the
VCFCLICK_DB_NAME env var inside a command propagates the choice to
ingest worker subprocesses and to the MCP server.
"""

from __future__ import annotations

import logging
import os
import sys

import click


def _set_db(name: str) -> None:
    """Make the named DB the active target for storage + subprocess workers."""
    os.environ["VCFCLICK_DB_NAME"] = name


def _setup_logging() -> None:
    """Surface library log records as plain stderr lines.

    Library modules (ingest, annotations, export) emit progress via the
    standard ``logging`` module. The CLI is the only intended consumer,
    so configure the root logger here with no timestamp / level prefix —
    matches the pre-logging UX where each step printed one line.
    """
    root = logging.getLogger()
    if root.handlers:
        return
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(message)s"))
    root.addHandler(handler)
    root.setLevel(logging.INFO)


@click.group()
@click.version_option(package_name="vcfclick")
def cli() -> None:
    """vcfclick — small VCF databases."""
    _setup_logging()


@cli.group()
def db() -> None:
    """Manage named VCF databases."""


@cli.group()
def annotations() -> None:
    """Manage the embedded annotation reference store (DuckDB).

    Annotations are reference data shared across all named databases:
    gene coordinates today (GENCODE), transcript / exon / CDS depth and
    ClinVar in Phase 2. They live at `annotations/annotations.duckdb`
    inside the installed package directory.
    """


# Importing the submodules registers their @db.command / @annotations.command
# functions against the groups defined above. importlib (not a `from … import`)
# keeps this a pure side-effect import — no unused symbol bound at module scope.
import importlib  # noqa: E402

importlib.import_module("cli.db")
importlib.import_module("cli.annotations")
importlib.import_module("cli.discover")


if __name__ == "__main__":
    cli()
