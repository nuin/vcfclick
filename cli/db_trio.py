"""`vcfclick db ped` + `db trio` — family-based analysis.

A PED file declares parents among already-ingested samples.
Load it, then `db trio` resolves a proband's parents and reports
candidate variants under each Mendelian inheritance model.

    vcfclick db ped  fam1 fam1.ped
    vcfclick db trio fam1 --proband CHILD

Inheritance models, with slivar-style genotype quality (GQ, depth,
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
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import click

from cli.main import _set_db, db
from cli.options import command_options
from cli.trio_analysis import Gates as Gates
from cli.trio_analysis import Trio as Trio
from cli.trio_analysis import (
    TrioAnalysis,
)
from cli.trio_analysis import _comphet_genes as _comphet_genes
from cli.trio_analysis import (
    _has_reference_rows,
    _resolve_parents,
    _sole_ingest_id,
)
from cli.trio_analysis import trio_sql as trio_sql
from cli.trio_output import _trio_json, trio_text


@dataclass(frozen=True)
class PedigreeOptions:
    name: str
    ped_path: str | None
    ingest_id: str | None
    proband: str | None
    father: str | None
    mother: str | None
    proband_sex: str


@dataclass(frozen=True)
class TrioOptions:
    name: str
    proband: str
    category: str
    min_gq: int
    min_dp: int
    max_af: float
    min_ab: float
    max_ab: float
    gnomad_max_af: float | None
    out_format: str
    limit: int

    @property
    def gates(self) -> Gates:
        return Gates(self.min_gq, self.min_dp, self.max_af, self.min_ab, self.max_ab)


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
    default=None,
    help="Without a PED file: the affected child's sample id.",
)
@click.option(
    "--father", default=None, help="Without a PED file: the father's sample id."
)
@click.option(
    "--mother", default=None, help="Without a PED file: the mother's sample id."
)
@click.option(
    "--proband-sex",
    type=click.Choice(["male", "female", "unknown"]),
    default="unknown",
    show_default=True,
    help="Without a PED file: the proband's sex.",
)
@command_options(PedigreeOptions)
def db_ped(options: PedigreeOptions) -> None:
    """Load family relationships into NAME, from a PED/FAM file or as a
    trio given by sample ids (--proband/--father/--mother).

    The individual ids must match sample ids already ingested under the
    target ingest_id (v1 assumes a joint-called trio, so all members share
    one ingest_id). Re-loading replaces the prior pedigree for that
    ingest_id.
    """
    name, ped_path, ingest_id = options.name, options.ped_path, options.ingest_id
    proband, father, mother = options.proband, options.father, options.mother
    proband_sex = options.proband_sex
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

    from ingest.pedigree import load_pedigree

    names = {"--proband": proband, "--father": father, "--mother": mother}
    if ped_path is None:
        missing = [k for k, v in names.items() if not v]
        if missing:
            raise click.ClickException(
                "give a PED file, or all of --proband, --father and --mother "
                f"(missing: {', '.join(missing)})"
            )
    elif any(names.values()):
        raise click.ClickException(
            "give either a PED file or --proband/--father/--mother, not both"
        )

    import tempfile

    tmp = None
    try:
        if ped_path is None:
            sex = {"male": "1", "female": "2", "unknown": "0"}[proband_sex]
            tmp = tempfile.NamedTemporaryFile("w", suffix=".ped", delete=False)
            tmp.write(
                f"fam1\t{father}\t0\t0\t1\t1\n"
                f"fam1\t{mother}\t0\t0\t2\t1\n"
                f"fam1\t{proband}\t{father}\t{mother}\t{sex}\t2\n"
            )
            tmp.close()
            ped_path = tmp.name
        n = load_pedigree(ingest_id, ped_path)
    except ValueError as e:
        raise click.ClickException(str(e)) from e
    finally:
        if tmp is not None:
            Path(tmp.name).unlink(missing_ok=True)
    click.echo(f"loaded pedigree: {n} individuals under ingest_id={ingest_id}")


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
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
    help="json: every model's count and candidates (with gene, gnomAD and "
    "ClinVar when loaded) in one document.",
)
@click.option(
    "--limit",
    type=int,
    default=1000,
    show_default=True,
    help="json: at most this many candidates per model.",
)
@command_options(TrioOptions)
def db_trio(options: TrioOptions) -> None:
    """Report candidate variants under Mendelian inheritance models for
    a trio, with genotype quality gates and an AF rarity filter."""
    from storage import db_path, get_session

    if not db_path(options.name).exists():
        raise click.ClickException(f"db {options.name!r} does not exist.")
    _set_db(options.name)
    sess = get_session(options.name)

    ingest_id = _sole_ingest_id(options.name)
    if ingest_id is None:
        raise click.ClickException(
            "database has multiple ingestions; trio analysis assumes one "
            "joint-called cohort. Re-ingest the trio as a single VCF "
            "(see `vcfclick merge`)."
        )
    father, mother = _resolve_parents(sess, ingest_id, options.proband)
    analysis = TrioAnalysis(
        session=sess,
        trio=Trio(ingest_id, options.proband, father, mother),
        gates=options.gates,
        has_reference=_has_reference_rows(sess),
        gnomad_max_af=options.gnomad_max_af,
    )
    if options.out_format == "json":
        result = _trio_json(analysis, options.category, options.limit)
        click.echo(json.dumps(result, indent=2))
    else:
        trio_text(analysis, options.category)
