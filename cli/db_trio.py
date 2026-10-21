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
from click.core import ParameterSource

from cli.family_checks import trio_af_gate
from cli.main import _set_db, db
from cli.options import command_options
from cli.trio_analysis import CATEGORY_CHOICES, DEFAULT_MAX_AF
from cli.trio_analysis import Gates as Gates
from cli.trio_analysis import Trio as Trio
from cli.trio_analysis import (
    TrioAnalysis,
)
from cli.trio_analysis import _comphet_genes as _comphet_genes
from cli.trio_analysis import (
    _children,
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
    proband: tuple[str, ...]
    father: str | None
    mother: str | None
    proband_sex: str


@dataclass(frozen=True)
class TrioOptions:
    name: str
    proband: str | None
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
    multiple=True,
    help="Without a PED file: an affected child's sample id. Repeat it for "
    "siblings who share the same parents.",
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
    help="Without a PED file: the proband's sex (a single --proband only; "
    "give each sibling's sex in a PED file).",
)
@command_options(PedigreeOptions)
def db_ped(options: PedigreeOptions) -> None:
    """Load family relationships into NAME, from a PED/FAM file or as a
    family given by sample ids (--proband, repeatable for siblings,
    with --father and --mother).

    The individual ids must match sample ids already ingested under the
    target ingest_id (v1 assumes a joint-called trio, so all members share
    one ingest_id). Re-loading replaces the prior pedigree for that
    ingest_id.
    """
    name, ped_path, ingest_id = options.name, options.ped_path, options.ingest_id
    probands, father, mother = options.proband, options.father, options.mother
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

    names = {"--proband": probands, "--father": father, "--mother": mother}
    if ped_path is None:
        missing = [k for k, v in names.items() if not v]
        if missing:
            raise click.ClickException(
                "give a PED file, or all of --proband, --father and --mother "
                f"(missing: {', '.join(missing)})"
            )
        if len(set(probands)) != len(probands):
            raise click.ClickException("--proband lists the same sample twice")
        if father in probands or mother in probands or father == mother:
            raise click.ClickException(
                "the probands, father and mother must be different samples"
            )
        if proband_sex != "unknown" and len(probands) > 1:
            raise click.ClickException(
                "--proband-sex applies to a single --proband; give each "
                "sibling's sex in a PED file"
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
            tmp.write(f"fam1\t{father}\t0\t0\t1\t1\nfam1\t{mother}\t0\t0\t2\t1\n")
            for child in probands:
                tmp.write(f"fam1\t{child}\t{father}\t{mother}\t{sex}\t2\n")
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
@click.option(
    "--proband",
    default=None,
    help="Sample id of the affected child. Omitted: every child in the "
    "pedigree with both parents declared, one report each.",
)
@click.option(
    "--category",
    type=click.Choice(CATEGORY_CHOICES),
    default="all",
    show_default=True,
    help="Inheritance model. 'all' prints per-model candidate counts.",
)
@click.option("--min-gq", default=20, show_default=True, type=int)
@click.option("--min-dp", default=10, show_default=True, type=int)
@click.option(
    "--max-af",
    default=DEFAULT_MAX_AF,
    show_default=True,
    type=float,
    help="Keep variants with population info_AF <= this (rarity filter). "
    "If INFO/AF turns out to be this cohort's own frequency (GATK-style "
    "family VCFs), the filter is skipped with a note unless this option is "
    "given explicitly.",
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
    "ClinVar when loaded) in one document; without --proband, a `trios` "
    "list with one document per child.",
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
    a trio, with genotype quality gates and an AF rarity filter. Without
    --proband, every child in the pedigree is analysed in turn."""
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
    if options.proband is not None:
        father, mother = _resolve_parents(sess, ingest_id, options.proband)
        families = [(options.proband, father, mother)]
    else:
        families = _children(sess, ingest_id)
        if not families:
            raise click.ClickException(
                "no child with both parents in the pedigree. Load one with "
                "`vcfclick db ped`, or name the proband with --proband."
            )
    has_ref = _has_reference_rows(sess)
    # Is INFO/AF this cohort's own frequency? Then --max-af would remove
    # every private candidate: say so, and skip the gate unless --max-af was
    # given on the command line (an explicit threshold is always honoured).
    explicit_af = (
        click.get_current_context().get_parameter_source("max_af")
        != ParameterSource.DEFAULT
    )
    af_gate, af_message = trio_af_gate(
        sess, ingest_id, options.max_af, explicit=explicit_af
    )
    if af_message:
        click.echo(f"note: {af_message}", err=True)
    gates = options.gates._replace(af_gate=af_gate)
    analyses = [
        TrioAnalysis(
            session=sess,
            trio=Trio(ingest_id, child, father, mother),
            gates=gates,
            has_reference=has_ref,
            gnomad_max_af=options.gnomad_max_af,
        )
        for child, father, mother in families
    ]
    if options.out_format == "json":
        docs = [_trio_json(a, options.category, options.limit) for a in analyses]
        result = docs[0] if options.proband is not None else {"trios": docs}
        click.echo(json.dumps(result, indent=2))
        return
    for i, analysis in enumerate(analyses):
        if i:
            click.echo()
        trio_text(analysis, options.category)
