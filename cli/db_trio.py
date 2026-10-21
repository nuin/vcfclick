"""`vcfclick db trio` — family-based analysis.

A PED file (or `db ped` flags, see cli/db_ped.py) declares parents among
already-ingested samples. Load it, then `db trio` resolves a proband's
parents and reports candidate variants under each Mendelian inheritance
model.

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
