"""The presets that ship with vcfclick, compiled in.

These are the AUTHOR'S STARTING POINTS, not validated clinical standards:
every one carries `status: starting_point`, and a clinical geneticist is
expected to review the numbers before anyone relies on them. They are
read-only (copy one with `vcfclick preset save NAME --from-preset BUILTIN`).

Governance. A built-in changes only when a step adds a model or a filter that
belongs in it; anything else is a NEW name, so one name never silently means
different filters. A change bumps `revision`, appends to `history` (with the
new fingerprint) and is a CHANGELOG entry; `tests/test_builtin_presets.py`
pins every fingerprint, so an accidental edit fails CI, and checks that a
built-in's category list omits a registered model only when the model is in
`excludes` with a reason. A step that registers a model (xlinked, mosaic, ...)
therefore fails that test until it has decided what each built-in does about
it.

No built-in carries genes, regions, HPO terms or a knob that does nothing.
Common to all of them before consequence filtering exists: there is NO
consequence or impact filter, so non-coding variants are included.
"""

from __future__ import annotations

from dataclasses import dataclass, field

STARTING_POINT = "starting_point"
_NO_CONSEQUENCE = (
    "There is no consequence or impact filter, so non-coding variants are included."
)
_NO_X = (
    "chrX gets no special treatment (no hemizygous male logic), and UPD and "
    "mosaic variants are not looked for."
)


@dataclass(frozen=True)
class BuiltinPreset:
    name: str
    title: str
    description: str
    revision: int
    covers: str
    does_not_cover: str
    filters: dict
    # models this preset deliberately does not run, with the reason
    excludes: dict[str, str] = field(default_factory=dict)
    # [{"revision", "fingerprint", "note"}], oldest first; the last entry is
    # the current revision
    history: tuple[dict, ...] = ()
    applies_to: tuple[str, ...] = ("trio",)
    status: str = STARTING_POINT

    def document(self) -> dict:
        """The preset as a file envelope (the shape a user preset has)."""
        return {
            "schema": "vcfclick.preset",
            "schema_version": 1,
            "name": self.name,
            "title": self.title,
            "description": self.description,
            "applies_to": list(self.applies_to),
            "filters": dict(self.filters),
            "meta": {},
        }


_ALL_MODELS = ["all"]

BUILTINS: dict[str, BuiltinPreset] = {
    p.name: p
    for p in (
        BuiltinPreset(
            name="broad-first-pass",
            title="Broad first pass",
            description=(
                "Every model, loose quality gates and a 1% rarity ceiling: a "
                "first look at a proband. Warn-only: a gate that cannot be "
                "enforced in this data is reported, not fatal."
            ),
            revision=1,
            covers=(
                "Every inheritance model vcfclick runs on a trio (de novo, "
                "recessive, dominant, compound het) at GQ >= 15, DP >= 8, "
                "allele balance 0.20-0.80 for proband hets, and a population "
                "frequency of at most 1% on INFO/AF and on gnomAD popmax."
            ),
            does_not_cover=(
                f"{_NO_CONSEQUENCE} {_NO_X} A variant absent from the loaded "
                "gnomAD slice is kept as rare. Alleles just under 1% are "
                "kept; common alleles are not."
            ),
            filters={
                "category": _ALL_MODELS,
                "min_gq": 15,
                "min_dp": 8,
                "min_ab": 0.20,
                "max_ab": 0.80,
                "max_af": 0.01,
                "gnomad_max_af": 0.01,
            },
            history=(
                {
                    "revision": 1,
                    "fingerprint": "cd3726cc3853",
                    "note": "first release",
                },
            ),
        ),
        BuiltinPreset(
            name="strict-denovo",
            title="Strict de novo",
            description=(
                "De novo calls in the proband with quality and frequency "
                "gates that must bind, FILTER=PASS, and parents whose stored "
                "0/0 is checked against their alt reads (a parent 0/0 with "
                ">= 2 alt reads and >= 5% alt fraction is excluded). Refuses "
                "to run when a gate cannot be enforced in the data "
                "(--no-require-bound-gates overrides for one run)."
            ),
            revision=1,
            covers=(
                "De novo single-nucleotide and indel calls in the proband "
                "where both parents have a stored 0/0 with under 5% alt "
                "reads, FILTER=PASS, GQ >= 30, DP >= 20, allele balance "
                "0.30-0.70, INFO/AF <= 0.001 and gnomAD popmax <= 0.0001. A "
                "missing GQ, DP or AD counts as a fail."
            ),
            does_not_cover=(
                "Inherited variants, and de novo + inherited compound-het "
                f"pairs. {_NO_X} Mapping artefacts beyond the FILTER column "
                f"are not excluded. {_NO_CONSEQUENCE}"
            ),
            filters={
                "category": ["denovo"],
                "min_gq": 30,
                "min_dp": 20,
                "min_ab": 0.30,
                "max_ab": 0.70,
                "max_af": 0.001,
                "gnomad_max_af": 0.0001,
                "strict_missing": True,
                "require_bound": True,
                "pass_only": True,
                "parent_alt_screen": True,
            },
            excludes={
                "recessive": "this preset is about new variants in the proband",
                "dominant": "this preset is about new variants in the proband",
                "comphet": "this preset is about new variants in the proband",
            },
            history=(
                {
                    "revision": 1,
                    "fingerprint": "30ffc72499a2",
                    "note": "first release",
                },
            ),
        ),
        BuiltinPreset(
            name="rare-recessive",
            title="Rare recessive",
            description=(
                "Homozygous and compound-heterozygous candidates with one 1% "
                "frequency ceiling for INFO/AF and gnomAD. 1% is the usual "
                "carrier-frequency ceiling for a recessive disorder, but "
                "alleles near 1% popmax (some CFTR and GJB2 alleles; verify "
                "against the gnomAD release in use) can still be hidden: for "
                "genes with a known carrier frequency start from "
                "broad-first-pass with a higher --gnomad-max-af."
            ),
            revision=1,
            covers=(
                "A homozygous-alt proband with both parents het (recessive) "
                "and genes with a rare het inherited from each parent "
                "(compound het), at GQ >= 20, DP >= 10, allele balance "
                "0.25-0.75 and a population frequency of at most 1%."
            ),
            does_not_cover=(
                "chrX in males, de novo + inherited compound-het pairs (the "
                "compound-het model takes only hets inherited from exactly "
                "one parent), UPD, mosaic; the allele balance of a "
                f"hom-alt call is not gated. {_NO_CONSEQUENCE}"
            ),
            filters={
                "category": ["recessive", "comphet"],
                "min_gq": 20,
                "min_dp": 10,
                "min_ab": 0.25,
                "max_ab": 0.75,
                "max_af": 0.01,
                "gnomad_max_af": 0.01,
            },
            excludes={
                "denovo": "a recessive search; use strict-denovo for new variants",
                "dominant": "a recessive search; use rare-dominant for dominant",
            },
            history=(
                {
                    "revision": 1,
                    "fingerprint": "c19bb0a98ece",
                    "note": "first release",
                },
            ),
        ),
        BuiltinPreset(
            name="rare-dominant",
            title="Rare dominant (inherited or de novo; parent phenotype not considered)",
            description=(
                "Heterozygous candidates that are new in the proband or "
                "inherited from exactly one parent, at a rarity of 0.1% "
                "(INFO/AF) and 0.01% (gnomAD)."
            ),
            revision=1,
            covers=(
                "A heterozygous proband call that is de novo, or that exactly "
                "one parent carries while the other has a stored 0/0, at GQ "
                ">= 25, DP >= 12, allele balance 0.25-0.75, INFO/AF <= 0.001 "
                "and gnomAD popmax <= 0.0001."
            ),
            does_not_cover=(
                "The dominant model keeps a variant inherited from EITHER "
                "parent: whether the carrying parent is affected is not "
                "considered, so it also lists inherited variants a healthy "
                f"parent carries. {_NO_X} {_NO_CONSEQUENCE}"
            ),
            filters={
                "category": ["dominant", "denovo"],
                "min_gq": 25,
                "min_dp": 12,
                "min_ab": 0.25,
                "max_ab": 0.75,
                "max_af": 0.001,
                "gnomad_max_af": 0.0001,
            },
            excludes={
                "recessive": "a dominant search; use rare-recessive for recessive",
                "comphet": "a dominant search; use rare-recessive for compound het",
            },
            history=(
                {
                    "revision": 1,
                    "fingerprint": "466b94313b5c",
                    "note": "first release",
                },
            ),
        ),
    )
}
