"""Which trios `db trio` analyses when no `--proband` is named.

Pure (a `Pedigree` in, a selection and notes out), so each rule has a unit
test. The rule is the one `db families` and the web UI use: a proband is an
analysable child that is not explicitly unaffected and either has no
children of its own or is affected (storage/pedigree_graph.py). In a
one-generation family whose children are not marked unaffected that is every
child, as it always was. `--all-trios` keeps the older rule: every child
with both parents declared, whatever its generation or phenotype.
"""

from __future__ import annotations

from dataclasses import dataclass

import click

from storage.pedigree_graph import Pedigree


@dataclass(frozen=True)
class Pick:
    """One trio to analyse."""

    family_id: str | None
    child: str
    father: str
    mother: str
    siblings: tuple[tuple[str, str], ...] = ()


@dataclass(frozen=True)
class Selection:
    picks: list[Pick]
    notes: list[str]


def make_pick(graph: Pedigree, child: str, father: str, mother: str) -> Pick:
    """`child` with the given parents, labelled with its family and its
    siblings from the pedigree (None and none when it is not in it)."""
    return Pick(
        graph.family_of(child),
        child,
        father,
        mother,
        tuple((s.sample, s.kind) for s in graph.siblings(child)),
    )


def pick_for(graph: Pedigree, child: str) -> Pick | None:
    """The trio of `child` as the pedigree declares it, or None if `child`
    is not in the graph or lacks a parent."""
    parents = graph.trio_of(child)
    return None if parents is None else make_pick(graph, child, *parents)


def _names(samples: list[str]) -> str:
    return ", ".join(samples)


def _why_skipped(graph: Pedigree, child: str) -> str:
    if graph.kind_of(child) != "trio":
        return "a parent is not a sample in the data"
    if graph.person(child).affected == "unaffected":
        return "unaffected"
    return "has children and is not affected"


def select_trios(
    graph: Pedigree, *, family: str | None = None, all_trios: bool = False
) -> Selection:
    """The trios to analyse, ordered by (family, sample id), plus the notes
    for stderr. Raises ClickException when there is nothing to analyse or
    the family is unknown."""
    if family is not None and family not in graph.families():
        known = ", ".join(graph.families()) or "(none: no pedigree loaded)"
        raise click.ClickException(
            f"no family {family!r} in the pedigree. Families: {known}."
        )
    scope = f"family {family!r}" if family else "the pedigree"
    trios = graph.all_children(family)
    if not trios:
        raise click.ClickException(
            f"no child with both parents in {scope}. Load one with "
            "`vcfclick db ped`, or name the proband with --proband."
        )
    chosen = trios if all_trios else graph.probands(family)
    notes: list[str] = []

    duos = [c for c in graph.analysable_children(family) if graph.kind_of(c) == "duo"]
    if duos:
        notes.append(
            f"skipping {_names(duos)}: one declared parent (a duo); duo "
            "analysis is not available in this release."
        )
    if not all_trios:
        skipped: dict[str, list[str]] = {}
        for child in trios:
            if child not in chosen:
                skipped.setdefault(_why_skipped(graph, child), []).append(child)
        if not chosen:
            detail = "; ".join(f"{_names(c)} ({why})" for why, c in skipped.items())
            raise click.ClickException(
                f"no proband in {scope}: every child with both parents "
                f"declared is left out ({detail}). Name one with --proband, "
                "or pass --all-trios to analyse every child."
            )
        if skipped:
            detail = "; ".join(f"{_names(c)} ({why})" for why, c in skipped.items())
            notes.append(
                f"not analysed as probands: {detail}. Pass --all-trios to "
                "analyse every child."
            )
    picks = [p for p in (pick_for(graph, c) for c in chosen) if p is not None]
    picks.sort(key=lambda p: (p.family_id or "", p.child))
    return Selection(picks, notes)
