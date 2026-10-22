"""How a `Pedigree` is reported: the `db families` document and text, and the
`families` list the web UI serves. Pure (a `Pedigree` in, plain data out).

`db families --format json` is the one place other programs (the desktop
apps, scripts) get family structure from, so its shape is a contract:

    {"ingest_id": "i1",
     "families": [{"family_id": "famA",
                   "members": [{"sample_id", "family_id", "father", "mother",
                                "sex", "affected", "genotyped",
                                "generation"}],
                   "probands": ["C1", ...],
                   "trios": [{"child", "father", "mother", "kind",
                              "is_proband",
                              "siblings": [{"sample", "kind"}]}]}],
     "problems": [{"family_id", "level", "code", "message"}]}

father / mother / sex / affected / generation are null when unknown (a
cycle gives no generation). `kind` is `trio` or `duo`. Additive only, like
the trio document (docs/TRIO.md).
"""

from __future__ import annotations

from storage.pedigree_graph import Pedigree

KIND_LABEL = {
    "full": "full",
    "half_paternal": "half, paternal",
    "half_maternal": "half, maternal",
}


def _siblings(graph: Pedigree, child: str) -> list[dict]:
    return [{"sample": s.sample, "kind": s.kind} for s in graph.siblings(child)]


def _family(graph: Pedigree, family_id: str) -> dict:
    generations = graph.generations()
    members = []
    for sid in graph.members(family_id):
        p = graph.person(sid)
        members.append(
            {
                "sample_id": sid,
                "family_id": family_id,
                "father": p.father,
                "mother": p.mother,
                "sex": p.sex,
                "affected": p.affected,
                "genotyped": graph.is_genotyped(sid),
                "generation": generations[sid],
            }
        )
    # founders first; a member in a cycle (no generation) last
    members.sort(
        key=lambda m: (m["generation"] is None, m["generation"] or 0, m["sample_id"])
    )
    trios = []
    for child in graph.analysable_children(family_id):
        father, mother = graph.parents_of(child)
        trios.append(
            {
                "child": child,
                "father": father,
                "mother": mother,
                "kind": graph.kind_of(child),
                "is_proband": graph.is_proband(child),
                "siblings": _siblings(graph, child),
            }
        )
    return {
        "family_id": family_id,
        "members": members,
        "probands": graph.probands(family_id),
        "trios": trios,
    }


def families_document(graph: Pedigree, ingest_id: str) -> dict:
    """The `db families --format json` document."""
    return {
        "ingest_id": ingest_id,
        "families": [_family(graph, f) for f in graph.families()],
        "problems": [q.as_dict() for q in graph.validate()],
    }


def web_families(graph: Pedigree) -> list[dict]:
    """The `families` list of /api/meta: per family its probands and every
    analysable child (trio or duo) with its parents, its phenotype and its
    siblings. Same rule and same words as `db families`."""
    out = []
    for fam in graph.families():
        trios = []
        for child in graph.analysable_children(fam):
            father, mother = graph.parents_of(child)
            trios.append(
                {
                    "child": child,
                    "father": father,
                    "mother": mother,
                    "kind": graph.kind_of(child),
                    "affected": graph.person(child).affected,
                    "siblings": _siblings(graph, child),
                }
            )
        out.append({"family_id": fam, "probands": graph.probands(fam), "trios": trios})
    return out


def _why_not_proband(graph: Pedigree, child: str) -> str:
    if graph.kind_of(child) == "duo":
        return "one parent: a duo, analysable from the duo release"
    if graph.person(child).affected == "unaffected":
        return "unaffected"
    return "has children and is not affected"


def _count(n: int, one: str, many: str) -> str:
    return f"{n} {one if n == 1 else many}"


def families_text(graph: Pedigree, ingest_id: str) -> str:
    """The human-readable `db families`."""
    fams = graph.families()
    lines = [
        f"ingest_id {ingest_id}: {_count(len(fams), 'family', 'families')}, "
        f"{_count(len(graph), 'individual', 'individuals')}"
    ]
    for fam in fams:
        members = graph.members(fam)
        generations = graph.generations()
        lines += [
            "",
            f"family {fam}: {_count(len(members), 'individual', 'individuals')}",
        ]
        by_gen: dict[object, list[str]] = {}
        for sid in members:
            by_gen.setdefault(generations[sid], []).append(sid)
        for gen in sorted(by_gen, key=lambda g: (g is None, g or 0)):
            label = "in or below a cycle" if gen is None else f"generation {gen}"
            lines.append(f"  {label}: {', '.join(by_gen[gen])}")
        probands = graph.probands(fam)
        lines.append(f"  probands: {', '.join(probands) if probands else '(none)'}")
        for child in graph.analysable_children(fam):
            father, mother = graph.parents_of(child)
            parents = ", ".join(
                f"{role} {who}"
                for role, who in (("father", father), ("mother", mother))
                if who is not None
            )
            status = (
                "proband"
                if graph.is_proband(child)
                else f"not a proband: {_why_not_proband(graph, child)}"
            )
            sibs = graph.siblings(child)
            sib_text = (
                "; siblings "
                + ", ".join(f"{s.sample} ({KIND_LABEL[s.kind]})" for s in sibs)
                if sibs
                else ""
            )
            lines.append(
                f"  {graph.kind_of(child)} {child}: {parents} [{status}]{sib_text}"
            )
    problems = graph.validate()
    if problems:
        lines += ["", "problems:"]
        for q in problems:
            where = f"{q.family_id}: " if q.family_id else ""
            lines.append(f"  {q.level:7s} {where}{q.message}")
    return "\n".join(lines)
