"""The pedigree as a graph: families, trios and duos, probands, siblings,
grandparents, generations and expected kinship.

Pure Python: no SQL, no click, no storage session, so `ingest/` (the write
path) and `cli/` (every command that needs family structure) import the
same code, and its tests need no database. `cli/pedigree_load.py` is the
one place that reads the `pedigree` table and builds a `Pedigree`.

Reads are TOLERANT. `Pedigree.from_rows` never raises on a cycle, a
cross-family parent or a self-parent, so a pedigree that loaded under an
older vcfclick still answers every question that does not touch the broken
part. Strict checking is `validate()` / `assert_valid()`, which the write
path and `db families` call.

One definition of "proband" (used by `db families`, `db trio`, the web UI
and, later, the desktop pickers):

    analysable child   a child with both parents declared and genotyped (a
                       trio), or exactly one (a duo)
    proband            an analysable child that is NOT explicitly
                       `unaffected` and that either has no children of its
                       own or is `affected`

`probands()` leaves duos out unless asked (`include_duos=True`): duo
analysis is a separate feature. `all_children()` is the older, wider rule
("every child with both parents declared") that `db trio --all-trios`
keeps.

Individuals are identified by sample id. A parent id that has no row of its
own is allowed (its family and sex are simply unknown). When the same sample
id appears in two rows, the later row wins.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import NamedTuple

SIBLING_KINDS: tuple[str, ...] = ("full", "half_paternal", "half_maternal")
GRANDPARENT_SLOTS: tuple[str, ...] = (
    "paternal_grandfather",
    "paternal_grandmother",
    "maternal_grandfather",
    "maternal_grandmother",
)
LEVELS: tuple[str, ...] = ("error", "warning", "info")


class PedigreeError(ValueError):
    """The pedigree cannot answer: a cyclic component, an unknown id, or a
    strict check (`assert_valid`) that found an error."""


class Problem(NamedTuple):
    """One finding of `Pedigree.validate()`."""

    level: str  # error | warning | info
    family_id: str | None
    code: str
    message: str

    def as_dict(self) -> dict:
        return {
            "family_id": self.family_id,
            "level": self.level,
            "code": self.code,
            "message": self.message,
        }


class Sibling(NamedTuple):
    """A sibling and how it is related: `full`, or `half_paternal` /
    `half_maternal` for the one declared parent they share."""

    sample: str
    kind: str


@dataclass(frozen=True)
class Person:
    sample_id: str
    family_id: str
    father: str | None  # None: no parent declared on that side
    mother: str | None
    sex: str | None  # male | female | None (unknown)
    affected: str | None  # affected | unaffected | None (unknown)


def _parent(value) -> str | None:
    """A declared parent id, or None for the PED founder sentinel ('0'),
    an empty string or NULL."""
    if value is None:
        return None
    text = str(value).strip()
    return None if text in ("", "0") else text


def _choice(value, allowed: tuple[str, ...]) -> str | None:
    return value if value in allowed else None


class Pedigree:
    """An immutable pedigree. Build it with `Pedigree.from_rows`."""

    def __init__(
        self, people: dict[str, Person], genotyped: frozenset[str] | None = None
    ):
        self._people = people
        # None: every member counts as genotyped (nothing is known to be
        # absent from the data).
        self._genotyped = genotyped
        children: dict[str, list[str]] = {}
        for sid in sorted(people):
            p = people[sid]
            for parent in {p.father, p.mother} - {None, sid}:  # not its own child
                children.setdefault(parent, []).append(sid)
        self._children = children
        self._topo_cache: tuple[list[str], frozenset[str]] | None = None
        self._rank_of: dict[str, int] = {}
        self._generation_cache: dict[str, int] | None = None
        self._phi: dict[frozenset[str], float] = {}

    @classmethod
    def from_rows(
        cls, rows: Iterable[Mapping], genotyped: Iterable[str] | None = None
    ) -> Pedigree:
        """Build from pedigree-table rows (dicts with sample_id, family_id,
        father_id, mother_id, sex, affected). '0', '' and NULL parents mean
        "none". `genotyped` is the set of sample ids present in the data;
        None means all of them. No ordering and no cycle check here."""
        people: dict[str, Person] = {}
        for r in rows:
            sid = str(r["sample_id"])
            people[sid] = Person(
                sample_id=sid,
                family_id=str(r.get("family_id") or ""),
                father=_parent(r.get("father_id")),
                mother=_parent(r.get("mother_id")),
                sex=_choice(r.get("sex"), ("male", "female")),
                affected=_choice(r.get("affected"), ("affected", "unaffected")),
            )
        return cls(people, None if genotyped is None else frozenset(genotyped))

    # ───────────────────────────── membership ─────────────────────────────

    def __len__(self) -> int:
        return len(self._people)

    def __contains__(self, sample_id: object) -> bool:
        return sample_id in self._people

    def person(self, sample_id: str) -> Person | None:
        return self._people.get(sample_id)

    def families(self) -> list[str]:
        """Family ids, sorted."""
        return sorted({p.family_id for p in self._people.values()})

    def members(self, family_id: str | None = None) -> list[str]:
        """Sample ids of one family (or of everyone), sorted."""
        return sorted(
            s
            for s, p in self._people.items()
            if family_id is None or p.family_id == family_id
        )

    def family_of(self, sample_id: str) -> str | None:
        p = self._people.get(sample_id)
        return None if p is None else p.family_id

    def is_genotyped(self, sample_id: str) -> bool:
        return self._genotyped is None or sample_id in self._genotyped

    def children_of(self, sample_id: str) -> list[str]:
        """Individuals that name `sample_id` as a parent, sorted."""
        return list(self._children.get(sample_id, ()))

    # ──────────────────── parents, trios, duos, probands ────────────────────

    def parents_of(self, child: str) -> tuple[str | None, str | None]:
        """(father, mother) as declared; None where none is declared."""
        p = self._people.get(child)
        return (None, None) if p is None else (p.father, p.mother)

    def trio_of(self, child: str) -> tuple[str, str] | None:
        """(father, mother) when BOTH parents are declared, else None. This
        is the older rule: it does not look at what is genotyped."""
        father, mother = self.parents_of(child)
        return (father, mother) if father and mother else None

    def kind_of(self, child: str) -> str | None:
        """`trio` (both parents declared and genotyped), `duo` (exactly one
        such parent) or None. The child must be genotyped too."""
        if child not in self._people or not self.is_genotyped(child):
            return None
        usable = [
            p for p in self.parents_of(child) if p is not None and self.is_genotyped(p)
        ]
        return {2: "trio", 1: "duo"}.get(len(usable))

    def all_children(self, family_id: str | None = None) -> list[str]:
        """Every individual with both parents declared (the older rule that
        `db trio --all-trios` keeps), in sample-id order."""
        return [c for c in self.members(family_id) if self.trio_of(c)]

    def analysable_children(self, family_id: str | None = None) -> list[str]:
        """Children that can be analysed: a trio or a duo, sorted."""
        return [c for c in self.members(family_id) if self.kind_of(c)]

    def is_proband(self, child: str, *, include_duos: bool = False) -> bool:
        """Is `child` a proband: analysable (a trio, or a duo with
        `include_duos`), not explicitly unaffected, and either affected or
        without children of their own?"""
        kind = self.kind_of(child)
        if kind is None or (kind == "duo" and not include_duos):
            return False
        affected = self._people[child].affected
        if affected == "unaffected":
            return False
        return affected == "affected" or not self._children.get(child)

    def probands(
        self, family_id: str | None = None, *, include_duos: bool = False
    ) -> list[str]:
        """The probands of one family (or of every family), in sample-id
        order. See the module docstring for the rule."""
        return [
            c
            for c in self.members(family_id)
            if self.is_proband(c, include_duos=include_duos)
        ]

    # ───────────────────────── siblings, grandparents ─────────────────────────

    def siblings(self, sample_id: str) -> list[Sibling]:
        """Everyone who shares a declared parent with `sample_id`, sorted
        by sample id. `full`: both parents shared; otherwise the single
        shared parent decides `half_paternal` or `half_maternal`. Two
        children who each name only the same one parent are half siblings:
        that parent is all that is known to be shared."""
        me = self._people.get(sample_id)
        if me is None:
            return []
        candidates: set[str] = set()
        for parent in (me.father, me.mother):
            if parent is not None:
                candidates.update(self._children.get(parent, ()))
        candidates.discard(sample_id)
        out = []
        for c in sorted(candidates):
            other = self._people[c]
            shares_father = me.father is not None and other.father == me.father
            shares_mother = me.mother is not None and other.mother == me.mother
            if shares_father and shares_mother:
                kind = "full"
            elif shares_father:
                kind = "half_paternal"
            elif shares_mother:
                kind = "half_maternal"
            else:
                continue  # a parent in the other slot: not a sibling
            out.append(Sibling(c, kind))
        return out

    def informative_for(self, sibling: str, parent_id: str) -> bool:
        """Is `sibling` evidence about an allele inherited from `parent_id`?
        Only a child of that parent is: a half sibling who does not carry
        the shared parent's allele says nothing about an allele that came
        from the other parent."""
        return parent_id in self.parents_of(sibling)

    def grandparents(self, sample_id: str) -> dict[str, str | None]:
        """The four grandparent slots (GRANDPARENT_SLOTS), None where
        unknown, not in the pedigree or not genotyped."""
        father, mother = self.parents_of(sample_id)
        out: dict[str, str | None] = {}
        for side, parent in (("paternal", father), ("maternal", mother)):
            gf, gm = self.parents_of(parent) if parent else (None, None)
            for role, g in (("grandfather", gf), ("grandmother", gm)):
                out[f"{side}_{role}"] = (
                    g if g is not None and self.is_genotyped(g) else None
                )
        return out

    # ───────────────────────────── generations ─────────────────────────────

    def _kahn(self, *, ignore_self: bool) -> tuple[list[str], frozenset[str]]:
        """(individuals in an order that puts parents before children, the
        individuals that cannot be ordered). The second set is every member
        on a cycle or descended from one. A self-parent is a cycle of one
        unless `ignore_self`. Parents without a row are always done."""
        waiting = {
            s: {
                x
                for x in (p.father, p.mother)
                if x is not None and x in self._people and not (ignore_self and x == s)
            }
            for s, p in self._people.items()
        }
        ready = deque(sorted(s for s, w in waiting.items() if not w))
        queued = set(ready)
        order: list[str] = []
        while ready:
            node = ready.popleft()
            order.append(node)
            for child in self._children.get(node, ()):
                pending = waiting[child]
                pending.discard(node)
                if not pending and child not in queued:
                    queued.add(child)
                    ready.append(child)
        return order, frozenset(self._people) - frozenset(order)

    def _topological(self) -> tuple[list[str], frozenset[str]]:
        if self._topo_cache is None:
            order, stuck = self._kahn(ignore_self=False)
            self._topo_cache = (order, stuck)
            self._rank_of = {s: i for i, s in enumerate(order)}
        return self._topo_cache

    def generations(self) -> dict[str, int | None]:
        """Generation of every member: 0 for an individual with no declared
        parent, else 1 + the highest generation among its parents (a parent
        without a row counts as 0). None for a member that is on or below a
        cycle. This is depth from the founders of its lineage: a married-in
        founder is 0 even when the spouse is in generation 1."""
        if self._generation_cache is None:
            order, _ = self._topological()
            gen: dict[str, int] = {}
            for s in order:
                p = self._people[s]
                parents = [x for x in (p.father, p.mother) if x is not None]
                gen[s] = 0 if not parents else 1 + max(gen.get(x, 0) for x in parents)
            self._generation_cache = gen
        return {s: self._generation_cache.get(s) for s in sorted(self._people)}

    def generation(self, sample_id: str) -> int:
        """`generations()[sample_id]`; raises PedigreeError for an unknown
        individual or one in a cyclic component."""
        if sample_id not in self._people:
            raise PedigreeError(f"{sample_id!r} is not in the pedigree")
        g = self.generations()[sample_id]
        if g is None:
            raise PedigreeError(
                f"{sample_id!r} is part of, or descended from, a cycle of "
                "ancestors, so it has no generation"
            )
        return g

    # ─────────────────────────────── kinship ───────────────────────────────

    def _require_acyclic(self, *samples: str) -> None:
        _, stuck = self._topological()
        for s in samples:
            if s not in self._people:
                raise PedigreeError(f"{s!r} is not in the pedigree")
            if s in stuck:
                raise PedigreeError(
                    f"{s!r} is part of, or descended from, a cycle of ancestors; "
                    "expected kinship is undefined"
                )

    def _rank(self, sample_id: str) -> int:
        """Position in the parents-before-children order; -1 for a parent
        that has no row (it has no ancestors, so it sorts first)."""
        self._topological()
        return self._rank_of.get(sample_id, -1)

    def _kinship(self, a: str, b: str) -> float:
        key = frozenset((a, b))
        hit = self._phi.get(key)
        if hit is not None:
            return hit
        if a == b:
            father, mother = self.parents_of(a)
            inbreeding = self._kinship(father, mother) if father and mother else 0.0
            value = (1.0 + inbreeding) / 2.0
        else:
            # Expand the later individual: a descendant is always after its
            # ancestors, so b can never be a descendant of the one expanded.
            if self._rank(a) < self._rank(b):
                a, b = b, a
            father, mother = self.parents_of(a)
            value = (
                sum(self._kinship(p, b) for p in (father, mother) if p is not None)
                / 2.0
            )
        self._phi[key] = value
        return value

    def expected_kinship(self, a: str, b: str) -> float:
        """The kinship coefficient the declared pedigree predicts for a and
        b, founders unrelated (the standard recursion
        phi(a, b) = (phi(father(a), b) + phi(mother(a), b)) / 2 and
        phi(a, a) = (1 + phi(father(a), mother(a))) / 2). Parent-child and
        full siblings 0.25; half siblings, grandparent-grandchild and
        aunt/uncle-niece/nephew 0.125; first cousins 0.0625. Raises
        PedigreeError for an unknown individual or a cyclic component."""
        self._require_acyclic(a, b)
        return self._kinship(a, b)

    def inbreeding(self, sample_id: str) -> float:
        """phi(father, mother): 0 unless both are declared and related."""
        self._require_acyclic(sample_id)
        father, mother = self.parents_of(sample_id)
        return self._kinship(father, mother) if father and mother else 0.0

    def is_inbred_pair(self, a: str, b: str) -> bool:
        """Is either individual's own kinship above 0.5, i.e. are its
        parents related? Expected kinship of such a pair is not the
        textbook value for its relationship."""
        return self.inbreeding(a) > 0 or self.inbreeding(b) > 0

    def expected_label(self, a: str, b: str) -> str:
        """The declared relationship, for display: parent-child,
        full-siblings, half-siblings, grandparent-grandchild,
        aunt/uncle-niece/nephew, first-cousins, else `relative` (some
        common ancestor) or `unrelated`."""
        if a == b:
            raise PedigreeError("a relationship needs two different individuals")
        self._require_acyclic(a, b)
        pa, pb = self.parents_of(a), self.parents_of(b)
        if a in pb or b in pa:
            return "parent-child"
        shared = [x for x in pa if x is not None and x in pb]
        if len(shared) == 2:
            return "full-siblings"
        if shared:
            return "half-siblings"
        grand_a = self._ancestors(a, 2) - self._ancestors(a, 1)
        grand_b = self._ancestors(b, 2) - self._ancestors(b, 1)
        if a in grand_b or b in grand_a:
            return "grandparent-grandchild"
        if self._sibling_among(pa, b) or self._sibling_among(pb, a):
            return "aunt/uncle-niece/nephew"
        if any(
            x is not None and y is not None and self._full_siblings(x, y)
            for x in pa
            for y in pb
        ):
            return "first-cousins"
        return "relative" if self._kinship(a, b) > 0 else "unrelated"

    def _ancestors(self, sample_id: str, depth: int) -> set[str]:
        """Ancestors up to `depth` generations back (1: parents)."""
        found: set[str] = set()
        frontier = {sample_id}
        for _ in range(depth):
            nxt = {x for s in frontier for x in self.parents_of(s) if x is not None}
            found |= nxt
            frontier = nxt
        return found

    def _sibling_among(
        self, parents: tuple[str | None, str | None], other: str
    ) -> bool:
        """Is `other` a full or half sibling of one of `parents`?"""
        return any(p is not None and self._shares_parent(p, other) for p in parents)

    def _shares_parent(self, a: str, b: str) -> bool:
        if a == b:
            return False
        pa, pb = self.parents_of(a), self.parents_of(b)
        return any(x is not None and x in pb for x in pa)

    def _full_siblings(self, a: str, b: str) -> bool:
        if a == b:
            return False
        pa, pb = self.parents_of(a), self.parents_of(b)
        return pa[0] is not None and pa[1] is not None and pa == pb

    # ───────────────────────────── validation ─────────────────────────────

    def validate(self) -> list[Problem]:
        """Everything wrong or worth knowing, errors first. Errors: a
        self-parent, the same individual as both parents, a parent in a
        different family, a cycle of ancestors. Warnings: a father recorded
        as female or a mother as male, an individual or declared parent that
        is not a sample in the data. Info: a duo (one declared parent) and
        a parent with no row of its own."""
        out: list[Problem] = []
        for sid in sorted(self._people):
            p = self._people[sid]
            fam = p.family_id
            if p.father is not None and p.father == p.mother:
                out.append(
                    Problem(
                        "error",
                        fam,
                        "same_parent",
                        f"{sid} lists {p.father} as both father and mother",
                    )
                )
            for slot, parent in (("father", p.father), ("mother", p.mother)):
                if parent is None:
                    continue
                if parent == sid:
                    out.append(
                        Problem(
                            "error",
                            fam,
                            "self_parent",
                            f"{sid} is listed as its own {slot}",
                        )
                    )
                    continue
                pp = self._people.get(parent)
                if pp is None:
                    out.append(
                        Problem(
                            "info",
                            fam,
                            "parent_without_row",
                            f"{slot} {parent} of {sid} has no pedigree row of "
                            "its own, so its family and sex are unknown",
                        )
                    )
                    continue
                if pp.family_id != fam:
                    out.append(
                        Problem(
                            "error",
                            fam,
                            "cross_family_parent",
                            f"{slot} {parent} of {sid} is in family "
                            f"{pp.family_id!r}, not {fam!r}",
                        )
                    )
                wrong = {"father": "female", "mother": "male"}[slot]
                if pp.sex == wrong:
                    out.append(
                        Problem(
                            "warning",
                            fam,
                            "parent_sex",
                            f"{slot} {parent} of {sid} is recorded as {wrong}",
                        )
                    )
            declared = [x for x in (p.father, p.mother) if x is not None]
            if len(declared) == 1 and declared[0] != sid:
                out.append(
                    Problem(
                        "info",
                        fam,
                        "duo",
                        f"{sid} has one declared parent ({declared[0]}): a duo, "
                        "analysable from the duo release, so not offered as a "
                        "proband yet",
                    )
                )
            if not self.is_genotyped(sid):
                out.append(
                    Problem(
                        "warning",
                        fam,
                        "not_genotyped",
                        f"{sid} is in the pedigree but is not a sample in the data",
                    )
                )
            for parent in declared:
                if parent in self._people and not self.is_genotyped(parent):
                    out.append(
                        Problem(
                            "warning",
                            fam,
                            "not_genotyped",
                            f"{parent}, a parent of {sid}, is not a sample in the data",
                        )
                    )
        out.extend(self._cycle_problems())
        # one finding per (code, message): a parent of several children
        # would otherwise repeat its not_genotyped warning
        unique = dict.fromkeys(out)
        rank = {lvl: i for i, lvl in enumerate(LEVELS)}
        return sorted(
            unique, key=lambda q: (rank[q.level], q.family_id or "", q.code, q.message)
        )

    def _cycle_problems(self) -> list[Problem]:
        """One `cycle` error per family that has a cycle of ancestors. A
        self-parent is reported as `self_parent`, not as a cycle."""
        _, stuck_set = self._kahn(ignore_self=True)
        stuck = set(stuck_set)
        # keep the members on a cycle: peel off those that only descend
        # from one (nothing below them is stuck)
        changed = True
        while changed:
            changed = False
            for s in sorted(stuck):
                if not any(c in stuck and c != s for c in self._children.get(s, ())):
                    stuck.discard(s)
                    changed = True
        by_family: dict[str, list[str]] = {}
        for s in sorted(stuck):
            by_family.setdefault(self._people[s].family_id, []).append(s)
        return [
            Problem(
                "error",
                fam,
                "cycle",
                "ancestry loops back on itself among " + ", ".join(members),
            )
            for fam, members in sorted(by_family.items())
        ]

    def errors(self) -> list[Problem]:
        return [q for q in self.validate() if q.level == "error"]

    def assert_valid(self) -> None:
        """Raise PedigreeError listing every error-level problem."""
        errors = self.errors()
        if errors:
            raise PedigreeError("; ".join(q.message for q in errors))
