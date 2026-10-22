"""`db ped` without a PED file: the flags, turned into pedigree rows.

Pure (no database, no click context), so every rule has a unit test. The
command (cli/db_ped.py) reads the existing pedigree, hands the rows to
`ingest.pedigree.load_rows` and decides what a replace may delete.

    db ped NAME --proband C1 --proband C2 --father F --mother M
    db ped NAME --family famB --merge --proband K1 --father FB --mother MB \\
               --parents FB=GF,GM --unaffected K2

What the flags declare:

  --proband S      an affected child of --father and --mother (repeatable)
  --father F       the proband's father, or 0 when only the mother is known
  --mother M       the proband's mother, or 0 when only the father is known
  --parents C=F,M  the parents of any other member (repeatable); either side
                   may be 0, both may not
  --family FID     the family every listed sample joins (default fam1)
  --affected / --unaffected / --unknown S
                   set the phenotype of any listed sample (repeatable)
  --sex ID=male|female|unknown
                   set the sex of any listed sample (repeatable)
  --proband-sex    the sex of a single proband

Phenotypes: probands are affected unless told otherwise; EVERYONE ELSE is
unknown (stored as NULL), parents included. A parent's status is an
assertion (it decides reduced penetrance), so only a PED file or an
explicit --affected/--unaffected makes one. Sex: a father is male and a
mother female unless --sex says otherwise (which is reported as a warning).

--merge changes only what the flags state. `keep_stored` lays the rows the
flags build over the stored ones: a person the flags name keeps their
stored parents, sex and phenotype where the flags say nothing about them
(a father named only as --father keeps his own parents), and everybody the
flags do not name is not in the rows at all. Without --merge the rows
replace the pedigree, so nothing is carried over.
"""

from __future__ import annotations

from dataclasses import dataclass, field

DEFAULT_FAMILY = "fam1"
PARENTS, SEX, AFFECTED = "parents", "sex", "affected"
SEX_CHOICES = ("male", "female", "unknown")
NO_PARENT = "0"

MISSING_HINT = (
    "give a PED file, or all of --proband, --father and --mother "
    "(missing: {missing}). Use 0 for a parent that is not known, e.g. "
    "--father 0"
)


@dataclass(frozen=True)
class PedFlags:
    """The `db ped` flags that declare a family without a PED file."""

    probands: tuple[str, ...] = ()
    father: str | None = None
    mother: str | None = None
    family: str = DEFAULT_FAMILY
    parents: tuple[str, ...] = ()  # raw "CHILD=FATHER,MOTHER" specs
    affected: tuple[str, ...] = ()
    unaffected: tuple[str, ...] = ()
    unknown: tuple[str, ...] = ()
    sex: tuple[str, ...] = ()  # raw "ID=male|female|unknown" specs
    proband_sex: str = "unknown"


@dataclass
class _Member:
    father: str | None = None
    mother: str | None = None
    sex: str | None = None  # an explicit choice
    role_sex: str | None = None  # implied by being a father or a mother
    affected: str | None = None
    affected_stated: bool = False  # a proband, or named by a phenotype flag
    parents_from: str = ""  # which flag set the parents, for messages


@dataclass
class FlagResult:
    rows: list[dict]
    warnings: list[str] = field(default_factory=list)
    # What the flags actually say about each sample: a subset of PARENTS,
    # SEX, AFFECTED. The rest of its row is a default, not a statement.
    stated: dict[str, frozenset[str]] = field(default_factory=dict)

    @property
    def samples(self) -> list[str]:
        return [r["sample_id"] for r in self.rows]


def _parent_id(value: str, flag: str) -> str | None:
    """A parent given as a sample id, or None for the sentinel 0."""
    value = value.strip()
    if not value:
        raise ValueError(f"{flag} is empty; use 0 for a parent that is not known")
    return None if value == NO_PARENT else value


def parse_parents(spec: str) -> tuple[str, str | None, str | None]:
    """`CHILD=FATHER,MOTHER` -> (child, father|None, mother|None)."""
    child, sep, pair = spec.partition("=")
    father, comma, mother = pair.partition(",")
    if not (sep and comma) or not child.strip() or "," in mother:
        raise ValueError(
            f"--parents {spec!r}: expected CHILD=FATHER,MOTHER (0 for a "
            "parent that is not known)"
        )
    f = _parent_id(father, "--parents")
    m = _parent_id(mother, "--parents")
    if f is None and m is None:
        raise ValueError(
            f"--parents {spec!r}: a child needs at least one parent; "
            "the father and the mother cannot both be 0"
        )
    return child.strip(), f, m


def parse_sex(spec: str) -> tuple[str, str]:
    """`ID=male|female|unknown` -> (id, sex)."""
    sample, sep, sex = spec.partition("=")
    if not sep or not sample.strip() or sex not in SEX_CHOICES:
        raise ValueError(f"--sex {spec!r}: expected ID=male, ID=female or ID=unknown")
    return sample.strip(), sex


def check_core(flags: PedFlags) -> None:
    """The legacy trio declaration is all-or-nothing: --proband, --father and
    --mother together (a parent may be 0)."""
    names = {
        "--proband": flags.probands,
        "--father": flags.father,
        "--mother": flags.mother,
    }
    missing = [k for k, v in names.items() if not v]
    if missing:
        raise ValueError(MISSING_HINT.format(missing=", ".join(missing)))


def build_rows(flags: PedFlags, stored: list[dict] | None = None) -> FlagResult:
    """Pedigree rows (the shape of `ingest.pedigree.parse_ped`) for the
    family the flags declare. Raises ValueError on anything contradictory,
    before the database is touched.

    `stored` is the stored pedigree, given only for a --merge call: then
    --affected, --unaffected, --unknown and --sex may also name a person who
    is stored but not otherwise in the call, to change just that."""
    check_core(flags)
    probands = flags.probands
    if len(set(probands)) != len(probands):
        raise ValueError("--proband lists the same sample twice")
    father = _parent_id(flags.father or "", "--father")
    mother = _parent_id(flags.mother or "", "--mother")
    if father is None and mother is None:
        raise ValueError(
            "a proband needs at least one parent: --father and --mother "
            "cannot both be 0"
        )
    named = [p for p in (father, mother) if p is not None]
    if any(p in probands for p in named) or (father is not None and father == mother):
        raise ValueError("the probands, father and mother must be different samples")
    if not flags.family.strip():
        raise ValueError("--family is empty")
    if flags.proband_sex != "unknown" and len(probands) > 1:
        raise ValueError(
            "--proband-sex applies to a single --proband; give each "
            "sibling's sex with --sex ID=male|female or in a PED file"
        )

    members: dict[str, _Member] = {}

    def member(sample: str) -> _Member:
        return members.setdefault(sample, _Member())

    def set_parents(child: str, f: str | None, m: str | None, source: str) -> None:
        me = member(child)
        if me.parents_from and (me.father, me.mother) != (f, m):
            raise ValueError(
                f"{child} is given two different pairs of parents "
                f"({me.parents_from}, then {source})"
            )
        me.father, me.mother, me.parents_from = f, m, source
        for parent, role in ((f, "male"), (m, "female")):
            if parent is None:
                continue
            if parent == child:
                raise ValueError(f"{child} cannot be its own parent")
            pm = member(parent)
            if pm.role_sex not in (None, role):
                raise ValueError(f"{parent} is used as both a father and a mother")
            pm.role_sex = role

    for child in probands:
        set_parents(child, father, mother, "--father/--mother")
        members[child].affected = "affected"
        members[child].affected_stated = True
    for spec in flags.parents:
        child, f, m = parse_parents(spec)
        set_parents(child, f, m, f"--parents {spec}")

    # Phenotype and sex overrides apply to listed samples only.
    stored_ids = {r["sample_id"] for r in stored or ()}

    def listed(flag: str, samples: tuple[str, ...]) -> list[str]:
        for s in samples:
            if s in stored_ids:
                member(s)  # stored: the flag changes this one thing about them
            elif s not in members:
                raise ValueError(
                    f"{flag} names {s!r}, who is not in the family being "
                    "declared (give them as --proband, --father, --mother or "
                    "in --parents)"
                )
        return list(samples)

    seen: dict[str, str] = {}
    for flag, samples, value in (
        ("--affected", flags.affected, "affected"),
        ("--unaffected", flags.unaffected, "unaffected"),
        ("--unknown", flags.unknown, None),
    ):
        for s in listed(flag, samples):
            if s in seen:
                raise ValueError(f"{s} is named by both {seen[s]} and {flag}")
            seen[s] = flag
            members[s].affected = value
            members[s].affected_stated = True

    explicit_sex: dict[str, str] = {}
    for spec in flags.sex:
        sample, sex = parse_sex(spec)
        listed("--sex", (sample,))
        if explicit_sex.get(sample, sex) != sex:
            raise ValueError(f"--sex gives {sample} two different sexes")
        explicit_sex[sample] = sex
    if flags.proband_sex != "unknown":
        (only,) = probands
        if explicit_sex.get(only, flags.proband_sex) != flags.proband_sex:
            raise ValueError(
                f"--proband-sex {flags.proband_sex} and --sex {only}="
                f"{explicit_sex[only]} disagree"
            )
        explicit_sex[only] = flags.proband_sex
    for sample, sex in explicit_sex.items():
        members[sample].sex = None if sex == "unknown" else sex

    warnings = []
    rows = []
    stated: dict[str, frozenset[str]] = {}
    for sample, me in members.items():
        stated[sample] = frozenset(
            name
            for name, said in (
                (PARENTS, bool(me.parents_from)),
                (SEX, sample in explicit_sex),
                (AFFECTED, me.affected_stated),
            )
            if said
        )
        sex = me.sex if sample in explicit_sex else me.role_sex
        if me.role_sex and me.sex and me.sex != me.role_sex:
            warnings.append(
                f"{sample} is a {'father' if me.role_sex == 'male' else 'mother'} "
                f"but --sex records {me.sex}"
            )
        rows.append(
            {
                "sample_id": sample,
                "family_id": flags.family,
                "father_id": me.father or NO_PARENT,
                "mother_id": me.mother or NO_PARENT,
                "sex": sex,
                "affected": me.affected,
            }
        )
    return FlagResult(rows, warnings, stated)


def keep_stored(result: FlagResult, stored: list[dict]) -> list[dict]:
    """The rows a --merge call writes: `result.rows` with, for every person
    already stored, whatever the flags do not state about them kept as
    stored (their parents, sex and phenotype). A person who is not stored
    yet is written as built. `stored` rows are not changed."""
    old = {r["sample_id"]: r for r in stored}
    out = []
    for row in result.rows:
        was = old.get(row["sample_id"])
        said = result.stated.get(row["sample_id"], frozenset())
        if was is not None:
            row = dict(row)
            if PARENTS not in said:
                row["father_id"], row["mother_id"] = was["father_id"], was["mother_id"]
            # a stored sex wins over the one a role implies (a father the PED
            # file recorded as female stays female; the load warns about it)
            if SEX not in said and was["sex"] is not None:
                row["sex"] = was["sex"]
            if AFFECTED not in said:
                row["affected"] = was["affected"]
        out.append(row)
    return out
