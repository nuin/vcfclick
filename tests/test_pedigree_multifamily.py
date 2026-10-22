"""Several families in one pedigree, through the CLI, on both backends.

Fixture: tests/fixtures/related_ext.* (18 samples, famA with three
generations and famB; see make_related_ext.py). The store is ingested once
per backend for the module. A test that changes the pedigree asks for the
`mutable` fixture, which reloads the base PED afterwards (one cheap
command), so every other test can read the base pedigree.

`db ped` is the subject of the first half: the flag surface, `--merge` (a
partial delete by sample_id list on chDB: a mutation that removes only the
rows of the samples the flags name from a ReplacingMergeTree, where every
`db ped` so far deleted the whole ingestion's rows; with a PED file it
removes the rows of the families the file lists), and the rule that every
check happens before anything is deleted. The second half reads the
families back (`db families`, `db trio`).
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

from tests.conftest import FIXTURES, run_cli
from tests.fixtures import make_related_ext

BACKENDS = ("chdb", "duckdb")
EXT_VCF = FIXTURES / "related_ext.vcf"
EXT_PED = FIXTURES / "related_ext.ped"
TRUTH = json.loads((FIXTURES / "related_ext.truth.json").read_text())

STORED = (
    "SELECT sample_id, family_id, father_id, mother_id, sex, affected "
    "FROM pedigree ORDER BY sample_id, ingested_at"  # a later duplicate wins
)


@pytest.fixture(scope="module", params=BACKENDS)
def backend(request):
    return request.param


@pytest.fixture(scope="module")
def ext_store(backend, tmp_path_factory):
    """The ingested store, once per backend: its home and whether a test
    left the pedigree changed."""
    home = tmp_path_factory.mktemp("ext") / "vcfclick"
    home.mkdir()
    run_cli(home, backend, "db", "create", "t")
    run_cli(
        home, backend, "db", "ingest", "t", str(EXT_VCF),
        "--ingest-id", "i1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(home, backend, "db", "ped", "t", str(EXT_PED))
    return {"home": home, "dirty": False}


def db(home, backend, *args, ok=True):
    return run_cli(home, backend, "db", *args, ok=ok)


@pytest.fixture
def ext_home(ext_store, backend):
    """The store's home with the base pedigree loaded."""
    if ext_store["dirty"]:  # a reset did not finish (a command was killed)
        db(ext_store["home"], backend, "ped", "t", str(EXT_PED))
        ext_store["dirty"] = False
    return ext_store["home"]


@pytest.fixture
def mutable(ext_store, ext_home, backend):
    """For a test that changes the pedigree: put the base one back after.
    Marked dirty meanwhile, so if the reset itself is interrupted (a chDB
    hang, then a rerun) the next test starts by finishing it."""
    ext_store["dirty"] = True
    yield
    db(ext_home, backend, "ped", "t", str(EXT_PED))
    ext_store["dirty"] = False


def stored(home, backend) -> dict[str, list]:
    """sample -> [family, father, mother, sex, affected] as stored."""
    r = db(home, backend, "query", "t", STORED, "--format", "JSONCompact")
    return {row[0]: row[1:] for row in json.loads(r.stdout)["data"]}


def by_family(rows: dict[str, list]) -> dict[str, int]:
    out: dict[str, int] = {}
    for fam, *_ in rows.values():
        out[fam] = out.get(fam, 0) + 1
    return dict(sorted(out.items()))


def ped_file(tmp_path, text: str, name="x.ped"):
    path = tmp_path / name
    path.write_text(text)
    return path


# ───────────────────────────── the fixture itself ─────────────────────────────


def test_the_generator_reproduces_the_committed_fixture(tmp_path):
    make_related_ext.write(tmp_path)
    for name in ("related_ext.vcf", "related_ext.ped", "related_ext.truth.json"):
        assert (tmp_path / name).read_text() == (FIXTURES / name).read_text(), name


def test_the_original_related_fixture_was_not_touched():
    """related.vcf / related.ped are pinned by the relatedness, trio-json
    and desktop suites (7 samples, 3 planted de novo sites): extending the
    pedigree fixtures must never change them."""
    digest = {
        name: hashlib.sha256((FIXTURES / name).read_bytes()).hexdigest()
        for name in ("related.vcf", "related.ped")
    }
    assert digest == {
        "related.vcf": "4f950928cada35665d8cfb3b2ed15034627c0655d4df2e0910b7de55326da1a2",
        "related.ped": "baa557125e9d7a8c255bb246a979fbefc630ff0b46f4c6eda11561174088c6ea",
    }


# ───────────────────────── --merge: what it replaces and what it keeps ─────────────────────────


def changed(before: dict, after: dict) -> dict[str, tuple]:
    """sample -> (row before, row after) for every sample that differs."""
    return {
        s: (before.get(s), after.get(s))
        for s in sorted({*before, *after})
        if before.get(s) != after.get(s)
    }


def test_a_ped_file_with_merge_replaces_the_families_it_lists_as_wholes(
    ext_home, backend, mutable, tmp_path
):
    base = stored(ext_home, backend)
    assert by_family(base) == {"famA": 13, "famB": 5}
    famb = [line for line in EXT_PED.read_text().splitlines() if "famB" in line]

    # famB listed without K2 and XB: they are removed, and the call says so
    small = ped_file(tmp_path, "\n".join(famb[:3]) + "\n", "small.ped")  # FB MB K1
    r = db(ext_home, backend, "ped", "t", str(small), "--merge")
    assert "families: famA(13), famB(3)" in r.stdout
    assert "removed from famB: K2, XB" in r.stdout
    assert by_family(stored(ext_home, backend)) == {"famA": 13, "famB": 3}

    # the whole famB again: nothing removed, famA never touched
    whole = ped_file(tmp_path, "\n".join(famb) + "\n", "whole.ped")
    r = db(ext_home, backend, "ped", "t", str(whole), "--merge")
    assert "families: famA(13), famB(5)" in r.stdout and "removed" not in r.stdout
    assert stored(ext_home, backend) == base

    # without --merge the same file is the whole pedigree, as it always was
    r = db(ext_home, backend, "ped", "t", str(whole))
    assert "families: famB(5)" in r.stdout
    assert by_family(stored(ext_home, backend)) == {"famB": 5}


def test_flag_merge_changes_only_the_people_it_names_and_what_it_states(
    ext_home, backend, mutable
):
    """The first thing to try on a family of 13: `--merge` with flags edits
    the people the flags name. The grandparents, the half siblings, the
    other family and what is stored about the parents all stay."""
    base = stored(ext_home, backend)

    # restating a stored trio rewrites exactly what is stored: F keeps his
    # parents, M hers, C1 his sex, and C2, H1, the grandparents are not touched
    r = db(
        ext_home, backend, "ped", "t", "--family", "famA", "--merge",
        "--proband", "C1", "--father", "F", "--mother", "M",
    )  # fmt: skip
    assert "loaded pedigree: 3 individuals" in r.stdout
    assert "families: famA(13), famB(5)" in r.stdout and "removed" not in r.stdout
    assert stored(ext_home, backend) == base

    # C1 gets another mother, C2 is marked affected: those two rows change
    r = db(
        ext_home, backend, "ped", "t", "--family", "famA", "--merge",
        "--proband", "C1", "--father", "F", "--mother", "M2", "--affected", "C2",
    )  # fmt: skip
    assert "families: famA(13), famB(5)" in r.stdout and "removed" not in r.stdout
    assert "replacing the stored phenotype" not in r.stderr  # C2 was asked for
    now = stored(ext_home, backend)
    assert changed(base, now) == {
        "C1": (base["C1"], ["famA", "F", "M2", "male", "affected"]),
        "C2": (base["C2"], ["famA", "F", "M", "female", "affected"]),
    }
    assert base["C2"][4] == "unaffected"

    # the old-app call with one parent as 0 turns K1 into a duo; the rest of
    # famB (MB, K2, XB) stays where it was
    r = db(
        ext_home, backend, "ped", "t", "--proband", "K1", "--father", "FB",
        "--mother", "0", "--family", "famB", "--merge",
    )  # fmt: skip
    assert "families: famA(13), famB(5)" in r.stdout and "removed" not in r.stdout
    after = stored(ext_home, backend)
    assert changed(now, after) == {
        "K1": (now["K1"], ["famB", "FB", "0", "male", "affected"])
    }


def test_flag_merge_can_change_one_thing_about_a_stored_person(
    ext_home, backend, mutable
):
    base = stored(ext_home, backend)
    db(
        ext_home, backend, "ped", "t", "--family", "famB", "--merge",
        "--proband", "K1", "--father", "FB", "--mother", "MB",
        "--unaffected", "K2", "--sex", "XB=female",
    )  # fmt: skip
    got = changed(base, stored(ext_home, backend))
    assert got == {
        "K2": (base["K2"], ["famB", "FB", "MB", "female", "unaffected"]),
        "XB": (base["XB"], ["famB", "FB", "MB", "female", None]),
    }


def test_a_flag_call_that_would_delete_others_fails_and_names_merge(
    ext_home, backend, mutable
):
    """The exact argument shape the desktop apps send: no --family, no
    --merge. They show nothing from a call that succeeds, so a replace that
    would drop a family has to fail, with a message they can show. The
    remedy it names must keep the people it lists."""
    base = stored(ext_home, backend)
    r = db(
        ext_home, backend, "ped", "t",
        "--proband", "C1", "--father", "F", "--mother", "M", ok=False,
    )  # fmt: skip
    assert "--merge" in r.stderr
    assert "famB" in r.stderr and "K1" in r.stderr  # what would have been lost
    assert "`--family famA --merge`" in r.stderr  # where C1, F and M are stored
    assert "keep everyone else as stored" in r.stderr
    assert "loaded pedigree" not in r.stdout
    assert stored(ext_home, backend) == base

    # the remedy, as printed: nobody is lost, not even the 10 other members
    # of famA that the message listed
    r = db(
        ext_home, backend, "ped", "t", "--family", "famA", "--merge",
        "--proband", "C1", "--father", "F", "--mother", "M",
    )  # fmt: skip
    assert "families: famA(13), famB(5)" in r.stdout
    assert stored(ext_home, backend) == base


def test_a_flag_call_without_merge_never_deletes_a_person_it_does_not_name(
    ext_home, backend, mutable, tmp_path
):
    """No exception for a small family. The apps send the trio flags with no
    --merge and show nothing from a call that succeeds, so a replace that
    drops an unaffected sibling (evidence the flags cannot restate), or
    swaps one proband for another, has to fail and leave the pedigree as
    it was. A call that names everyone stored still replaces."""
    fam1 = ped_file(
        tmp_path,
        "fam1\tF\t0\t0\t1\t0\nfam1\tM\t0\t0\t2\t0\n"
        "fam1\tC1\tF\tM\t1\t2\nfam1\tC2\tF\tM\t2\t1\n",
    )
    db(ext_home, backend, "ped", "t", str(fam1))
    base = stored(ext_home, backend)
    assert by_family(base) == {"fam1": 4} and base["C2"][4] == "unaffected"

    # the trio as the apps save it, over a family with an unaffected sibling
    r = db(
        ext_home, backend, "ped", "t", "--proband", "C1", "--father", "F",
        "--mother", "M", ok=False,
    )  # fmt: skip
    assert "--merge" in r.stderr and "fam1: C2" in r.stderr
    assert "loaded pedigree" not in r.stdout
    # another proband over the stored one: C1 would go
    r = db(
        ext_home, backend, "ped", "t", "--proband", "C2", "--father", "F",
        "--mother", "M", ok=False,
    )  # fmt: skip
    assert "--merge" in r.stderr and "fam1: C1" in r.stderr
    # another mother: M would go
    r = db(
        ext_home, backend, "ped", "t", "--proband", "C1", "--proband", "C2",
        "--father", "F", "--mother", "M2", ok=False,
    )  # fmt: skip
    assert "--merge" in r.stderr and "fam1: M" in r.stderr
    assert stored(ext_home, backend) == base

    # naming everybody stored deletes nobody and keeps what was asserted
    r = db(
        ext_home, backend, "ped", "t", "--proband", "C1", "--proband", "C2",
        "--father", "F", "--mother", "M", "--unaffected", "C2",
    )  # fmt: skip
    assert "families: fam1(4)" in r.stdout and "removed" not in r.stdout
    replaced = stored(ext_home, backend)
    assert by_family(replaced) == {"fam1": 4} and replaced["C2"][4] == "unaffected"

    # the remedy the message names keeps the sibling exactly as stored
    r = db(
        ext_home, backend, "ped", "t", "--merge", "--proband", "C1",
        "--father", "F", "--mother", "M",
    )  # fmt: skip
    assert "families: fam1(4)" in r.stdout and "removed" not in r.stdout
    assert stored(ext_home, backend) == replaced


def test_a_flag_call_that_shares_nobody_with_the_stored_family_fails(
    ext_home, backend, mutable, tmp_path
):
    """A stored family the call never names is not replaced just because it
    is the only one: C1, F and M are not in famB."""
    famb = [line for line in EXT_PED.read_text().splitlines() if "famB" in line]
    db(ext_home, backend, "ped", "t", str(ped_file(tmp_path, "\n".join(famb) + "\n")))
    base = stored(ext_home, backend)
    assert by_family(base) == {"famB": 5}

    r = db(
        ext_home, backend, "ped", "t", "--proband", "C1", "--father", "F",
        "--mother", "M", ok=False,
    )  # fmt: skip
    assert (
        "--merge" in r.stderr and "famB: FB, K1, K2, MB" in r.stderr
    )  # first four, then ...
    assert "loaded pedigree" not in r.stdout
    assert stored(ext_home, backend) == base


# ───────────────────────── checks happen before any delete ─────────────────────────

CYCLE_PED = "famA\tPGF\tF\t0\t1\t0\nfamA\tF\tPGF\tPGM\t1\t0\nfamA\tPGM\t0\t0\t2\t0\n"
SELF_PED = "famA\tC1\tC1\tM\t1\t2\nfamA\tM\t0\t0\t2\t0\n"
CROSS_PED = "famA\tF\t0\t0\t1\t0\nfamB\tMB\t0\t0\t2\t0\nfamB\tK1\tF\tMB\t1\t2\n"


def test_nothing_is_deleted_when_a_check_fails(ext_home, backend, tmp_path):
    before = stored(ext_home, backend)
    assert len(before) == 18

    cases = [
        (["ped", "t", str(ped_file(tmp_path, CYCLE_PED, "cycle.ped"))], "loops back"),
        (["ped", "t", str(ped_file(tmp_path, SELF_PED, "self.ped"))], "its own father"),
        (["ped", "t", str(ped_file(tmp_path, CROSS_PED, "cross.ped"))], "family 'famA', not 'famB'"),
        (
            # a cycle made with the flags: F's father is its own child C1
            ["ped", "t", "--family", "famA", "--merge", "--proband", "C1",
             "--father", "F", "--mother", "M", "--parents", "F=C1,PGM"],
            "loops back",
        ),
        (
            # C1 is stored in famA; --merge keeps people in their family, it does not move them
            ["ped", "t", "--family", "famB", "--merge", "--proband", "C1",
             "--father", "FB", "--mother", "MB"],
            "already in family 'famA'",
        ),
        (
            ["ped", "t", "--family", "famA", "--merge", "--proband", "GHOST",
             "--father", "F", "--mother", "M"],
            "GHOST",
        ),
        (
            # a PED file does not mix with the flags
            ["ped", "t", str(EXT_PED), "--proband", "C1"],
            "not both",
        ),
    ]  # fmt: skip
    for args, expected in cases:
        r = db(ext_home, backend, *args, ok=False)
        assert expected in r.stderr, (args, r.stderr)
        assert "loaded pedigree" not in r.stdout
    assert stored(ext_home, backend) == before


def test_a_problem_in_another_family_does_not_block_an_edit(ext_home, backend, mutable):
    """A cycle in famB (an older vcfclick let such rows in) is reported when
    famA is edited, not made the reason that famA cannot be edited. A cycle
    the call itself would make is still refused. (The inserted row carries a
    later ingested_at: that column has one-second resolution, so a row
    inserted in the same second as the stored one would tie with it.)"""
    db(
        ext_home, backend, "query", "t",
        "INSERT INTO pedigree (ingest_id, sample_id, family_id, father_id, "
        "mother_id, sex, affected, ingested_at) VALUES "
        "('i1', 'FB', 'famB', 'K1', '0', 'male', NULL, '2030-01-01 00:00:00')",
    )  # fmt: skip
    before = stored(ext_home, backend)
    edit = ["ped", "t", "--family", "famA", "--merge", "--proband", "C1",
            "--father", "F", "--mother", "M"]  # fmt: skip
    r = db(ext_home, backend, *edit, "--affected", "C2")
    assert "already has a problem" in r.stderr
    assert "loops back on itself among FB, K1 (family famB)" in r.stderr
    assert changed(before, stored(ext_home, backend)).keys() == {"C2"}

    r = db(ext_home, backend, *edit, "--parents", "F=C1,PGM", ok=False)
    assert "nothing was changed" in r.stderr and "(family famA)" in r.stderr
    assert "famB" not in r.stderr


def test_a_father_recorded_as_female_loads_with_a_warning(
    ext_home, backend, mutable, tmp_path
):
    lines = EXT_PED.read_text().splitlines()
    swapped = "\n".join(
        line.replace("famA\tF\tPGF\tPGM\t1", "famA\tF\tPGF\tPGM\t2") for line in lines
    )
    assert swapped != "\n".join(lines)
    r = db(ext_home, backend, "ped", "t", str(ped_file(tmp_path, swapped + "\n")))
    assert "warning: father F of C1 is recorded as female" in r.stderr
    assert "loaded pedigree: 18 individuals" in r.stdout


# ─────────────────────────── the flag surface ───────────────────────────


def test_flags_declare_a_pedigree_with_unknown_parents_and_explicit_status(
    ext_home, backend, mutable, tmp_path
):
    # start from one stored person, so a flag call is allowed to replace it
    one = ped_file(tmp_path, "famA\tPGF\t0\t0\t1\t0\n")
    db(ext_home, backend, "ped", "t", str(one))
    flags = [
        "ped", "t", "--family", "famA",
        "--proband", "C1", "--proband", "C2", "--father", "F", "--mother", "M",
        "--parents", "F=PGF,PGM", "--parents", "H1=F,M2",
    ]  # fmt: skip
    r = db(
        ext_home, backend, *flags,
        "--unaffected", "C2", "--affected", "F", "--sex", "C2=female",
    )  # fmt: skip
    assert "loaded pedigree: 8 individuals" in r.stdout
    assert "families: famA(8)" in r.stdout
    rows = stored(ext_home, backend)
    assert rows["C1"] == ["famA", "F", "M", None, "affected"]  # sex unknown
    assert rows["C2"] == ["famA", "F", "M", "female", "unaffected"]
    assert rows["H1"] == ["famA", "F", "M2", None, None]  # not a --proband
    # parents are male/female by role and have NO phenotype unless told
    assert rows["F"] == ["famA", "PGF", "PGM", "male", "affected"]
    assert rows["M"] == ["famA", "0", "0", "female", None]
    assert rows["PGF"][3:] == ["male", None] and rows["PGM"][3:] == ["female", None]

    # the same call without the phenotype flags changes what was asserted:
    # that is said on stderr (the call is still allowed, every person is named)
    r = db(ext_home, backend, *flags)
    assert "replacing the stored phenotype of" in r.stderr
    assert "C2 (unaffected -> affected)" in r.stderr
    assert "F (affected -> unknown)" in r.stderr
    assert stored(ext_home, backend)["C2"][4] == "affected"


def test_zero_declares_a_duo(ext_home, backend, mutable):
    r = db(
        ext_home, backend, "ped", "t", "--family", "famA", "--merge",
        "--proband", "H1", "--father", "0", "--mother", "M2",
    )  # fmt: skip
    assert "loaded pedigree: 2 individuals" in r.stdout
    rows = stored(ext_home, backend)
    assert rows["H1"][1:3] == ["0", "M2"] and by_family(rows) == {"famA": 13, "famB": 5}
    r = db(
        ext_home, backend, "ped", "t", "--proband", "H1", "--father", "0",
        "--mother", "0", ok=False,
    )  # fmt: skip
    assert "cannot both be 0" in r.stderr

    # A duo is listed (with an info), is not a proband, and `db trio` skips
    # it with a note; naming it is still the old error until duos are analysed.
    doc = families_json(ext_home, backend)
    fam_a = next(f for f in doc["families"] if f["family_id"] == "famA")
    duos = [t for t in fam_a["trios"] if t["kind"] == "duo"]
    assert [(t["child"], t["is_proband"]) for t in duos] == [("H1", False)]
    assert fam_a["probands"] == ["C1", "C2T", "H2"]  # base: H1 was one too
    assert ("info", "duo", "famA") in {
        (q["level"], q["code"], q["family_id"]) for q in doc["problems"]
    }
    r = db(ext_home, backend, "trio", "t", "--category", "denovo", "--format", "json")
    assert "skipping H1: one declared parent (a duo)" in r.stderr
    assert [t["trio"]["proband"] for t in json.loads(r.stdout)["trios"]] == [
        "C1", "C2T", "H2", "K1", "K2", "XB",
    ]  # fmt: skip
    r = db(ext_home, backend, "trio", "t", "--proband", "H1", ok=False)
    assert "missing a parent" in r.stderr


# ───────────────────────────────── db families ─────────────────────────────────

GOLDEN_FAMILIES = Path(__file__).parent / "goldens" / "families" / "related_ext.json"


def families_json(home, backend) -> dict:
    return json.loads(db(home, backend, "families", "t", "--format", "json").stdout)


def test_families_lists_the_two_families_as_the_fixture_says(ext_home, backend):
    doc = families_json(ext_home, backend)
    assert doc["ingest_id"] == "i1" and doc["problems"] == []
    fams = {f["family_id"]: f for f in doc["families"]}
    assert list(fams) == ["famA", "famB"]
    for fam, probands in TRUTH["probands"].items():
        assert fams[fam]["probands"] == probands
    assert {m["sample_id"] for m in fams["famB"]["members"]} == set(
        TRUTH["families"]["famB"]
    )

    trios = {t["child"]: t for f in fams.values() for t in f["trios"]}
    # F and M are trios with children, C2 and UF are unaffected: none is a proband
    assert {c for c, t in trios.items() if not t["is_proband"]} == set(
        TRUTH["trios_not_probands"]
    )
    assert all(t["kind"] == "trio" for t in trios.values())
    for child, sibs in TRUTH["siblings"].items():
        got = {s["sample"]: s["kind"] for s in trios[child]["siblings"]}
        assert got == sibs, child
    assert trios["C1"]["father"] == "F" and trios["C1"]["mother"] == "M"

    members = {m["sample_id"]: m for f in fams.values() for m in f["members"]}
    assert {s: m["generation"] for s, m in members.items()} == TRUTH["generation"]
    assert members["F"]["affected"] == "unaffected" and members["F"]["sex"] == "male"
    assert members["PGF"]["affected"] is None and members["PGF"]["father"] is None
    assert all(m["genotyped"] for m in members.values())


def test_families_json_is_pinned_and_the_same_on_both_backends(ext_home, backend):
    text = json.dumps(families_json(ext_home, backend), indent=2) + "\n"
    if os.environ.get("VCFCLICK_UPDATE_GOLDENS") == "1":
        GOLDEN_FAMILIES.parent.mkdir(parents=True, exist_ok=True)
        GOLDEN_FAMILIES.write_text(text)
    assert text == GOLDEN_FAMILIES.read_text()


def test_families_text_reads_like_a_pedigree(ext_home, backend):
    out = db(ext_home, backend, "families", "t").stdout
    assert out.startswith("ingest_id i1: 2 families, 18 individuals")
    assert "family famA: 13 individuals" in out and "family famB: 5 individuals" in out
    assert "  probands: C1, C2T, H1, H2" in out
    assert "trio C1: father F, mother M [proband]; siblings C2 (full)" in out
    assert "[not a proband: unaffected]" in out
    assert "problems:" not in out


def test_families_reads_a_pedigree_that_would_not_load(ext_home, backend, mutable):
    """A cycle and a cross-family parent that got into the table (an older
    vcfclick, a script) are reported, not fatal: `db families` exits 0 and
    `db trio --proband` still resolves its trio."""
    db(
        ext_home, backend, "query", "t",
        "INSERT INTO pedigree (ingest_id, sample_id, family_id, father_id, "
        "mother_id, sex, affected, ingested_at) VALUES "
        "('i1', 'PGF', 'famA', 'F', '0', 'male', NULL, '2030-01-01 00:00:00'), "
        "('i1', 'K2', 'famB', 'F', 'MB', 'female', NULL, '2030-01-01 00:00:00')",
    )  # fmt: skip
    doc = families_json(ext_home, backend)
    problems = {(q["code"], q["family_id"]) for q in doc["problems"]}
    assert ("cycle", "famA") in problems and ("cross_family_parent", "famB") in problems
    members = {m["sample_id"]: m for f in doc["families"] for m in f["members"]}
    assert members["PGF"]["generation"] is None and members["F"]["generation"] is None
    assert members["FB"]["generation"] == 0  # the rest still has one
    trio = json.loads(
        db(
            ext_home, backend, "trio", "t", "--proband", "C1",
            "--category", "denovo", "--format", "json",
        ).stdout
    )  # fmt: skip
    assert trio["trio"]["father"] == "F" and trio["models"]["denovo"]["count"] == 1


# ─────────────────────────────────── db trio ───────────────────────────────────

GOLDEN_TRIO = Path(__file__).parent / "goldens" / "trio_full" / "related_ext_C1.json"
PROBANDS = [p for fam in ("famA", "famB") for p in TRUTH["probands"][fam]]


def trio_json(home, backend, *args) -> dict:
    r = db(home, backend, "trio", "t", "--format", "json", *args)
    return json.loads(r.stdout)


def planted_counts(model: str, *, include_non_probands=False) -> dict[str, int]:
    """Candidates per sample under `model`, from the fixture's truth."""
    out: dict[str, int] = {}
    sites = [s["candidate_for"] for s in TRUTH["planted"].values()]
    if include_non_probands:
        sites += list(TRUTH["candidates_for_non_probands"].values())
    for candidate_for in sites:
        for sample in candidate_for.get(model, []):
            out[sample] = out.get(sample, 0) + 1
    return out


def test_trio_without_a_proband_runs_the_probands_of_every_family(ext_home, backend):
    doc = trio_json(ext_home, backend)
    assert doc["families"] == ["famA", "famB"]
    got = [t["trio"] for t in doc["trios"]]
    assert [t["proband"] for t in got] == PROBANDS
    assert [t["family_id"] for t in got] == ["famA"] * 4 + ["famB"] * 3
    assert all(t["kind"] == "trio" for t in got)
    for t in got:
        want = TRUTH["siblings"].get(t["proband"])
        if want is not None:
            assert {s["sample"]: s["kind"] for s in t["siblings"]} == want

    # the planted sites land on exactly the probands the truth names
    for model in ("denovo", "dominant", "recessive"):
        want = planted_counts(model)
        counts = {
            d["trio"]["proband"]: d["models"][model]["count"] for d in doc["trios"]
        }
        assert counts == {p: want.get(p, 0) for p in PROBANDS}, model
    # F and M are trios too, but not probands: no report for them
    assert not {"F", "M", "C2", "UF"} & {t["proband"] for t in got}


def test_all_trios_also_runs_unaffected_children_and_parents_with_children(
    ext_home, backend
):
    doc = trio_json(ext_home, backend, "--all-trios", "--category", "dominant")
    got = [t["trio"]["proband"] for t in doc["trios"]]
    assert got == ["C1", "C2", "C2T", "F", "H1", "H2", "M", "UF", "K1", "K2", "XB"]
    counts = {
        d["trio"]["proband"]: d["models"]["dominant"]["count"] for d in doc["trios"]
    }
    want = planted_counts("dominant", include_non_probands=True)
    assert counts == {p: want.get(p, 0) for p in got}
    assert counts["F"] == 1  # PGF's variant, passed to F and on to H1


def test_trio_family_and_its_errors(ext_home, backend):
    doc = trio_json(ext_home, backend, "--family", "famB", "--category", "denovo")
    assert doc["families"] == ["famB"]
    assert [t["trio"]["proband"] for t in doc["trios"]] == ["K1", "K2", "XB"]
    r = db(ext_home, backend, "trio", "t", "--family", "nope", ok=False)
    assert "famA, famB" in r.stderr
    for extra in ("--family", "famA"), ("--all-trios",):
        r = db(ext_home, backend, "trio", "t", "--proband", "C1", *extra, ok=False)
        assert "cannot be combined with --proband" in r.stderr


def test_trio_text_labels_families_and_says_what_it_left_out(ext_home, backend):
    r = db(ext_home, backend, "trio", "t", "--category", "denovo")
    lines = r.stdout.splitlines()
    assert lines.count("family: famA") == 4 and lines.count("family: famB") == 3
    assert lines.index("family: famA") < lines.index(
        "trio: proband=C1 father=F mother=M"
    )
    assert "trio: proband=K1 father=FB mother=MB" in lines
    assert "note: not analysed as probands: " in r.stderr
    assert "C2, F, UF (unaffected)" in r.stderr  # F is marked unaffected in the PED
    assert "M (has children and is not affected)" in r.stderr
    # one family: no label line
    one = db(ext_home, backend, "trio", "t", "--family", "famB", "--category", "denovo")
    assert "family:" not in one.stdout


def _normalised(doc):
    """Floats in `af` rounded to 6 places: Float32 on chDB, double on DuckDB
    (the one known backend difference, see tests/test_trio_goldens.py)."""
    if isinstance(doc, dict):
        return {
            k: round(v, 6) if k == "af" and isinstance(v, float) else _normalised(v)
            for k, v in doc.items()
        }
    if isinstance(doc, list):
        return [_normalised(v) for v in doc]
    return doc


def test_the_whole_trio_document_is_pinned(ext_home, backend):
    """The full document of one proband, new keys included: the golden of
    the OLD keys (tests/test_trio_goldens.py) strips them, this one does not."""
    doc = _normalised(trio_json(ext_home, backend, "--proband", "C1"))
    if os.environ.get("VCFCLICK_UPDATE_GOLDENS") == "1":
        GOLDEN_TRIO.parent.mkdir(parents=True, exist_ok=True)
        GOLDEN_TRIO.write_text(json.dumps(doc, indent=2) + "\n")
    assert doc == json.loads(GOLDEN_TRIO.read_text())
    assert list(doc["trio"]) == [
        "ingest_id", "proband", "father", "mother", "family_id", "kind", "siblings",
    ]  # fmt: skip


def test_probands_follow_the_phenotypes_in_the_ped(
    ext_home, backend, mutable, tmp_path
):
    lines = EXT_PED.read_text().splitlines()

    def with_phenotypes(**by_sample: str) -> str:
        out = []
        for line in lines:
            f = line.split("\t")
            if f[1] in by_sample:
                f[5] = by_sample[f[1]]
            out.append("\t".join(f))
        return "\n".join(out) + "\n"

    # every child unaffected: no proband, and the message says what to do
    kids = ("C1", "C2", "C2T", "H1", "H2", "K1", "K2", "XB", "UF", "F", "M")
    db(
        ext_home,
        backend,
        "ped",
        "t",
        str(ped_file(tmp_path, with_phenotypes(**dict.fromkeys(kids, "1")))),
    )
    r = db(ext_home, backend, "trio", "t", ok=False)
    assert "no proband in the pedigree" in r.stderr and "--all-trios" in r.stderr
    r = db(ext_home, backend, "trio", "t", "--family", "famB", ok=False)
    assert "no proband in family 'famB'" in r.stderr
    doc = trio_json(ext_home, backend, "--all-trios", "--category", "denovo")
    assert len(doc["trios"]) == len(kids)

    # F affected: a parent with children is a proband once it is affected
    db(
        ext_home,
        backend,
        "ped",
        "t",
        str(ped_file(tmp_path, with_phenotypes(F="2"), "f.ped")),
    )
    fam_a = next(
        f
        for f in families_json(ext_home, backend)["families"]
        if f["family_id"] == "famA"
    )
    assert "F" in fam_a["probands"]
    assert next(t for t in fam_a["trios"] if t["child"] == "F")["is_proband"] is True
