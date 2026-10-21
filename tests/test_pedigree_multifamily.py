"""Several families in one pedigree, through the CLI, on both backends.

Fixture: tests/fixtures/related_ext.* (18 samples, famA with three
generations and famB; see make_related_ext.py). The store is ingested once
per backend for the module; a test that changes the pedigree starts by
reloading the base PED, which is one cheap command.

`db ped` is the subject of this part: the flag surface, `--merge` (a partial
delete of the pedigree table, which chDB runs as a mutation on a
ReplacingMergeTree), and the rule that every check happens before anything
is deleted.
"""

from __future__ import annotations

import hashlib
import json

import pytest

from tests.conftest import FIXTURES, run_cli
from tests.fixtures import make_related_ext

BACKENDS = ("chdb", "duckdb")
EXT_VCF = FIXTURES / "related_ext.vcf"
EXT_PED = FIXTURES / "related_ext.ped"
TRUTH = json.loads((FIXTURES / "related_ext.truth.json").read_text())

STORED = (
    "SELECT sample_id, family_id, father_id, mother_id, sex, affected "
    "FROM pedigree ORDER BY sample_id"
)


@pytest.fixture(scope="module", params=BACKENDS)
def backend(request):
    return request.param


@pytest.fixture(scope="module")
def ext_home(backend, tmp_path_factory):
    home = tmp_path_factory.mktemp("ext") / "vcfclick"
    home.mkdir()
    run_cli(home, backend, "db", "create", "t")
    run_cli(
        home, backend, "db", "ingest", "t", str(EXT_VCF),
        "--ingest-id", "i1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(home, backend, "db", "ped", "t", str(EXT_PED))
    return home


def db(home, backend, *args, ok=True):
    return run_cli(home, backend, "db", *args, ok=ok)


def reset(home, backend):
    """Back to the base pedigree: the whole PED replaces whatever is stored."""
    db(home, backend, "ped", "t", str(EXT_PED))


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


# ─────────────────────────── --merge replaces whole families ───────────────────────────


def test_merge_replaces_only_the_family_it_names(ext_home, backend, tmp_path):
    reset(ext_home, backend)
    assert by_family(stored(ext_home, backend)) == {"famA": 13, "famB": 5}

    # flags: famB shrinks to K1 and its parents; famA is untouched
    r = db(
        ext_home, backend, "ped", "t", "--family", "famB", "--merge",
        "--proband", "K1", "--father", "FB", "--mother", "MB",
    )  # fmt: skip
    assert "loaded pedigree: 3 individuals under ingest_id=i1" in r.stdout
    assert "families: famA(13), famB(3)" in r.stdout
    rows = stored(ext_home, backend)
    assert by_family(rows) == {"famA": 13, "famB": 3}
    assert "K2" not in rows and "XB" not in rows
    assert rows["C1"] == ["famA", "F", "M", "male", "affected"]  # famA as loaded

    # a PED file with --merge puts famB back and still leaves famA alone
    famb = "".join(
        line + "\n" for line in EXT_PED.read_text().splitlines() if "famB" in line
    )
    r = db(ext_home, backend, "ped", "t", str(ped_file(tmp_path, famb)), "--merge")
    assert "families: famA(13), famB(5)" in r.stdout
    assert by_family(stored(ext_home, backend)) == {"famA": 13, "famB": 5}

    # without --merge the same PED file is the whole pedigree, as it always was
    r = db(ext_home, backend, "ped", "t", str(ped_file(tmp_path, famb)))
    assert "families: famB(5)" in r.stdout
    assert by_family(stored(ext_home, backend)) == {"famB": 5}


def test_a_flag_call_that_would_delete_others_fails_and_names_merge(ext_home, backend):
    """The exact argument shape the desktop apps send: no --family, no
    --merge. They show nothing from a call that succeeds, so a replace that
    would drop a family has to fail, with a message they can show."""
    reset(ext_home, backend)
    r = db(
        ext_home, backend, "ped", "t",
        "--proband", "C1", "--father", "F", "--mother", "M", ok=False,
    )  # fmt: skip
    assert "--merge" in r.stderr
    assert "famB" in r.stderr and "K1" in r.stderr  # what would have been lost
    assert "`--family famA --merge`" in r.stderr  # where C1, F and M are stored
    assert "loaded pedigree" not in r.stdout
    assert by_family(stored(ext_home, backend)) == {"famA": 13, "famB": 5}

    # --merge is the way to do it on purpose: famA becomes the three, famB stays
    db(
        ext_home, backend, "ped", "t", "--family", "famA", "--merge",
        "--proband", "C1", "--father", "F", "--mother", "M",
    )  # fmt: skip
    assert by_family(stored(ext_home, backend)) == {"famA": 3, "famB": 5}


# ───────────────────────── checks happen before any delete ─────────────────────────

CYCLE_PED = "famA\tPGF\tF\t0\t1\t0\nfamA\tF\tPGF\tPGM\t1\t0\nfamA\tPGM\t0\t0\t2\t0\n"
SELF_PED = "famA\tC1\tC1\tM\t1\t2\nfamA\tM\t0\t0\t2\t0\n"
CROSS_PED = "famA\tF\t0\t0\t1\t0\nfamB\tMB\t0\t0\t2\t0\nfamB\tK1\tF\tMB\t1\t2\n"


def test_nothing_is_deleted_when_a_check_fails(ext_home, backend, tmp_path):
    reset(ext_home, backend)
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
            # C1 is stored in famA; --merge replaces families, it does not move people
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


def test_a_father_recorded_as_female_loads_with_a_warning(ext_home, backend, tmp_path):
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
    ext_home, backend, tmp_path
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


def test_zero_declares_a_duo(ext_home, backend):
    reset(ext_home, backend)
    r = db(
        ext_home, backend, "ped", "t", "--family", "famA", "--merge",
        "--proband", "H1", "--father", "0", "--mother", "M2",
    )  # fmt: skip
    assert "loaded pedigree: 2 individuals" in r.stdout
    rows = stored(ext_home, backend)
    assert rows["H1"][1:3] == ["0", "M2"] and "F" not in rows
    r = db(
        ext_home, backend, "ped", "t", "--proband", "H1", "--father", "0",
        "--mother", "0", ok=False,
    )  # fmt: skip
    assert "cannot both be 0" in r.stderr
