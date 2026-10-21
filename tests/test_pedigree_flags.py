"""`db ped` flags -> pedigree rows (cli/pedigree_flags.py): pure, no database.

The rules the flag surface promises: probands affected and everyone else
unknown, roles fix sex, 0 declares a duo, and every contradiction is an
error raised before anything could be written.
"""

from __future__ import annotations

import pytest

from cli.pedigree_flags import PedFlags, build_rows, parse_parents, parse_sex


def rows_by_id(**kw) -> dict[str, dict]:
    flags = PedFlags(**kw)
    return {r["sample_id"]: r for r in build_rows(flags).rows}


def trio(**kw) -> dict:
    base = {"probands": ("C1",), "father": "F", "mother": "M"}
    return {**base, **kw}


def error(**kw) -> str:
    with pytest.raises(ValueError) as e:
        build_rows(PedFlags(**kw))
    return str(e.value)


def test_a_trio_has_an_affected_proband_and_unknown_parents():
    rows = rows_by_id(**trio())
    assert set(rows) == {"C1", "F", "M"}
    assert rows["C1"] == {
        "sample_id": "C1",
        "family_id": "fam1",
        "father_id": "F",
        "mother_id": "M",
        "sex": None,
        "affected": "affected",
    }
    # parents: roles fix the sex, but nobody asserted a phenotype
    assert (rows["F"]["sex"], rows["F"]["affected"]) == ("male", None)
    assert (rows["M"]["sex"], rows["M"]["affected"]) == ("female", None)
    assert rows["F"]["father_id"] == rows["F"]["mother_id"] == "0"


def test_siblings_share_the_parents_and_the_family_flag_names_the_family():
    rows = rows_by_id(**trio(probands=("C1", "C2"), family="famB"))
    assert {r["family_id"] for r in rows.values()} == {"famB"}
    assert all(rows[c]["affected"] == "affected" for c in ("C1", "C2"))
    assert rows["C2"]["father_id"] == "F" and rows["C2"]["mother_id"] == "M"


def test_phenotype_flags_set_any_listed_sample():
    rows = rows_by_id(
        **trio(
            probands=("C1", "C2"), unaffected=("C2",), affected=("F",), unknown=("C1",)
        )
    )
    assert rows["C1"]["affected"] is None  # a proband can be unknown
    assert rows["C2"]["affected"] == "unaffected"
    assert rows["F"]["affected"] == "affected" and rows["M"]["affected"] is None


def test_zero_declares_a_duo_in_either_slot():
    rows = rows_by_id(**trio(father="0"))
    assert set(rows) == {"C1", "M"}  # no father row
    assert rows["C1"]["father_id"] == "0" and rows["C1"]["mother_id"] == "M"
    rows = rows_by_id(**trio(mother="0"))
    assert set(rows) == {"C1", "F"} and rows["C1"]["mother_id"] == "0"


def test_both_parents_zero_is_an_error():
    assert "cannot both be 0" in error(**trio(father="0", mother="0"))
    assert "cannot both be 0" in error(**trio(parents=("F=0,0",)))


def test_the_trio_flags_are_all_or_nothing_and_the_hint_names_zero():
    msg = error(probands=("C1",), father="F")
    assert "--mother" in msg and "missing" in msg and "0" in msg
    msg = error(probands=("C1",))
    assert "--father" in msg and "--mother" in msg
    assert "--proband" in error(father="F", mother="M")
    assert "give a PED file" in error()


def test_grandparents_and_half_siblings_come_from_parents():
    rows = rows_by_id(
        **trio(
            probands=("C1", "H1"),
            parents=("F=PGF,PGM", "M=MGF,MGM"),
            family="famA",
        )
    )
    assert set(rows) == {"C1", "H1", "F", "M", "PGF", "PGM", "MGF", "MGM"}
    assert rows["F"]["father_id"] == "PGF" and rows["F"]["mother_id"] == "PGM"
    assert (rows["PGF"]["sex"], rows["PGM"]["sex"]) == ("male", "female")
    assert rows["F"]["affected"] is None and rows["PGF"]["affected"] is None

    half = rows_by_id(**trio(parents=("H1=F,M2",)))
    assert half["H1"]["father_id"] == "F" and half["H1"]["mother_id"] == "M2"
    assert half["H1"]["affected"] is None  # not a --proband: unknown
    assert half["M2"]["sex"] == "female"


def test_the_documented_example_builds_what_the_docs_say():
    """docs/TRIO.md: grandparents, an unaffected sibling and an affected half
    sibling in a named family."""
    rows = rows_by_id(
        probands=("C1",),
        father="F",
        mother="M",
        family="famA",
        parents=("F=PGF,PGM", "C2=F,M", "H1=F,M2"),
        unaffected=("C2",),
        affected=("H1",),
        sex=("C1=male",),
    )
    assert set(rows) == {"C1", "C2", "H1", "F", "M", "M2", "PGF", "PGM"}
    assert {r["family_id"] for r in rows.values()} == {"famA"}
    assert (rows["C1"]["sex"], rows["C1"]["affected"]) == ("male", "affected")
    assert rows["C2"]["affected"] == "unaffected" and rows["C2"]["father_id"] == "F"
    assert rows["H1"]["affected"] == "affected" and rows["H1"]["mother_id"] == "M2"
    assert rows["F"]["father_id"] == "PGF" and rows["F"]["affected"] is None


def test_a_parent_may_be_zero_in_parents_too():
    rows = rows_by_id(**trio(parents=("D=0,M",)))
    assert rows["D"]["father_id"] == "0" and rows["D"]["mother_id"] == "M"


def test_a_child_given_the_same_parents_twice_is_fine_different_ones_are_not():
    rows_by_id(**trio(parents=("C1=F,M",)))
    assert "two different pairs" in error(**trio(parents=("C1=F,M2",)))
    assert "two different pairs" in error(
        **trio(parents=("K=F,M", "K=F,M2")),
    )


def test_one_sample_cannot_be_both_a_father_and_a_mother():
    assert "both a father and a mother" in error(**trio(parents=("K=M,X",)))


def test_self_parent_and_same_person_in_both_slots():
    assert "own parent" in error(**trio(parents=("F=F,X",)))
    assert "different samples" in error(probands=("C1",), father="X", mother="X")
    assert "different samples" in error(probands=("C1",), father="C1", mother="M")
    # the same sample twice
    assert "twice" in error(**trio(probands=("C1", "C1")))


def test_sex_follows_the_role_unless_told_otherwise():
    result = build_rows(PedFlags(**trio(sex=("C1=male", "F=female", "M=unknown"))))
    rows = {r["sample_id"]: r for r in result.rows}
    assert rows["C1"]["sex"] == "male"
    assert rows["F"]["sex"] == "female"  # the explicit choice wins ...
    assert rows["M"]["sex"] is None  # ... and unknown clears the role's
    assert result.warnings == ["F is a father but --sex records female"]


def test_proband_sex_for_one_proband_and_its_errors():
    assert rows_by_id(**trio(proband_sex="female"))["C1"]["sex"] == "female"
    msg = error(**trio(probands=("C1", "C2"), proband_sex="male"))
    assert "PED" in msg and "--sex" in msg
    assert "disagree" in error(**trio(proband_sex="male", sex=("C1=female",)))
    ok = rows_by_id(**trio(proband_sex="male", sex=("C1=male",)))
    assert ok["C1"]["sex"] == "male"


def test_listed_flags_must_name_members_once():
    assert "'Z'" in error(**trio(affected=("Z",)))
    assert "who is not in the family" in error(**trio(unaffected=("Z",)))
    assert "'Z'" in error(**trio(sex=("Z=male",)))
    assert "both --affected and --unaffected" in error(
        **trio(affected=("F",), unaffected=("F",))
    )
    assert "two different sexes" in error(**trio(sex=("F=male", "F=female")))


def test_spec_parsers():
    assert parse_parents("F=PGF,PGM") == ("F", "PGF", "PGM")
    assert parse_parents("K=0,M") == ("K", None, "M")
    assert parse_sex("C1=male") == ("C1", "male")
    for bad in ("F", "F=PGF", "F=PGF,PGM,X", "=A,B", "F=,M"):
        with pytest.raises(ValueError, match="--parents"):
            parse_parents(bad)
    for bad in ("C1", "C1=other", "=male"):
        with pytest.raises(ValueError, match="--sex"):
            parse_sex(bad)


def test_an_empty_family_id_is_refused():
    assert "--family" in error(**trio(family=" "))
