"""What a `db ped` write keeps and what it refuses: pure, no database.

`resulting_rows` says what a load leaves behind (by family for a PED file,
by sample for flags), `removed_rows` who is gone, and `check_result` which
problems refuse a load and which are only reported. The CLI behaviour on
both backends is in tests/test_pedigree_multifamily.py.
"""

from __future__ import annotations

import pytest

from ingest.pedigree import check_result, removed_rows, resulting_rows


def row(sample, family, father="0", mother="0", sex=None, affected=None):
    return {
        "sample_id": sample,
        "family_id": family,
        "father_id": father,
        "mother_id": mother,
        "sex": sex,
        "affected": affected,
    }


A = [  # famA: a trio, a sibling and a grandfather of C1
    row("GF", "famA", sex="male"),
    row("F", "famA", "GF", "0", "male"),
    row("M", "famA", sex="female"),
    row("C1", "famA", "F", "M", "male", "affected"),
    row("C2", "famA", "F", "M", "female", "unaffected"),
]
B = [
    row("FB", "famB", sex="male"),
    row("MB", "famB", sex="female"),
    row("K1", "famB", "FB", "MB", "male", "affected"),
]
STORED = A + B


def ids(rows):
    return sorted(r["sample_id"] for r in rows)


# ───────────────────────────── resulting_rows ─────────────────────────────


def test_a_replace_is_the_new_rows_alone():
    new = [row("C1", "famA", "F", "M")]
    assert resulting_rows(new, STORED, merge=False) == new


def test_by_sample_replaces_the_named_rows_and_nobody_else():
    """The flags: C1's parents change; the rest of famA, and famB, stay."""
    new = [row("C1", "famA", "F", "X"), row("F", "famA", "0", "0", "male")]
    out = resulting_rows(new, STORED, merge=True, by_sample=True)
    assert ids(out) == ids(STORED)  # nobody added or lost
    by_id = {r["sample_id"]: r for r in out}
    assert by_id["C1"]["mother_id"] == "X"
    assert by_id["F"]["father_id"] == "0"  # replaced as given
    assert by_id["C2"] == STORED[4] and by_id["GF"] == STORED[0]


def test_by_sample_cannot_move_a_person_to_another_family():
    new = [row("C1", "famB", "FB", "MB")]
    with pytest.raises(ValueError, match=r"already in family 'famA'.*--family famA"):
        resulting_rows(new, STORED, merge=True, by_sample=True)


def test_by_family_replaces_the_families_as_wholes():
    """A PED file: famA is exactly the file's famA, so C2 and GF are gone."""
    new = [row("F", "famA", sex="male"), row("C1", "famA", "F", "0")]
    out = resulting_rows(new, STORED, merge=True)
    assert ids(out) == ["C1", "F", "FB", "K1", "MB"]
    assert ids(removed_rows(STORED, out)) == ["C2", "GF", "M"]


def test_by_family_still_refuses_to_move_a_person_between_families():
    new = [row("C1", "famB", "FB", "MB"), row("K1", "famB")]
    with pytest.raises(ValueError, match="already in family 'famA'"):
        resulting_rows(new, STORED, merge=True)
    # unless the family they leave is replaced too
    both = [*new, row("F", "famA")]
    assert "C1" in ids(resulting_rows(both, STORED, merge=True))


def test_a_merge_needs_the_stored_rows():
    with pytest.raises(ValueError, match="stored pedigree rows"):
        resulting_rows(A, None, merge=True)


def test_removed_rows_is_the_stored_rows_that_are_gone():
    assert removed_rows(None, STORED) == [] and removed_rows([], STORED) == []
    assert ids(removed_rows(STORED, A)) == ["FB", "K1", "MB"]
    assert removed_rows(STORED, STORED) == []


# ───────────────────────────── check_result ─────────────────────────────

# a cycle only an older vcfclick could have stored: FB's father is K1
CYCLIC_B = [
    row("FB", "famB", "K1", "0", "male"),
    row("MB", "famB", sex="female"),
    row("K1", "famB", "FB", "MB", "male", "affected"),
]


def test_a_new_error_refuses_the_load_and_names_the_family():
    looped = [row("GF", "famA", "C1", "0", "male"), *A[1:]]  # C1 -> F -> GF -> C1
    with pytest.raises(ValueError, match=r"loops back.*\(family famA\)") as e:
        check_result(looped, STORED, {"famA"})
    assert "nothing was changed" in str(e.value)
    # a replace has nothing earlier to excuse it
    with pytest.raises(ValueError, match="loops back"):
        check_result(CYCLIC_B, None, {"famB"})


def test_a_problem_the_stored_pedigree_already_had_is_reported_not_fatal():
    """famB has a cycle that got in earlier. Editing famA is not blocked, and
    the warning says where the problem is."""
    stored = A + CYCLIC_B
    result = [
        *(r for r in A if r["sample_id"] != "C1"),
        row("C1", "famA", "F", "M"),
        *CYCLIC_B,
    ]
    warnings = check_result(result, stored, {"famA"})
    assert len(warnings) == 1
    assert "already has a problem" in warnings[0]
    assert "loops back on itself among FB, K1" in warnings[0]
    assert "(family famB)" in warnings[0]

    # also in the family being written, if this call did not cause it
    assert check_result(stored, stored, {"famB"}) == warnings

    # but a second, new cycle in famA is still refused, and the old one is
    # not what the message is about
    bad = [
        row("GF", "famA", "C2", "0", "male"),
        *[r for r in result if r["sample_id"] != "GF"],
    ]
    with pytest.raises(ValueError) as e:
        check_result(bad, stored, {"famA"})
    assert "(family famA)" in str(e.value) and "famB" not in str(e.value)


def test_a_parent_recorded_with_the_wrong_sex_is_a_warning_for_a_written_family():
    flipped = [
        r if r["sample_id"] != "F" else row("F", "famA", "GF", "0", "female") for r in A
    ]
    assert check_result(flipped, None, {"famA"}) == [
        "father F of C1 is recorded as female",
        "father F of C2 is recorded as female",
    ]
    assert check_result(flipped, None, {"famB"}) == []  # not a family written
    # already stored that way: said when it was stored, not at every merge
    assert check_result(flipped, flipped, {"famA"}) == []
    assert check_result(A + B, A + B, {"famA"}) == []
