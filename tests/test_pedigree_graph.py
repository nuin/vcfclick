"""The pedigree graph (storage/pedigree_graph.py): pure Python, no database.

Everything here builds a `Pedigree` from rows, so no chDB session is opened
and the file costs a fraction of a second. The expected-kinship table is
checked twice: against closed-form values, and against a gene-drop Monte
Carlo that never calls the recursion under test.
"""

from __future__ import annotations

import random

import pytest

from storage.pedigree_graph import (
    GRANDPARENT_SLOTS,
    Pedigree,
    PedigreeError,
    Problem,
    Sibling,
)


def row(sample, father="0", mother="0", fam="famA", sex=None, aff=None):
    return {
        "sample_id": sample,
        "family_id": fam,
        "father_id": father,
        "mother_id": mother,
        "sex": sex,
        "affected": aff,
    }


def codes(ped: Pedigree, level: str | None = None) -> list[str]:
    return [q.code for q in ped.validate() if level is None or q.level == level]


# A three-generation family with a half sibling, a cousin and an inbred line.
#
#   PGF PGM        MGF MGM       M2
#     |  \           |            \
#   F    UF          M              \
#   |\    \         /               |
#   | \    UC      /                |
#   |  C1  C2 ----+                 |
#   H1 (F,M2)                       |
FAMILY = [
    row("PGF", sex="male"),
    row("PGM", sex="female"),
    row("MGF", sex="male"),
    row("MGM", sex="female"),
    row("M2", sex="female"),
    row("UM", sex="female"),
    row("F", "PGF", "PGM", sex="male"),
    row("UF", "PGF", "PGM", sex="male"),
    row("M", "MGF", "MGM", sex="female"),
    row("C1", "F", "M", aff="affected"),
    row("C2", "F", "M", aff="unaffected"),
    row("H1", "F", "M2", aff="affected"),
    row("UC", "UF", "UM"),  # first cousin of C1 and C2
    # sibling mating: C1 and C2 have a child
    row("I", "C1", "C2"),
]


@pytest.fixture
def family() -> Pedigree:
    return Pedigree.from_rows(FAMILY)


# ───────────────────────────── building ─────────────────────────────


def test_sentinels_and_types_are_normalised():
    ped = Pedigree.from_rows(
        [
            row("A", "0", ""),
            {"sample_id": 7, "family_id": "x", "father_id": None, "mother_id": "A"},
            row("B", "A", "0", sex="other", aff="maybe"),
        ]
    )
    assert ped.parents_of("A") == (None, None)
    assert ped.parents_of("7") == (None, "A")
    assert ped.person("B").sex is None and ped.person("B").affected is None
    assert ped.members() == ["7", "A", "B"]


def test_a_repeated_sample_keeps_the_later_row():
    ped = Pedigree.from_rows(
        [row("C", "F", "M", aff="unaffected"), row("C", "F", "M", aff="affected")]
    )
    assert len(ped) == 1 and ped.person("C").affected == "affected"


def test_families_and_members_are_sorted():
    ped = Pedigree.from_rows([row("b", fam="z"), row("a", fam="z"), row("c", fam="y")])
    assert ped.families() == ["y", "z"]
    assert ped.members("z") == ["a", "b"] and ped.members() == ["a", "b", "c"]


# ─────────────────────── trios, duos and the proband rule ───────────────────────


def test_kinds_and_the_older_all_children_rule(family):
    assert family.kind_of("C1") == "trio" and family.kind_of("F") == "trio"
    assert family.kind_of("PGF") is None  # a founder is nobody's child
    assert family.trio_of("H1") == ("F", "M2") and family.trio_of("PGF") is None
    assert family.all_children() == ["C1", "C2", "F", "H1", "I", "M", "UC", "UF"]


def test_probands_leave_out_unaffected_children_and_unaffected_intermediates(family):
    # C1 and H1 are affected; C2 is explicitly unaffected; F, M and UF are
    # trios that have children and no explicit affected status; UC and I
    # have no children and no explicit phenotype.
    assert family.probands() == ["C1", "H1", "I", "UC"]


def test_an_affected_parent_is_a_proband_and_an_unaffected_leaf_is_not():
    rows = [
        row("G1"),
        row("G2"),
        row("P", "G1", "G2", aff="affected"),
        row("Q", "G1", "G2", aff="unaffected"),
        row("K", "P", "X"),
        row("X"),
    ]
    ped = Pedigree.from_rows(rows)
    assert ped.probands() == ["K", "P"]
    assert ped.is_proband("P") and not ped.is_proband("Q")


def test_probands_are_scoped_to_a_family():
    ped = Pedigree.from_rows(
        [
            row("F1"),
            row("M1"),
            row("A", "F1", "M1", fam="famA"),
            row("F2", fam="famB"),
            row("M2", fam="famB"),
            row("B", "F2", "M2", fam="famB"),
        ]
    )
    assert ped.probands("famA") == ["A"] and ped.probands("famB") == ["B"]
    assert ped.probands() == ["A", "B"]


def test_a_duo_is_analysable_but_not_a_proband_by_default():
    ped = Pedigree.from_rows([row("D"), row("K", "D", "0")])
    assert ped.kind_of("K") == "duo"
    assert ped.probands() == [] and ped.probands(include_duos=True) == ["K"]
    assert ped.all_children() == []  # the older rule needs both parents
    assert ped.trio_of("K") is None
    assert ped.analysable_children() == ["K"]


def test_a_parent_missing_from_the_data_makes_a_duo_of_a_trio():
    rows = [row("F"), row("M"), row("C", "F", "M")]
    ped = Pedigree.from_rows(rows, genotyped={"F", "C"})
    assert ped.kind_of("C") == "duo" and ped.probands() == []
    assert ped.trio_of("C") == ("F", "M")  # declared, whatever the data holds
    assert "not_genotyped" in codes(ped, "warning")
    assert Pedigree.from_rows(rows, genotyped={"F", "M", "C"}).probands() == ["C"]
    assert Pedigree.from_rows(rows, genotyped={"F", "M"}).kind_of("C") is None


# ───────────────────────── siblings, informativeness ─────────────────────────


def test_siblings_are_full_or_half_by_the_shared_parent(family):
    assert family.siblings("C1") == [
        Sibling("C2", "full"),
        Sibling("H1", "half_paternal"),
    ]
    assert family.siblings("H1") == [
        Sibling("C1", "half_paternal"),
        Sibling("C2", "half_paternal"),
    ]
    assert family.siblings("F") == [Sibling("UF", "full")]
    assert family.siblings("PGF") == [] and family.siblings("nobody") == []


def test_a_maternal_half_sibling_and_a_strangers_child():
    ped = Pedigree.from_rows(
        [
            row("F1"),
            row("F2"),
            row("M"),
            row("A", "F1", "M"),
            row("B", "F2", "M"),
            row("C", "F2", "M3"),
        ]
    )
    assert ped.siblings("A") == [Sibling("B", "half_maternal")]
    assert ped.siblings("B") == [
        Sibling("A", "half_maternal"),
        Sibling("C", "half_paternal"),
    ]


def test_a_child_with_an_unknown_parent_shares_only_the_known_one():
    ped = Pedigree.from_rows([row("F"), row("A", "F", "M"), row("B", "F", "0")])
    assert ped.siblings("A") == [Sibling("B", "half_paternal")]
    assert ped.siblings("B") == [Sibling("A", "half_paternal")]


def test_informative_for_is_true_only_for_children_of_that_parent(family):
    # C1's variant came from M: C2 (M's child) is informative, H1 (F and M2's
    # child) is not -- its lacking the allele says nothing about M's.
    assert family.informative_for("C2", "M") is True
    assert family.informative_for("H1", "M") is False
    assert family.informative_for("H1", "F") is True
    assert family.informative_for("nobody", "F") is False


# ───────────────────────────── grandparents, generations ─────────────────────────────


def test_grandparents_fill_the_four_slots(family):
    assert tuple(family.grandparents("C1")) == GRANDPARENT_SLOTS
    assert family.grandparents("C1") == {
        "paternal_grandfather": "PGF",
        "paternal_grandmother": "PGM",
        "maternal_grandfather": "MGF",
        "maternal_grandmother": "MGM",
    }


def test_grandparents_are_none_where_unknown_or_not_genotyped(family):
    partial = family.grandparents("H1")  # M2 is a founder: no maternal side
    assert partial["paternal_grandfather"] == "PGF"
    assert partial["maternal_grandfather"] is None
    assert partial["maternal_grandmother"] is None
    assert all(v is None for v in family.grandparents("F").values())
    thin = Pedigree.from_rows(FAMILY, genotyped=set(family.members()) - {"PGM"})
    assert thin.grandparents("C1")["paternal_grandmother"] is None
    assert thin.grandparents("C1")["paternal_grandfather"] == "PGF"


def test_generation_counts_from_the_founders(family):
    gen = {s: family.generation(s) for s in family.members()}
    assert gen["PGF"] == gen["M2"] == gen["UM"] == 0
    assert gen["F"] == gen["M"] == gen["UF"] == 1
    assert gen["C1"] == gen["H1"] == gen["UC"] == 2
    assert gen["I"] == 3
    with pytest.raises(PedigreeError, match="not in the pedigree"):
        family.generation("nobody")


# ────────────────────── expected kinship and labels ──────────────────────


# Closed forms, founders unrelated. I is the child of full siblings C1, C2.
KINSHIP = [
    ("F", "C1", 0.25, "parent-child"),
    ("C1", "M", 0.25, "parent-child"),
    ("C1", "C2", 0.25, "full-siblings"),
    ("C1", "H1", 0.125, "half-siblings"),
    ("PGF", "C1", 0.125, "grandparent-grandchild"),
    ("MGM", "C2", 0.125, "grandparent-grandchild"),
    ("UF", "C1", 0.125, "aunt/uncle-niece/nephew"),
    ("C1", "UC", 0.0625, "first-cousins"),
    ("C2", "UC", 0.0625, "first-cousins"),
    ("F", "M", 0.0, "unrelated"),
    ("PGF", "MGF", 0.0, "unrelated"),
    ("H1", "UC", 0.0625, "first-cousins"),
    # the child of full siblings is inbred: phi(I, I) = (1 + 0.25) / 2, and
    # phi(I, C1) = (phi(C1, C1) + phi(C2, C1)) / 2 = (0.5 + 0.25) / 2
    ("I", "I", 0.625, None),
    ("I", "C1", 0.375, "parent-child"),
    ("I", "F", 0.25, "grandparent-grandchild"),
]


@pytest.mark.parametrize("a, b, phi, label", KINSHIP)
def test_expected_kinship_matches_the_closed_forms(family, a, b, phi, label):
    assert family.expected_kinship(a, b) == pytest.approx(phi)
    assert family.expected_kinship(b, a) == pytest.approx(phi)
    if label is not None:
        assert family.expected_label(a, b) == label
        assert family.expected_label(b, a) == label


def test_inbreeding_is_the_kinship_of_the_parents(family):
    assert family.inbreeding("I") == pytest.approx(0.25)
    assert family.inbreeding("C1") == 0.0 and family.inbreeding("PGF") == 0.0
    assert family.is_inbred_pair("I", "H1") and not family.is_inbred_pair("C1", "H1")


def test_expected_kinship_agrees_with_a_gene_drop_simulation(family):
    """The independent oracle: drop labelled haplotypes down the pedigree
    and count how often two individuals carry the same label. It does not
    use the recursion, so a mistake in the recursion cannot hide here."""
    order = [r["sample_id"] for r in FAMILY]  # parents are listed first
    parents = {
        r["sample_id"]: (r["father_id"], r["mother_id"])
        for r in FAMILY
        if r["father_id"] != "0"
    }
    pairs = [(a, b) for a, b, *_ in KINSHIP]
    trials = 40_000
    rng = random.Random(20261022)
    same = dict.fromkeys(pairs, 0.0)
    for _ in range(trials):
        hap: dict[str, tuple[str, str]] = {}
        for s in order:
            if s in parents:
                f, m = parents[s]
                hap[s] = (rng.choice(hap[f]), rng.choice(hap[m]))
            else:
                hap[s] = (f"{s}.a", f"{s}.b")
        for a, b in pairs:
            same[(a, b)] += sum(x == y for x in hap[a] for y in hap[b]) / 4
    for a, b, phi, _ in KINSHIP:
        assert same[(a, b)] / trials == pytest.approx(phi, abs=0.01), (a, b)


def test_a_pair_from_different_families_is_unrelated():
    ped = Pedigree.from_rows([row("A", fam="x"), row("B", fam="y")])
    assert ped.expected_kinship("A", "B") == 0.0
    assert ped.expected_kinship("A", "A") == 0.5
    assert ped.expected_label("A", "B") == "unrelated"


def test_expected_kinship_refuses_unknown_individuals(family):
    with pytest.raises(PedigreeError, match="not in the pedigree"):
        family.expected_kinship("C1", "nobody")
    with pytest.raises(PedigreeError):
        family.expected_label("C1", "C1")


# ─────────────────────────────── validation ───────────────────────────────


def test_a_valid_family_has_no_errors_and_one_info_per_duo(family):
    assert family.errors() == []
    family.assert_valid()
    ped = Pedigree.from_rows([row("D"), row("K", "D", "0")])
    problems = ped.validate()
    assert [(q.level, q.code) for q in problems] == [("info", "duo")]
    assert "duo release" in problems[0].message
    ped.assert_valid()  # a duo is not an error


def test_hard_errors_self_parent_same_parent_cross_family_cycle():
    self_parent = Pedigree.from_rows([row("X", "X", "M"), row("M")])
    assert codes(self_parent, "error") == ["self_parent"]

    same = Pedigree.from_rows([row("P"), row("K", "P", "P")])
    assert codes(same, "error") == ["same_parent"]

    cross = Pedigree.from_rows(
        [row("F", fam="b"), row("M", fam="a"), row("K", "F", "M", fam="a")]
    )
    problems = cross.errors()
    assert [q.code for q in problems] == ["cross_family_parent"]
    assert problems[0].family_id == "a" and "'b'" in problems[0].message

    loop = Pedigree.from_rows(
        [row("A", "B", "0"), row("B", "C", "0"), row("C", "A", "0"), row("D", "C", "0")]
    )
    errors = loop.errors()
    assert [q.code for q in errors] == ["cycle"]
    assert "A, B, C" in errors[0].message and "D" not in errors[0].message
    with pytest.raises(PedigreeError, match="loops back"):
        loop.assert_valid()


def test_sex_slips_and_missing_data_are_warnings_not_errors():
    ped = Pedigree.from_rows(
        [
            row("F", sex="female"),
            row("M", sex="male"),
            row("K", "F", "M"),
        ],
        genotyped={"F", "K"},
    )
    assert ped.errors() == []
    messages = " ".join(q.message for q in ped.validate() if q.level == "warning")
    assert "father F" in messages and "mother M" in messages
    assert "M, a parent of K, is not a sample" in messages


def test_a_parent_without_a_row_is_info_and_does_not_break_anything():
    ped = Pedigree.from_rows([row("K", "GHOST", "M"), row("M")])
    assert ("info", "parent_without_row") in [(q.level, q.code) for q in ped.validate()]
    assert ped.generation("K") == 1 and ped.expected_kinship("K", "M") == 0.25


def test_problems_are_sorted_errors_first_and_serialise():
    ped = Pedigree.from_rows([row("D"), row("K", "D", "0"), row("S", "S", "0")])
    levels = [q.level for q in ped.validate()]
    assert levels == sorted(levels, key=("error", "warning", "info").index)
    first = ped.validate()[0]
    assert isinstance(first, Problem)
    assert set(first.as_dict()) == {"family_id", "level", "code", "message"}


# ───────────────────────── tolerance of a broken pedigree ─────────────────────────


def test_a_cycle_does_not_stop_the_graph_from_loading_or_answering():
    """Reads are tolerant: only a question about the cyclic component raises."""
    ped = Pedigree.from_rows(
        [
            row("A", "B", "0"),
            row("B", "A", "0"),
            row("K", "A", "M"),  # below the cycle
            row("M"),
            row("OK", "F", "M"),  # unrelated to the cycle
            row("F"),
        ]
    )
    assert ped.parents_of("K") == ("A", "M")
    assert ped.kind_of("K") == "trio"
    assert ped.siblings("K") == [
        Sibling("B", "half_paternal"),
        Sibling("OK", "half_maternal"),
    ]
    assert ped.generation("OK") == 1 and ped.expected_kinship("OK", "F") == 0.25
    for s in ("A", "B", "K"):
        with pytest.raises(PedigreeError, match="cycle"):
            ped.generation(s)
    with pytest.raises(PedigreeError, match="cycle"):
        ped.expected_kinship("K", "M")
    assert ped.generations() == {
        "A": None,
        "B": None,
        "F": 0,
        "K": None,
        "M": 0,
        "OK": 1,
    }


def test_a_self_parent_and_a_cross_family_parent_still_answer_reads():
    ped = Pedigree.from_rows(
        [
            row("F", fam="b"),
            row("M", fam="a"),
            row("K", "F", "M", fam="a", aff="affected"),
            row("X", "X", "M", fam="a"),
        ]
    )
    assert ped.probands() == ["K", "X"]
    assert ped.trio_of("K") == ("F", "M")
    assert ped.expected_kinship("K", "F") == 0.25
    with pytest.raises(PedigreeError):
        ped.generation("X")
