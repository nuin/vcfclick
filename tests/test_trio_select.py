"""Which trios `db trio` analyses without --proband (cli/trio_select.py).

Pure: a Pedigree built from rows, no database.
"""

from __future__ import annotations

import click
import pytest

from cli.trio_select import make_pick, select_trios
from storage.pedigree_graph import Pedigree


def row(sample, father="0", mother="0", fam="famA", aff=None):
    return {
        "sample_id": sample,
        "family_id": fam,
        "father_id": father,
        "mother_id": mother,
        "sex": None,
        "affected": aff,
    }


# famA: three generations (G1, G2 -> F; F, M -> C1 affected, C2 unaffected);
# famB: a plain trio with two children. A duo in famB.
ROWS = [
    row("G1"),
    row("G2"),
    row("F", "G1", "G2"),
    row("M"),
    row("C1", "F", "M", aff="affected"),
    row("C2", "F", "M", aff="unaffected"),
    row("FB", fam="famB"),
    row("MB", fam="famB"),
    row("K2", "FB", "MB", fam="famB"),
    row("K1", "FB", "MB", fam="famB", aff="affected"),
    row("D", "FB", "0", fam="famB"),
]


@pytest.fixture
def graph() -> Pedigree:
    return Pedigree.from_rows(ROWS)


def test_without_options_the_probands_of_every_family_in_family_order(graph):
    sel = select_trios(graph)
    assert [(p.family_id, p.child) for p in sel.picks] == [
        ("famA", "C1"),
        ("famB", "K1"),
        ("famB", "K2"),
    ]
    c1 = sel.picks[0]
    assert (c1.father, c1.mother) == ("F", "M")
    assert c1.siblings == (("C2", "full"),)


def test_notes_name_what_was_left_out_and_the_duo(graph):
    notes = select_trios(graph).notes
    assert any(
        n.startswith("not analysed as probands: ")
        and "C2 (unaffected)" in n
        and "F (has children and is not affected)" in n
        and "--all-trios" in n
        for n in notes
    ), notes
    assert any(
        n.startswith("skipping D: one declared parent (a duo)") for n in notes
    ), notes


def test_all_trios_is_every_child_with_both_parents_declared(graph):
    sel = select_trios(graph, all_trios=True)
    assert [p.child for p in sel.picks] == ["C1", "C2", "F", "K1", "K2"]
    assert not any(n.startswith("not analysed") for n in sel.notes)
    assert any("duo" in n for n in sel.notes)  # a duo is skipped either way


def test_a_family_selects_its_own_probands(graph):
    assert [p.child for p in select_trios(graph, family="famA").picks] == ["C1"]
    sel = select_trios(graph, family="famB", all_trios=True)
    assert [p.child for p in sel.picks] == ["K1", "K2"]


def test_an_unknown_family_lists_the_known_ones(graph):
    with pytest.raises(click.ClickException) as e:
        select_trios(graph, family="nope")
    assert "nope" in e.value.message and "famA, famB" in e.value.message


def test_nothing_to_analyse_says_why_and_points_at_all_trios():
    unaffected = Pedigree.from_rows(
        [
            row("F"),
            row("M"),
            row("A", "F", "M", aff="unaffected"),
            row("B", "F", "M", aff="unaffected"),
        ]
    )
    with pytest.raises(click.ClickException) as e:
        select_trios(unaffected)
    msg = e.value.message
    assert "no proband in the pedigree" in msg and "A, B (unaffected)" in msg
    assert "--all-trios" in msg and "--proband" in msg
    assert [p.child for p in select_trios(unaffected, all_trios=True).picks] == [
        "A",
        "B",
    ]
    with pytest.raises(click.ClickException) as e:
        select_trios(unaffected, family="famA")
    assert "family 'famA'" in e.value.message


def test_no_trio_at_all_keeps_the_original_message():
    with pytest.raises(click.ClickException) as e:
        select_trios(Pedigree.from_rows([row("F"), row("K", "F", "0")]))
    assert e.value.message == (
        "no child with both parents in the pedigree. Load one with "
        "`vcfclick db ped`, or name the proband with --proband."
    )
    with pytest.raises(click.ClickException):
        select_trios(Pedigree.from_rows([]))


def test_a_child_whose_parent_is_not_in_the_data_is_left_out_not_dropped_silently():
    ped = Pedigree.from_rows(
        [row("F"), row("M"), row("C", "F", "M"), row("D", "F", "M")],
        genotyped={"F", "C", "D"},
    )
    with pytest.raises(click.ClickException) as e:
        select_trios(ped)
    assert "a parent is not a sample in the data" in e.value.message
    assert [p.child for p in select_trios(ped, all_trios=True).picks] == ["C", "D"]


def test_make_pick_labels_with_the_given_parents(graph):
    pick = make_pick(graph, "C1", "F", "M")
    assert (pick.family_id, pick.child) == ("famA", "C1")
    assert make_pick(graph, "nobody", "X", "Y").family_id is None
