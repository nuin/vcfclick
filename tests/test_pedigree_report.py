"""The shapes `db families` and /api/meta serve (cli/pedigree_report.py).

Pure: a Pedigree built from rows in, plain data out. The keys are a contract
with the desktop apps and scripts, so they are listed here in full.
"""

from __future__ import annotations

from cli.pedigree_report import families_document, families_text, web_families
from storage.pedigree_graph import Pedigree


def row(sample, father="0", mother="0", fam="famA", sex=None, aff=None):
    return {
        "sample_id": sample,
        "family_id": fam,
        "father_id": father,
        "mother_id": mother,
        "sex": sex,
        "affected": aff,
    }


ROWS = [
    row("F", sex="male"),
    row("M", sex="female"),
    row("C1", "F", "M", sex="male", aff="affected"),
    row("C2", "F", "M", aff="unaffected"),
    row("D", "F", "0"),  # a duo
    row("GHOST", "F", "M"),  # in the pedigree, not in the data
    row("X1", fam="famB"),
    row("X2", "X1", "0", fam="famB", aff="affected"),
]
GENOTYPED = {"F", "M", "C1", "C2", "D", "X1", "X2"}


def doc():
    return families_document(Pedigree.from_rows(ROWS, genotyped=GENOTYPED), "i1")


def test_the_document_has_exactly_the_documented_keys():
    d = doc()
    assert list(d) == ["ingest_id", "families", "problems"]
    fam = d["families"][0]
    assert list(fam) == ["family_id", "members", "probands", "trios"]
    assert list(fam["members"][0]) == [
        "sample_id", "father", "mother", "sex", "affected", "genotyped", "generation",
    ]  # fmt: skip
    assert list(fam["trios"][0]) == [
        "child", "father", "mother", "kind", "is_proband", "siblings",
    ]  # fmt: skip
    assert list(d["problems"][0]) == ["family_id", "level", "code", "message"]
    assert fam["trios"][0]["siblings"][0].keys() == {"sample", "kind"}


def test_members_probands_trios_and_duos():
    fam = {f["family_id"]: f for f in doc()["families"]}["famA"]
    members = {m["sample_id"]: m for m in fam["members"]}
    assert [m["sample_id"] for m in fam["members"]][:2] == ["F", "M"]  # founders first
    assert members["C1"] == {
        "sample_id": "C1", "father": "F", "mother": "M", "sex": "male",
        "affected": "affected", "genotyped": True, "generation": 1,
    }  # fmt: skip
    assert members["GHOST"]["genotyped"] is False
    # C2 is unaffected, GHOST is not in the data, D is a duo: only C1 is a proband
    assert fam["probands"] == ["C1"]
    trios = {t["child"]: t for t in fam["trios"]}
    assert trios["D"]["kind"] == "duo" and trios["D"]["mother"] is None
    assert trios["D"]["is_proband"] is False
    assert trios["C1"]["is_proband"] is True and trios["C2"]["is_proband"] is False
    assert {s["sample"]: s["kind"] for s in trios["C1"]["siblings"]} == {
        "C2": "full",
        "D": "half_paternal",
        "GHOST": "full",
    }
    assert "GHOST" not in trios  # its parents are fine, but it is not genotyped


def test_problems_say_what_is_wrong_and_what_is_merely_worth_knowing():
    problems = {(q["level"], q["code"]) for q in doc()["problems"]}
    assert ("info", "duo") in problems
    assert ("warning", "not_genotyped") in problems
    assert not [q for q in doc()["problems"] if q["level"] == "error"]


def test_a_cycle_has_no_generation_and_does_not_stop_the_report():
    loop = Pedigree.from_rows([row("A", "B"), row("B", "A"), row("C")])
    d = families_document(loop, "i1")
    gens = {m["sample_id"]: m["generation"] for m in d["families"][0]["members"]}
    assert gens == {"C": 0, "A": None, "B": None}
    assert [m["sample_id"] for m in d["families"][0]["members"]] == ["C", "A", "B"]
    assert d["problems"][0]["level"] == "error"
    assert "in or below a cycle" in families_text(loop, "i1")


def test_an_empty_pedigree_is_an_empty_document():
    assert families_document(Pedigree.from_rows([]), "i1") == {
        "ingest_id": "i1",
        "families": [],
        "problems": [],
    }


def test_the_web_shape_is_probands_and_trios_per_family():
    fams = web_families(Pedigree.from_rows(ROWS, genotyped=GENOTYPED))
    assert [f["family_id"] for f in fams] == ["famA", "famB"]
    assert list(fams[0]) == ["family_id", "probands", "trios"]
    assert fams[0]["probands"] == ["C1"] and fams[1]["probands"] == []
    c1 = next(t for t in fams[0]["trios"] if t["child"] == "C1")
    assert c1 == {
        "child": "C1",
        "father": "F",
        "mother": "M",
        "kind": "trio",
        "affected": "affected",
        "siblings": [
            {"sample": "C2", "kind": "full"},
            {"sample": "D", "kind": "half_paternal"},
            {"sample": "GHOST", "kind": "full"},
        ],
    }
    x2 = next(t for t in fams[1]["trios"] if t["child"] == "X2")
    assert x2["kind"] == "duo" and x2["mother"] is None


def test_the_text_names_each_family_its_probands_and_why_not():
    text = families_text(Pedigree.from_rows(ROWS, genotyped=GENOTYPED), "i1")
    assert text.startswith("ingest_id i1: 2 families, 8 individuals")
    assert "  probands: C1" in text
    assert "  probands: (none)" in text  # famB's only child is a duo
    assert "trio C2: father F, mother M [not a proband: unaffected]" in text
    assert "duo D: father F [not a proband: one parent: a duo" in text
    assert "problems:" in text
