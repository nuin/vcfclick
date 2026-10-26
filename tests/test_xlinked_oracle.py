"""The X pattern table against first principles, and the pure-Python units.

An allele-transmission oracle: a son's single X is one of his mother's two;
a daughter has her father's single X and one of her mother's. Over every
stored genotype combination, each de novo pattern must be Mendel-inconsistent,
each recessive / female_het pattern consistent, and every other combination
one of the documented classes. In-process, no database.
"""

from __future__ import annotations

import itertools

import pytest

from cli.sex import (
    SexEvidence,
    infer_sex,
    resolve_trio_sexes,
    transmission_call,
)
from cli.xchrom import GRCH38
from cli.xlinked import (
    ANY,
    PATTERNS,
    _in,
    classify,
    patterns_for,
    x_guard_sql,
)

GT = (0, 1, 2, None)


def _mother_alleles(m: int) -> list[tuple[int, int]]:
    return {0: [(0, 0)], 1: [(0, 1)], 2: [(1, 1)]}[m]


def consistent(sex: str, p: int, f: int, m: int) -> bool:
    """Can a child with stored gt `p` come from these parents? Male calls are
    haploid: stored 0 is REF, 2 is ALT (1 is not a male genotype)."""
    mats = {a for pair in _mother_alleles(m) for a in pair}
    if sex == "male":
        return (1 if p == 2 else 0) in mats
    fa = 1 if f == 2 else 0
    return any(fa + a == p for a in mats)


def test_every_pattern_is_on_the_right_side_of_mendel():
    for sex in ("male", "female"):
        for p, f, m in itertools.product((1, 2), (0, 2), (0, 1, 2)):
            if sex == "male" and p == 1:
                continue
            cls = classify(sex, p, f, m)
            ok = consistent(sex, p, f, m)
            mode = next((x.mode for x in PATTERNS if x.name == cls), None)
            if mode == "denovo":
                assert not ok, (sex, p, f, m)
            elif mode in ("recessive", "female_het"):
                assert ok, (sex, p, f, m)
            elif cls == "consistent_unpatterned":
                assert ok, (sex, p, f, m)
            else:
                assert cls == "mendel_error_other" and not ok, (sex, p, f, m, cls)


def test_artefacts_and_missing_calls():
    assert classify("male", 1, 0, 1) == "male_het"
    assert classify("female", 1, 1, 0) == "father_het"
    assert classify("female", 2, None, 1) is None
    assert classify("male", 2, None, 1) == "xlr_son"  # the father is not read
    assert classify("male", 2, 2, None) is None
    assert classify(None, 2, 0, 1) is None
    assert classify("female", 0, 2, 1) is None
    # one new ALT allele in a daughter: maternal, or paternal (the mirror)
    assert classify("female", 2, 2, 0) == "x_denovo_daughter_mat"
    assert classify("female", 2, 0, 1) == "x_denovo_daughter_pat"
    assert classify("female", 2, 0, 2) == "x_denovo_daughter_pat"
    # two new alleles, or an ALT lost: an error, not a de novo
    assert classify("female", 2, 0, 0) == "mendel_error_other"
    assert classify("female", 1, 2, 2) == "mendel_error_other"


@pytest.mark.parametrize("sex", ["male", "female"])
def test_classify_is_total_and_patterns_do_not_overlap(sex):
    for p, f, m in itertools.product(GT, GT, GT):
        classify(sex, p, f, m)  # never raises
        matches = [
            x.name
            for x in patterns_for(sex)
            if all(
                cell is ANY or value in cell
                for cell, value in ((x.p, p), (x.f, f), (x.m, m))
            )
        ]
        assert len(matches) <= 1, (p, f, m, matches)


def test_in_renders_null_and_any():
    assert _in("t.p_gt", (2,)) == "(t.p_gt IN (2))"
    assert _in("t.m_gt", (1, 2)) == "(t.m_gt IN (1, 2))"
    assert _in("t.f_gt", (0, None)) == "(t.f_gt IN (0) OR t.f_gt IS NULL)"
    assert _in("t.f_gt", ANY) == "1 = 1"


def test_guard_sql_text_and_the_mcp_briefing_say_the_same():
    from vcfclick_mcp.server import SCHEMA_DESCRIPTION

    male = x_guard_sql(GRCH38, father=True, proband_male=True)
    assert male == (
        "NOT ((g.chrom IN ('chrX', 'X', 'chr23', '23') AND "
        "(g.pos BETWEEN 2781480 AND 155701382)) AND (f.gt = 1 OR g.gt = 1))"
    )
    assert x_guard_sql(GRCH38, father=True, proband_male=False).endswith(
        "AND (f.gt = 1))"
    )
    assert x_guard_sql(GRCH38, father=False, proband_male=False) is None
    for phrase in (
        "g.chrom IN ('chrX', 'X', 'chr23', '23')",
        "(f.gt = 1 OR g.gt = 1)",
        "hemizygous",
        "PAR",
        "xlinked",
        "2,781,480-155,701,382",
    ):
        assert phrase in SCHEMA_DESCRIPTION, phrase


# ───────────────────────────── sex ─────────────────────────────


def test_infer_sex_thresholds():
    assert infer_sex(1, 19) == ("unknown", None)
    assert infer_sex(2, 20)[0] == "male"
    assert infer_sex(3, 20)[0] == "ambiguous"
    assert infer_sex(4, 20)[0] == "female"


def test_transmission_veto():
    assert transmission_call(4, 4) == (None, None)  # too few sites
    assert transmission_call(19, 1)[0] == "male"
    assert transmission_call(21, 17)[0] == "female"
    assert transmission_call(10, 5)[0] is None


def _ev(declared=None, inferred="unknown", frac=None, transmitted=None):
    return SexEvidence(declared, inferred, frac, transmitted)


def test_resolution_precedence():
    ev = {"P": _ev(), "F": _ev(), "M": _ev()}
    r = resolve_trio_sexes("P", "F", "M", ev)
    assert (r.father_sex, r.father_sex_source) == ("male", "role")
    assert r.blocked_reason == "needs_sex" and r.proband_sex is None

    ev["P"] = _ev(inferred="female", frac=0.6)
    r = resolve_trio_sexes("P", "F", "M", ev)
    assert (r.proband_sex, r.proband_sex_source, r.blocked_reason) == (
        "female", "inferred", None,
    )  # fmt: skip

    ev["P"] = _ev("male", "female", 0.6)  # declared wins, with a warning
    r = resolve_trio_sexes("P", "F", "M", ev)
    assert (r.proband_sex, r.proband_sex_source) == ("male", "pedigree")
    assert [w["code"] for w in r.warnings] == ["sex_conflict"]

    ev["P"] = _ev(inferred="male", frac=0.05, transmitted="female")
    r = resolve_trio_sexes("P", "F", "M", ev)
    assert r.blocked_reason == "needs_sex"
    assert [w["code"] for w in r.warnings] == ["sex_ambiguous"]


def test_a_parent_whose_sex_contradicts_the_role_blocks():
    ev = {"P": _ev("male"), "F": _ev(inferred="female", frac=0.5), "M": _ev()}
    r = resolve_trio_sexes("P", "F", "M", ev)
    assert r.blocked_reason == "parent_sex_conflict"
    assert r.conflict == "F is the father, chrX says female"
    ev["F"] = _ev("male", "female", 0.5)
    r = resolve_trio_sexes("P", "F", "M", ev)
    assert r.blocked_reason == "parent_sex_conflict"
    assert [w["code"] for w in r.warnings] == ["sex_conflict"]
    ev["F"] = _ev(declared="female", inferred="male", frac=0.02)
    assert resolve_trio_sexes("P", "F", "M", ev).conflict == (
        "F is the father, declared female"
    )
