"""The SQL of the opt-in trio gates, as text.

`strict_missing`, `pass_only` and `parent_alt_screen` change the WHERE clause of
`db trio` only when they are on; with all three off the text is the classic
text (tests/test_trio_sql_text.py pins that). In-process, no database.
"""

from __future__ import annotations

from cli.trio_analysis import Gates, Trio, _where_parts, trio_sql

GATES = Gates(20, 10, 0.01, 0.25, 0.75)


TRIO = Trio("i1", "C1", "F", "M")


def parts(gates=GATES, model="MODEL", **kw):
    return dict(_where_parts(TRIO, gates, model, **kw))


def test_the_gate_parts_are_named_and_ordered_and_absent_when_off():
    assert list(parts()) == [
        "model", "quality_proband", "quality_father", "quality_mother",
        "allele_balance", "population_af",
    ]  # fmt: skip
    on = GATES._replace(pass_only=True)
    assert list(parts(on, parent_alt_screen=True))[-2:] == [
        "pass_only",
        "parent_alt_screen",
    ]
    skipped = GATES._replace(af_gate="skipped_cohort_derived")
    assert "population_af" not in parts(skipped)


def test_strict_missing_turns_a_null_into_a_fail():
    tolerant = parts()
    strict = parts(GATES._replace(strict_missing=True))
    assert "g.gq IS NULL OR" in tolerant["quality_proband"]
    for who in ("proband", "father", "mother"):
        alias = {"proband": "g", "father": "f", "mother": "m"}[who]
        sql = strict[f"quality_{who}"]
        assert f"({alias}.gq IS NOT NULL AND {alias}.gq >= 20)" in sql
        assert f"({alias}.dp IS NOT NULL AND {alias}.dp >= 10)" in sql
        assert (
            f"{alias}.gt != 1 OR ({alias}.ad_ref IS NOT NULL"
            in strict["allele_balance"]
        )
    assert strict["allele_balance"].startswith(tolerant["allele_balance"])
    assert "IS NULL" not in strict["quality_proband"]


def test_pass_only_accepts_pass_unset_and_dot_only():
    sql = parts(GATES._replace(pass_only=True))["pass_only"]
    assert sql == "(v.filter IS NULL OR v.filter IN ('PASS', '.', ''))"


def test_the_parent_alt_screen_is_the_parent_ref_ok_guard_on_both_parents_denovo_only():
    gates = GATES._replace(parent_alt_screen=True, max_parent_alt_frac=0.07)
    denovo = trio_sql("denovo", TRIO, gates, count_only=True)
    assert "f.ad_alt IS NULL OR f.ad_alt < 2" in denovo
    assert "m.ad_alt IS NULL OR m.ad_alt < 2" in denovo
    assert "< 0.07" in denovo
    for other in ("recessive", "dominant"):
        assert "ad_alt < 2" not in trio_sql(other, TRIO, gates, count_only=True)
    # off: the text is the classic text
    assert "ad_alt < 2" not in trio_sql("denovo", TRIO, GATES, count_only=True)
