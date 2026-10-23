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


OFF = GATES._replace(min_gq=0, min_dp=0, min_ab=0.0, max_ab=1.0)


def test_strict_missing_holds_a_call_only_to_the_gates_that_are_asked_for():
    """A minimum of 0 and a band of 0..1 ask nothing of the field, so a call
    without it is not a fail: strict on top of them is the text without strict."""
    assert not OFF.ab_requested and GATES.ab_requested
    strict = parts(OFF._replace(strict_missing=True))
    assert strict == parts(OFF)
    assert "IS NOT NULL" not in strict["quality_proband"] + strict["allele_balance"]
    # one asked for, the other not
    gq_only = parts(OFF._replace(strict_missing=True, min_gq=20))
    assert "(g.gq IS NOT NULL AND g.gq >= 20)" in gq_only["quality_proband"]
    assert "(g.dp IS NULL OR g.dp >= 0)" in gq_only["quality_proband"]
    assert "ad_ref IS NOT NULL" not in gq_only["allele_balance"]
    # a band narrower than 0..1, on either side, is a gate
    for band in ({"min_ab": 0.1}, {"max_ab": 0.9}):
        with_band = parts(OFF._replace(strict_missing=True, **band))
        assert "ad_ref IS NOT NULL" in with_band["allele_balance"]


def test_a_gate_exempted_for_one_run_is_tolerant_again():
    gates = GATES._replace(strict_missing=True, strict_exempt=("min_gq", "min_ab"))
    assert [gates.strict_for(k) for k in ("min_gq", "min_dp", "min_ab")] == [
        False, True, False,
    ]  # fmt: skip
    p = parts(gates)
    assert "(g.gq IS NULL OR g.gq >= 20)" in p["quality_proband"]
    assert "(g.dp IS NOT NULL AND g.dp >= 10)" in p["quality_proband"]
    assert p["allele_balance"] == parts()["allele_balance"]
    # no strict_missing: nothing to exempt from
    assert not GATES._replace(strict_exempt=("min_gq",)).strict_for("min_dp")


def test_strict_missing_makes_the_parent_screen_ask_for_allele_depths():
    gates = GATES._replace(parent_alt_screen=True)
    plain = parts(gates, parent_alt_screen=True)["parent_alt_screen"]
    strict_gates = gates._replace(strict_missing=True)
    strict = parts(strict_gates, parent_alt_screen=True)["parent_alt_screen"]
    assert "f.ad_alt IS NULL OR" in plain  # a parent 0/0 with no depths is trusted
    assert not plain.startswith("(f.ad_ref IS NOT NULL AND f.ad_alt IS NOT NULL)")
    assert strict.startswith(
        "(f.ad_ref IS NOT NULL AND f.ad_alt IS NOT NULL) AND "
        "(m.ad_ref IS NOT NULL AND m.ad_alt IS NOT NULL) AND "
    )
    assert strict.endswith(plain)
    exempt = strict_gates._replace(strict_exempt=("parent_alt_screen",))
    assert parts(exempt, parent_alt_screen=True)["parent_alt_screen"] == plain


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
