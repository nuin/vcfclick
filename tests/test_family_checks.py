"""Shared family-analysis helpers, exercised on real stores on both backends.

  * parent_ref_ok: when is a parent's stored 0/0 row trusted?
  * the null-safe genotype pivot (an absent parent is NULL, never 0)
  * family_has_reference on sparse vs --keep-reference stores
  * the cohort-AF probe and the `gates.af_gate` states it drives

Nothing opens a chDB session in the pytest process: stores are built with
the CLI and read in a fresh interpreter (`run_python`).
"""

from __future__ import annotations

import json
import textwrap

import pytest

from tests.conftest import run_cli
from tests.fixtures._vcfgen import (
    TRIO_PED,
    TRIO_SAMPLES,
    Site,
    trio_sites,
    write_vcf,
)

BACKENDS = ("chdb", "duckdb")


def _ingest(home, backend, vcf, *, keep_reference=True, name="t", ingest_id="i1"):
    run_cli(home, backend, "db", "create", name)
    args = ["db", "ingest", name, str(vcf), "--ingest-id", ingest_id, "--serial"]
    if keep_reference:
        args.append("--keep-reference")
    run_cli(home, backend, *args)


def _py(run_python, home, backend, code: str):
    out = run_python(home, textwrap.dedent(code), VCFCLICK_BACKEND=backend)
    return json.loads(out.strip().splitlines()[-1])


# ───────────────────────────── parent_ref_ok ─────────────────────────────

# (case id, FORMAT AD value or None for no AD, trusted at the default 0.05)
# One row per (ad_ref, ad_alt) shape the rule has to get right.
REF_OK_TABLE = [
    ("no_ad_field", None, True),
    ("ad_missing", ".", True),
    ("alt_missing_ref_5", "5,.", True),
    ("ref_missing_1_alt_read", ".,1", True),  # alt < 2: noise
    ("ref_missing_2_alt_reads", ".,2", False),  # NULL ref, >= 2 alt: not trusted
    ("ref_missing_5_alt_reads", ".,5", False),
    ("zero_reads", "0,0", True),
    ("0_alt_reads", "30,0", True),
    ("1_alt_read", "29,1", True),
    ("1_alt_read_in_20_is_exactly_0.05", "19,1", True),  # alt < 2 wins
    ("2_alt_reads_of_100", "98,2", True),  # 0.02
    ("2_alt_reads_of_40_is_exactly_0.05", "38,2", False),
    ("fraction_0.049", "951,49", True),
    ("fraction_0.05", "950,50", False),  # strict: not below the cut
    ("fraction_0.051", "949,51", False),
    ("2_alt_reads_no_ref", "0,2", False),  # fraction 1.0
    ("max_ref_2_alt", "65535,2", True),  # 65537 total: no USMALLINT overflow
    ("max_ref_1_alt", "65535,1", True),
    ("max_both", "65535,65535", False),  # 131070 total, fraction 0.5
    ("max_alt_zero_ref", "0,65535", False),
]

# Expected at other thresholds: {frac: {case id: trusted}} for the cases
# that change. 1.0 trusts everything whose ref count is present.
OTHER_FRACTIONS = {
    0.10: {"fraction_0.05": True, "fraction_0.051": True, "max_both": False},
    1.0: {
        "fraction_0.05": True,
        "2_alt_reads_no_ref": False,  # 1.0 is not < 1.0
        "max_both": True,
        "ref_missing_2_alt_reads": False,  # NULL ref is never trusted
    },
    0.0: {"0_alt_reads": True, "1_alt_read": True, "fraction_0.049": False},
}


def _ref_ok_vcf(tmp_path):
    sites = []
    for i, (_, ad, _ok) in enumerate(REF_OK_TABLE, start=1):
        gt = "0/0" if ad is None else f"0/0:{ad}"
        sites.append(
            Site("chr1", i, "A", "T", (gt,), fmt="GT" if ad is None else "GT:AD")
        )
    return write_vcf(tmp_path / "ref_ok.vcf", ["S"], sites)


@pytest.mark.parametrize("backend", BACKENDS)
def test_parent_ref_ok_table(backend, vcfclick_home, tmp_path, run_python):
    _ingest(vcfclick_home, backend, _ref_ok_vcf(tmp_path))
    fracs = [0.05, *OTHER_FRACTIONS]
    got = _py(
        run_python,
        vcfclick_home,
        backend,
        f"""
        import json
        from cli.trio_analysis import Gates, parent_ref_ok
        from storage import get_session

        sess = get_session("t")
        out = {{}}
        for frac in {fracs!r}:
            gates = Gates(20, 10, 0.01, 0.25, 0.75, max_parent_alt_frac=frac)
            expr = parent_ref_ok("g", gates)
            raw = sess.query(
                f"SELECT g.pos, CASE WHEN {{expr}} IS NULL THEN -1 "
                f"WHEN {{expr}} THEN 1 ELSE 0 END AS ok FROM genotypes g "
                "WHERE g.sample_id = 'S' AND g.gt = 0 ORDER BY g.pos",
                "JSONCompact",
            ).bytes().decode()
            out[str(frac)] = {{int(p): int(v) for p, v in json.loads(raw)["data"]}}
        print(json.dumps(out))
        """,
    )
    ids = {i: case for i, (case, _, _) in enumerate(REF_OK_TABLE, start=1)}

    def verdicts(frac):
        return {ids[int(p)]: v for p, v in got[str(frac)].items()}

    default = verdicts(0.05)
    assert len(default) == len(REF_OK_TABLE), "every row evaluated"
    # never NULL: a NULL would silently drop the row from the model
    assert -1 not in default.values()
    assert default == {case: int(ok) for case, _, ok in REF_OK_TABLE}
    for frac, changes in OTHER_FRACTIONS.items():
        v = verdicts(frac)
        assert -1 not in v.values(), frac
        for case, ok in changes.items():
            assert v[case] == int(ok), (frac, case)


def test_parent_ref_ok_rejects_a_fraction_outside_0_1():
    from cli.trio_analysis import Gates, parent_ref_ok

    for bad in (-0.01, 1.01, float("nan")):
        with pytest.raises(ValueError):
            parent_ref_ok("g", Gates(20, 10, 0.01, 0.25, 0.75, max_parent_alt_frac=bad))


def test_classic_models_do_not_read_the_parent_fraction():
    """Changing max_parent_alt_frac must not move the legacy SQL."""
    from cli.trio_analysis import Gates, Trio, trio_sql

    trio = Trio("i1", "C", "F", "M")
    a = Gates(20, 10, 0.01, 0.25, 0.75)
    b = a._replace(max_parent_alt_frac=0.5)
    for cat in ("denovo", "recessive", "dominant"):
        assert trio_sql(cat, trio, a, count_only=False) == trio_sql(
            cat, trio, b, count_only=False
        )


# ───────────────────────── null-safe genotype pivot ─────────────────────────


@pytest.mark.parametrize("backend", BACKENDS)
def test_pivot_reads_an_absent_parent_as_null_not_zero(
    backend, vcfclick_home, tmp_path, run_python
):
    """Site 100: father has a stored 0/0. Site 200: father is a no-call, so
    he has no row at all. A LEFT JOIN on genotypes reads the missing father
    as gt = 0 on chDB, which is how an absent parent looks like a proven
    reference; the pivot must give NULL there on both backends."""
    het, ref, nocall = "0/1:15,15:99:30", "0/0:30,0:99:30", "./."
    sites = [
        Site("chr1", 100, "A", "G", (het, ref, ref), fmt="GT:AD:GQ:DP"),
        Site("chr1", 200, "A", "G", (het, nocall, ref), fmt="GT:AD:GQ:DP"),
        Site("chr1", 300, "A", "G", (het, het, nocall), fmt="GT:AD:GQ:DP"),
    ]
    vcf = write_vcf(tmp_path / "pivot.vcf", TRIO_SAMPLES, sites)
    _ingest(vcfclick_home, backend, vcf)
    got = _py(
        run_python,
        vcfclick_home,
        backend,
        """
        import json
        from storage import get_session
        from storage.sql import genotype_pivot_sql

        pivot = genotype_pivot_sql(
            "i1", {"p": "CHILD", "f": "FATHER", "m": "MOTHER"}, ("gt",)
        )
        raw = get_session("t").query(
            "SELECT x.pos AS pos, x.p_gt AS p, x.f_gt AS f, x.m_gt AS m "
            f"FROM ({pivot}) x ORDER BY x.pos", "JSONCompact"
        ).bytes().decode()
        print(json.dumps(json.loads(raw)["data"]))
        """,
    )
    assert [[int(r[0]), *r[1:]] for r in got] == [
        [100, 1, 0, 0],  # stored reference calls read 0
        [200, 1, None, 0],  # father no-call: NULL, not 0
        [300, 1, 1, None],  # mother no-call
    ]


# ─────────────────────────── family_has_reference ───────────────────────────


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("keep_reference", [True, False])
def test_family_has_reference_on_sparse_vs_keep_reference(
    backend, keep_reference, vcfclick_home, tmp_path, run_python
):
    vcf = write_vcf(tmp_path / "t.vcf", TRIO_SAMPLES, trio_sites(4, af="none"))
    _ingest(vcfclick_home, backend, vcf, keep_reference=keep_reference)
    got = _py(
        run_python,
        vcfclick_home,
        backend,
        """
        import json
        from cli.family_checks import family_has_reference, samples_with_reference
        from storage import get_session

        sess = get_session("t")
        print(json.dumps({
            "per_sample": samples_with_reference(sess, "i1", ["CHILD", "FATHER", "NOBODY"]),
            "family": family_has_reference(sess, "i1", ["CHILD", "FATHER", "MOTHER"]),
            "empty": family_has_reference(sess, "i1", []),
            "other_ingest": family_has_reference(sess, "zzz", ["FATHER"]),
        }))
        """,
    )
    # CHILD is het at every site, so it never has a 0/0 row of its own; the
    # parents do, but only in a --keep-reference store.
    assert got["per_sample"] == {
        "CHILD": False,
        "FATHER": keep_reference,
        "NOBODY": False,
    }
    assert got["family"] is False  # CHILD has no reference row in either store
    assert got["empty"] is False  # nothing listed proves nothing
    assert got["other_ingest"] is False


@pytest.mark.parametrize("backend", BACKENDS)
def test_family_has_reference_true_when_everyone_has_a_reference_row(
    backend, vcfclick_home, tmp_path, run_python
):
    ref, het = "0/0:30,0:99:30", "0/1:15,15:99:30"
    sites = [
        Site("chr1", 100, "A", "G", (het, ref, ref), fmt="GT:AD:GQ:DP"),
        Site("chr1", 200, "A", "G", (ref, het, ref), fmt="GT:AD:GQ:DP"),
    ]
    _ingest(vcfclick_home, backend, write_vcf(tmp_path / "t.vcf", TRIO_SAMPLES, sites))
    got = _py(
        run_python,
        vcfclick_home,
        backend,
        """
        import json
        from cli.family_checks import family_has_reference
        from storage import get_session

        print(json.dumps(family_has_reference(
            get_session("t"), "i1", ["CHILD", "FATHER", "MOTHER"]
        )))
        """,
    )
    assert got is True


# ─────────────────────── cohort-derived allele frequency ───────────────────────


def _trio_store(home, backend, tmp_path, *, af, n=60, keep_reference=True):
    vcf = write_vcf(tmp_path / f"trio_{af}.vcf", TRIO_SAMPLES, trio_sites(n, af=af))
    ped = tmp_path / "trio.ped"
    ped.write_text(TRIO_PED)
    _ingest(home, backend, vcf, keep_reference=keep_reference)
    run_cli(home, backend, "db", "ped", "t", str(ped))


def _trio_json(home, backend, *extra):
    r = run_cli(
        home,
        backend,
        "db",
        "trio",
        "t",
        "--proband",
        "CHILD",
        "--format",
        "json",
        *extra,
    )
    return json.loads(r.stdout), r.stderr


@pytest.mark.parametrize("backend", BACKENDS)
def test_cohort_derived_af_is_detected_and_skipped(backend, vcfclick_home, tmp_path):
    """A callset-AF trio: every private het has AF 1/6 or 1/3, so the default
    --max-af 0.01 would drop all of them. The probe sees AF == ac/an at every
    site, skips the gate, and says so."""
    _trio_store(vcfclick_home, backend, tmp_path, af="cohort")
    doc, err = _trio_json(vcfclick_home, backend)
    assert doc["gates"]["af_gate"] == "skipped_cohort_derived"
    assert doc["gates"]["max_af"] == 0.01  # the number is still reported
    assert doc["models"]["denovo"]["count"] == 30  # every de novo site survives
    assert "cohort" in err and "--max-af" in err and "skipped" in err
    assert "--gnomad-max-af" in err

    # Text mode prints the same note on stderr and keeps stdout clean.
    r = run_cli(vcfclick_home, backend, "db", "trio", "t", "--proband", "CHILD")
    assert "skipped" in r.stderr and "note:" not in r.stdout
    assert "denovo" in r.stdout


@pytest.mark.parametrize("backend", BACKENDS)
def test_explicit_max_af_is_always_honoured(backend, vcfclick_home, tmp_path):
    _trio_store(vcfclick_home, backend, tmp_path, af="cohort")
    doc, err = _trio_json(vcfclick_home, backend, "--max-af", "0.01")
    assert doc["gates"]["af_gate"] == "cohort_derived_suspected"
    assert doc["models"]["denovo"]["count"] == 0  # the gate ran, as asked
    assert "applied as requested" in err
    # a threshold above the cohort values keeps them, still reported
    doc, _ = _trio_json(vcfclick_home, backend, "--max-af", "0.5")
    assert doc["gates"]["af_gate"] == "cohort_derived_suspected"
    assert doc["models"]["denovo"]["count"] == 30


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_population_style_af_is_applied(backend, vcfclick_home, tmp_path):
    """Hand-set, population-like frequencies (0.0005 / 0.002): not the
    cohort's own, so the default gate runs and passes them."""
    _trio_store(vcfclick_home, backend, tmp_path, af="popular")
    doc, err = _trio_json(vcfclick_home, backend)
    assert doc["gates"]["af_gate"] == "applied"
    assert doc["models"]["denovo"]["count"] == 30
    assert "cohort" not in err
    doc, _ = _trio_json(vcfclick_home, backend, "--max-af", "0.001")
    assert doc["gates"]["af_gate"] == "applied"
    assert doc["models"]["denovo"]["count"] == 30  # 0.0005 <= 0.001
    assert doc["models"]["dominant"]["count"] == 0  # 0.002 > 0.001: gate bites


@pytest.mark.parametrize("backend", BACKENDS)
def test_probe_unavailable_without_af_or_with_too_few_sites(
    backend, vcfclick_home, tmp_path
):
    """No INFO/AF at all, or fewer than 20 sites to judge: the probe cannot
    say, the gate runs as it always did."""
    _trio_store(vcfclick_home, backend, tmp_path, af="none")
    doc, err = _trio_json(vcfclick_home, backend)
    assert doc["gates"]["af_gate"] == "unavailable" and "cohort" not in err


@pytest.mark.parametrize("backend", BACKENDS)
def test_probe_unavailable_with_few_sites(backend, vcfclick_home, tmp_path):
    _trio_store(vcfclick_home, backend, tmp_path, af="cohort", n=10)
    doc, err = _trio_json(vcfclick_home, backend)
    assert doc["gates"]["af_gate"] == "unavailable"
    assert "cohort" not in err
    # the gate ran: cohort AFs 1/6 and 1/3 are above 0.01
    assert doc["models"]["denovo"]["count"] == 0


@pytest.mark.parametrize("backend", BACKENDS)
def test_probe_unavailable_on_a_database_without_called_counts(
    backend, vcfclick_home, tmp_path, run_python
):
    """A store ingested before an_called / ac_called existed."""
    _trio_store(vcfclick_home, backend, tmp_path, af="cohort")
    run_python(
        vcfclick_home,
        textwrap.dedent(
            """
            from storage import get_session
            s = get_session("t")
            for col in ("an_called", "ac_called", "n_called"):
                s.query(f"ALTER TABLE variants DROP COLUMN {col}")
            """
        ),
        VCFCLICK_BACKEND=backend,
    )
    doc, err = _trio_json(vcfclick_home, backend)
    assert doc["gates"]["af_gate"] == "unavailable"
    assert doc["models"]["denovo"]["count"] == 0  # gate applied as before


def test_resolve_af_gate_covers_the_whole_enum():
    from cli.family_checks import CohortAfProbe, resolve_af_gate
    from cli.trio_analysis import AF_GATE_VALUES

    off = CohortAfProbe(False, reason="x")
    like = CohortAfProbe(True, 100, 95)
    unlike = CohortAfProbe(True, 100, 20)
    edge = CohortAfProbe(True, 100, 90)  # exactly 90% is suspected
    assert resolve_af_gate(off, explicit=False) == "unavailable"
    assert resolve_af_gate(off, explicit=True) == "unavailable"
    assert resolve_af_gate(unlike, explicit=False) == "applied"
    assert resolve_af_gate(like, explicit=False) == "skipped_cohort_derived"
    assert resolve_af_gate(like, explicit=True) == "cohort_derived_suspected"
    assert resolve_af_gate(edge, explicit=False) == "skipped_cohort_derived"
    assert AF_GATE_VALUES == (
        "applied",
        "cohort_derived_suspected",
        "skipped_cohort_derived",
        "unavailable",
    )
