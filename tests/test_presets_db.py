"""Presets against real data, on both backends.

One ingest per backend and fixture (module-scoped), then `db trio` runs against
it, about 40 per backend; everything that does not need a database is in the
in-process preset tests. The fixture is tests/fixtures/make_preset_vcf.py: sites
planted exactly on the boundaries of the built-in presets' gates, with the
expected rows written by hand in presets_truth.json. Assertions are on planted
POSITIONS, never whole-file totals.

related.vcf cannot test a gate (it has no GQ, DP or INFO/AF, and the gates are
NULL-tolerant), which `test_related_vcf_cannot_tell_a_strict_gate_from_a_loose_one`
shows, so a test written against it would pass vacuously.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from unittest import mock

import pytest
from click.testing import CliRunner

from cli import filter_registry as fr
from cli.builtin_presets import BUILTINS
from cli.main import cli
from tests.conftest import FIXTURES, run_cli
from tests.fixtures._vcfgen import TRIO_PED, TRIO_SAMPLES, trio_sites, write_vcf
from tests.test_trio_contract import assert_trio_document

BACKENDS = ("chdb", "duckdb")
TRUTH = json.loads((FIXTURES / "presets_truth.json").read_text())
POS: dict[str, int] = TRUTH["positions"]
NAME_AT = {v: k for k, v in POS.items()}
DE_NOVO = {name for name in POS if name.startswith("dn_")}
KEPT = TRUTH["kept"]
ENGINE_NOISE = re.compile(r"^Default dirty pages decay period: \d+ms\n", re.M)


@dataclass
class World:
    home: Path
    backend: str
    ann: dict[str, Path]


@pytest.fixture(scope="module", params=BACKENDS)
def backend(request):
    return request.param


@pytest.fixture(scope="module")
def world(backend, tmp_path_factory):
    from annotations.db import get_connection

    base = tmp_path_factory.mktemp(f"presets-{backend}")
    home = base / "home"
    home.mkdir()
    ann = {
        "loaded": base / "loaded.duckdb",
        "empty": base / "empty.duckdb",
        # one gnomAD record that reaches no candidate: on a chromosome the
        # cohort does not have, and on its chromosome at a position it does not
        "wrong": base / "wrong.duckdb",
        "far": base / "far.duckdb",
        # a gnomAD that holds the cohort's inherited variants and none of its
        # de novo ones: it judges the cohort, and the candidates are novel
        "novel": base / "novel.duckdb",
    }
    with mock.patch.dict(os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(ann["loaded"])}):
        conn = get_connection()  # in-process DuckDB, never a chDB session
        sym, chrom, start, end = (
            TRUTH["gene"][k] for k in ("symbol", "chrom", "start", "end")
        )
        conn.execute(
            "INSERT INTO refseq_genes VALUES (?, ?, ?, ?, '+', 'id', 'protein_coding')",
            [sym, chrom, start, end],
        )
        for pos, popmax in TRUTH["gnomad_popmax"].items():
            conn.execute(
                "INSERT INTO gnomad_af VALUES ('chr1', ?, 'A', 'G', ?, ?)",
                [int(pos), popmax / 2, popmax],
            )
        # the 30 common filler sites are in gnomAD too, as common variants are: a
        # real load has a record for most of a cohort's variants, which is what
        # the pre-run probe of the cohort's own variants asks for
        for i in range(30):
            conn.execute(
                "INSERT INTO gnomad_af VALUES ('chr1', ?, 'A', 'G', 0.1, 0.2)",
                [20000 + 100 * i],
            )
        conn.close()
    with mock.patch.dict(os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(ann["empty"])}):
        get_connection().close()
    for key, chrom in (("wrong", "chr9"), ("far", "chr1")):
        with mock.patch.dict(os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(ann[key])}):
            conn = get_connection()
            conn.execute(
                "INSERT INTO gnomad_af VALUES (?, 5, 'A', 'G', 0.0001, 0.0002)", [chrom]
            )
            conn.close()

    with mock.patch.dict(os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(ann["novel"])}):
        conn = get_connection()
        for i in range(1, 30, 2):  # trio_sites: the odd sites are the inherited ones
            conn.execute(
                "INSERT INTO gnomad_af VALUES ('chr1', ?, 'A', 'G', 0.0001, 0.0002)",
                [1000 + i * 1000],
            )
        conn.close()

    ped = FIXTURES / "presets.ped"
    cohort = write_vcf(
        base / "cohort.vcf", TRIO_SAMPLES, trio_sites(30, af="cohort")
    )  # INFO/AF is the callset's own frequency, as GATK writes it
    cohort_ped = base / "cohort.ped"
    cohort_ped.write_text(TRIO_PED)
    stores = (
        ("p", FIXTURES / "presets.vcf.gz", ped),
        ("n", FIXTURES / "presets_nogq.vcf.gz", ped),
        ("c", cohort, cohort_ped),
    )
    for name, vcf, pedigree in stores:
        run_cli(home, backend, "db", "create", name)
        run_cli(
            home, backend, "db", "ingest", name, str(vcf),
            "--cohort", "fam", "--ingest-id", "fam1", "--serial", "--keep-reference",
        )  # fmt: skip
        run_cli(home, backend, "db", "ped", name, str(pedigree))
    return World(home, backend, ann)


def trio(w: World, db: str, *args: str, ann: str = "loaded", ok: bool = True):
    env = {"VCFCLICK_ANNOTATIONS_DB": str(w.ann[ann])}
    with mock.patch.dict(os.environ, env):
        return run_cli(
            w.home, w.backend, "db", "trio", db, "--proband", "CHILD", *args, ok=ok
        )


def jtrio(w: World, db: str, *args: str, ann: str = "loaded") -> dict:
    r = trio(w, db, *args, "--format", "json", "--limit", "200", ann=ann)
    doc = json.loads(r.stdout)
    assert_trio_document(doc)  # the structural contract the apps decode, every time
    return doc


def names(doc: dict, model: str) -> set[str]:
    """The planted sites (or comphet genes) a model reports."""
    m = doc["models"][model]
    if "genes" in m:
        return {g["gene"] for g in m["genes"]}
    return {NAME_AT[c["pos"]] for c in m["candidates"]}


def quiet(stderr: str) -> str:
    return ENGINE_NOISE.sub("", stderr)


def flags_for(filters: dict) -> list[str]:
    """The `db trio` flags that say what a preset's filters say."""
    out = []
    for key, value in filters.items():
        spec = fr.REGISTRY[key]
        if spec.type == "category":
            out += [spec.flag, ",".join(value)]
        elif spec.is_bool:
            on, off = spec.flag.split("/")
            out.append(on if value else off)
        else:
            out += [spec.flag, repr(value)]
    return out


@pytest.fixture(scope="module")
def preset_runs(world):
    """`db trio --preset NAME` once per built-in, with the annotations loaded."""
    return {n: jtrio(world, "p", "--preset", n) for n in BUILTINS}


# ───────────────────── what each built-in keeps, by hand ─────────────────────


@pytest.mark.parametrize("name", list(BUILTINS))
def test_each_builtin_keeps_exactly_the_planted_rows(name, preset_runs):
    doc = preset_runs[name]
    expected = KEPT[name]
    assert set(doc["models"]) == set(expected), (
        "the preset's category is the model list"
    )
    for model, kept in expected.items():
        assert names(doc, model) == set(kept), (name, model)
    assert doc["preset"]["name"] == name
    assert doc["preset"]["fingerprint"] == BUILTINS[name].history[-1]["fingerprint"]


def test_strict_denovo_is_stricter_than_the_others_on_every_planted_boundary(
    preset_runs,
):
    strict = names(preset_runs["strict-denovo"], "denovo")
    broad = names(preset_runs["broad-first-pass"], "denovo")
    dominant = names(preset_runs["rare-dominant"], "denovo")
    assert strict < dominant < broad  # nested, and strictly so
    # the sites each strict gate exists for are in broad and not in strict
    for site in ("dn_gq29", "dn_dp19", "dn_ab029", "dn_ab071", "dn_af0011",
                 "dn_parent_alt", "dn_nonpass", "dn_null_dp_parent", "dn_gnomad_common"):  # fmt: skip
        assert site in broad and site not in strict, site
    # ... and the site just inside each boundary survives
    for site in ("dn_gq30", "dn_dp20", "dn_ab031", "dn_ab069", "dn_af0009",
                 "dn_parent_alt_low", "dn_gnomad_rare"):  # fmt: skip
        assert site in strict, site


@pytest.mark.parametrize("name", list(BUILTINS))
def test_a_preset_gives_the_same_answer_as_its_values_typed_as_flags(
    name, world, preset_runs
):
    """The equivalence oracle: on a store whose INFO/AF is a population
    frequency, `--preset P` and the same values as flags differ only in the
    `preset` and `effective` objects."""
    via_flags = jtrio(world, "p", *flags_for(BUILTINS[name].filters))
    via_preset = json.loads(json.dumps(preset_runs[name]))
    for key in ("preset", "effective"):
        via_preset.pop(key)
        via_flags.pop(key, None)  # a flag run that insists computes `effective` too
    # a preset (unlike a bare flag) also shows the parents' reads on de novo rows
    for doc in (via_preset, via_flags):
        for model in doc["models"].values():
            for candidate in model.get("candidates", ()):
                candidate.pop("parent_alt_reads", None)
    assert via_preset == via_flags


def test_related_vcf_cannot_tell_a_strict_gate_from_a_loose_one(
    vcfclick_home, monkeypatch, tmp_path
):
    """Why the fixture above exists: related.vcf has no GQ, DP or INFO/AF, and
    the gates are NULL-tolerant, so the same assertions would pass vacuously
    there: its counts do not move however strict the gates are."""
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    h, backend = vcfclick_home, "duckdb"  # the gates are the same SQL on both
    run_cli(h, backend, "db", "create", "t")
    run_cli(
        h, backend, "db", "ingest", "t", str(FIXTURES / "related.vcf"),
        "--ingest-id", "g1", "--serial", "--keep-reference",
    )  # fmt: skip
    run_cli(
        h,
        backend,
        "db",
        "ped",
        "t",
        "--proband",
        "C1",
        "--father",
        "F",
        "--mother",
        "M",
    )

    def counts(*gates):
        r = run_cli(
            h, backend, "db", "trio", "t", "--proband", "C1", "--format", "json", *gates
        )
        return {m: v["count"] for m, v in json.loads(r.stdout)["models"].items()}

    loose = counts()
    assert loose["denovo"] == 3 and loose["recessive"] > 0
    assert counts("--min-gq", "99", "--min-dp", "60", "--max-af", "0.0001") == loose


# ───────────────────────────── the echo ─────────────────────────────


def test_a_plain_run_carries_no_preset_block_and_prints_no_new_stderr(world):
    r = trio(world, "p", "--format", "json")
    doc = json.loads(r.stdout)
    assert list(doc) == ["trio", "keep_reference", "gates", "annotations", "models"]
    assert list(doc["gates"]) == [
        "min_gq", "min_dp", "max_af", "min_ab", "max_ab", "gnomad_max_af", "af_gate"
    ]  # fmt: skip
    assert quiet(r.stderr) == ""
    assert "parent_alt_reads" not in doc["models"]["denovo"]["candidates"][0]


def test_the_preset_object_says_what_ran_and_where_each_value_came_from(world):
    doc = jtrio(world, "p", "--preset", "rare-recessive", "--max-af", "0.005")
    p = doc["preset"]
    assert list(p) == [
        "name", "source", "revision", "status", "schema_version", "fingerprint",
        "pinned", "meta", "path", "requested", "resolved", "overridden", "skipped",
        "models_not_run", "covers", "does_not_cover", "notes",
    ]  # fmt: skip
    assert (p["source"], p["revision"], p["status"]) == ("builtin", 1, "starting_point")
    assert p["path"] is None and p["pinned"] is False
    assert p["meta"] == {"approved_by": None, "approved_on": None, "locked": False}
    assert p["resolved"]["max_af"] == {"value": 0.005, "from": "flag", "bound": True}
    assert p["resolved"]["min_gq"] == {"value": 20, "from": "preset", "bound": True}
    assert p["resolved"]["min_dp"]["from"] == "preset"
    assert p["overridden"] == [{"key": "max_af", "preset": 0.01, "flag": 0.005}]
    assert [m["model"] for m in p["models_not_run"]] == ["denovo", "dominant"]
    assert (
        doc["gates"]["max_af"] == 0.005
    )  # the gates block carries the effective values
    assert list(doc["effective"]) == [
        "rows", "gq", "dp", "ad", "median_dp", "info_af", "filter_column", "gnomad",
        "population_filter", "warnings",
    ]  # fmt: skip
    assert doc["effective"]["gnomad"]["candidates_without_record"] >= 1


def test_text_has_the_preset_and_gates_lines_after_the_trio_line(world):
    out = trio(world, "p", "--preset", "strict-denovo", "--max-af", "0.0005").stdout
    lines = out.splitlines()
    assert lines[0].startswith("trio: proband=CHILD")
    assert lines[1].startswith("preset: strict-denovo (builtin, rev 1, ")
    assert lines[1].endswith("starting point)")
    assert lines[2].startswith("  gates: ")
    assert "max_af=0.0005 [flag; preset 0.001]" in lines[2]
    assert "min_gq=30 [preset]" in lines[2]
    assert "strict_missing=on [preset]" in lines[2]
    # the de novo listing shows the parents' reads at each site
    assert "chr1:7100 A>G" in out
    assert "parents' alt reads: father 1/40, mother 0/40" in out
    assert "chr1:6000" not in out  # AF 0.0009 is over the typed 0.0005
    # the gnomAD slice holds two of the candidates: the rest it cannot judge, and
    # the line says how much of the cohort it does reach
    assert lines[-1] == (
        "  gnomAD: no record for 6 of 8 candidates (kept as rare); it has a "
        "record for 62.7% of a sample of the cohort's variants"
    )


def test_parent_allele_depths_are_shown_on_de_novo_rows(preset_runs):
    for c in preset_runs["strict-denovo"]["models"]["denovo"]["candidates"]:
        assert set(c["parent_alt_reads"]) == {"father", "mother"}
    reads = {
        NAME_AT[c["pos"]]: c["parent_alt_reads"]["father"]
        for c in preset_runs["strict-denovo"]["models"]["denovo"]["candidates"]
    }
    assert reads["dn_parent_alt_low"] == {"ref": 39, "alt": 1}
    assert reads["dn_ok"] == {"ref": 40, "alt": 0}
    # a broad pass shows the screen's reason for dropping dn_parent_alt: 4 of 40
    broad = {
        NAME_AT[c["pos"]]: c["parent_alt_reads"]
        for c in preset_runs["broad-first-pass"]["models"]["denovo"]["candidates"]
    }
    assert broad["dn_parent_alt"]["father"] == {"ref": 36, "alt": 4}


# ───────────────── the new gates, alone, without a preset ─────────────────


def test_the_opt_in_gates_remove_exactly_their_planted_site(world):
    base = names(jtrio(world, "p", "--category", "denovo"), "denovo")
    assert {"dn_parent_alt", "dn_nonpass", "dn_null_dp_parent"} <= base
    for flag, site in (
        ("--pass-only", "dn_nonpass"),
        ("--parent-alt-screen", "dn_parent_alt"),
        ("--strict-missing", "dn_null_dp_parent"),
    ):
        doc = jtrio(world, "p", "--category", "denovo", flag)
        assert names(doc, "denovo") == base - {site}, flag
        assert doc["gates"][flag[2:].replace("-", "_")] is True
        if flag == "--parent-alt-screen":
            assert doc["gates"]["max_parent_alt_frac"] == 0.05
            assert doc["effective"]["warnings"] == []  # it binds here: AD is in the VCF
            assert "preset" not in doc
            assert "parent_alt_reads" in doc["models"]["denovo"]["candidates"][0]
        elif flag == "--pass-only":
            assert "effective" not in doc  # no preset, nothing insisted on


def test_a_parent_alt_screen_that_screens_nothing_says_so_and_says_only_that(world):
    """The store has no allele depths: every parent 0/0 would be trusted. A typed
    opt-in gate that judges nothing must not look like a screen that worked."""
    doc = jtrio(world, "n", "--category", "denovo", "--parent-alt-screen")
    assert names(doc, "denovo") == DE_NOVO  # nothing was screened out
    (w,) = doc["effective"]["warnings"]  # not the GQ / DP gates it did not ask for
    assert (w["code"], w["key"]) == ("gate_unbound", "parent_alt_screen")
    assert "every parent 0/0 would be trusted" in w["message"]
    text = trio(world, "n", "--category", "denovo", "--parent-alt-screen").stdout
    assert [x for x in text.splitlines() if x.startswith("  WARNING:")] == [
        f"  WARNING: {w['message']}"
    ]  # on stdout, and it did not refuse: nothing insisted on the gate


# ─────────────────────────── gates that cannot bind ───────────────────────────


def test_a_warn_only_preset_runs_on_data_without_gq_dp_ad_and_says_so(world):
    r = trio(world, "n", "--preset", "broad-first-pass")
    warnings = [x for x in r.stdout.splitlines() if x.startswith("  WARNING:")]
    assert len(warnings) == 3  # on stdout, not only stderr
    for field in ("GQ", "DP", "AD"):
        assert any(f"carry {field}" in x for x in warnings), field
    counts = {
        x.split()[0]: x.split()[1]
        for x in r.stdout.splitlines()
        if x.split()[:1] in (["denovo"], ["recessive"], ["dominant"], ["comphet"])
    }
    assert counts["denovo"] == "17"  # tolerant gates: it ran, and removed nothing


def test_a_strict_preset_refuses_to_run_when_a_gate_cannot_bind(world):
    r = trio(world, "n", "--preset", "strict-denovo", ok=False)
    assert r.stdout == ""
    for text in ("preset 'strict-denovo' needs every gate to bind", "min_gq=30",
                 "min_dp=20", "parent_alt_screen", "Remedy:", "broad-first-pass",
                 "--no-require-bound-gates"):  # fmt: skip
        assert text in r.stderr, text
    # ... and says what the gate does with a call that lacks the field
    assert "strict_missing counts a call without GQ as a fail" in r.stderr
    # strict_missing alone insists on the GQ / DP / AD gates, with no preset
    alone = trio(world, "n", "--strict-missing", ok=False)
    assert "this run needs every gate to bind" in alone.stderr
    assert "NO population frequency filter" not in alone.stderr


def test_no_require_bound_gates_overrides_the_refusal_for_one_run(world):
    doc = jtrio(world, "n", "--preset", "strict-denovo", "--no-require-bound-gates")
    assert doc["preset"]["resolved"]["require_bound"] == {
        "value": False,
        "from": "flag",
    }
    assert {"key": "require_bound", "preset": True, "flag": False} in doc["preset"][
        "overridden"
    ]
    # "run anyway" is not "fail every call": the gates that cannot bind are
    # exempt from strict_missing, so the other gates decide, as they would on
    # data that carries the fields
    assert names(doc, "denovo") == DE_NOVO - {
        "dn_nonpass",  # pass_only still binds
        "dn_af0011",  # so does max_af
        "dn_gnomad_common",  # and gnomAD
    }
    assert doc["gates"]["strict_missing"] is True
    assert doc["gates"]["strict_exempt"] == [
        "min_gq", "min_dp", "min_ab", "parent_alt_screen"
    ]  # fmt: skip
    eff = doc["effective"]
    assert [eff[k]["bound"] for k in ("gq", "dp", "ad")] == [False, False, False]
    assert (
        eff["info_af"]["bound"] is True
    )  # INFO/AF is there; the format fields are not
    assert {w["code"] for w in eff["warnings"]} == {"gate_unbound"}
    for w in eff["warnings"]:  # and each says the gate lets the call through
        assert "strict_missing is not applied to it for this run" in w["message"]
    text = trio(world, "n", "--preset", "strict-denovo", "--no-require-bound-gates")
    assert "a call without GQ passes" in text.stdout
    assert "denovo candidates: 14" in text.stdout


def test_strict_missing_is_only_held_against_the_gates_that_are_asked_for(world):
    """The store has no GQ, DP or AD at all. Asking for none of them, strict_missing
    must change nothing, and must not claim a gate it cannot judge."""
    none = ["--min-gq", "0", "--min-dp", "0", "--min-ab", "0", "--max-ab", "1"]
    plain = jtrio(world, "n", "--category", "denovo", *none)
    strict = jtrio(world, "n", "--category", "denovo", "--strict-missing", *none)
    assert names(strict, "denovo") == names(plain, "denovo") == DE_NOVO
    assert strict["effective"]["warnings"] == []
    assert "strict_exempt" not in strict["gates"]
    # the way out the refusal gives (flags that do not need the fields) works
    via_flags = jtrio(
        world, "n", "--preset", "strict-denovo", *none, "--no-parent-alt-screen",
        "--preset-ignore", "gnomad_max_af",
    )  # fmt: skip
    assert names(via_flags, "denovo") == DE_NOVO - {"dn_nonpass", "dn_af0011"}
    assert via_flags["effective"]["warnings"] == []


def test_depth_missing_on_a_few_rows_binds_and_is_reported_as_partial(preset_runs):
    eff = preset_runs["strict-denovo"]["effective"]
    assert eff["dp"]["bound"] is True and 0.98 < eff["dp"]["father"] < 1.0
    (partial,) = [w for w in eff["warnings"] if w["code"] == "gate_partial"]
    assert partial["key"] == "min_dp" and "counts as a fail" in partial["message"]
    assert eff["rows"] == {"proband": 51, "father": 51, "mother": 51}
    assert eff["median_dp"]["proband"] == 40


def test_a_median_depth_near_min_dp_is_a_warning(world):
    doc = jtrio(world, "p", "--preset", "broad-first-pass", "--min-dp", "30")
    (w,) = [x for x in doc["effective"]["warnings"] if x["code"] == "thin_coverage"]
    assert "median proband DP is 40 but min_dp is 30" in w["message"]


def test_gnomad_asked_for_and_not_loaded(world):
    doc = jtrio(world, "p", "--preset", "broad-first-pass", ann="empty")
    codes = [w["code"] for w in doc["effective"]["warnings"]]
    assert "gnomad_not_loaded" in codes
    assert doc["effective"]["gnomad"]["loaded"] is False
    assert doc["effective"]["population_filter"] is True  # INFO/AF still narrows
    r = trio(world, "p", "--preset", "strict-denovo", ann="empty", ok=False)
    assert "gnomad_max_af=0.0001 cannot be enforced" in r.stderr
    assert "annotations load-gnomad" in r.stderr


def test_a_gnomad_slice_on_other_chromosomes_is_not_a_gate_that_binds(world):
    r = trio(world, "c", "--preset", "strict-denovo", ann="wrong", ok=False)
    assert "no record on any chromosome this cohort has variants on" in r.stderr
    assert "(cohort: chr1; gnomAD: chr9)" in r.stderr
    assert "Remedy: load gnomAD for this cohort's chromosomes" in r.stderr
    doc = jtrio(world, "c", "--preset", "broad-first-pass", ann="wrong")  # warns
    eff = doc["effective"]
    assert (eff["gnomad"]["loaded"], eff["gnomad"]["bound"]) == (True, False)
    assert "gnomad_no_coverage" in [w["code"] for w in eff["warnings"]]
    assert eff["population_filter"] is False


def test_a_gnomad_with_one_row_on_the_right_chromosome_is_not_a_gate_that_binds(world):
    """The chromosome matches and a record exists, but gnomAD has a record for none
    of the cohort's variants: the filter judges nothing, and a run that insists on
    it refuses before it starts (and a broad one says so, in the same document)."""
    r = trio(world, "c", "--preset", "strict-denovo", ann="far", ok=False)
    assert r.stdout == ""
    assert "gnomad_max_af=0.0001 cannot be enforced" in r.stderr
    assert (
        "a record for only 0% of a sample of 30 of this cohort's variants" in r.stderr
    )
    assert "Remedy: load gnomAD for this cohort's chromosomes" in r.stderr
    # the same single row does not satisfy it on the GQ-bearing fixture either
    refused = trio(world, "p", "--preset", "strict-denovo", ann="far", ok=False)
    assert "a record for only 0% of a sample of 51" in refused.stderr
    # --no-require-bound-gates runs it anyway, and the document agrees with itself
    doc = jtrio(
        world, "c", "--preset", "strict-denovo", "--no-require-bound-gates", ann="far"
    )
    eff = doc["effective"]
    assert (eff["gnomad"]["loaded"], eff["gnomad"]["bound"]) == (True, False)
    assert (
        eff["gnomad"]["cohort_variants_sampled"],
        eff["gnomad"]["cohort_variants_with_record"],
    ) == (30, 0)
    assert "gnomad_no_coverage" in [w["code"] for w in eff["warnings"]]
    assert eff["population_filter"] is False
    # a broad preset only warns
    broad = jtrio(world, "c", "--preset", "broad-first-pass", ann="far")
    assert broad["effective"]["gnomad"]["bound"] is False
    assert "gnomad_no_coverage" in [w["code"] for w in broad["effective"]["warnings"]]


def test_a_gnomad_that_judges_the_cohort_but_not_its_novel_candidates_says_so_after_the_run(
    world,
):
    """gnomAD holds half of the cohort's variants (the inherited ones), so the filter
    binds, and none of the de novo candidates is in it: what a de novo list looks like.
    That is a fact about the candidates, said after the run, and not a refusal."""
    doc = jtrio(world, "c", "--preset", "strict-denovo", ann="novel")
    eff = doc["effective"]
    assert eff["gnomad"]["bound"] is True
    assert eff["gnomad"]["cohort_variants_sampled"] == 30
    assert eff["gnomad"]["cohort_variants_with_record"] == 15
    assert eff["gnomad"]["candidates_checked"] == 15
    assert eff["gnomad"]["candidates_without_record"] == 15
    (late,) = [w for w in eff["warnings"] if w["code"] == "gnomad_no_record"]
    assert "judged none of the 15 candidates" in late["message"]
    assert "a record for 50% of a sample of the cohort's variants" in late["message"]
    assert "probably novel" in late["message"]
    assert eff["population_filter"] is False  # INFO/AF is the cohort's own
    out = trio(world, "c", "--preset", "strict-denovo", ann="novel").stdout
    assert out.splitlines()[-1].startswith(
        "  WARNING: gnomad_max_af=0.0001 judged none of the 15"
    )
    # a store that does hold some candidates is not that warning
    ok = jtrio(world, "p", "--preset", "strict-denovo")
    assert "gnomad_no_record" not in [w["code"] for w in ok["effective"]["warnings"]]
    assert ok["effective"]["gnomad"]["candidates_checked"] == 8 + 1  # the common one
    assert ok["effective"]["gnomad"]["cohort_variants_sampled"] == 51
    assert ok["effective"]["gnomad"]["cohort_variants_with_record"] == 32


def test_preset_ignore_switches_one_key_off_for_one_run(world):
    doc = jtrio(
        world, "p", "--preset", "strict-denovo", "--preset-ignore", "gnomad_max_af"
    )
    assert names(doc, "denovo") == set(KEPT["strict-denovo"]["denovo"]) | {
        "dn_gnomad_common"
    }
    assert doc["gates"]["gnomad_max_af"] is None
    assert {"key": "gnomad_max_af", "preset": 0.0001, "flag": "off"} in doc["preset"][
        "overridden"
    ]


# ───────────────────── cohort-derived AF and a preset's max_af ─────────────────────


def test_a_presets_max_af_does_not_defeat_the_cohort_af_skip(world):
    doc = jtrio(world, "c", "--preset", "broad-first-pass", ann="empty")
    assert doc["gates"]["af_gate"] == "skipped_cohort_derived"
    assert doc["gates"]["max_af"] == 0.01  # the preset's value, shown, not applied
    assert doc["preset"]["resolved"]["max_af"]["from"] == "preset"
    assert doc["models"]["denovo"]["count"] == 15  # every private het survives
    codes = [w["code"] for w in doc["effective"]["warnings"]]
    assert "no_population_frequency_filter" in codes
    assert "preset_max_af_not_applied" in codes
    assert any("was not applied" in n for n in doc["preset"]["notes"])
    assert doc["effective"]["info_af"] == {
        "state": "skipped_cohort_derived",
        "bound": False,
    }
    text = trio(world, "c", "--preset", "broad-first-pass", ann="empty").stdout
    assert "NO population frequency filter applied" in text
    assert "max_af=0.01 [preset; skipped: INFO/AF is cohort-derived]" in text


def test_typing_max_af_forces_the_gate_even_with_a_preset(world):
    doc = jtrio(
        world, "c", "--preset", "broad-first-pass", "--max-af", "0.01", ann="empty"
    )
    assert doc["gates"]["af_gate"] == "cohort_derived_suspected"
    assert doc["models"]["denovo"]["count"] == 0  # AF 1/6 > 0.01: removed, as asked
    assert doc["preset"]["resolved"]["max_af"]["from"] == "flag"


def test_strict_denovo_without_a_population_filter_refuses_and_can_be_overridden(world):
    r = trio(world, "c", "--preset", "strict-denovo", ann="empty", ok=False)
    assert "NO population frequency filter applied" in r.stderr
    assert "annotations load-gnomad" in r.stderr
    doc = jtrio(
        world, "c", "--preset", "strict-denovo", "--no-require-bound-gates", ann="empty"
    )
    assert doc["models"]["denovo"]["count"] == 15
    assert doc["effective"]["population_filter"] is False


# ───────────────────────────── --category lists ─────────────────────────────


def test_a_comma_list_runs_those_models_and_the_text_has_a_line_each(world):
    default = jtrio(world, "p")
    listed = jtrio(world, "p", "--category", "recessive,comphet")
    assert list(listed["models"]) == ["recessive", "comphet"]
    assert listed["models"] == {
        m: default["models"][m] for m in ("recessive", "comphet")
    }
    text = trio(world, "p", "--category", "recessive,comphet").stdout.splitlines()
    assert [ln.split()[0] for ln in text[1:]] == ["recessive", "comphet"]


def test_a_bad_category_has_todays_wording():
    import click

    from cli.trio_analysis import CATEGORY_CHOICES

    with pytest.raises(click.BadParameter) as plain:
        click.Choice(CATEGORY_CHOICES).convert("nope", None, None)
    for value in ("nope", "recessive,nope"):
        r = CliRunner().invoke(cli, ["db", "trio", "t", "--category", value])
        assert r.exit_code == 2
        assert "Invalid value for '--category'" in r.output
        assert plain.value.message in r.output
