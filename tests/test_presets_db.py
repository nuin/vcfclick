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

import dataclasses
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
from tests.fixtures import make_preset_vcf
from tests.fixtures._vcfgen import (
    TRIO_PED,
    TRIO_SAMPLES,
    Site,
    trio_sites,
    write_vcf,
)
from tests.test_trio_contract import assert_trio_document

BACKENDS = ("chdb", "duckdb")
TRUTH = json.loads((FIXTURES / "presets_truth.json").read_text())
POS: dict[str, int] = TRUTH["positions"]
NAME_AT = {v: k for k, v in POS.items()}
DE_NOVO = {name for name in POS if name.startswith("dn_")}
KEPT = TRUTH["kept"]
ENGINE_NOISE = re.compile(r"^Default dirty pages decay period: \d+ms\n", re.M)


def mismatch_sites() -> list[Site]:
    """10 de novo hets whose parents' 0/0 rows have GQ and DP but no AD, then 200
    hom-alt sites where both parents are het, with AD. INFO/AF is a population
    frequency, so no other gate gets in the way."""
    de_novo = [
        Site("chr1", 1000 * (i + 1), "A", "G",
             ("0/1:15,15:99:30", "0/0:.:99:30", "0/0:.:99:30"),
             info="AF=0.0005", fmt="GT:AD:GQ:DP")
        for i in range(10)
    ]  # fmt: skip
    inherited = [
        Site("chr1", 20000 + 100 * i, "A", "G",
             ("1/1:0,30:99:30", "0/1:15,15:99:30", "0/1:15,15:99:30"),
             info="AF=0.002", fmt="GT:AD:GQ:DP")
        for i in range(200)
    ]  # fmt: skip
    return de_novo + inherited


def gq_shape_sites() -> list[Site]:
    """10 de novo hets whose proband GQ is '.', then 200 sites where the proband is
    hom-alt and both parents het, all with GQ. GQ is on 200 of the proband's 210
    non-reference rows (95.2%) and on none of the rows the de novo model reads."""
    de_novo = [
        Site("chr1", 1000 * (i + 1), "A", "G",
             ("0/1:15,15:.:30", "0/0:30,0:99:30", "0/0:30,0:99:30"),
             info="AF=0.0005", fmt="GT:AD:GQ:DP")
        for i in range(10)
    ]  # fmt: skip
    inherited = [
        Site("chr1", 20000 + 100 * i, "A", "G",
             ("1/1:0,30:99:30", "0/1:15,15:99:30", "0/1:15,15:99:30"),
             info="AF=0.002", fmt="GT:AD:GQ:DP")
        for i in range(200)
    ]  # fmt: skip
    return de_novo + inherited


def ad_shape_sites() -> list[Site]:
    """10 de novo hets whose parents' 0/0 rows carry no AD, then 400 hets inherited
    from one parent (200 from each), where the other parent's 0/0 row has AD. Each
    parent has 210 0/0 rows, 200 with AD (95.2%), and none at the de novo sites."""
    de_novo = [
        Site("chr1", 1000 * (i + 1), "A", "G",
             ("0/1:15,15:99:30", "0/0:.:99:30", "0/0:.:99:30"),
             info="AF=0.0005", fmt="GT:AD:GQ:DP")
        for i in range(10)
    ]  # fmt: skip
    from_father = [
        Site("chr1", 20000 + 100 * i, "A", "G",
             ("0/1:15,15:99:30", "0/1:15,15:99:30", "0/0:30,0:99:30"),
             info="AF=0.002", fmt="GT:AD:GQ:DP")
        for i in range(200)
    ]  # fmt: skip
    from_mother = [
        Site("chr1", 50000 + 100 * i, "A", "G",
             ("0/1:15,15:99:30", "0/0:30,0:99:30", "0/1:15,15:99:30"),
             info="AF=0.002", fmt="GT:AD:GQ:DP")
        for i in range(200)
    ]  # fmt: skip
    return de_novo + from_father + from_mother


def af_gap_sites() -> list[Site]:
    """Three de novo hets (INFO/AF 0.0005, 0.5 and none) and 100 inherited hets, 5 of
    them with no INFO/AF: the AF gate binds (97 of 103 variants carry it) but not on
    every variant, and a variant without INFO/AF passes it."""
    denovo = [
        Site("chr1", pos, "A", "G",
             ("0/1:15,15:99:30", "0/0:30,0:99:30", "0/0:30,0:99:30"),
             info=info, fmt="GT:AD:GQ:DP")
        for pos, info in ((1000, "AF=0.0005"), (2000, "AF=0.5"), (3000, ".")) 
    ]  # fmt: skip
    inherited = [
        Site("chr1", 20000 + 100 * i, "A", "G",
             ("0/1:15,15:99:30", "0/1:15,15:99:30", "0/0:30,0:99:30"),
             info="." if i < 5 else "AF=0.002", fmt="GT:AD:GQ:DP")
        for i in range(100)
    ]  # fmt: skip
    return denovo + inherited


def reference_depth_sites() -> list[Site]:
    """Ten de novo hets at DP 24, then 100 sites where the proband is 0/0 at DP 4
    (and the father het): the proband's stored rows are mostly shallow reference
    calls, and its calls are at DP 24."""
    denovo = [
        Site("chr1", 1000 * (i + 1), "A", "G",
             ("0/1:12,12:99:24", "0/0:30,0:99:30", "0/0:30,0:99:30"),
             info="AF=0.0005", fmt="GT:AD:GQ:DP")
        for i in range(10)
    ]  # fmt: skip
    shallow = [
        Site("chr1", 20000 + 100 * i, "A", "G",
             ("0/0:4,0:99:4", "0/1:15,15:99:30", "0/0:30,0:99:30"),
             info="AF=0.002", fmt="GT:AD:GQ:DP")
        for i in range(100)
    ]  # fmt: skip
    return denovo + shallow


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
    # the planted sites again, with the FILTER column '.' everywhere (the caller
    # wrote no labels): the LowQual site is no longer marked
    unlabelled = write_vcf(
        base / "unlabelled.vcf",
        make_preset_vcf.SAMPLES,
        [dataclasses.replace(x, filter=".") for x in make_preset_vcf.sites()],
        filter_header=make_preset_vcf.FILTERS,
    )
    # a caller whose reference calls carry no allele depths while its het calls do:
    # over ALL the parents' rows AD is on 95.2%, and on the 0/0 rows a de novo call
    # rests on it is on none
    mismatch = write_vcf(base / "mismatch.vcf", TRIO_SAMPLES, mismatch_sites())
    # the same shape, judged on the rows the candidates are made of: the field is on
    # 95.2% of the class the gate reads and on none of the rows at the candidates
    gq_shape = write_vcf(base / "gq_shape.vcf", TRIO_SAMPLES, gq_shape_sites())
    ad_shape = write_vcf(base / "ad_shape.vcf", TRIO_SAMPLES, ad_shape_sites())
    af_gap = write_vcf(base / "af_gap.vcf", TRIO_SAMPLES, af_gap_sites())
    ref_depth = write_vcf(base / "ref_depth.vcf", TRIO_SAMPLES, reference_depth_sites())
    stores = (
        ("p", FIXTURES / "presets.vcf.gz", ped),
        ("n", FIXTURES / "presets_nogq.vcf.gz", ped),
        ("c", cohort, cohort_ped),
        ("f", unlabelled, ped),
        ("m", mismatch, cohort_ped),
        ("g", gq_shape, cohort_ped),
        ("a", ad_shape, cohort_ped),
        ("u", af_gap, cohort_ped),
        ("d", ref_depth, cohort_ped),
    )
    for name, vcf, pedigree in stores:
        run_cli(home, backend, "db", "create", name)
        run_cli(
            home, backend, "db", "ingest", name, str(vcf),
            "--cohort", "fam", "--ingest-id", "fam1", "--serial", "--keep-reference",
        )  # fmt: skip
        run_cli(home, backend, "db", "ped", name, str(pedigree))
    # the planted sites ingested WITHOUT --keep-reference: no stored 0/0 calls, so a
    # parent cannot be shown to be 0/0 and the models that need one cannot run
    run_cli(home, backend, "db", "create", "r")
    run_cli(
        home, backend, "db", "ingest", "r", str(FIXTURES / "presets.vcf.gz"),
        "--cohort", "fam", "--ingest-id", "fam1", "--serial",
    )  # fmt: skip
    run_cli(home, backend, "db", "ped", "r", str(ped))
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
    assert lines[1].startswith("preset: strict-denovo (builtin, rev 2, ")
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
    assert lines[-2] == (
        "  gnomAD: no record for 6 of 8 candidates (kept as rare); it has a "
        "record for 62.7% of a sample of the cohort's variants"
    )
    # and last, what strict_missing dropped: dn_null_dp_parent, whose parent's depth
    # is missing, passes every other gate
    assert lines[-1].startswith(
        "  WARNING: strict_missing dropped 1 denovo candidates that pass every other "
        "gate"
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
            # told about, as the parent-allele screen is, but only about that gate
            assert "preset" not in doc
            assert doc["effective"]["warnings"] == []  # the LowQual site is a label
            assert doc["effective"]["filter_column"]["non_pass_variants"] == 1
            assert doc["effective"]["gq"]["bound"] is None  # not asked: not examined


def test_a_parent_alt_screen_that_screens_nothing_says_so_and_says_only_that(world):
    """The store has no allele depths: every parent 0/0 would be trusted. A typed
    opt-in gate that judges nothing must not look like a screen that worked."""
    # a gnomAD that judges nothing here (one row at chr1:5) is not probed by this
    # run, so the report must not call it bound nor a population filter in effect
    doc = jtrio(
        world, "n", "--category", "denovo", "--parent-alt-screen",
        "--gnomad-max-af", "0.0001", ann="far",
    )  # fmt: skip
    assert names(doc, "denovo") == DE_NOVO  # nothing was screened out
    eff = doc["effective"]
    (w,) = eff["warnings"]  # not the GQ / DP gates it did not ask for
    assert (w["code"], w["key"]) == ("gate_unbound", "parent_alt_screen")
    assert "every parent 0/0 would be trusted" in w["message"]
    assert (eff["gnomad"]["loaded"], eff["gnomad"]["bound"]) == (True, None)
    assert eff["gnomad"]["cohort_variants_sampled"] is None  # never looked at
    assert (eff["population_filter"], eff["ad"]["bound"]) == (None, False)
    assert eff["gq"]["bound"] is None and eff["info_af"]["bound"] is None
    text = trio(world, "n", "--category", "denovo", "--parent-alt-screen").stdout
    assert [x for x in text.splitlines() if x.startswith("  WARNING:")] == [
        f"  WARNING: {w['message']}"
    ]  # on stdout, and it did not refuse: nothing insisted on the gate


# ─────────────────────────── gates that cannot bind ───────────────────────────


def warning_lines(stdout: str) -> list[str]:
    return [x for x in stdout.splitlines() if x.startswith("  WARNING:")]


def test_pass_only_on_a_vcf_without_filter_labels_says_it_removes_nothing(world):
    """FILTER '.' everywhere: the planted LowQual site is not marked, so the gate
    cannot remove it. strict-denovo lists it as a strict candidate, and must say
    that its FILTER gate did nothing (a warning: an all-PASS VCF is legitimate)."""
    doc = jtrio(world, "f", "--preset", "strict-denovo")
    assert "dn_nonpass" in names(doc, "denovo")
    eff = doc["effective"]
    assert eff["filter_column"] == {"present": True, "non_pass_variants": 0}
    (w,) = [x for x in eff["warnings"] if x["code"] == "pass_only_no_labels"]
    assert (w["level"], w["key"]) == ("warning", "pass_only")
    assert doc["preset"]["resolved"]["pass_only"]["bound"] is False
    text = trio(world, "f", "--preset", "strict-denovo").stdout
    assert f"  WARNING: {w['message']}" in warning_lines(text)
    # the same preset on a VCF that does label a site: bound, silent, and it works
    labelled = jtrio(world, "p", "--preset", "strict-denovo")
    assert labelled["preset"]["resolved"]["pass_only"]["bound"] is True
    assert "dn_nonpass" not in names(labelled, "denovo")
    assert "pass_only_no_labels" not in [
        x["code"] for x in labelled["effective"]["warnings"]
    ]


def test_a_typed_pass_only_is_told_about_and_only_about_that_gate(world):
    doc = jtrio(world, "f", "--category", "denovo", "--pass-only")
    assert names(doc, "denovo") == DE_NOVO  # nothing was removed
    assert "preset" not in doc
    eff = doc["effective"]
    (w,) = eff["warnings"]  # not the GQ / DP gates it did not ask for
    assert (w["code"], w["key"]) == ("pass_only_no_labels", "pass_only")
    assert eff["population_filter"] is None and eff["gnomad"]["bound"] is None
    out = trio(world, "f", "--category", "denovo", "--pass-only").stdout
    assert warning_lines(out) == [f"  WARNING: {w['message']}"]
    # without the flag the run is the plain one: no effective object, no line
    plain = jtrio(world, "f", "--category", "denovo")
    assert "effective" not in plain
    assert warning_lines(trio(world, "f", "--category", "denovo").stdout) == []


def test_a_gate_is_judged_on_the_rows_it_reads_so_a_strict_preset_refuses(world):
    """The parents' 0/0 rows (10) carry no allele depths and their het rows (200)
    do: AD is on 95.2% of the parents' rows and on none of the rows the parent-allele
    screen reads, the 0/0 rows at the de novo sites. strict_missing then fails every
    de novo candidate, so a strict preset used to run, say every gate bound, and
    print `denovo candidates: 0`."""
    ignore = ("--preset-ignore", "gnomad_max_af")
    r = trio(world, "m", "--preset", "strict-denovo", *ignore, ok=False)
    assert r.stdout == ""
    for text in ("parent_alt_screen cannot be enforced",
                 "only 0% of father/mother 0/0 rows at denovo candidate sites carry AD",
                 "no de novo candidate can pass", "Remedy: AD is missing",
                 "--no-require-bound-gates"):  # fmt: skip
        assert text in r.stderr, text
    # the gates that do bind are not blamed
    assert "min_gq" not in r.stderr and "min_ab" not in r.stderr
    # run anyway: the gate that cannot bind is exempt, and the others decide
    doc = jtrio(
        world, "m", "--preset", "strict-denovo", *ignore, "--no-require-bound-gates"
    )
    assert doc["models"]["denovo"]["count"] == 10
    eff = doc["effective"]
    assert eff["ad"]["father"] == eff["ad"]["mother"] == 0.9524  # the misleading share
    assert eff["ad"]["bound"] is False  # ... which no longer decides
    assert (eff["gq"]["bound"], eff["dp"]["bound"]) == (True, True)
    assert doc["gates"]["strict_exempt"] == ["parent_alt_screen"]
    assert doc["preset"]["resolved"]["parent_alt_screen"]["bound"] is False
    assert doc["preset"]["resolved"]["min_ab"]["bound"] is True


def test_a_typed_parent_alt_screen_that_screens_nothing_says_so_whatever_the_other_rows(
    world,
):
    doc = jtrio(world, "m", "--category", "denovo", "--parent-alt-screen")
    assert doc["models"]["denovo"]["count"] == 10  # nothing was screened out
    (w,) = doc["effective"]["warnings"]
    assert (w["code"], w["key"]) == ("gate_unbound", "parent_alt_screen")
    assert (
        "only 0% of father/mother 0/0 rows at denovo candidate sites carry AD"
        in w["message"]
    )
    assert "every parent 0/0 would be trusted" in w["message"]
    out = trio(world, "m", "--category", "denovo", "--parent-alt-screen").stdout
    assert warning_lines(out) == [f"  WARNING: {w['message']}"]
    # a model that never reads the parents' 0/0 rows is not blamed for them: the
    # screen is simply not applicable there, which is said, and is not unbound
    rec = jtrio(world, "m", "--category", "recessive", "--parent-alt-screen")
    (n,) = rec["effective"]["warnings"]
    assert (n["code"], n["key"]) == ("gate_not_applicable", "parent_alt_screen")
    assert "the screen is read by the denovo model only" in n["message"]
    assert rec["effective"]["ad"]["bound"] is False


def test_a_gate_no_selected_model_reads_is_not_bound_and_says_so(world):
    """`--category recessive --min-ab 0.4`: the proband is hom-alt, which the band
    does not judge, so the gate removes nothing. It used to read as bound."""
    args = ("--preset", "rare-recessive", "--category", "recessive", "--min-ab", "0.4")
    doc = jtrio(world, "p", *args)
    assert names(doc, "recessive") == {"rec_ok"}
    assert doc["preset"]["resolved"]["min_ab"] == {
        "value": 0.4, "from": "flag", "bound": False,
    }  # fmt: skip
    assert doc["preset"]["resolved"]["min_gq"]["bound"] is True
    (w,) = [
        x for x in doc["effective"]["warnings"] if x["code"] == "gate_not_applicable"
    ]
    assert (
        w["key"] == "min_ab"
        and "min_ab=0.4/max_ab=0.75 is not applied by" in w["message"]
    )
    assert "(recessive)" in w["message"]
    out = trio(world, "p", *args)  # a warning, never a refusal
    assert f"  WARNING: {w['message']}" in warning_lines(out.stdout)
    # the same band on a run that includes a model that reads het calls binds
    both = jtrio(
        world, "p", *args[:2], "--category", "recessive,dominant", "--min-ab", "0.4"
    )
    assert both["preset"]["resolved"]["min_ab"]["bound"] is True


DN_POSITIONS = {1000 * (i + 1) for i in range(10)}


def positions(doc: dict, model: str = "denovo") -> set[int]:
    return {c["pos"] for c in doc["models"][model]["candidates"]}


IGNORE_GNOMAD = ("--preset-ignore", "gnomad_max_af")


def test_gq_on_the_inherited_calls_and_not_on_the_de_novo_ones_refuses_strict_denovo(
    world,
):
    """The proband's non-reference rows carry GQ on 95.2% of them (10 de novo hets
    without, 200 inherited hom-alts with), and none of the rows a de novo call is made
    on does. strict_missing then fails every candidate: it used to run, say every gate
    bound, and print `denovo candidates: 0`. It refuses, with a remedy that names GQ."""
    r = trio(world, "g", "--preset", "strict-denovo", *IGNORE_GNOMAD, ok=False)
    assert r.stdout == ""
    for text in ("min_gq=30 cannot be enforced",
                 "only 0% of proband 0/1 or 1/1 rows at denovo candidate sites carry GQ",
                 "strict_missing counts a call without GQ as a fail",
                 "Remedy: GQ is missing where the gates read it",
                 "--no-require-bound-gates"):  # fmt: skip
        assert text in r.stderr, text
    assert "min_dp" not in r.stderr and "min_ab" not in r.stderr  # the others bind
    # run anyway: the gate is exempt, the 10 are listed, and the document says why
    doc = jtrio(
        world, "g", "--preset", "strict-denovo", *IGNORE_GNOMAD,
        "--no-require-bound-gates",
    )  # fmt: skip
    assert positions(doc) == DN_POSITIONS and doc["models"]["denovo"]["count"] == 10
    assert doc["gates"]["strict_exempt"] == ["min_gq"]
    eff = doc["effective"]
    assert (
        eff["gq"]["proband"] == 0.9524
    )  # the share that misled, and no longer decides
    assert (eff["gq"]["bound"], eff["dp"]["bound"], eff["ad"]["bound"]) == (
        False, True, True,
    )  # fmt: skip
    (w,) = [x for x in eff["warnings"] if x["code"] == "gate_unbound"]
    assert w["key"] == "min_gq"
    assert (
        "only 0% of proband 0/1 or 1/1 rows at denovo candidate sites" in w["message"]
    )
    assert "strict_missing is not applied to it for this run" in w["message"]
    out = trio(
        world, "g", "--preset", "strict-denovo", *IGNORE_GNOMAD,
        "--no-require-bound-gates",
    ).stdout  # fmt: skip
    assert f"  WARNING: {w['message']}" in warning_lines(out)
    assert "denovo candidates: 10" in out


def test_gq_missing_on_the_candidates_is_a_warning_for_a_warn_only_preset(world):
    out = trio(world, "g", "--preset", "broad-first-pass", *IGNORE_GNOMAD).stdout
    assert (
        "  WARNING: min_gq=15 cannot be enforced: only 0% of proband 0/1 or 1/1 rows "
        "at denovo candidate sites carry GQ"
    ) in warning_lines(out)
    doc = jtrio(world, "g", "--preset", "broad-first-pass", *IGNORE_GNOMAD)
    assert doc["effective"]["gq"]["bound"] is False
    assert doc["preset"]["resolved"]["min_gq"]["bound"] is False
    assert doc["preset"]["resolved"]["min_dp"]["bound"] is True
    assert positions(doc) == DN_POSITIONS  # tolerant: it ran, and removed nothing
    # a model whose candidates carry GQ is not blamed: the recessive ones do
    rec = jtrio(
        world, "g", "--preset", "rare-recessive", "--category", "recessive",
        *IGNORE_GNOMAD,
    )  # fmt: skip
    assert rec["effective"]["gq"]["bound"] is True
    assert not [w for w in rec["effective"]["warnings"] if w["code"] == "gate_unbound"]


def test_ad_on_the_inherited_reference_calls_and_not_on_the_de_novo_ones_refuses(world):
    """Each parent's 0/0 rows carry AD on 95.2% of them and on none of the 10 at the
    de novo sites, which is where the parent-allele screen reads them."""
    r = trio(world, "a", "--preset", "strict-denovo", *IGNORE_GNOMAD, ok=False)
    assert r.stdout == ""
    for text in ("parent_alt_screen cannot be enforced",
                 "only 0% of father/mother 0/0 rows at denovo candidate sites carry AD",
                 "no de novo candidate can pass",
                 "Remedy: AD is missing where the gates read it"):  # fmt: skip
        assert text in r.stderr, text
    assert "min_gq" not in r.stderr and "min_ab" not in r.stderr
    doc = jtrio(
        world, "a", "--preset", "strict-denovo", *IGNORE_GNOMAD,
        "--no-require-bound-gates",
    )  # fmt: skip
    assert positions(doc) == DN_POSITIONS
    assert doc["gates"]["strict_exempt"] == ["parent_alt_screen"]
    eff = doc["effective"]
    assert eff["ad"]["father"] == eff["ad"]["mother"] == round(400 / 410, 4)
    assert eff["ad"]["bound"] is False and eff["gq"]["bound"] is True
    # a warn-only run with the screen on: the warning is printed, and nothing refuses
    out = trio(
        world, "a", "--preset", "broad-first-pass", *IGNORE_GNOMAD,
        "--parent-alt-screen",
    ).stdout  # fmt: skip
    assert any(
        "parent_alt_screen cannot be enforced: only 0% of father/mother 0/0 rows "
        "at denovo candidate sites carry AD" in x
        for x in warning_lines(out)
    )
    typed = jtrio(world, "a", "--category", "denovo", "--parent-alt-screen")
    (w,) = typed["effective"]["warnings"]
    assert (w["code"], w["key"]) == ("gate_unbound", "parent_alt_screen")


NO_REFERENCE = (
    "denovo cannot run: this database has no stored hom-reference calls; "
    "re-ingest with --keep-reference"
)


def test_a_strict_preset_refuses_a_database_its_model_cannot_run_on(world):
    """No stored 0/0 calls: strict-denovo would print `denovo candidates: 0`, exit 0
    and look like a clean result. It refuses, with the remedy."""
    r = trio(world, "r", "--preset", "strict-denovo", ok=False)
    assert r.stdout == ""
    for text in ("preset 'strict-denovo' needs every gate to bind", NO_REFERENCE,
                 "Remedy:", "--keep-reference", "--no-require-bound-gates"):  # fmt: skip
        assert text in r.stderr, text
    # --no-require-bound-gates runs it, and the run says in stdout and in the document
    out = trio(world, "r", "--preset", "strict-denovo", "--no-require-bound-gates")
    assert f"  WARNING: {NO_REFERENCE}" in warning_lines(out.stdout)
    assert "denovo candidates: 0" in out.stdout
    doc = jtrio(world, "r", "--preset", "strict-denovo", "--no-require-bound-gates")
    (w,) = [x for x in doc["effective"]["warnings"] if x["code"] == "model_blocked"]
    assert (w["level"], w["key"], w["message"]) == ("warning", "category", NO_REFERENCE)
    assert doc["models"]["denovo"]["blocked"] is True
    assert doc["models"]["denovo"]["blocked_reason"] == "needs_keep_reference"
    # the same preset on a database that does store them is not touched
    ok = jtrio(world, "p", "--preset", "strict-denovo")
    assert not [x for x in ok["effective"]["warnings"] if x["code"] == "model_blocked"]


def test_a_typed_parent_alt_screen_on_a_database_without_reference_calls_says_so(world):
    """The screen is about the de novo model, which has no stored 0/0 to screen on
    this database: it would print `denovo candidates: 0` and look like a screen that
    found nothing."""
    out = trio(world, "r", "--category", "denovo", "--parent-alt-screen")
    assert warning_lines(out.stdout) == [f"  WARNING: {NO_REFERENCE}"]
    assert "denovo candidates: 0" in out.stdout
    doc = jtrio(world, "r", "--category", "denovo", "--parent-alt-screen")
    (w,) = doc["effective"]["warnings"]
    assert (w["code"], w["key"], w["message"]) == (
        "model_blocked",
        "category",
        NO_REFERENCE,
    )
    assert "preset" not in doc and doc["models"]["denovo"]["blocked"] is True
    # a database that stores them, or a model the screen is not about: nothing blocked
    ok = jtrio(world, "p", "--category", "denovo", "--parent-alt-screen")
    assert ok["effective"]["warnings"] == []
    rec = jtrio(world, "r", "--category", "recessive", "--parent-alt-screen")
    assert [w["code"] for w in rec["effective"]["warnings"]] == ["gate_not_applicable"]
    # --pass-only typed alone is not about the de novo model
    flt = jtrio(world, "r", "--category", "denovo", "--pass-only")
    assert not [w for w in flt["effective"]["warnings"] if w["code"] == "model_blocked"]


def test_a_warn_only_preset_names_each_selected_model_that_cannot_run(world):
    """rare-recessive selects recessive (which needs no reference) and comphet
    (which needs the reference and the gene table). Nothing is refused; the
    models that did not run are named, and the one that did is not."""
    r = trio(world, "r", "--preset", "rare-recessive")
    assert [x for x in warning_lines(r.stdout) if "cannot run" in x] == [
        "  WARNING: comphet cannot run: this database has no stored hom-reference "
        "calls; re-ingest with --keep-reference"
    ]
    assert "denovo cannot run" not in r.stdout  # not selected: not mentioned
    doc = jtrio(world, "r", "--preset", "rare-recessive")
    assert names(doc, "recessive") == {"rec_ok"}
    # no gene table: comphet cannot run on a database that does store the 0/0 calls
    genes = trio(world, "p", "--preset", "rare-recessive", ann="empty")
    assert any(
        x.startswith("  WARNING: comphet cannot run: no gene annotations")
        for x in warning_lines(genes.stdout)
    )
    gdoc = jtrio(world, "p", "--preset", "rare-recessive", ann="empty")
    assert gdoc["models"]["comphet"]["blocked_reason"] == "needs_genes"
    assert [
        x["message"]
        for x in gdoc["effective"]["warnings"]
        if x["code"] == "model_blocked"
    ] == [
        "comphet cannot run: no gene annotations are loaded "
        "(`vcfclick annotations load`)"
    ]
    # a run that insists on its gates refuses for it too
    refused = trio(
        world, "p", "--preset", "rare-recessive", "--require-bound-gates", ann="empty",
        ok=False,
    )  # fmt: skip
    assert "comphet cannot run: no gene annotations are loaded" in refused.stderr


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


def test_candidates_dropped_only_for_a_missing_field_are_counted_after_the_run(world):
    """DP binds (it is on 94% of the rows at the candidate sites), so strict-denovo
    runs, and the one candidate whose parent has no depth is dropped for that alone:
    dn_null_dp_parent passes every other gate. Said on stdout and in the document."""
    ignore = (*IGNORE_GNOMAD, "--max-af", "0.001")
    doc = jtrio(world, "p", "--preset", "strict-denovo", *ignore)
    assert "dn_null_dp_parent" not in names(doc, "denovo")
    (w,) = [
        x for x in doc["effective"]["warnings"] if x["code"] == "strict_missing_dropped"
    ]
    assert (w["level"], w["key"]) == ("warning", "strict_missing")
    assert w["message"].startswith(
        "strict_missing dropped 1 denovo candidates that pass every other gate, only "
        "because a GQ, DP or allele depth"
    )
    out = trio(world, "p", "--preset", "strict-denovo", *ignore).stdout
    assert warning_lines(out)[-1] == f"  WARNING: {w['message']}"
    # the same candidates without strict_missing: the one it dropped is back, and
    # nothing is said
    loose = jtrio(
        world, "p", "--preset", "strict-denovo", *ignore, "--no-strict-missing"
    )
    assert names(loose, "denovo") == names(doc, "denovo") | {"dn_null_dp_parent"}
    assert "strict_missing_dropped" not in [
        x["code"] for x in loose["effective"]["warnings"]
    ]
    # a run that asks for no strictness is not told
    broad = jtrio(world, "p", "--preset", "broad-first-pass")
    assert "strict_missing_dropped" not in [
        x["code"] for x in broad["effective"]["warnings"]
    ]
    # --strict-missing typed on its own is the same check, in the same two places
    alone = jtrio(world, "p", "--category", "denovo", "--strict-missing")
    assert [x["code"] for x in alone["effective"]["warnings"]].count(
        "strict_missing_dropped"
    ) == 1
    text = trio(world, "p", "--category", "denovo", "--strict-missing").stdout
    assert any("strict_missing dropped" in x for x in warning_lines(text))


def test_an_info_af_gate_that_binds_on_some_variants_says_how_many(world):
    """INFO/AF is on 97 of 103 variants, so the gate binds, and on the other six a
    variant with no INFO/AF passes it. It used to bind in silence."""
    args = ("--preset", "broad-first-pass", "--category", "denovo", *IGNORE_GNOMAD)
    doc = jtrio(world, "u", *args)
    (w,) = [x for x in doc["effective"]["warnings"] if x["key"] == "max_af"]
    assert (w["code"], w["level"]) == ("gate_partial", "info")
    assert w["message"] == (
        "max_af=0.01 enforced on 94.1% of variants (a variant with no INFO/AF passes)"
    )
    assert doc["effective"]["info_af"]["bound"] is True
    # the gate does what the message says: 0.5 is removed, no INFO/AF passes
    assert positions(doc) == {1000, 3000}
    r = trio(world, "u", *args)
    assert f"note: {w['message']}" in r.stderr  # an info note, not a WARNING line
    assert not [x for x in warning_lines(r.stdout) if "max_af" in x]
    # a store that carries INFO/AF on every variant says nothing of the kind
    full = jtrio(world, "p", *args)
    assert not [x for x in full["effective"]["warnings"] if x["key"] == "max_af"]


def test_strict_missing_that_drops_nothing_says_nothing(world):
    """No candidate lacks a field that a gate needs: nothing was lost to it."""
    doc = jtrio(world, "n", "--preset", "strict-denovo", "--no-require-bound-gates")
    assert "strict_missing_dropped" not in [
        x["code"] for x in doc["effective"]["warnings"]
    ]  # the fields are missing, the gates are exempt, so nothing is dropped
    ok = jtrio(world, "g", "--preset", "strict-denovo", *IGNORE_GNOMAD,
               "--no-require-bound-gates")  # fmt: skip
    assert (
        doc["models"]["denovo"]["count"] > 0 and ok["models"]["denovo"]["count"] == 10
    )
    assert "strict_missing_dropped" not in [
        x["code"] for x in ok["effective"]["warnings"]
    ]


def test_a_median_depth_near_min_dp_is_a_warning(world):
    doc = jtrio(world, "p", "--preset", "broad-first-pass", "--min-dp", "30")
    (w,) = [x for x in doc["effective"]["warnings"] if x["code"] == "thin_coverage"]
    assert "median proband DP is 40 but min_dp is 30" in w["message"]


def test_thin_coverage_is_judged_on_the_depth_of_the_proband_s_calls(world):
    """Its stored rows are 100 hom-ref calls at DP 4 and 10 calls at DP 24. The median
    over all of them is 4, which made every assay look thin; the calls are at 24."""
    base = ("--preset", "broad-first-pass", *IGNORE_GNOMAD)
    deep = jtrio(world, "d", *base, "--min-dp", "10")  # 1.5 x 10 = 15 < 24
    assert deep["effective"]["median_dp"]["proband"] == 4  # the data, all rows
    assert not [
        x for x in deep["effective"]["warnings"] if x["code"] == "thin_coverage"
    ]
    thin = jtrio(world, "d", *base, "--min-dp", "20")  # 1.5 x 20 = 30 > 24
    (w,) = [x for x in thin["effective"]["warnings"] if x["code"] == "thin_coverage"]
    assert w["message"].startswith(
        "median proband DP is 24 but min_dp is 20 (the median"
    )
    out = trio(world, "d", *base, "--min-dp", "20").stdout
    assert f"  WARNING: {w['message']}" in warning_lines(out)
    assert not [x for x in warning_lines(trio(world, "d", *base, "--min-dp", "10").stdout)
                if "thin" in x]  # fmt: skip


def test_gnomad_asked_for_and_not_loaded(world):
    doc = jtrio(world, "p", "--preset", "broad-first-pass", ann="empty")
    codes = [w["code"] for w in doc["effective"]["warnings"]]
    assert "gnomad_not_loaded" in codes
    assert doc["effective"]["gnomad"]["loaded"] is False
    # the filter ran over an empty store: no count of candidates it found no
    # record for, because it could not have found one
    gn = doc["effective"]["gnomad"]
    assert (gn["bound"], gn["candidates_checked"], gn["candidates_without_record"]) == (
        False, None, None,
    )  # fmt: skip
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
    # no candidate has a record either, and the probe has said gnomAD judges
    # nothing here: the candidates are not "probably novel", the build is the suspect
    (late,) = [w for w in eff["warnings"] if w["code"] == "gnomad_no_record"]
    assert "probably novel" not in late["message"]
    assert "Check that its genome build and chromosome naming" in late["message"]
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


def test_a_gnomad_that_was_never_asked_for_is_neither_bound_nor_probed(world):
    """The preset's gnomad_max_af is switched off, so the loaded store (here one row
    that judges nothing in this cohort) is not looked at: `bound` is null, as for a
    run that never asked, and it says nothing about candidates."""
    doc = jtrio(
        world, "p", "--preset", "strict-denovo", "--preset-ignore", "gnomad_max_af",
        ann="far",
    )  # fmt: skip
    gn = doc["effective"]["gnomad"]
    assert gn["loaded"] is True and gn["bound"] is None
    assert gn["cohort_variants_sampled"] is None
    assert (gn["candidates_checked"], gn["candidates_without_record"]) == (None, None)
    assert not [
        w for w in doc["effective"]["warnings"] if w["code"].startswith("gnomad_")
    ]


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
