"""`vcfclick db segregate` end to end, on both backends.

The family (tests/fixtures/make_segregation_vcf.py): FA MO, S1 S2 affected,
S3 S4 unaffected, S5 unknown, a hand-planted site for each rule and 90 random
background sites. Every expectation is by planted label (the position names a
site), never by a whole-file total, so the fixture's RNG is not part of the
contract. The CLI runs in subprocesses (chDB allows one server per process).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import pytest

from tests.conftest import FIXTURES, bgzip_vcf, run_cli

sys.path.insert(0, str(FIXTURES))
import make_segregation_vcf as plant  # noqa: E402

BACKENDS = ("chdb", "duckdb")
NAME = {pos: name for name, (pos, _, _) in plant.SITES.items()}
NAME.update(
    {
        pos: f"bg{pos}"
        for pos, _ in plant.background(__import__("random").Random(20261020))
    }
)


@pytest.fixture(scope="module")
def annotations(tmp_path_factory):
    """A gene store with G1..G5 for the whole module (a plain env var: the CLI
    subprocesses inherit it)."""
    from annotations.loaders import gencode_genes

    path = tmp_path_factory.mktemp("segann") / "ann.duckdb"
    previous = os.environ.get("VCFCLICK_ANNOTATIONS_DB")
    os.environ["VCFCLICK_ANNOTATIONS_DB"] = str(path)
    gencode_genes.load(FIXTURES / "segregation_genes.gff3")
    yield path
    if previous is None:
        os.environ.pop("VCFCLICK_ANNOTATIONS_DB", None)
    else:
        os.environ["VCFCLICK_ANNOTATIONS_DB"] = previous


@pytest.fixture(scope="module")
def gz(tmp_path_factory):
    base = tmp_path_factory.mktemp("segvcf")
    return {
        "hand": bgzip_vcf(FIXTURES / "segregation.vcf", base / "seg.vcf.gz"),
        "cohort": bgzip_vcf(
            FIXTURES / "segregation_cohortaf.vcf", base / "segco.vcf.gz"
        ),
    }


class Store:
    """One database on one backend."""

    def __init__(self, home: Path, backend: str):
        self.home, self.backend = home, backend

    def cli(self, *args, ok=True):
        return run_cli(self.home, self.backend, *args, ok=ok)

    @lru_cache(maxsize=None)  # noqa: B019 - a module-lifetime cache is the point
    def doc(self, *args) -> dict:
        r = self.cli("db", "segregate", "t", "--format", "json", *args)
        return json.loads(r.stdout)

    def fam(self, *args) -> dict:
        """The one sibship's document (the fixture has one family)."""
        d = self.doc(*args)
        return d["families"][0] if "families" in d else d

    def text(self, *args) -> str:
        return self.cli("db", "segregate", "t", *args).stdout


def _build(tmp_path_factory, backend, vcf, *, keep_reference, ped=None, ingest_id="i1"):
    home = tmp_path_factory.mktemp(f"seg-{backend}")
    run_cli(home, backend, "db", "create", "t")
    args = ["db", "ingest", "t", str(vcf), "--ingest-id", ingest_id, "--serial"]
    if keep_reference:
        args.append("--keep-reference")
    run_cli(home, backend, *args)
    run_cli(home, backend, "db", "ped", "t", str(ped or FIXTURES / "segregation.ped"))
    return Store(home, backend)


@pytest.fixture(scope="module", params=BACKENDS)
def kr(request, tmp_path_factory, gz, annotations):
    """--keep-reference ingest, hand-set AF."""
    return _build(tmp_path_factory, request.param, gz["hand"], keep_reference=True)


@pytest.fixture(scope="module", params=BACKENDS)
def sparse(request, tmp_path_factory, gz, annotations):
    return _build(tmp_path_factory, request.param, gz["hand"], keep_reference=False)


@pytest.fixture(scope="module", params=BACKENDS)
def cohort(request, tmp_path_factory, gz, annotations):
    return _build(tmp_path_factory, request.param, gz["cohort"], keep_reference=True)


def by_label(doc: dict, model: str) -> dict[str, dict]:
    """The listed candidates of a row model by planted name."""
    return {NAME[c["pos"]]: c for c in doc["models"][model]["candidates"]}


def statuses(doc: dict) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for model in ("denovo", "recessive", "dominant"):
        for label, c in by_label(doc, model).items():
            out.setdefault(label, {})[model] = c["status"]
    return out


def comphet(doc: dict) -> dict[str, list[str]]:
    return {
        g["gene"]: [p["status"] for p in g["pairs"]]
        for g in doc["models"]["comphet"]["genes"]
    }


# ───────────────────────────── planted labels ─────────────────────────────


def test_every_planted_site_has_its_expected_status(kr):
    got = statuses(kr.fam("--allow-unresolved"))
    background = {k for k in got if k.startswith("bg")}
    seen = {k: v for k, v in got.items() if k in plant.EXPECTED}
    assert seen == plant.EXPECTED
    for label in plant.ABSENT_EVERYWHERE:
        assert label not in got, (label, got[label])
    # a gene-pair site is a plain dominant site where one parent carries
    assert {k for k in got if k not in plant.EXPECTED and k not in background} <= {
        "g3_mat", "g4_pat", "g5_pat", "g5_mat",
    }  # fmt: skip
    assert not background  # AF 0.30 is dropped by the default --max-af


def test_the_default_listing_leaves_out_compatible_but_counts_it(kr):
    doc = kr.fam()
    dom = doc["models"]["dominant"]
    labels = set(by_label(doc, "dominant"))
    assert {"dom_pos", "dom_mat", "dom_homalt_parent"} <= labels
    assert "dom_unaff_nocall" not in labels
    assert dom["count"] == len(dom["candidates"]) == dom["by_status"]["consistent"]
    assert dom["n_compatible"] == dom["by_status"]["compatible"] > 0
    assert doc["models"]["denovo"]["count"] == 1
    assert doc["models"]["denovo"]["n_compatible"] == 2  # anc_nocall_dn, dn_parent_half


def test_allow_unresolved_lists_what_the_counts_already_said(kr):
    plain, loose = kr.fam(), kr.fam("--allow-unresolved")
    for model in ("denovo", "recessive", "dominant"):
        assert (
            loose["models"][model]["by_status"] == plain["models"][model]["by_status"]
        )
        assert loose["models"][model]["count"] == (
            plain["models"][model]["count"] + plain["models"][model]["n_compatible"]
        )


def test_the_unresolved_child_is_named(kr):
    c = by_label(kr.fam("--allow-unresolved"), "dominant")["dom_unaff_nocall"]
    assert c["unresolved"] == ["S4"] and c["contradicting"] == []
    s4 = next(k for k in c["calls"] if k["sample"] == "S4")
    assert s4["state"] == "no_call" and s4["gt"] is None and s4["evidence"] is None
    low = by_label(kr.fam("--allow-unresolved"), "dominant")["dom_aff_lowab"]
    s2 = next(k for k in low["calls"] if k["sample"] == "S2")
    assert s2["state"] == "low_quality" and s2["gt"] == 1 and s2["evidence"] == "stored"


def test_a_candidate_does_not_depend_on_which_affected_child_is_fine(kr):
    """The anchor sites: the first (or second) affected child is a no-call, low
    AB or low GQ, and the site is still found; with S1 at 0/0 it is not."""
    got = statuses(kr.fam("--allow-unresolved"))
    for label in ("anc_nocall_d", "anc_lowab_d", "anc_lowgq_d", "anc_nocall_dm"):
        assert got[label] == {"dominant": "compatible"}, label
    assert "anc_ref_d" not in got


def test_de_novo_flags_parental_mosaicism_and_drops_alt_read_parents(kr):
    denovo = by_label(kr.fam("--allow-unresolved"), "denovo")
    pos = denovo["dn_pos"]
    assert pos["parental_mosaic_suspect"] is True  # two affected children share it
    assert pos["parent_alt_reads"] == {"father": 0, "mother": 0}
    assert pos["n_informative_meioses"] is None and pos["p_random"] is None
    assert "dn_parent_altreads" not in denovo  # FA 0/0 with 3 alt reads of 30
    assert "dn_unaff_altreads" not in denovo  # the same in an unaffected child
    half = denovo["dn_parent_half"]  # a partial call in S5: the 0/0s are untrusted
    assert half["status"] == "compatible" and half["partial_call_site"] is True
    assert half["father_state"] == "ref_partial_site"
    assert kr.fam()["counts"]["partial_call_sites"] >= 2


def test_a_hom_alt_carrier_parent_is_reported_and_uninformative(kr):
    c = by_label(kr.fam(), "dominant")["dom_homalt_parent"]
    assert c["status"] == "consistent"
    assert c["carrier_parent"] == {"sample": "FA", "affected": None}
    assert c["father_state"] == "hom_alt" and c["father_gt"] == 2
    assert c["n_informative_meioses"] == 0 and c["p_random"] == 1.0


def test_informativeness_numbers(kr):
    doc = kr.fam()
    dom = by_label(doc, "dominant")["dom_pos"]
    assert (dom["n_informative_meioses"], dom["p_random"]) == (4, 0.0625)
    assert (
        dom["n_affected_carrying"],
        dom["n_unaffected_carrying"],
        dom["n_unaffected_ref"],
    ) == (2, 0, 2)
    rec = by_label(doc, "recessive")["rec_pos"]
    assert (
        rec["p_random"] == round(0.25**2 * 0.75**2, 10)
        and rec["n_informative_meioses"] == 4
    )


def test_a_low_ab_parent_het_is_not_gated(kr):
    """Parents' hets are not allele-balance gated, as in `db trio`."""
    c = by_label(kr.fam(), "recessive")["rec_parent_lowab"]
    assert c["status"] == "consistent" and c["father_state"] == "het"


def test_compound_het_pairs(kr):
    doc = kr.fam("--allow-unresolved")
    assert comphet(doc) == {"G1": ["consistent"], "G4": ["compatible"]}
    g1 = next(g for g in doc["models"]["comphet"]["genes"] if g["gene"] == "G1")
    assert g1["pairs"][0]["paternal"] == "chr1:5000100:A:G"
    assert g1["pairs"][0]["maternal"] == "chr1:5000900:A:G"
    assert g1["pairs"][0]["n_informative_meioses"] == 4
    assert [v["calls"][0]["sample"] for v in g1["paternal"]] == ["S1"]
    default = kr.fam("--category", "comphet")
    assert comphet(default) == {"G1": ["consistent"]}  # G4 is compatible
    assert default["models"]["comphet"]["n_compatible"] == 1
    assert default["models"]["comphet"]["by_status"]["consistent"] == 1


def test_sex_and_mitochondrial_contigs_are_skipped_and_counted(kr):
    doc = kr.fam("--allow-unresolved")
    assert doc["counts"]["skipped_sites"] == 3  # chrX, chrM, chrY
    assert any("sex or mitochondrial" in n for n in doc["notes"])
    for model in ("denovo", "recessive", "dominant"):
        assert all(c["chrom"] == "chr1" for c in doc["models"][model]["candidates"])


# ───────────────────────────── reduced penetrance ─────────────────────────────


def test_max_unaffected_carriers_labels_tolerated_candidates(kr):
    doc = kr.fam("--max-unaffected-carriers", "1")
    dom, denovo, rec = (by_label(doc, m) for m in ("dominant", "denovo", "recessive"))
    for table, label in (
        (dom, "dom_unaff_carrier"),
        (denovo, "dn_unaff_carrier"),
        (rec, "rec_unaff_hom"),
    ):
        c = table[label]
        assert c["status"] == "reduced_penetrance" and c["reduced_penetrance"] is True
        assert c["contradicting"] and c["n_unaffected_carrying"] >= 1
    assert "dom_two_unaff" not in dom  # two carriers > N
    assert dom["dom_pos"]["status"] == "consistent"
    assert doc["models"]["dominant"]["by_status"]["reduced_penetrance"] >= 1
    # never labelled consistent
    for model in ("denovo", "recessive", "dominant"):
        for c in doc["models"][model]["candidates"]:
            assert (c["status"] == "consistent") == (
                not c["contradicting"] and not c["unresolved"]
            )


def test_max_unaffected_carriers_must_leave_the_filter_in_force(kr):
    r = kr.cli("db", "segregate", "t", "--max-unaffected-carriers", "2", ok=False)
    assert "would exclude nothing" in r.stderr
    r = kr.cli("db", "segregate", "t", "--max-unaffected-carriers", "-1", ok=False)
    assert "must be >= 0" in r.stderr


def test_max_unaffected_carriers_is_ignored_for_comphet_with_a_note(kr):
    doc = kr.fam("--max-unaffected-carriers", "1")
    assert comphet(doc) == {"G1": ["consistent"]}
    assert any("ignored for comphet" in n for n in doc["notes"])
    assert doc["options"]["max_unaffected_carriers"] == 1


# ───────────────────────────── sparse ingests ─────────────────────────────


def test_sparse_strict_blocks_what_needs_reference_and_keeps_recessive(sparse, kr):
    doc = sparse.fam("--allow-unresolved")
    for model, reason in (
        ("denovo", "needs_keep_reference"),
        ("dominant", "needs_keep_reference"),
        ("comphet", "needs_keep_reference"),
    ):
        m = doc["models"][model]
        assert (
            m["blocked"] is True and m["blocked_reason"] == reason and m["count"] == 0
        )
    assert doc["keep_reference"] is False
    rec = statuses(doc)
    keep = statuses(kr.fam("--allow-unresolved"))
    for label in (
        "rec_pos",
        "rec_unaff_het_pair",
        "rec_unaff_nocall",
        "anc_nocall_r",
        "rec_aff_homalt_lowab",
        "rec_parent_lowab",
    ):
        assert rec[label]["recessive"] == keep[label]["recessive"], label
    text = sparse.text()
    assert "(parents must be stored reference: needs --keep-reference)" in text
    assert "(needs --keep-reference or --absent called)" in text


def test_absent_called_infers_reference_and_is_never_consistent(sparse):
    doc = sparse.fam("--absent", "called", "--allow-unresolved")
    dom = by_label(doc, "dominant")
    assert dom["dom_pos"]["status"] == "consistent_inferred"
    assert dom["dom_pos"]["inferred_reference"] is True
    s3 = next(k for k in dom["dom_pos"]["calls"] if k["sample"] == "S3")
    assert s3["state"] == "absent_called" and s3["evidence"] == "inferred"
    assert all(c["status"] != "consistent" for c in dom.values())
    # de novo needs stored parental 0/0 under every policy
    assert doc["models"]["denovo"]["blocked"] is True
    assert comphet(doc) == {"G1": ["consistent_inferred"], "G4": ["compatible"]}
    assert doc["absent"] == "called"


def test_absent_called_changes_nothing_on_a_keep_reference_ingest_but_its_label(kr):
    strict = statuses(kr.fam("--allow-unresolved"))
    called = statuses(kr.fam("--allow-unresolved", "--absent", "called"))
    assert strict == called  # every reference proof is a stored 0/0


# ───────────────────────────── the AF guard ─────────────────────────────


def test_hand_set_af_is_gated(kr):
    assert kr.fam()["gates"]["af_gate"] == "applied"
    assert "dom_common" not in by_label(kr.fam("--allow-unresolved"), "dominant")
    assert kr.fam()["models"]["dominant"]["n_af_filtered"] >= 1


def test_cohort_derived_af_skips_the_gate_unless_it_is_typed(cohort):
    doc = cohort.fam()
    assert doc["gates"]["af_gate"] == "skipped_cohort_derived"
    assert any("cohort" in n.lower() for n in doc["notes"])
    assert doc["models"]["dominant"]["count"] > 6  # background sites now listed
    typed = cohort.fam("--max-af", "0.01")
    assert typed["gates"]["af_gate"] == "cohort_derived_suspected"
    assert typed["models"]["dominant"]["n_af_filtered"] > 0
    assert typed["models"]["dominant"]["count"] < doc["models"]["dominant"]["count"]


# ───────────────────────────── selecting sibships ─────────────────────────────


def test_selectors_decide_the_shape_of_the_json(kr):
    assert list(kr.doc()) == ["families"] and len(kr.doc()["families"]) == 1
    one = kr.doc("--family", "fam1")
    assert "families" not in one and one["family"]["family_id"] == "fam1"
    assert kr.doc("--father", "FA", "--mother", "MO")["family"]["mother"] == "MO"


def test_unknown_selection_lists_what_exists(kr):
    r = kr.cli("db", "segregate", "t", "--family", "nope", ok=False)
    assert "no sibship matches" in r.stderr and "fam1" in r.stderr


def _ped(tmp_path, text: str) -> Path:
    p = tmp_path / "x.ped"
    p.write_text(text)
    return p


def _reped(store: Store, ped: Path):
    store.cli("db", "ped", "t", str(ped))


@pytest.fixture
def fresh(tmp_path_factory, gz, annotations):
    """A DuckDB copy whose pedigree a test may replace."""
    return _build(tmp_path_factory, "duckdb", gz["hand"], keep_reference=True)


def test_half_siblings_are_separate_sibships(fresh, tmp_path):
    # S3 and S4 share FA with S1 and S2 but have S5 as their mother
    ped = _ped(
        tmp_path,
        "fam1 FA 0 0 1 0\nfam1 MO 0 0 2 0\nfam1 S5 0 0 2 0\n"
        "fam1 S1 FA MO 1 2\nfam1 S2 FA MO 2 1\n"
        "fam1 S3 FA S5 1 2\nfam1 S4 FA S5 2 1\n",
    )
    _reped(fresh, ped)
    doc = fresh.doc()
    assert [
        (f["family"]["father"], f["family"]["mother"]) for f in doc["families"]
    ] == [("FA", "MO"), ("FA", "S5")]
    r = fresh.cli("db", "segregate", "t", "--family", "fam1", ok=False)
    assert "half siblings" in r.stderr
    assert "FA x MO: S1, S2" in r.stderr and "FA x S5: S3, S4" in r.stderr
    one = fresh.doc("--family", "fam1", "--mother", "S5")
    assert one["family"]["mother"] == "S5"


def test_a_single_child_is_a_trio_and_the_message_says_so(fresh, tmp_path):
    _reped(
        fresh, _ped(tmp_path, "fam1 FA 0 0 1 0\nfam1 MO 0 0 2 0\nfam1 S1 FA MO 1 2\n")
    )
    r = fresh.cli("db", "segregate", "t", ok=False)
    assert (
        "no sibship with >=2 children, both parents declared and >=1 affected child"
        in r.stderr
    )
    r = fresh.cli("db", "segregate", "t", "--family", "fam1", ok=False)
    assert "db trio --proband S1" in r.stderr


def test_no_affected_child_is_refused(fresh, tmp_path):
    _reped(
        fresh,
        _ped(
            tmp_path,
            "fam1 FA 0 0 1 0\nfam1 MO 0 0 2 0\nfam1 S1 FA MO 1 1\nfam1 S2 FA MO 2 1\n",
        ),
    )
    r = fresh.cli("db", "segregate", "t", "--family", "fam1", ok=False)
    assert "no affected child" in r.stderr


def test_no_unaffected_child_filters_by_shared_by_all_affected_with_a_note(
    fresh, tmp_path
):
    _reped(
        fresh,
        _ped(
            tmp_path,
            "fam1 FA 0 0 1 0\nfam1 MO 0 0 2 0\nfam1 S1 FA MO 1 2\nfam1 S2 FA MO 2 2\n",
        ),
    )
    doc = fresh.fam("--allow-unresolved")
    assert doc["counts"]["unaffected"] == 0
    assert any("no unaffected child" in n for n in doc["notes"])
    labels = by_label(doc, "dominant")
    # an unaffected carrier no longer excludes: nothing is unaffected
    assert "dom_unaff_carrier" in labels and "dom_two_unaff" in labels


def test_a_single_affected_child_with_unaffected_siblings(fresh, tmp_path):
    _reped(
        fresh,
        _ped(
            tmp_path,
            "fam1 FA 0 0 1 0\nfam1 MO 0 0 2 0\nfam1 S1 FA MO 1 2\nfam1 S2 FA MO 2 1\nfam1 S3 FA MO 1 1\n",
        ),
    )
    doc = fresh.fam()
    assert doc["counts"] == {
        **doc["counts"],
        "affected": 1,
        "unaffected": 2,
        "children": 3,
    }
    # S2 is now unaffected: dom_aff_missing (S2 0/0) has its one affected child
    assert "dom_aff_missing" in by_label(doc, "dominant")


def test_unknown_children_never_constrain_and_are_listed(kr):
    doc = kr.fam()
    kids = {c["sample"]: c["affected"] for c in doc["family"]["children"]}
    assert kids == {
        "S1": "affected",
        "S2": "affected",
        "S3": "unaffected",
        "S4": "unaffected",
        "S5": None,
    }
    # S5 is 0/1 at dom_pos and dom_pos is still consistent
    assert by_label(doc, "dominant")["dom_pos"]["status"] == "consistent"


def test_mz_pair_collapses_informativeness_and_is_validated(kr):
    mz = kr.fam("--mz-pair", "S1,S2")
    c = by_label(mz, "dominant")["dom_pos"]
    assert c["n_informative_meioses"] == 3 and c["p_random"] == 0.125
    assert mz["family"]["mz_pairs"] == [["S1", "S2"]]
    assert by_label(mz, "denovo")["dn_pos"]["parental_mosaic_suspect"] is False
    r = kr.cli("db", "segregate", "t", "--mz-pair", "S1,S3", ok=False)
    assert "must share a phenotype" in r.stderr
    r = kr.cli("db", "segregate", "t", "--mz-pair", "S1,NOPE", ok=False)
    assert "not two children of one sibship" in r.stderr
    r = kr.cli("db", "segregate", "t", "--mz-pair", "S1", ok=False)
    assert "two different sample ids" in r.stderr


def test_parent_phenotype_use_downgrades_an_unaffected_carrier_parent(fresh, tmp_path):
    ped = plant.PED.replace("fam1 FA 0 0 1 0", "fam1 FA 0 0 1 1")
    _reped(fresh, _ped(tmp_path, ped))
    ignore = by_label(fresh.fam(), "dominant")
    assert ignore["dom_pos"]["status"] == "consistent"
    assert ignore["dom_pos"]["carrier_parent"] == {
        "sample": "FA",
        "affected": "unaffected",
    }
    use = fresh.fam("--parent-phenotype", "use")
    dom = by_label(use, "dominant")
    assert dom["dom_pos"]["status"] == "reduced_penetrance"
    assert dom["dom_pos"]["contradicting"] == ["FA"]
    assert dom["dom_mat"]["status"] == "consistent"  # the mother carries there


# ───────────────────────────── the JSON shape ─────────────────────────────

DOC_KEYS = [
    "family", "counts", "keep_reference", "absent", "partial_calls",
    "parent_phenotype", "gates", "options", "annotations", "notes", "models",
]  # fmt: skip
CANDIDATE_KEYS = [
    "chrom", "pos", "ref", "alt", "af", "genes", "gnomad_popmax", "clinvar",
    "status", "father_gt", "mother_gt", "father_state", "mother_state", "calls",
    "n_affected_carrying", "n_unaffected_carrying", "n_unaffected_ref",
    "unresolved", "contradicting", "inferred_reference", "reduced_penetrance",
    "partial_call_site", "parental_mosaic_suspect", "parent_alt_reads",
    "carrier_parent", "n_informative_meioses", "p_random",
]  # fmt: skip


def test_the_json_shape_is_pinned_key_by_key(kr):
    """Any rename fails here: the desktop apps decode this document."""
    doc = kr.fam()
    assert list(doc) == DOC_KEYS
    assert list(doc["family"]) == [
        "ingest_id",
        "family_id",
        "father",
        "mother",
        "children",
        "mz_pairs",
    ]
    assert list(doc["family"]["children"][0]) == ["sample", "affected", "sex"]
    assert doc["family"]["children"][0] == {
        "sample": "S1",
        "affected": "affected",
        "sex": "male",
    }
    assert list(doc["counts"]) == [
        "children",
        "affected",
        "unaffected",
        "unknown",
        "skipped_sites",
        "partial_call_sites",
    ]
    assert list(doc["gates"]) == [
        "min_gq", "min_dp", "max_af", "af_gate", "min_ab", "max_ab",
        "min_homalt_ab", "max_parent_alt_frac", "gnomad_max_af",
    ]  # fmt: skip
    assert list(doc["options"]) == [
        "allow_unresolved",
        "max_unaffected_carriers",
        "limit",
    ]
    assert list(doc["models"]) == ["denovo", "recessive", "dominant", "comphet"]
    for model in ("denovo", "recessive", "dominant"):
        m = doc["models"][model]
        assert list(m) == [
            "count", "n_compatible", "by_status", "n_af_filtered", "blocked",
            "blocked_reason", "truncated", "candidates",
        ]  # fmt: skip
        assert list(m["by_status"]) == [
            "consistent",
            "consistent_inferred",
            "compatible",
            "reduced_penetrance",
        ]
    cand = doc["models"]["dominant"]["candidates"][0]
    assert list(cand) == CANDIDATE_KEYS
    assert list(cand["calls"][0]) == ["sample", "gt", "state", "evidence"]
    assert [c["sample"] for c in cand["calls"]] == ["S1", "S2", "S3", "S4", "S5"]
    com = kr.fam("--category", "comphet")["models"]["comphet"]
    assert list(com) == [
        "count", "n_compatible", "by_status", "n_af_filtered", "blocked",
        "blocked_reason", "needs_genes", "truncated", "genes",
    ]  # fmt: skip
    gene = com["genes"][0]
    assert list(gene) == ["gene", "paternal", "maternal", "pairs"]
    assert list(gene["paternal"][0]) == [
        "chrom",
        "pos",
        "ref",
        "alt",
        "af",
        "calls",
        "genes",
        "gnomad_popmax",
        "clinvar",
    ]
    assert list(gene["pairs"][0]) == [
        "paternal", "maternal", "status", "contradicting", "unresolved",
        "n_informative_meioses", "p_random",
    ]  # fmt: skip


def test_the_state_vocabulary_is_closed(kr):
    from cli.segregation_verdict import STATE_NAMES

    seen = set()
    for model in ("denovo", "recessive", "dominant"):
        for c in kr.fam("--allow-unresolved")["models"][model]["candidates"]:
            seen |= {c["father_state"], c["mother_state"]} | {
                k["state"] for k in c["calls"]
            }
    assert (
        seen <= set(STATE_NAMES)
        and {"hom_ref", "het", "no_call", "low_quality"} <= seen
    )


def test_limit_truncates_and_says_so(kr):
    doc = kr.fam("--limit", "2", "--category", "dominant")
    m = doc["models"]["dominant"]
    assert len(m["candidates"]) == 2 and m["truncated"] is True and m["count"] == 6


def test_text_marks_a_listing_cut_by_limit(kr):
    one = kr.text("--limit", "2", "--category", "dominant")
    assert "dominant candidates: 6" in one and "(showing 2, raise --limit)" in one
    many = kr.text("--limit", "2")
    assert "(truncated, raise --limit)" in many
    assert "raise --limit" not in kr.text("--category", "dominant")


def test_a_pair_cap_is_reported_even_when_no_pair_is_listed(monkeypatch):
    from types import SimpleNamespace

    from cli import segregation as seg

    monkeypatch.setattr(seg, "MAX_PAIRS_PER_GENE", 3)
    monkeypatch.setattr(seg, "children_of", lambda sib, site: ())
    monkeypatch.setattr(seg.sv, "judge_pair", lambda **kw: None)

    def site(pos):
        return SimpleNamespace(
            chrom="chr1",
            pos=pos,
            ref="A",
            alt="G",
            pc=False,
            states=("het", "hom_ref"),
            key=("chr1", pos, "A", "G"),
        )

    sides = {"paternal": [site(1), site(2)], "maternal": [site(3), site(4)]}
    ctx = SimpleNamespace(sib=None, opts=SimpleNamespace())
    entry, cut = seg._gene_entry("G", sides, ctx, None, {})
    assert entry is None and cut is True


def test_a_preset_category_with_xlinked_is_dropped_with_a_note(kr):
    r = kr.cli("db", "segregate", "t", "--format", "json", "--preset", "rare-recessive")
    doc = json.loads(r.stdout)
    doc = doc["families"][0] if "families" in doc else doc
    assert set(doc["models"]) == {"recessive", "comphet"}
    assert "xlinked" in r.stderr
    typed = kr.cli("db", "segregate", "t", "--category", "xlinked", ok=False)
    assert typed.returncode != 0 and "no xlinked model" in typed.stderr


# ───────────────────────────── text ─────────────────────────────


def test_text_header_summary_and_listing(kr):
    out = kr.text()
    lines = out.splitlines()
    assert (
        lines[0]
        == "segregation: family=fam1 father=FA mother=MO affected=S1,S2 unaffected=S3,S4 unknown=S5"
    )
    assert "  denovo          1  (+2 compatible)" in lines
    assert any(x.startswith("  recessive") and "(+3 compatible)" in x for x in lines)
    assert any(
        x.startswith("  comphet") and x.endswith("genes  (+1 compatible)")
        for x in lines
    )
    listing = kr.text("--category", "dominant", "--allow-unresolved")
    assert (
        "  chr1:1000000 A>G  FA=0/1 MO=0/0  S1=0/1* S2=0/1* S3=0/0 S4=0/0 S5=0/1  consistent"
        in listing
    )
    assert "S4=./." in listing and "S4=0/0?" in listing  # no-call, low-quality 0/0
    assert "S2=0/1?" in listing  # the low-allele-balance het: a stored gt with a ?
    # ./1 is stored as a het, and at a partial-call site every stored 0/0 is untrusted
    assert "MO=0/0!  S1=0/1* S2=0/1* S3=0/1 S4=0/0!" in listing
    first_partial = kr.text("--category", "denovo", "--allow-unresolved")
    assert "FA=0/0! MO=0/0!" in first_partial  # untrusted stored reference


def test_a_blocked_model_alone_says_so(sparse):
    out = sparse.text("--category", "dominant")
    assert "dominant: blocked" in out


# ───────────────────────────── gates and presets ─────────────────────────────


def test_gate_flags_take_effect(kr):
    strict = kr.fam("--min-gq", "60", "--allow-unresolved")
    assert statuses(strict) == {}  # every call is GQ 50 (or 99 at the alt-read sites)
    loose_ab = kr.fam("--min-ab", "0.05", "--max-ab", "0.95", "--allow-unresolved")
    assert by_label(loose_ab, "dominant")["dom_aff_lowab"]["status"] == "consistent"
    homalt = kr.fam("--min-homalt-ab", "0.4", "--allow-unresolved")
    assert (
        by_label(homalt, "recessive")["rec_aff_homalt_lowab"]["status"] == "consistent"
    )
    frac = kr.fam("--max-parent-alt-frac", "0.2")
    assert "dn_parent_altreads" in by_label(
        frac, "denovo"
    )  # 3/30 = 0.1 < 0.2 is trusted
    assert frac["gates"]["max_parent_alt_frac"] == 0.2


def test_pass_only_runs(kr):
    doc = kr.fam("--pass-only")
    assert doc["gates"]["pass_only"] is True
    assert doc["models"]["dominant"]["count"] == 6  # every record is PASS


def test_a_preset_is_echoed_and_its_keys_resolve(kr):
    doc = kr.fam("--preset", "rare-dominant")
    assert doc["preset"]["name"] == "rare-dominant"
    resolved = doc["preset"]["resolved"]
    assert resolved["min_gq"]["from"] in ("preset", "default")
    assert "effective" in doc and "warnings" in doc["effective"]
    text = kr.text("--preset", "rare-dominant")
    assert "preset: rare-dominant" in text


def test_a_strict_run_reports_whether_the_gates_bind(kr):
    doc = kr.fam("--strict-missing", "--require-bound-gates")
    assert doc["gates"]["strict_missing"] is True
    assert "effective" in doc
    assert doc["models"]["dominant"]["count"] > 0


def test_scope_regions_limit_every_model(kr):
    region = "chr1:1000000-1000900"
    doc = kr.fam("--region", region, "--allow-unresolved")
    assert doc["scope"]["applied"] == "sql"
    for model in ("denovo", "recessive", "dominant"):
        for c in doc["models"][model]["candidates"]:
            assert 1_000_000 <= c["pos"] <= 1_000_900
    assert {NAME[c["pos"]] for c in doc["models"]["dominant"]["candidates"]} >= {
        "dom_pos",
        "dom_mat",
    }
    assert doc["models"]["denovo"]["count"] == 0
    # compound het: inclusions act after the grouping
    g = kr.fam("--category", "comphet", "--region", "chr1:5000000-5001000")
    assert comphet(g) == {"G1": ["consistent"]}
    excluded = kr.fam(
        "--exclude-region", "chr1:5000100-5000100", "--category", "comphet"
    )
    assert comphet(excluded) == {}


def test_xlinked_is_not_a_category_here(kr):
    r = kr.cli("db", "segregate", "t", "--category", "xlinked", ok=False)
    assert "no xlinked model" in r.stderr


def test_the_new_options_are_registered_with_their_reasons():
    from cli import filter_registry as fr

    for key in ("absent", "max_unaffected_carriers"):
        assert key in fr.NEVER_PRESET and "weaken" in fr.NEVER_PRESET[key]
    segregate_only = {k.file_key for k in fr.keys_for("segregate")} - {
        k.file_key for k in fr.keys_for("trio")
    }
    assert segregate_only == {
        "min_homalt_ab",
        "partial_calls",
        "parent_phenotype",
        "allow_unresolved",
    }


# ───────────────────── backends agree, order does not matter ─────────────────────

CONFIGS = [
    (),
    ("--allow-unresolved",),
    ("--allow-unresolved", "--max-unaffected-carriers", "1"),
    ("--allow-unresolved", "--absent", "called"),
    ("--allow-unresolved", "--partial-calls", "ignore"),
    ("--allow-unresolved", "--mz-pair", "S1,S2"),
]


@pytest.fixture(scope="module")
def both(tmp_path_factory, gz, annotations):
    return {
        b: _build(tmp_path_factory, b, gz["hand"], keep_reference=True)
        for b in BACKENDS
    }


@pytest.fixture(scope="module")
def both_sparse(tmp_path_factory, gz, annotations):
    return {
        b: _build(tmp_path_factory, b, gz["hand"], keep_reference=False)
        for b in BACKENDS
    }


@pytest.mark.parametrize("config", CONFIGS, ids=lambda c: " ".join(c) or "default")
def test_chdb_and_duckdb_give_the_same_document(both, config):
    docs = {b: both[b].doc(*config) for b in BACKENDS}
    assert docs["chdb"] == docs["duckdb"]


@pytest.mark.parametrize("config", CONFIGS[:4], ids=lambda c: " ".join(c) or "default")
def test_backends_agree_on_a_sparse_ingest(both_sparse, config):
    docs = {b: both_sparse[b].doc(*config) for b in BACKENDS}
    assert docs["chdb"] == docs["duckdb"]


_REVERSED = """
import json, sys
from dataclasses import replace
from cli import segregation as seg
from cli.pedigree_load import load_pedigree_graph
from cli.segregation import RunContext, SegOptions
from cli.trio_analysis import Gates, _has_reference_rows
from cli.trio_output import open_render
from storage import get_session

sess = get_session("t")
graph = load_pedigree_graph(sess, "i1")
sib = seg.with_ingest(seg.find_sibships(graph)[0][0], "i1")
gates = Gates(20, 10, 0.01, 0.25, 0.75)
opts = SegOptions(models=seg.MODELS, gates=gates, allow_unresolved=True)
out = {}
with open_render(None, 1000) as render:
    for label, s in (("forward", sib), ("reversed", replace(sib, children=tuple(reversed(sib.children))))):
        ctx = RunContext(sess=sess, sib=s, opts=opts, has_reference=_has_reference_rows(sess),
                         caps=seg.probe_caps(), n_samples=seg.ingest_n_samples(sess, "i1"))
        doc = seg.segregation_document(ctx, render, [])
        out[label] = {
            m: sorted((c["pos"], c["status"], tuple(sorted(c["unresolved"]))) for c in doc["models"][m]["candidates"])
            for m in ("denovo", "recessive", "dominant")
        }
print(json.dumps(out))
"""


def test_the_result_is_invariant_to_child_order(both, run_python):
    for backend, store in both.items():
        out = run_python(
            store.home, _REVERSED, VCFCLICK_BACKEND=backend, VCFCLICK_DB_NAME="t"
        )
        d = json.loads(out[out.index("{") :])
        assert d["forward"] == d["reversed"], backend
        assert d["forward"]["dominant"], backend


_API = """
import json
from cli.segregation_api import segregate_api
print(json.dumps(segregate_api(family="fam1", allow_unresolved=True, limit=1000)))
"""


def test_the_read_only_api_path_gives_the_cli_document(both, run_python):
    """The web endpoint and the MCP tool run on the engine-locked read-only
    session (chDB `readonly = 1`, DuckDB external access off): the same
    document the command prints, on both backends."""
    for backend, store in both.items():
        out = run_python(
            store.home, _API, VCFCLICK_BACKEND=backend, VCFCLICK_DB_NAME="t"
        )
        api = json.loads(out[out.index("{") :])
        cli = store.doc("--family", "fam1", "--allow-unresolved")
        assert api["models"] == cli["models"], backend
        assert api["counts"] == cli["counts"], backend


def test_a_scope_line_is_printed_without_a_preset(kr):
    text = kr.text("--category", "dominant", "--region", "chr1:1000000-1000900")
    assert "  scope: 1 region" in text


# ───────────────────────────── old databases ─────────────────────────────

_DROP = """
from storage import get_session
sess = get_session("t")
sess.query("DROP TABLE missing_genotypes", "CSV")
for col in ("n_called", "an_called", "ac_called"):
    sess.query(f"ALTER TABLE variants DROP COLUMN {col}", "CSV")
"""


def test_an_old_database_gets_the_reduced_builder_not_an_error(
    tmp_path_factory, gz, annotations, run_python
):
    """No missing_genotypes table and no variants.n_called: state 3 and 7 never
    arise, every absent member is absent_unknown, stored 0/0 rows still prove
    reference, and a note says to re-ingest."""
    kr = _build(tmp_path_factory, "duckdb", gz["hand"], keep_reference=True)
    run_python(kr.home, _DROP, VCFCLICK_BACKEND="duckdb", VCFCLICK_DB_NAME="t")
    doc = kr.fam("--allow-unresolved")
    assert any("no complete missing-call accounting" in n for n in doc["notes"])
    assert doc["counts"]["partial_call_sites"] == 0
    dom = by_label(doc, "dominant")
    assert dom["dom_pos"]["status"] == "consistent"  # stored reference proves it
    assert dom["dom_unaff_nocall"]["status"] == "compatible"
    # without the partial-call accounting a '0/.' elsewhere no longer demotes
    assert by_label(doc, "denovo")["dn_pos"]["status"] == "consistent"
    sp = _build(tmp_path_factory, "duckdb", gz["hand"], keep_reference=False)
    run_python(sp.home, _DROP, VCFCLICK_BACKEND="duckdb", VCFCLICK_DB_NAME="t")
    rec = by_label(sp.fam("--allow-unresolved"), "recessive")
    assert rec["rec_pos"]["status"] == "compatible"
    states = {
        k["state"] for k in rec["rec_pos"]["calls"] if k["sample"] in ("S3", "S4")
    }
    assert states == {"absent_unknown"}


# ───────────────────────────── independent oracles ─────────────────────────────


def _vcf_calls() -> dict[int, list[dict]]:
    """pos -> per-sample {gt, gq, dp, ad} read from the VCF text alone."""
    out: dict[int, list[dict]] = {}
    for line in (FIXTURES / "segregation.vcf").read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[0] != "chr1":
            continue
        calls = []
        for cell in f[9:]:
            gt, gq, dp, ad = cell.split(":")
            ref, alt = (None, None) if ad == "." else map(int, ad.split(","))
            calls.append(
                {
                    "gt": gt,
                    "gq": None if gq == "." else int(gq),
                    "dp": None if dp == "." else int(dp),
                    "ref": ref,
                    "alt": alt,
                }
            )
        out[int(f[1])] = calls
    return out


def _oracle(pos: int, calls: list[dict], af: float, model: str) -> bool:
    """`consistent` for one site, derived from the VCF text with no SQL and no
    verdict function: the rules of cli/segregation_verdict.py written again."""
    if af > 0.01 + 1e-7:
        return False
    if any("." in c["gt"] and c["gt"] != "./." for c in calls):
        return False  # a partial call anywhere: reference is not proven

    def ok(c, role):
        if c["gt"] == "./.":
            return False
        if c["gq"] < 20 or c["dp"] < 10:
            return False
        frac = c["alt"] / (c["alt"] + c["ref"]) if c["alt"] + c["ref"] else None
        if (
            role in ("child", "affected")
            and c["gt"] == "0/1"
            and not 0.25 <= frac <= 0.75
        ):
            return False
        if role == "affected" and c["gt"] == "1/1" and frac < 0.85:
            return False
        return True

    def ref(c, strict_alt):
        if c["gt"] != "0/0" or not ok(c, "any"):
            return False
        return not (
            strict_alt and c["alt"] >= 2 and c["alt"] / (c["alt"] + c["ref"]) >= 0.05
        )

    fa, mo, s1, s2, s3, s4, _s5 = calls
    het = lambda c, role: c["gt"] == "0/1" and ok(c, role)  # noqa: E731
    hom = lambda c, role: c["gt"] == "1/1" and ok(c, role)  # noqa: E731
    carries = lambda c, role: (c["gt"] in ("0/1", "1/1")) and ok(c, role)  # noqa: E731
    if model == "dominant":
        return (
            het(s1, "affected") and het(s2, "affected") and ref(s3, False) and ref(s4, False)
            and ((carries(fa, "p") and ref(mo, False)) or (ref(fa, False) and carries(mo, "p")))
        )  # fmt: skip
    if model == "denovo":
        return (
            carries(s1, "affected") and carries(s2, "affected") and ref(fa, True) and ref(mo, True)
            and ref(s3, True) and ref(s4, True)
        )  # fmt: skip
    return (
        hom(s1, "affected") and hom(s2, "affected") and het(fa, "p") and het(mo, "p")
        and all(ref(c, False) or het(c, "child") for c in (s3, s4))
    )  # fmt: skip


@pytest.mark.parametrize("model", ["dominant", "denovo", "recessive"])
def test_consistent_equals_a_reading_of_the_vcf_text(kr, model):
    calls = _vcf_calls()
    afs = {}
    for line in (FIXTURES / "segregation.vcf").read_text().splitlines():
        if not line.startswith("#") and line.startswith("chr1\t"):
            f = line.split("\t")
            afs[int(f[1])] = float(f[7].split("=")[1])
    want = {p for p, c in calls.items() if _oracle(p, c, afs[p], model)}
    got = {
        c["pos"]
        for c in kr.fam()["models"][model]["candidates"]
        if c["status"] == "consistent"
    }
    assert got == want and want


def _bcftools_positions(expr: str) -> set[int]:
    r = subprocess.run(
        ["bcftools", "query", "-i", expr, "-f", "%POS\n", str(FIXTURES / "segregation.vcf")],
        capture_output=True, text=True, check=True,
    )  # fmt: skip
    return {int(x) for x in r.stdout.split()}


def _q(i: int) -> str:
    return f"FMT/GQ[{i}]>=20 && FMT/DP[{i}]>=10"


def _ab(i: int, lo: float, hi: float) -> str:
    return f"(FMT/AD[{i}:1]/(FMT/AD[{i}:0]+FMT/AD[{i}:1]))>={lo} && (FMT/AD[{i}:1]/(FMT/AD[{i}:0]+FMT/AD[{i}:1]))<={hi}"


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not on PATH")
@pytest.mark.parametrize("model", ["dominant", "denovo", "recessive"])
def test_consistent_equals_independent_bcftools_expressions(kr, model):
    """No vcfclick code on one side: bcftools evaluates the segregation rule on
    the VCF text. Sites with a partial call are removed from its answer, as
    their reference is deliberately not proven."""
    quality = " && ".join(_q(i) for i in range(6))
    partial = {
        int(line.split("\t")[1])
        for line in (FIXTURES / "segregation.vcf").read_text().splitlines()
        if not line.startswith("#")
        and any(c.split(":")[0] in ("0/.", "./1", "./0") for c in line.split("\t")[9:])
    }
    if model == "dominant":
        expr = (
            'INFO/AF<=0.01 && GT[2]="het" && GT[3]="het" && GT[4]="RR" && GT[5]="RR" '
            '&& ((GT[0]="alt" && GT[1]="RR") || (GT[0]="RR" && GT[1]="alt")) '
            f"&& {quality} && {_ab(2, 0.25, 0.75)} && {_ab(3, 0.25, 0.75)}"
        )
    elif model == "denovo":
        ceiling = " && ".join(
            f"(FMT/AD[{i}:1]<2 || FMT/AD[{i}:1]/(FMT/AD[{i}:0]+FMT/AD[{i}:1])<0.05)"
            for i in (0, 1, 4, 5)
        )
        expr = (
            'INFO/AF<=0.01 && GT[2]="alt" && GT[3]="alt" && GT[0]="RR" && GT[1]="RR" '
            f'&& GT[4]="RR" && GT[5]="RR" && {quality} && {ceiling} '
            f"&& {_ab(2, 0.25, 0.75)} && {_ab(3, 0.25, 0.75)}"
        )
    else:
        unaff = " && ".join(
            f'(GT[{i}]="RR" || (GT[{i}]="het" && {_ab(i, 0.25, 0.75)}))' for i in (4, 5)
        )
        expr = (
            'INFO/AF<=0.01 && GT[2]="AA" && GT[3]="AA" && GT[0]="het" && GT[1]="het" '
            f"&& {unaff} && {quality} && {_ab(2, 0.85, 1)} && {_ab(3, 0.85, 1)}"
        )
    want = _bcftools_positions(f'CHROM="chr1" && {expr}') - partial
    got = {
        c["pos"]
        for c in kr.fam()["models"][model]["candidates"]
        if c["status"] == "consistent"
    }
    assert got == want and want


# ───────────────────────────── in-process unit pieces ─────────────────────────────


def test_select_models_orders_and_refuses():
    import click

    from cli.segregation import select_models

    assert select_models(("all",)) == ("denovo", "recessive", "dominant", "comphet")
    assert select_models(("comphet", "denovo")) == ("denovo", "comphet")
    with pytest.raises(click.UsageError, match="no xlinked model"):
        select_models(("xlinked",))


def test_mz_groups_join_overlapping_pairs():
    from cli.segregation import mz_groups
    from cli.segregation_sql import ChildInfo, Sibship

    sib = Sibship(
        "i",
        "f",
        "FA",
        "MO",
        tuple(ChildInfo(s, "affected") for s in ("A", "B", "C", "D")),
    )
    g = mz_groups([("A", "B"), ("B", "C")], sib)
    assert g["A"] == g["B"] == g["C"] and "D" not in g


def test_the_sql_has_no_dialect_specific_text(kr):
    """build_sql is plain SQL both engines parse: no ClickHouse-only function
    and no DuckDB-only syntax appears outside the storage helpers."""
    from cli.segregation_sql import Caps, ChildInfo, Sibship, SqlOptions, build_sql
    from cli.trio_analysis import Gates

    sib = Sibship(
        "i1",
        "f",
        "FA",
        "MO",
        (ChildInfo("A", "affected"), ChildInfo("U", "unaffected")),
    )
    sql = build_sql(
        "dominant",
        sib,
        Gates(20, 10, 0.01, 0.25, 0.75),
        SqlOptions(n_samples=4),
        Caps(),
    )
    for banned in ("countIf", "arrayJoin", "FILTER (", "toUInt", "QUALIFY", "::"):
        assert banned not in sql
    assert "n_samples - " not in sql  # UInt32 would wrap on chDB
    assert "n_called + coalesce(x.n_miss, 0)" in sql
