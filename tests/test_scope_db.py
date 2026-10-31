"""`db trio` with a gene list or regions, against real data on both backends.

One ingest per backend and contig naming (module-scoped): tests/fixtures/
make_scope_vcf.py plants candidates of every model inside and outside a small
panel (C9orf72, C19orf12, C12orf57, whose HGNC symbols have lower-case
letters), on chr7/chrX (`c`) and on bare 7/X (`b`). Every site passes every
default gate, so only the scope removes one. Assertions are on planted
POSITIONS and gene names (scope_truth.json), never totals. About 30 CLI runs
per backend; the pure checks are in tests/test_scope_filter.py.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from unittest import mock

import pytest

from tests.conftest import FIXTURES, _run_python, run_cli

BACKENDS = ("chdb", "duckdb")
TRUTH = json.loads((FIXTURES / "scope_truth.json").read_text())
SITES = TRUTH["sites"]
NAME_AT = {(s["chrom"].removeprefix("chr"), s["pos"]): n for n, s in SITES.items()}
PANEL = ",".join(TRUTH["panel"])

# `db trio` in a fresh interpreter with the scope forced into one SQL mode:
# `python` (the backstop alone), `sql_buckets` (the bucket superset) or `sql`
# (exact BETWEENs), whatever the number of intervals
_FORCED = """
import sys
import cli.main
import cli.scope
limits = {"python": (0, 0), "sql_buckets": (10_000, 0), "sql": (10_000, 10_000)}
cli.scope.SCOPE_SQL_MAX_INTERVALS, cli.scope.SCOPE_SQL_EXACT_MAX = limits[sys.argv[1]]
from click.testing import CliRunner
r = CliRunner().invoke(cli.main.cli, sys.argv[2:])
sys.stdout.write(r.output)
sys.exit(r.exit_code)
"""


@dataclass(frozen=True)
class World:
    home: Path
    backend: str
    ann: dict[str, Path]


@pytest.fixture(scope="module", params=BACKENDS)
def backend(request):
    return request.param


def _load_annotations(path: Path, *, hgnc: bool) -> None:
    from annotations.loaders import gencode_genes, hgnc_aliases

    with mock.patch.dict(os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(path)}):
        gencode_genes.load(FIXTURES / "scope_genes.gff3")
        if hgnc:
            hgnc_aliases.load(FIXTURES / "hgnc_aliases.tsv")


@pytest.fixture(scope="module")
def world(backend, tmp_path_factory):
    base = tmp_path_factory.mktemp(f"scope-{backend}")
    home = base / "home"
    home.mkdir()
    ann = {k: base / f"{k}.duckdb" for k in ("genes", "aliases", "empty")}
    _load_annotations(ann["genes"], hgnc=False)
    _load_annotations(ann["aliases"], hgnc=True)
    from annotations.db import get_connection

    with mock.patch.dict(os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(ann["empty"])}):
        get_connection().close()
    for db, vcf in (("c", "scope.vcf.gz"), ("b", "scope_bare.vcf.gz")):
        run_cli(home, backend, "db", "create", db)
        run_cli(
            home, backend, "db", "ingest", db, str(FIXTURES / vcf),
            "--ingest-id", "g1", "--serial", "--keep-reference",
        )  # fmt: skip
        run_cli(home, backend, "db", "ped", db, str(FIXTURES / "scope.ped"))
    return World(home, backend, ann)


def trio(w: World, db: str, *args: str, ann: str = "genes", ok: bool = True):
    with mock.patch.dict(os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(w.ann[ann])}):
        return run_cli(
            w.home, w.backend, "db", "trio", db, "--proband", "CHILD", *args, ok=ok
        )


def jtrio(w: World, db: str, *args: str, ann: str = "genes") -> dict:
    return json.loads(trio(w, db, *args, "--format", "json", ann=ann).stdout)


def forced(w: World, db: str, mode: str, *args: str) -> dict:
    out = _run_python(
        w.home, _FORCED, mode,
        "db", "trio", db, "--proband", "CHILD", *args, "--format", "json",
        VCFCLICK_BACKEND=w.backend, VCFCLICK_ANNOTATIONS_DB=str(w.ann["genes"]),
    )  # fmt: skip
    return json.loads(out[out.index("{") :])


def names(doc: dict, model: str) -> list[str]:
    """The planted names a model reports: rows by name, genes by symbol."""
    m = doc["models"][model]
    if model == "comphet":
        return [g["gene"] for g in m["genes"]]
    return [
        NAME_AT[(c["chrom"].removeprefix("chr"), c["pos"])] for c in m["candidates"]
    ]


def models_by_name(doc: dict) -> dict[str, list[str]]:
    return {m: names(doc, m) for m in TRUTH["in_panel"]}


@pytest.fixture(scope="module")
def panel_runs(world):
    return {db: jtrio(world, db, "--genes", PANEL) for db in ("c", "b")}


# ───────────────────────────── genes ─────────────────────────────


@pytest.mark.parametrize("db", ["c", "b"])
def test_the_panel_keeps_exactly_the_planted_rows(panel_runs, db):
    """chr7 and bare 7 alike (the gene store says chr7), and the lower-case
    HGNC symbols resolve."""
    assert models_by_name(panel_runs[db]) == TRUTH["in_panel"]


def test_the_panel_echo(panel_runs):
    scope = panel_runs["c"]["scope"]
    assert scope["genes"]["requested"] == 3 and scope["genes"]["genes"] == 3
    assert scope["genes"]["resolved_by"] == {"symbol": 3}
    assert scope["applied"] == "sql" and scope["intervals"] == 3
    # OVL1 overlaps C19orf12 and is not listed: dropped, with OTHER1
    assert scope["comphet_genes_dropped_out_of_scope"] == 2
    assert scope["n_outside"] == 0  # gene spans: every reported partner is inside
    assert "scope" not in jtrio_plain(panel_runs)


_PLAIN: dict = {}


def jtrio_plain(panel_runs):  # the run with no scope, shared by the tests below
    return _PLAIN["doc"]


@pytest.fixture(scope="module", autouse=True)
def _plain(world):
    _PLAIN["doc"] = jtrio(world, "c")
    yield
    _PLAIN.clear()


def test_without_a_scope_every_planted_row_is_there(world):
    doc = _PLAIN["doc"]
    assert "scope" not in doc
    every = {n for n, s in SITES.items() if s["model"] == "denovo"}
    assert set(names(doc, "denovo")) == every
    assert sorted(names(doc, "comphet")) == ["C12orf57", "C19orf12", "OTHER1", "OVL1"]


def test_sql_and_python_paths_agree(world, panel_runs):
    """The backstop alone, the bucket superset and the exact SQL give the same
    rows and counts as the default, with genes and an exclusion."""
    args = ("--genes", PANEL, "--exclude-region", "chr7:35000-36000")
    want = jtrio(world, "c", *args)
    for mode in ("python", "sql_buckets", "sql"):
        doc = forced(world, "c", mode, *args)
        assert doc["models"] == want["models"], mode
        assert doc["scope"]["applied"] == mode


@pytest.mark.parametrize("token", ["c9orf72", "ENSG00000147894.99"])
def test_one_gene_by_lower_case_symbol_or_ensembl_id(world, token):
    doc = jtrio(world, "b", "--genes", token, "--category", "denovo")
    assert names(doc, "denovo") == ["dn_c9"]


def test_flank_reaches_a_promoter(world):
    doc = jtrio(world, "c", "--genes", "OTHER1", "--gene-flank", "5000")
    assert names(doc, "dominant") == ["ch_other_pat", "ch_other_mat", "dom_out"]
    assert names(doc, "recessive") == []  # rec_out at 90000 is 10 kb past the end


def test_exclusions_act_before_the_comphet_grouping(world):
    doc = jtrio(
        world, "c", "--exclude-genes", "C9orf72",
        "--exclude-region", "chr7:59000-59000", "--exclude-region", "chr7:29999-29999",
    )  # fmt: skip
    assert "dn_c9" not in names(doc, "denovo") and "dn_before_edge" not in names(
        doc, "denovo"
    )
    assert "dn_at_start" in names(doc, "denovo")
    # C12orf57 lost its maternal partner before grouping, so it is no compound het
    assert sorted(names(doc, "comphet")) == ["C19orf12", "OTHER1", "OVL1"]
    assert doc["scope"]["exclude_genes"]["genes"] == 1


def test_unresolved_genes_stop_a_short_list(world):
    r = trio(world, "c", "--genes", "C9orf72,FAM123X", ok=False)
    assert "1 of 2 --genes tokens did not resolve: FAM123X" in r.stderr
    r = trio(world, "c", "--genes", "C9orf72,FAM123X", "--allow-unresolved-genes")
    assert "  scope: 1 gene (1 not found: FAM123X)" in r.stdout
    assert "  WARNING: 1 of 2 --genes tokens did not resolve: FAM123X" in r.stdout


def test_a_gene_list_needs_the_gene_table(world):
    r = trio(world, "c", "--genes", "C9orf72", ann="empty", ok=False)
    assert "vcfclick annotations load" in r.stderr
    # regions need no gene table
    doc = jtrio(world, "c", "--region", "chr7:15000-15000", ann="empty")
    assert names(doc, "denovo") == ["dn_c9"]


def test_previous_symbols_need_the_hgnc_table(world):
    r = trio(world, "c", "--genes", "FTDALS1", ok=False)
    assert "load-hgnc" in r.stderr
    doc = jtrio(world, "c", "--genes", "FTDALS1,HGNC:25443", ann="aliases")
    assert doc["scope"]["genes"]["resolved_by"] == {"previous": 1, "hgnc_id": 1}
    assert names(doc, "denovo") == ["dn_c9", "dn_at_start", "dn_at_end"]


def test_a_grch37_database_refuses_a_gene_list(world):
    r = trio(world, "c", "--genes", "C9orf72", "--build", "GRCh37", ok=False)
    assert "GRCh37" in r.stderr


# ───────────────────────────── regions ─────────────────────────────


@pytest.mark.parametrize("db, region", [("c", "chr7:30000-40000"), ("b", "chr7:30000-40000"), ("c", "7:30000-40000")])  # fmt: skip
def test_region_ends_are_inclusive_whatever_the_contig_spelling(world, db, region):
    doc = jtrio(world, db, "--region", region, "--category", "denovo")
    assert names(doc, "denovo") == ["dn_at_start", "dn_at_end"]


def test_a_bed_file_is_zero_based(world, tmp_path):
    bed = tmp_path / "edge.bed"
    bed.write_text("chr7\t29999\t40000\n")  # = chr7:30000-40000
    doc = jtrio(world, "c", "--region-file", str(bed), "--category", "denovo")
    assert names(doc, "denovo") == ["dn_at_start", "dn_at_end"]


def test_comphet_keeps_a_gene_with_one_partner_in_a_region(world):
    """A capture edge: the paternal het at 55000 is inside chr7:50000-56000, the
    maternal one at 59000 is not. The gene is kept and the partner flagged."""
    doc = jtrio(world, "c", "--region", "chr7:50000-56000")
    [gene] = doc["models"]["comphet"]["genes"]
    assert gene["gene"] == "C12orf57"
    assert [v["in_scope"] for v in gene["paternal"]] == [True]
    assert [v["in_scope"] for v in gene["maternal"]] == [False]
    assert doc["scope"]["n_outside"] == 1
    assert doc["scope"]["comphet_genes_dropped_out_of_scope"] == 3
    text = trio(world, "c", "--region", "chr7:50000-56000", "--category", "comphet")
    assert "maternal chr7:59000 A>G  AF=0.0005" in text.stdout
    assert "(outside the scope)" in text.stdout


def test_a_listed_and_an_unlisted_gene_on_one_variant(world):
    doc = jtrio(world, "c", "--genes", "C19orf12", "--category", "comphet")
    assert names(doc, "comphet") == ["C19orf12"]  # not OVL1


def _bcftools_positions(vcf: Path, *args: str) -> set[int]:
    out = subprocess.run(
        ["bcftools", "view", "-H", *args, str(vcf)],
        capture_output=True, text=True, check=True,
    ).stdout  # fmt: skip
    return {int(line.split("\t")[1]) for line in out.splitlines() if line}


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not on PATH")
def test_regions_and_genes_agree_with_bcftools(world, tmp_path):
    """The scope's rows are the unscoped rows that bcftools finds in the same
    intervals, by region (-r) and by gene expressed as a BED (-R)."""
    vcf = FIXTURES / "scope.vcf.gz"
    plain = _PLAIN["doc"]
    every = {
        m: {c["pos"] for c in plain["models"][m]["candidates"]}
        for m in ("denovo", "dominant", "recessive")
    }
    region = "chr7:30000-60000"
    inside = _bcftools_positions(vcf, "-r", region)
    doc = jtrio(world, "c", "--region", region)
    for m, positions in every.items():
        assert {
            c["pos"] for c in doc["models"][m]["candidates"]
        } == positions & inside, m
    bed = tmp_path / "panel.bed"
    genes = {"C9orf72": (10_000, 20_000), "C12orf57": (50_000, 60_000)}
    bed.write_text("".join(f"chr7\t{lo - 1}\t{hi}\n" for lo, hi in genes.values()))
    inside = _bcftools_positions(vcf, "-R", str(bed))
    doc = jtrio(world, "c", "--genes", ",".join(genes))
    for m, positions in every.items():
        assert {
            c["pos"] for c in doc["models"][m]["candidates"]
        } == positions & inside, m


# ───────────────────────────── counts, presets, text ─────────────────────────────


def test_the_count_lines_equal_the_listings(world, panel_runs):
    """`--category all` text counts each model with analysis.count(); the JSON
    count is the length of its rows. Both under the scope, on both paths."""
    text = trio(world, "c", "--genes", PANEL).stdout.splitlines()
    counts = {ln.split()[0]: int(ln.split()[1]) for ln in text[2:] if ln.strip()}
    doc = panel_runs["c"]
    assert counts == {m: doc["models"][m]["count"] for m in doc["models"]}
    assert text[1].startswith("  scope: 3 genes; gene spans only")


def test_a_preset_carries_its_gene_list(world):
    with mock.patch.dict(
        os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(world.ann["genes"])}
    ):
        run_cli(
            world.home, world.backend, "preset", "save", "id-panel", "--force",
            "--genes", PANEL, "--genes-label", "test panel v1", "--min-gq", "20",
        )  # fmt: skip
    doc = jtrio(world, "c", "--preset", "id-panel")
    assert models_by_name(doc) == TRUTH["in_panel"]
    assert doc["preset"]["requested"]["genes"] == TRUTH["panel"]
    assert doc["preset"]["resolved"]["genes"] == {
        "value": TRUTH["panel"],
        "from": "preset",
    }
    assert doc["scope"]["genes"]["label"] == "test panel v1"
    assert "effective" in doc
    text = trio(world, "c", "--preset", "id-panel").stdout.splitlines()
    assert text[1].startswith("preset: id-panel (user, ")
    assert "genes" not in text[2]  # the gates line; the list has its own line
    assert any(ln.startswith("  scope: 3 genes [test panel v1]") for ln in text)
    # look outside the panel for one run
    off = jtrio(world, "c", "--preset", "id-panel", "--preset-ignore", "genes")
    assert "scope" not in off
    assert off["preset"]["overridden"][0]["flag"] == "off"
    # a typed list replaces the preset's, and the preset's label goes with it
    typed = jtrio(world, "c", "--preset", "id-panel", "--genes", "C9orf72")
    assert names(typed, "denovo") == ["dn_c9"]
    assert typed["scope"]["genes"]["label"] is None


def test_unresolved_genes_are_an_effective_warning_with_a_preset(world):
    doc = jtrio(
        world, "c", "--preset", "broad-first-pass", "--genes", "C9orf72,FAM123X",
        "--allow-unresolved-genes",
    )  # fmt: skip
    codes = [w["code"] for w in doc["effective"]["warnings"]]
    assert "genes_unresolved" in codes
