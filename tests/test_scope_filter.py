"""The scope (cli/scope.py) without a genotype database.

Parsing, intervals, the SQL predicate, gene resolution against an isolated
annotation store (DuckDB in process, never chDB), the unresolved-genes rule,
the echo, and how a gene list travels through the registry, Click and presets.
The database-backed checks are in tests/test_scope_db.py and the model drift
test in tests/test_scope_models_drift.py.
"""

from __future__ import annotations

import json

import click
import pytest
from click.testing import CliRunner

from cli import filter_registry as fr
from cli import scope as sc
from cli.filter_options import file_option_tokens, gate_options, typed_list_values
from cli.preset_options import PresetCommand, preset_option, resolve_filters
from cli.preset_schema import fingerprint, validate_document
from tests.conftest import FIXTURES

GFF = FIXTURES / "scope_genes.gff3"
HGNC = FIXTURES / "hgnc_aliases.tsv"
PANEL = ["C9orf72", "C19orf12", "C12orf57"]


@pytest.fixture
def genes_loaded(isolated_annotation_db):
    from annotations.loaders.gencode_genes import load

    load(GFF)
    return isolated_annotation_db


def _scope(**kw) -> sc.Scope | None:
    build = kw.pop("build", "GRCh38")
    return sc.build_scope(sc.ScopeRequest.from_values(kw), build=build)


# ───────────────────────────── values and intervals ─────────────────────────────


def test_region_tokens_and_their_checks():
    key = fr.REGISTRY["regions"]
    assert fr.check_value(key, ["chr7:100-200", "7:5-5", "X:1-10"])[1] is None
    for bad in (
        "chr7:200-100",
        "chr7:0-10",
        "chr7:100",
        "chr7_alt:1-2",
        "chr7:1,000-2",
    ):
        assert fr.check_value(key, [bad])[1] is not None, bad
    assert fr.parse_region("chr7:30000-40000") == ("chr7", 30000, 40000)


def test_gene_tokens_split_dedupe_and_check():
    key = fr.REGISTRY["genes"]
    value, why = fr.check_value(key, ["C9orf72,brca1", "BRCA1", " c9ORF72 "])
    assert why is None and value == ["C9orf72", "brca1"]
    assert fr.check_value(key, ["BAD TOKEN"])[1] is not None
    assert fr.check_value(key, ["x" * 41])[1] is not None
    too_many = [f"G{i}" for i in range(fr.MAX_GENES + 1)]
    assert "at most" in fr.check_value(key, too_many)[1]


def test_interval_set_merges_and_matches_every_contig_spelling():
    ivs = sc.IntervalSet([("7", 100, 200), ("chr7", 150, 300), ("chr7", 301, 400)])
    assert ivs.items() == [("chr7", 100, 400)]
    assert ivs.contains("7", 100) and ivs.contains("chr7", 400)
    assert not ivs.contains("7", 99) and not ivs.contains("7", 401)
    x = sc.IntervalSet([("chrX", 10, 20)])
    assert x.contains("X", 15) and x.contains("23", 15) and not x.contains("chr7", 15)


def test_scope_sql_text_and_threshold(monkeypatch):
    scope = _scope(
        regions=["chr7:100-200", "7:500-600"], exclude_regions=["chr7:150-160"]
    )
    sql = sc.scope_sql("g", scope)
    assert sql == (
        "((g.chrom IN ('7', 'chr7') AND (g.pos BETWEEN 100 AND 200 OR "
        "g.pos BETWEEN 500 AND 600))) AND NOT ((g.chrom IN ('7', 'chr7') AND "
        "(g.pos BETWEEN 150 AND 160)))"
    )
    assert sc.scope_sql("g", scope, part="exclude") == (
        "NOT ((g.chrom IN ('7', 'chr7') AND (g.pos BETWEEN 150 AND 160)))"
    )
    assert sc.scope_sql("g", None) is None
    assert sc.scope_sql("g", _scope(regions=["chr7:1-2"]), part="exclude") is None
    monkeypatch.setattr(sc, "SCOPE_SQL_MAX_INTERVALS", 2)
    assert sc.scope_sql("g", scope) is None  # 3 intervals: Python only
    assert not scope.pushdown


@pytest.mark.parametrize("backend", ["chdb", "duckdb"])
def test_the_bucket_superset_and_its_whole_bucket_subset(backend, monkeypatch):
    """Above SCOPE_SQL_EXACT_MAX the SQL keeps every bucket an inclusion
    touches (a superset) and removes only the buckets an exclusion covers
    whole (so it never removes a row outside the exclusion)."""
    monkeypatch.setenv("VCFCLICK_BACKEND", backend)
    monkeypatch.setattr(sc, "SCOPE_SQL_EXACT_MAX", 0)
    scope = _scope(regions=["chr7:300-700"], exclude_regions=["chr7:256-1023"])
    assert scope.sql_mode == "sql_buckets" and scope.include.bucket_size() == 256
    sql = sc.scope_sql("g", scope)
    div = "intDiv(g.pos, 256)" if backend == "chdb" else "(g.pos // 256)"
    # 300-700 touches buckets 1 and 2 (256-767); 256-1023 covers 1, 2 and 3 whole
    assert f"{div} IN (1, 2))" in sql
    assert f"NOT ((g.chrom IN ('7', 'chr7') AND {div} IN (1, 2, 3)))" in sql
    partial = _scope(exclude_regions=["chr7:300-700"])
    assert f"{div} IN ()" not in sc.scope_sql("g", partial)  # no whole bucket
    assert "(1 = 0)" in sc.scope_sql("g", partial)


def test_bucket_size_grows_to_keep_the_set_small(monkeypatch):
    monkeypatch.setattr(sc, "MAX_BUCKETS", 10)
    ivs = sc.IntervalSet([("chr1", 1, 100_000)])
    assert ivs.bucket_size() == 16384  # 7 buckets; 8192 would give 13


def test_scope_none_leaves_the_model_sql_unchanged_and_a_scope_is_in_it():
    from cli.trio_analysis import Gates, Trio, _comphet_sql, trio_sql

    trio, gates = Trio("i1", "C1", "F", "M"), Gates(20, 10, 0.01, 0.25, 0.75)
    scope = _scope(regions=["chr7:100-200"], exclude_regions=["chr7:150-160"])
    for cat in ("denovo", "recessive", "dominant"):
        plain = trio_sql(cat, trio, gates, count_only=False)
        assert trio_sql(cat, trio, gates, count_only=False, scope=None) == plain
        scoped = trio_sql(cat, trio, gates, count_only=False, scope=scope)
        assert sc.scope_sql("g", scope) in scoped and "ORDER BY" in scoped
    comphet = _comphet_sql(trio, gates, scope=scope)
    assert sc.scope_sql("g", scope, part="exclude") in comphet
    assert scope.include.sql("g") not in comphet  # inclusions act after grouping


def test_filter_rows_is_the_backstop():
    scope = _scope(regions=["chr7:100-200"], exclude_regions=["7:150-150"])
    rows = [("7", 99), ("chr7", 100), ("7", 150), ("7", 200), ("chr8", 150)]
    assert sc.filter_rows(rows, scope) == [("chr7", 100), ("7", 200)]
    assert sc.drop_excluded(rows, scope) == [r for r in rows if r != ("7", 150)]


def test_compound_het_rule_after_grouping():
    """Kept: a listed gene, or a gene with a variant in a listed region.
    Exclusion alone keeps every gene it did not remove."""
    v = lambda c, p: (c, p, "A", "G", None)  # noqa: E731
    genes = {
        "IN": {"paternal": [v("chr7", 55000)], "maternal": [v("chr7", 59000)]},
        "OUT": {"paternal": [v("chr7", 72000)], "maternal": [v("chr7", 78000)]},
    }
    region = _scope(regions=["chr7:50000-56000"])
    kept, dropped = sc.genes_in_scope(genes, region)
    assert list(kept) == ["IN"] and dropped == 1
    assert region.contains("chr7", 55000) and not region.contains("chr7", 59000)
    assert sc.genes_in_scope(genes, _scope(exclude_regions=["chr1:1-2"]))[1] == 0


def test_bed_file_is_zero_based_half_open(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("track name=x\n# a comment\nchr7\t99\t200\textra\n7\t0\t1\n")
    assert file_option_tokens(fr.REGISTRY["regions"], str(bed)) == [
        "chr7:100-200",
        "7:1-1",
    ]
    bad = tmp_path / "bad.bed"
    bad.write_text("chr7\t200\t100\n")
    with pytest.raises(click.BadParameter):
        file_option_tokens(fr.REGISTRY["regions"], str(bad))


def test_gene_file_reads_one_token_per_line_with_comments(tmp_path):
    f = tmp_path / "panel.txt"
    f.write_text("# ID panel v4\nC9orf72\nC19orf12  # NBIA\n\nENSG00000111678, BRCA1\n")
    assert file_option_tokens(fr.REGISTRY["genes"], str(f)) == [
        "C9orf72",
        "C19orf12",
        "ENSG00000111678",
        "BRCA1",
    ]


# ───────────────────────────── gene resolution ─────────────────────────────


@pytest.mark.parametrize(
    "token", ["c9orf72", "C9ORF72", "C9orf72", "ENSG00000147894.3"]
)
def test_genes_resolve_case_insensitively_and_by_ensembl_id(genes_loaded, token):
    scope = _scope(genes=[token])
    assert scope.echo["genes"]["genes"] == 1
    assert scope.contains("chr7", 15000) and scope.contains("7", 10000)
    assert not scope.contains("chr7", 9999) and not scope.contains("chr7", 20001)


def test_echo_shape_is_pinned(genes_loaded):
    scope = _scope(
        genes=PANEL, genes_label="test panel v1", gene_flank=0,
        exclude_genes=["OVL1"], regions=["chr7:1-100"], exclude_regions=["chr7:5-6"],
        regions_build="GRCh38",
    )  # fmt: skip
    doc = sc.scope_document(scope, {})
    assert set(doc) == {
        "genes", "exclude_genes", "note", "regions", "exclude_regions",
        "regions_build", "intervals", "applied", "comphet_genes_dropped_out_of_scope",
        "n_outside", "notes",
    }  # fmt: skip
    assert set(doc["genes"]) == {
        "requested", "resolved", "genes", "resolved_by", "unresolved", "ambiguous",
        "sha256_12", "build", "build_verified", "aliases_loaded", "flank", "label",
        "excluded",
    }  # fmt: skip
    g = doc["genes"]
    assert (g["requested"], g["resolved"], g["genes"], g["excluded"]) == (3, 3, 3, 1)
    assert g["resolved_by"] == {"symbol": 3} and g["unresolved"] == []
    assert g["build"] == "GRCh38" and g["build_verified"] is True
    assert g["label"] == "test panel v1" and len(g["sha256_12"]) == 12
    assert doc["applied"] == "sql" and doc["intervals"] == 6  # 4 included, 2 excluded
    json.dumps(doc)  # plain JSON


def test_the_list_hash_ignores_order_and_case(genes_loaded):
    a = _scope(genes=PANEL).echo["genes"]["sha256_12"]
    b = _scope(genes=[g.upper() for g in reversed(PANEL)]).echo["genes"]["sha256_12"]
    assert a == b != _scope(genes=PANEL[:2]).echo["genes"]["sha256_12"]


def test_a_short_list_may_miss_no_gene(genes_loaded):
    with pytest.raises(
        sc.ScopeError, match="1 of 4 --genes tokens did not resolve: FAM123X"
    ):
        _scope(genes=[*PANEL, "FAM123X"])
    scope = sc.build_scope(
        sc.ScopeRequest.from_values(
            {"genes": [*PANEL, "FAM123X"]}, allow_unresolved=True
        ),
        build="GRCh38",
    )
    [w] = scope.warnings
    assert w["code"] == "genes_unresolved" and w["level"] == "warning"
    assert scope.echo["genes"]["unresolved"] == ["FAM123X"]
    assert "1 not found: FAM123X" in sc.scope_line(scope)


@pytest.mark.parametrize("n, missed, error", [(20, 1, True), (21, 1, False), (21, 2, True), (40, 2, False), (40, 3, True)])  # fmt: skip
def test_the_unresolved_threshold(n, missed, error):
    assert sc.unresolved_is_error(n, missed) is error
    assert sc.unresolved_is_error(n, 0) is False


def test_a_long_list_with_few_misses_warns(genes_loaded):
    tokens = PANEL * 7 + ["NOPE1"]  # 22 tokens (repeats count), one miss: 4.5%
    tokens = [f"{t}" for t in tokens]
    scope = sc.build_scope(sc.ScopeRequest(genes=tuple(tokens)), build="GRCh38")
    assert scope.warnings[0]["code"] == "genes_unresolved"


def test_aliases_need_the_hgnc_table(genes_loaded):
    with pytest.raises(sc.ScopeError, match="load-hgnc"):
        _scope(genes=["FTDALS1"])
    from annotations.loaders.hgnc_aliases import load

    load(HGNC)
    scope = _scope(genes=["FTDALS1", "HGNC:25443", "113246"])
    g = scope.echo["genes"]
    assert g["resolved_by"] == {"previous": 1, "hgnc_id": 1, "entrez": 1}
    assert g["aliases_loaded"] is True and g["genes"] == 3


def test_an_alias_of_two_genes_is_unresolved(genes_loaded):
    from annotations.loaders.hgnc_aliases import load

    load(HGNC)
    with pytest.raises(sc.ScopeError, match="name more than one gene"):
        _scope(genes=["SHARED1"])


def test_no_gene_table_is_an_error(isolated_annotation_db):
    with pytest.raises(sc.ScopeError, match="vcfclick annotations load"):
        _scope(genes=["C9orf72"])
    assert _scope(regions=["chr7:1-2"]) is not None  # regions need no genes


def test_a_grch37_database_cannot_take_a_gene_list(genes_loaded):
    with pytest.raises(sc.ScopeError, match="GRCh37"):
        _scope(genes=["C9orf72"], build="GRCh37")
    unknown = _scope(genes=["C9orf72"], build="unknown")
    assert unknown.echo["genes"]["build_verified"] is False
    assert any("not verified" in n for n in unknown.notes)


def test_regions_build_mismatch_warns(genes_loaded):
    scope = _scope(regions=["chr7:1-2"], regions_build="GRCh37", build="GRCh38")
    assert [w["code"] for w in scope.warnings] == ["scope_build_mismatch"]
    assert (
        _scope(regions=["chr7:1-2"], regions_build="GRCh37", build="unknown").warnings
        == ()
    )


def test_flank_widens_each_gene(genes_loaded):
    scope = _scope(genes=["C9orf72"], gene_flank=500)
    assert scope.contains("chr7", 9500) and not scope.contains("chr7", 9499)
    assert scope.contains("chr7", 20500) and scope.echo["genes"]["flank"] == 500


def test_a_par_gene_says_it_has_one_locus(genes_loaded):
    scope = _scope(genes=["XG"])
    assert scope.contains("chrX", 2_750_000) and not scope.contains("chrY", 2_750_000)
    assert any("pseudo-autosomal" in n and "chrX only" in n for n in scope.notes)


def test_two_thousand_symbols_open_one_connection(genes_loaded, monkeypatch):
    import annotations.db as adb

    opened = []
    real = adb.get_connection
    monkeypatch.setattr(adb, "get_connection", lambda: opened.append(1) or real())
    tokens = [f"NOPE{i}" for i in range(1997)] + PANEL
    scope = sc.build_scope(
        sc.ScopeRequest(genes=tuple(tokens), allow_unresolved=True), build="GRCh38"
    )
    assert len(opened) == 1
    assert scope.echo["genes"]["genes"] == 3
    assert scope.echo["genes"]["unresolved"] == {
        "count": 1997,
        "sha256_12": sc._hash([f"NOPE{i}" for i in range(1997)], upper=False),
    }


# ───────────────────────────── registry, Click, presets ─────────────────────────────


@click.command(cls=PresetCommand, name="trio")
@gate_options("trio")
@preset_option()
def _dummy(**kw):
    ctx = click.get_current_context()
    resolved, _ = resolve_filters(ctx, "trio", None)
    click.echo(json.dumps({k: resolved[k].value for k in ("genes", "regions")}))
    click.echo(json.dumps(typed_list_values(ctx, "trio")))


def test_genes_and_gene_file_typed_together_are_unioned(tmp_path):
    f = tmp_path / "p.txt"
    f.write_text("C19orf12\nc9orf72\n")
    r = CliRunner().invoke(_dummy, ["--genes", "C9orf72,BRCA1", "--gene-file", str(f), "--genes", "TP53"])  # fmt: skip
    assert r.exit_code == 0, r.output
    values = json.loads(r.output.splitlines()[0])
    assert values["genes"] == ["C9orf72", "BRCA1", "TP53", "C19orf12"]


def test_a_region_file_alone_is_typed(tmp_path):
    bed = tmp_path / "r.bed"
    bed.write_text("chr7\t29999\t40000\n")
    r = CliRunner().invoke(_dummy, ["--region-file", str(bed)])
    assert r.exit_code == 0, r.output
    assert json.loads(r.output.splitlines()[0])["regions"] == ["chr7:30000-40000"]


def test_a_preset_gene_list_is_replaced_by_a_typed_one(tmp_path, monkeypatch):
    from cli.presets import new_document, save_user_preset

    monkeypatch.setenv("VCFCLICK_HOME", str(tmp_path))
    save_user_preset(new_document("panel", {"genes": PANEL, "genes_label": "v1"}))

    @click.command(cls=PresetCommand, name="trio")
    @gate_options("trio")
    @preset_option()
    def cmd(**kw):
        from cli.preset_options import preset_state

        ctx = click.get_current_context()
        resolved, over = resolve_filters(ctx, "trio", preset_state(ctx))
        click.echo(
            json.dumps([resolved["genes"].value, resolved["genes"].source, over])
        )

    plain = json.loads(CliRunner().invoke(cmd, ["--preset", "panel"]).output)
    assert plain[:2] == [PANEL, "preset"]
    typed = json.loads(
        CliRunner().invoke(cmd, ["--preset", "panel", "--genes", "BRCA1"]).output
    )
    assert typed[:2] == [["BRCA1"], "flag"]
    assert typed[2] == [{"key": "genes", "preset": PANEL, "flag": ["BRCA1"]}]
    off = json.loads(
        CliRunner()
        .invoke(cmd, ["--preset", "panel", "--preset-ignore", "genes"])
        .output
    )
    assert off[:2] == [[], "default"] and off[2][0]["flag"] == "off"


def test_preset_save_embeds_the_gene_file_and_label(tmp_path, monkeypatch):
    from cli.main import cli

    monkeypatch.setenv("VCFCLICK_HOME", str(tmp_path))
    f = tmp_path / "panel.txt"
    f.write_text("C9orf72\nC19orf12\n")
    r = CliRunner().invoke(
        cli,
        ["preset", "save", "id-panel", "--gene-file", str(f), "--genes", "C12orf57",
         "--genes-label", "ID v4.2", "--region", "chr7:1-100", "--regions-build", "GRCh38"],
    )  # fmt: skip
    assert r.exit_code == 0, r.output
    doc = json.loads((tmp_path / "presets" / "id-panel.json").read_text())
    assert doc["filters"]["genes"] == ["C12orf57", "C9orf72", "C19orf12"]
    assert doc["filters"]["genes_label"] == "ID v4.2"
    assert doc["filters"]["regions"] == ["chr7:1-100"]
    assert doc["filters"]["regions_build"] == "GRCh38"
    assert "gene_file" not in doc["filters"]  # the symbols, never the path


def test_preset_files_check_scope_keys():
    def problems(filters):
        doc = {"schema": "vcfclick.preset", "schema_version": 1, "filters": filters}
        return [(p.code, p.key) for p in validate_document(doc).problems]

    assert problems({"genes": ["C9orf72"], "regions": ["chr7:1-2"]}) == []
    assert problems({"genes": ["BAD TOKEN"]}) == [("bad_value", "genes")]
    assert problems({"regions": ["chr7:2-1"]}) == [("bad_value", "regions")]
    assert problems({"gene_file": "/etc/x"}) == [("never_presettable", "gene_file")]
    assert problems({"allow_unresolved_genes": True}) == [
        ("never_presettable", "allow_unresolved_genes")
    ]
    assert problems({"regions_build": "hg19"}) == [("bad_value", "regions_build")]


def test_the_fingerprint_ignores_gene_order_and_case():
    def fp(filters):
        return fingerprint({"schema_version": 1, "filters": filters})

    assert fp({"genes": ["C9orf72", "BRCA1"]}) == fp({"genes": ["brca1", "C9ORF72"]})
    assert fp({"genes": ["C9orf72"]}) != fp({"genes": ["C9orf72", "BRCA1"]})
    assert fp({"regions": ["chr7:1-2", "chr7:5-6"]}) == fp(
        {"regions": ["chr7:5-6", "chr7:1-2"]}
    )


def test_relabelling_keeps_the_fingerprint_and_the_builtins_are_unchanged():
    """genes_label and regions_build filter nothing: a relabelled panel keeps
    its fingerprint (and every NAME@FINGERPRINT pin on it). Dropping them from
    the fingerprint moves no built-in's."""
    from cli.presets import builtin_preset
    from tests.test_builtin_presets import PINNED

    def fp(filters):
        return fingerprint({"schema_version": 1, "filters": filters})

    panel = {"genes": ["C9orf72"], "regions": ["chr7:1-2"]}
    assert fp(panel) == fp({**panel, "genes_label": "v1"})
    assert fp({**panel, "genes_label": "v1"}) == fp({**panel, "genes_label": "v2"})
    assert fp(panel) == fp({**panel, "regions_build": "GRCh38"})
    for name, (_, pinned) in PINNED.items():
        assert builtin_preset(name).fingerprint == pinned, name


def test_a_long_list_is_echoed_as_count_and_hash():
    from cli.preset_echo import compact

    genes = [f"G{i}" for i in range(501)]
    out = compact("genes", genes)
    assert out["count"] == 501 and len(out["sha256_12"]) == 12
    assert compact("genes", genes[:500]) == genes[:500]
    assert compact("min_gq", 20) == 20


def test_from_result_copies_a_long_list_back_from_its_preset(tmp_path, monkeypatch):
    """A 501-gene list is echoed as count and hash; `preset save --from-result`
    copies it from the preset only while that preset is unchanged."""
    from cli.preset_echo import compact
    from cli.preset_exchange import filters_from_result
    from cli.preset_schema import PresetError
    from cli.presets import find_preset, new_document, save_user_preset

    monkeypatch.setenv("VCFCLICK_HOME", str(tmp_path))
    genes = [f"G{i}" for i in range(501)]
    save_user_preset(new_document("big", {"genes": genes, "genes_label": "L"}))
    fp = find_preset("big")[0].fingerprint
    result = {
        "preset": {
            "name": "big",
            "source": "user",
            "fingerprint": fp,
            "requested": {"genes": compact("genes", genes), "genes_label": "L"},
            "resolved": {"genes": {"value": compact("genes", genes), "from": "preset"}},
        }
    }
    filters, _ = filters_from_result(json.dumps(result))
    assert filters == {"genes": genes, "genes_label": "L"}
    result["preset"]["fingerprint"] = "0" * 12
    with pytest.raises(PresetError, match="no longer has fingerprint"):
        filters_from_result(json.dumps(result))
