"""The pedigree loader (cli/pedigree_load.py), on both backends.

One store per backend (trio.vcf.gz plus its PED) is built with the CLI and
read in a fresh interpreter, which also changes the store behind the
loader's back: a duplicated row (what a ReplacingMergeTree can show before
a background merge) and a database with no pedigree table at all.
"""

from __future__ import annotations

import json
import textwrap

import pytest

from tests.conftest import FIXTURES, run_cli

BACKENDS = ("chdb", "duckdb")

SCRIPT = textwrap.dedent(
    """
    import json, os
    from storage import get_session
    from cli.main import _set_db
    from cli.pedigree_load import (
        NO_PEDIGREE, NoPedigree, load_pedigree_graph, pedigree_rows,
    )
    from cli.trio_analysis import _resolve_parents

    _set_db("t")
    sess = get_session("t")
    out = {}
    rows = pedigree_rows(sess, "fam1")
    graph = load_pedigree_graph(sess, "fam1")
    out["rows"] = rows
    out["probands"] = graph.probands()
    out["kinds"] = {s: graph.kind_of(s) for s in graph.members()}
    out["other_ingest"] = pedigree_rows(sess, "nope")

    # The same sample inserted twice: the later ingested_at wins.
    ins = "INSERT INTO pedigree VALUES ('fam1', 'CHILD', 'FAM1', 'FATHER', " \\
          "'MOTHER', 'male', '%s', '%s')"
    sess.query(ins % ("unaffected", "2030-01-01 00:00:00"))
    dup = [r for r in pedigree_rows(sess, "fam1") if r["sample_id"] == "CHILD"]
    out["after_duplicate"] = dup
    sess.query(ins % ("affected", "2031-01-01 00:00:00"))
    dup = [r for r in pedigree_rows(sess, "fam1") if r["sample_id"] == "CHILD"]
    out["after_newer"] = dup

    # An empty table is "no rows"; a missing one is NoPedigree.
    from storage import delete_where_sql
    sess.query(delete_where_sql("pedigree", "1 = 1"))
    out["empty"] = pedigree_rows(sess, "fam1")
    try:
        load_pedigree_graph(sess, "fam1", require=True)
        out["empty_required"] = None
    except NoPedigree as e:
        out["empty_required"] = e.message

    # The commands say so too, in one sentence, and `db families` exits 0.
    from click.testing import CliRunner
    from cli.main import cli

    def run(*args):
        r = CliRunner().invoke(cli, list(args))
        return {"rc": r.exit_code, "out": r.stdout, "err": r.stderr}

    out["families_empty_text"] = run("db", "families", "t")
    out["families_empty_json"] = run("db", "families", "t", "--format", "json")
    sess.query("DROP TABLE pedigree")
    out["families_missing_text"] = run("db", "families", "t")
    out["families_missing_json"] = run("db", "families", "t", "--format", "json")
    out["trio_missing"] = run("db", "trio", "t", "--proband", "CHILD")
    for key, call in (
        ("missing_rows", lambda: pedigree_rows(sess, "fam1")),
        ("missing_graph", lambda: load_pedigree_graph(sess, "fam1")),
        ("missing_resolve", lambda: _resolve_parents(sess, "fam1", "CHILD")),
    ):
        try:
            call()
            out[key] = None
        except NoPedigree as e:
            out[key] = e.message
    out["constant"] = NO_PEDIGREE
    print(json.dumps(out))
    """
)


@pytest.mark.parametrize("backend", BACKENDS)
def test_loader_reads_collapses_duplicates_and_reports_no_pedigree(
    backend, vcfclick_home, run_python
):
    run_cli(vcfclick_home, backend, "db", "create", "t")
    run_cli(
        vcfclick_home, backend, "db", "ingest", "t", str(FIXTURES / "trio.vcf.gz"),
        "--ingest-id", "fam1", "--serial",
    )  # fmt: skip
    run_cli(vcfclick_home, backend, "db", "ped", "t", str(FIXTURES / "trio.ped"))
    out = json.loads(
        run_python(vcfclick_home, SCRIPT, VCFCLICK_BACKEND=backend)
        .strip()
        .splitlines()[-1]
    )

    by_id = {r["sample_id"]: r for r in out["rows"]}
    assert sorted(by_id) == ["CHILD", "FATHER", "MOTHER"]
    assert set(by_id["CHILD"]) == {
        "sample_id", "family_id", "father_id", "mother_id", "sex", "affected",
    }  # fmt: skip
    assert by_id["CHILD"]["father_id"] == "FATHER"
    assert out["probands"] == ["CHILD"]
    assert out["kinds"] == {"CHILD": "trio", "FATHER": None, "MOTHER": None}
    assert out["other_ingest"] == []

    # one row per sample: the newest insert is the one that counts
    assert [r["affected"] for r in out["after_duplicate"]] == ["unaffected"]
    assert [r["affected"] for r in out["after_newer"]] == ["affected"]

    assert out["empty"] == []
    for key in ("empty_required", "missing_rows", "missing_graph", "missing_resolve"):
        assert out[key] == out["constant"], key
    assert "no pedigree loaded" in out["constant"]

    for key in ("families_empty_text", "families_missing_text"):
        assert out[key]["rc"] == 0 and "no pedigree loaded" in out[key]["out"], key
    empty = {"ingest_id": "fam1", "families": [], "problems": []}
    for key in ("families_empty_json", "families_missing_json"):
        assert out[key]["rc"] == 0, key
        assert json.loads(out[key]["out"]) == empty
        assert "no pedigree loaded" in out[key]["err"]
    assert out["trio_missing"]["rc"] != 0
    assert "no pedigree loaded" in out["trio_missing"]["err"]
