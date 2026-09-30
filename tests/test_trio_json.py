"""Tests for the machine-readable trio workflow the desktop apps use:
`db ped --proband/--father/--mother` (no PED file) and `db trio --format json`.

Fixture: tests/fixtures/related.vcf. C1 is the child of F and M; everything
is Mendelian except three planted de novo sites in C1.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
REL_VCF = Path(__file__).parent / "fixtures" / "related.vcf"
DE_NOVO = {151_000_000, 151_500_000, 152_000_000}


def _vc(home: Path, *args: str, ok: bool = True):
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    if ok:
        assert r.returncode == 0, f"{' '.join(args)} failed:\n{r.stderr}"
    return r


def _setup(home: Path, *, keep_reference: bool = True) -> None:
    _vc(home, "db", "create", "t")
    args = ["db", "ingest", "t", str(REL_VCF), "--ingest-id", "i1", "--serial"]
    if keep_reference:
        args.append("--keep-reference")
    _vc(home, *args)
    _vc(home, "db", "ped", "t", "--proband", "C1", "--father", "F", "--mother", "M")


def _trio(home: Path, *extra: str) -> dict:
    return json.loads(
        _vc(
            home, "db", "trio", "t", "--proband", "C1", "--format", "json", *extra
        ).stdout
    )


def test_ped_from_sample_names(vcfclick_home):
    _setup(vcfclick_home)
    rows = json.loads(
        _vc(
            vcfclick_home,
            "db",
            "query",
            "t",
            "SELECT sample_id, father_id, mother_id, sex, affected FROM pedigree ORDER BY sample_id",
            "--format",
            "JSONCompact",
        ).stdout
    )["data"]
    by_id = {r[0]: r for r in rows}
    assert by_id["C1"][1:3] == ["F", "M"]
    assert by_id["F"][3] == "male" and by_id["M"][3] == "female"
    assert by_id["C1"][4] == "affected"


def test_ped_needs_a_file_or_all_three_names(vcfclick_home):
    _vc(vcfclick_home, "db", "create", "t")
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "t",
        str(REL_VCF),
        "--ingest-id",
        "i1",
        "--serial",
    )
    r = _vc(
        vcfclick_home, "db", "ped", "t", "--proband", "C1", "--father", "F", ok=False
    )
    assert r.returncode != 0 and "--mother" in r.stderr
    r = _vc(
        vcfclick_home,
        "db",
        "ped",
        "t",
        "--proband",
        "C1",
        "--father",
        "F",
        "--mother",
        "NOPE",
        ok=False,
    )
    assert r.returncode != 0 and "NOPE" in r.stderr


def test_json_reports_every_model(vcfclick_home):
    _setup(vcfclick_home)
    res = _trio(vcfclick_home)
    assert res["trio"] == {
        "ingest_id": "i1",
        "proband": "C1",
        "father": "F",
        "mother": "M",
    }
    assert res["keep_reference"] is True
    assert set(res["models"]) == {"denovo", "recessive", "dominant", "comphet"}

    dn = res["models"]["denovo"]
    assert dn["blocked"] is False
    assert {c["pos"] for c in dn["candidates"]} == DE_NOVO
    assert dn["count"] == 3
    for c in dn["candidates"]:
        assert (c["proband_gt"], c["father_gt"], c["mother_gt"]) == (1, 0, 0)

    rec = res["models"]["recessive"]
    assert rec["count"] == len(rec["candidates"]) > 0
    assert all(
        (c["proband_gt"], c["father_gt"], c["mother_gt"]) == (2, 1, 1)
        for c in rec["candidates"]
    )

    dom = res["models"]["dominant"]
    assert all(
        c["proband_gt"] == 1 and (c["father_gt"] > 0) != (c["mother_gt"] > 0)
        for c in dom["candidates"]
    )


def test_json_gates_and_limit(vcfclick_home):
    _setup(vcfclick_home)
    res = _trio(vcfclick_home, "--category", "recessive", "--limit", "2")
    assert set(res["models"]) == {"recessive"}
    rec = res["models"]["recessive"]
    assert len(rec["candidates"]) == 2 and rec["count"] > 2 and rec["truncated"] is True
    assert res["gates"]["min_gq"] == 20 and res["gates"]["max_af"] == 0.01


def test_sparse_ingest_marks_models_that_need_reference(vcfclick_home):
    _setup(vcfclick_home, keep_reference=False)
    res = _trio(vcfclick_home)
    assert res["keep_reference"] is False
    for m in ("denovo", "dominant", "comphet"):
        assert res["models"][m]["blocked"] is True, m
    assert res["models"]["recessive"]["blocked"] is False


def test_text_output_unchanged(vcfclick_home):
    _setup(vcfclick_home)
    out = _vc(vcfclick_home, "db", "trio", "t", "--proband", "C1").stdout
    assert "trio: proband=C1 father=F mother=M" in out
    assert "denovo" in out and "recessive" in out
