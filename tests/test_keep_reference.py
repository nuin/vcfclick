"""Tests for `db ingest --keep-reference`.

Default ingest is sparse: only non-reference genotype calls are stored,
so a parent absent at a site is indistinguishable between 0/0 and ./.
--keep-reference additionally stores confident hom-reference calls
(gt=0) so trio de-novo analysis can PROVE a parent is 0/0 — while still
dropping no-calls (./.), which assert nothing.

The trio fixture (tests/fixtures/trio.vcf.gz):
  chr1:100  CHILD het, FATHER 0/0,  MOTHER 0/0   (de novo)
  chr1:200  CHILD homalt, FATHER het, MOTHER het (recessive)
  chr1:300  CHILD het, FATHER het,  MOTHER 0/0   (dominant from father)
  chr1:400  CHILD het, FATHER 0/0,  MOTHER het   (dominant from mother)
  chr1:500  CHILD het, FATHER ./.,  MOTHER 0/0   (NOT de novo: father no-call)
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
TRIO_VCF = Path(__file__).parent / "fixtures" / "trio.vcf.gz"


def _vc(home: Path, *args: str):
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    assert r.returncode == 0, (
        f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
        f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    )
    return r


def _rows(home: Path, db: str, sql: str) -> list[list]:
    r = _vc(home, "db", "query", db, sql, "--format", "JSONCompact")
    return json.loads(r.stdout)["data"]


def _ingest(home: Path, *extra: str) -> None:
    _vc(home, "db", "create", "trio")
    _vc(
        home,
        "db",
        "ingest",
        "trio",
        str(TRIO_VCF),
        "--cohort",
        "fam",
        "--ingest-id",
        "fam1",
        "--serial",
        *extra,
    )


def test_default_ingest_is_sparse_no_homref(vcfclick_home):
    """Without --keep-reference, hom-ref calls are NOT stored. At the
    de-novo site (chr1:100) only the CHILD (het) row exists."""
    _ingest(vcfclick_home)
    rows = _rows(
        vcfclick_home,
        "trio",
        "SELECT sample_id, gt FROM genotypes WHERE pos = 100 ORDER BY sample_id",
    )
    assert rows == [["CHILD", 1]], f"expected only CHILD het, got {rows}"
    # No gt=0 rows anywhere in a sparse ingest.
    zero = _rows(vcfclick_home, "trio", "SELECT count(*) FROM genotypes WHERE gt = 0")
    assert zero == [[0]], "sparse ingest must not store hom-ref (gt=0) rows"


def test_keep_reference_stores_homref_at_variant_sites(vcfclick_home):
    """With --keep-reference, the de-novo site stores all three members:
    CHILD het (1) + both parents confident hom-ref (0). This is what
    makes de novo provable."""
    _ingest(vcfclick_home, "--keep-reference")
    rows = _rows(
        vcfclick_home,
        "trio",
        "SELECT sample_id, gt FROM genotypes WHERE pos = 100 ORDER BY sample_id",
    )
    assert rows == [
        ["CHILD", 1],
        ["FATHER", 0],
        ["MOTHER", 0],
    ], f"de-novo site should have child het + both parents hom-ref: {rows}"


def test_keep_reference_still_drops_no_calls(vcfclick_home):
    """A no-call (./.) must NOT be stored even with --keep-reference —
    it asserts nothing and must stay distinguishable from a stored 0/0.
    At chr1:500 the FATHER is ./. , so only CHILD + MOTHER have rows."""
    _ingest(vcfclick_home, "--keep-reference")
    rows = _rows(
        vcfclick_home,
        "trio",
        "SELECT sample_id FROM genotypes WHERE pos = 500 ORDER BY sample_id",
    )
    assert rows == [
        ["CHILD"],
        ["MOTHER"],
    ], f"FATHER no-call must be absent, not stored: {rows}"
    # This is the crux: de novo at 500 is NOT provable (father unknown),
    # whereas at 100 it IS (both parents have a stored gt=0).
    father_100 = _rows(
        vcfclick_home,
        "trio",
        "SELECT gt FROM genotypes WHERE pos = 100 AND sample_id = 'FATHER'",
    )
    father_500 = _rows(
        vcfclick_home,
        "trio",
        "SELECT gt FROM genotypes WHERE pos = 500 AND sample_id = 'FATHER'",
    )
    assert father_100 == [[0]], "father provably hom-ref at de-novo site"
    assert father_500 == [], "father no-call at chr1:500 is unprovable (absent)"


def test_keep_reference_parallel_path(vcfclick_home):
    """The parallel ingester must honour --keep-reference too (it threads
    the flag through the worker args tuple)."""
    _vc(vcfclick_home, "db", "create", "trio")
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "trio",
        str(TRIO_VCF),
        "--cohort",
        "fam",
        "--ingest-id",
        "fam1",
        "--workers",
        "2",
        "--keep-reference",
    )
    rows = _rows(
        vcfclick_home,
        "trio",
        "SELECT sample_id, gt FROM genotypes WHERE pos = 100 ORDER BY sample_id",
    )
    assert rows == [["CHILD", 1], ["FATHER", 0], ["MOTHER", 0]]
