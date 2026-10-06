"""Shared pytest fixtures.

Every test gets its own VCFCLICK_HOME under tmp_path so nothing touches
the user's real ~/.vcfclick/dbs/. The fixture VCF is a tiny committed
file under tests/fixtures/ — 5 variants, 3 samples, bgzip+tabix indexed.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def vcfclick_home(tmp_path, monkeypatch) -> Path:
    """Isolated VCFCLICK_HOME for one test.

    Tests invoke the CLI via subprocess (see test_cli._vc), so they pick
    this env var up at process startup — no module reload needed.
    """
    home = tmp_path / "vcfclick"
    home.mkdir()
    monkeypatch.setenv("VCFCLICK_HOME", str(home))
    monkeypatch.delenv("VCFCLICK_DB_NAME", raising=False)
    return home


@pytest.fixture(autouse=True)
def _no_chdb_session_in_pytest_process():
    """Fail the test that opens a chDB session in the pytest process.

    chDB keeps one embedded server per process, bound to the first path
    it opened, so a session left behind here makes a later, unrelated
    test fail with "EmbeddedServer already initialized with path ...".
    Use the `run_python` fixture or the CLI instead.
    """
    yield
    sdb = sys.modules.get("storage.db")
    if sdb is None:
        return
    leaked = [k for k in sdb._sessions if k.startswith("chdb::")]
    for k in leaked:
        session = sdb._sessions.pop(k)
        if hasattr(session, "close"):
            session.close()
    assert not leaked, (
        f"chDB session(s) opened inside the pytest process: {leaked}. "
        "Run that work in a subprocess (the run_python fixture or the CLI)."
    )


def _run_python(home: Path, code: str, *args: str, **env: str) -> str:
    """Run `code` (argv = `args`) in a fresh interpreter against `home`;
    return its stdout.

    chDB allows one embedded server per process, bound to the first path it
    opens. A test that opens a chDB session inside the pytest process either
    inherits another test's server or leaves its own behind, so tests that
    need the storage layer in-process (to monkeypatch it, say) do that work
    in a subprocess instead.
    """
    full = {**os.environ, "VCFCLICK_HOME": str(home), **env}
    r = subprocess.run(
        [sys.executable, "-c", code, *args],
        cwd=Path(__file__).resolve().parent.parent,
        env=full,
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, (
        f"python subprocess failed (rc={r.returncode}):\n"
        f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    )
    return r.stdout


@pytest.fixture
def run_python():
    """`run_python(home, code, *args, **env) -> stdout`, in a subprocess."""
    return _run_python


@pytest.fixture
def tiny_vcf() -> Path:
    """Path to the committed 5-variant / 3-sample fixture VCF."""
    p = FIXTURES / "tiny.vcf.gz"
    assert p.exists(), f"missing fixture: {p}"
    assert (FIXTURES / "tiny.vcf.gz.tbi").exists(), "missing tabix index"
    return p


@pytest.fixture
def isolated_annotation_db(tmp_path, monkeypatch) -> Path:
    """Redirect annotations.db.DUCKDB_PATH to a temp file per test so
    the ClinVar loader / position_for_gene / clinvar_lookup tests
    cannot touch the user's real annotations/annotations.duckdb.

    get_connection() reads DUCKDB_PATH at call time, so patching the
    module attribute is enough — no reload needed.
    """
    import annotations.db as adb

    path = tmp_path / "test_annotations.duckdb"
    monkeypatch.setattr(adb, "DUCKDB_PATH", path)
    return path


def bgzip_vcf(src: Path, out: Path) -> Path:
    """bgzip + tabix a plain-text VCF fixture to `out` (skip without htslib)."""
    import shutil

    if not (shutil.which("bgzip") and shutil.which("tabix")):
        pytest.skip("bgzip/tabix not on PATH")
    with open(out, "wb") as fh:
        subprocess.run(["bgzip", "-c", str(src)], stdout=fh, check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(out)], check=True)
    return out


@pytest.fixture
def popgen_vcf(tmp_path) -> Path:
    """The population-genetics fixture (tests/fixtures/popgen.vcf), bgzipped
    and indexed: 11 samples in three panel populations plus one unlabelled
    sample, with missing (./.), partial (./1), haploid and phased calls,
    a split multi-allelic site, an indel, a filtered site and chrX."""
    return bgzip_vcf(FIXTURES / "popgen.vcf", tmp_path / "popgen.vcf.gz")


REPO = Path(__file__).resolve().parent.parent
CLI_TIMEOUT_S = 120


def run_cli(home: Path, backend: str, *args: str, ok: bool = True):
    """Run `vcfclick *args` against `home` on `backend` in a subprocess.

    A command that exceeds CLI_TIMEOUT_S is killed (subprocess.run kills
    the child, so nothing is orphaned) and reported by name. The message
    starts with "Timeout" so the suite's rerun policy for the known
    intermittent chDB subprocess hang (pyproject addopts) applies.
    """
    import shutil

    exe = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
    env = {**os.environ, "VCFCLICK_HOME": str(home), "VCFCLICK_BACKEND": backend}
    env.pop("VCFCLICK_DB_NAME", None)
    try:
        r = subprocess.run(
            [exe, *args],
            cwd=REPO,
            env=env,
            capture_output=True,
            text=True,
            timeout=CLI_TIMEOUT_S,
        )
    except subprocess.TimeoutExpired:
        pytest.fail(
            f"Timeout: `vcfclick {' '.join(args)}` ({backend}) did not finish "
            f"in {CLI_TIMEOUT_S}s and was killed"
        )
    if ok:
        assert r.returncode == 0, (
            f"`vcfclick {' '.join(args)}` ({backend}) failed (rc={r.returncode}):\n"
            f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
        )
    else:
        assert r.returncode != 0, (
            f"`vcfclick {' '.join(args)}` should fail:\n{r.stdout}"
        )
    return r
