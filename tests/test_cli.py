"""End-to-end smoke tests for the `vcfclick db` CLI.

These exercise the full lifecycle a real user goes through:
create → ingest → query → info → dump → push → pull → rm.

Each command is invoked as a subprocess rather than via Click's
CliRunner. The reason is chDB: its embedded ClickHouse server can only
initialise once per process, so a multi-DB test (push from one DB,
pull into another) cannot run in-process. Subprocesses also more
faithfully mirror how users actually invoke the CLI.

Each test gets its own VCFCLICK_HOME under tmp_path (see conftest), so
nothing touches the user's real ~/.vcfclick/dbs/.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tarfile
from pathlib import Path

import pyarrow.parquet as pq


REPO = Path(__file__).resolve().parent.parent
# Prefer the installed entry-point script (handles the importlib shim in
# cli/main.py correctly — `python -m cli.main` would run the file as
# `__main__`, which is a *different* module from `cli.main`, so the
# subcommand registrations in cli/db.py and cli/annotations.py would
# attach to the wrong group object).
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")

# Expected counts in the committed fixture VCF.
FIXTURE_VARIANTS = 5
FIXTURE_SAMPLES = 3
# Genotypes table only stores rows where GT differs from 0/0 (sparse).
# The fixture has, per variant: (S1,S2,S3) = (0/0,0/1,1/1), (0/0,0/0,0/1),
# (0/1,0/1,0/1), (1/1,0/1,1/1), (0/0,0/0,0/1) → 2+1+3+3+1 = 10 non-ref calls.
FIXTURE_GENOTYPES = 10


def _vc(home: Path, *args: str) -> str:
    """Run `vcfclick <args>` against the isolated home. Return stdout."""
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, (
        f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
        f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    )
    return r.stdout


def _ingest(home: Path, db: str, vcf: Path, ingest_id: str = "batch_a") -> None:
    _vc(
        home,
        "db",
        "ingest",
        db,
        str(vcf),
        "--cohort",
        "demo",
        "--ingest-id",
        ingest_id,
        "--serial",
    )


def test_create_list_rm(vcfclick_home):
    out = _vc(vcfclick_home, "db", "list")
    assert "no databases yet" in out.lower()

    out = _vc(vcfclick_home, "db", "create", "smoke")
    assert "smoke" in out

    out = _vc(vcfclick_home, "db", "list")
    assert "smoke" in out

    _vc(vcfclick_home, "db", "rm", "smoke", "--yes")

    out = _vc(vcfclick_home, "db", "list")
    assert "smoke" not in out


def test_ingest_query_info(vcfclick_home, tiny_vcf):
    _vc(vcfclick_home, "db", "create", "smoke")
    _ingest(vcfclick_home, "smoke", tiny_vcf)

    out = _vc(vcfclick_home, "db", "query", "smoke", "SELECT count() FROM variants")
    assert str(FIXTURE_VARIANTS) in out

    out = _vc(vcfclick_home, "db", "query", "smoke", "SELECT count() FROM samples")
    assert str(FIXTURE_SAMPLES) in out

    out = _vc(vcfclick_home, "db", "query", "smoke", "SELECT count() FROM genotypes")
    assert str(FIXTURE_GENOTYPES) in out

    out = _vc(vcfclick_home, "db", "info", "smoke")
    assert f"variants:  {FIXTURE_VARIANTS}" in out
    assert f"samples:   {FIXTURE_SAMPLES}" in out
    assert f"genotypes: {FIXTURE_GENOTYPES}" in out


def test_idempotent_reingest(vcfclick_home, tiny_vcf):
    """Same ingest_id twice must not double-count rows (ReplacingMergeTree)."""
    _vc(vcfclick_home, "db", "create", "smoke")
    _ingest(vcfclick_home, "smoke", tiny_vcf)
    _ingest(vcfclick_home, "smoke", tiny_vcf)
    out = _vc(
        vcfclick_home, "db", "query", "smoke", "SELECT count() FROM variants FINAL"
    )
    assert str(FIXTURE_VARIANTS) in out


def test_dump_produces_parquet(vcfclick_home, tiny_vcf, tmp_path):
    _vc(vcfclick_home, "db", "create", "smoke")
    _ingest(vcfclick_home, "smoke", tiny_vcf)

    out_dir = tmp_path / "dump"
    _vc(vcfclick_home, "db", "dump", "smoke", "--out", str(out_dir))

    expected = {
        "variants.parquet",
        "genotypes.parquet",
        "samples.parquet",
        "ingestions.parquet",
    }
    actual = {p.name for p in out_dir.iterdir()}
    assert expected.issubset(actual), f"missing files: {expected - actual}"
    for p in out_dir.iterdir():
        assert p.stat().st_size > 0, f"empty parquet: {p}"


def test_push_pull_roundtrip(vcfclick_home, tiny_vcf, tmp_path):
    _vc(vcfclick_home, "db", "create", "src")
    _ingest(vcfclick_home, "src", tiny_vcf)

    bundle = tmp_path / "src.tar.gz"
    _vc(vcfclick_home, "db", "push", "src", str(bundle))
    assert bundle.exists() and bundle.stat().st_size > 0

    _vc(vcfclick_home, "db", "pull", "dst", str(bundle))

    src_out = _vc(vcfclick_home, "db", "query", "src", "SELECT count() FROM variants")
    dst_out = _vc(vcfclick_home, "db", "query", "dst", "SELECT count() FROM variants")
    # Both should report the same variant count.
    assert str(FIXTURE_VARIANTS) in src_out
    assert str(FIXTURE_VARIANTS) in dst_out


def test_pull_accepts_bundle_with_old_variants_schema(vcfclick_home, tiny_vcf, tmp_path):
    """Older demo bundles lack newer nullable/defaulted variants columns."""
    _vc(vcfclick_home, "db", "create", "src")
    _ingest(vcfclick_home, "src", tiny_vcf)

    bundle = tmp_path / "src.tar.gz"
    _vc(vcfclick_home, "db", "push", "src", str(bundle))

    old_bundle = tmp_path / "old-src.tar.gz"
    with tarfile.open(bundle, "r:gz") as tar:
        tar.extractall(tmp_path / "bundle", filter="data")

    variants_path = tmp_path / "bundle" / "variants.parquet"
    table = pq.read_table(variants_path)
    old_table = table.drop(
        [
            "info_FractionInformativeReads",
            "info_HAPCOMP",
            "info_HAPDOM",
            "info_DragenSnvHardQUAL",
            "info_DragenIndelHardQUAL",
        ]
    )
    pq.write_table(old_table, variants_path)

    with tarfile.open(old_bundle, "w:gz") as tar:
        for p in sorted((tmp_path / "bundle").iterdir()):
            tar.add(p, arcname=p.name)

    _vc(vcfclick_home, "db", "pull", "dst", str(old_bundle))

    out = _vc(vcfclick_home, "db", "query", "dst", "SELECT count() FROM variants")
    assert str(FIXTURE_VARIANTS) in out


def test_query_genomic_region(vcfclick_home, tiny_vcf):
    """Region predicate hits chr1 fixture rows. Variants at 100, 250, 500, 750, 900."""
    _vc(vcfclick_home, "db", "create", "smoke")
    _ingest(vcfclick_home, "smoke", tiny_vcf)

    out = _vc(
        vcfclick_home,
        "db",
        "query",
        "smoke",
        "SELECT count() FROM variants WHERE chrom='chr1' AND pos BETWEEN 200 AND 800",
    )
    # 250, 500, 750 → three matches.
    assert "3" in out
