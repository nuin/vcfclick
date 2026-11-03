"""Parallel ingest lands exactly the rows the serial ingest lands.

Each test builds a VCF whose shape once made the parallel path drop or
duplicate rows while exiting 0, then compares `db ingest` (parallel)
against `db ingest --serial` on the active backend:

  * 150 SNVs 20 kb apart (3 Mb, one BGZF block): every 16 kb tabix
    window had a byte delta of 0 and was skipped, so the planned region
    was the last 100 kb bucket alone and 5 of 150 variants landed.
  * 20,000 SNVs at POS 1000..200990: the last tabix window
    (196608..212991) was credited to the 100 kb bucket of its first
    base only, so the 99 variants past 200000 were never read.
  * The same 20,000 SNVs plus 9 bp deletions at POS 99998 and 199998:
    a region query returns every record overlapping the region, so a
    deletion crossing a region boundary landed in two workers' output.
  * A header-only VCF re-ingested under an id that held data: the
    commit INSERTed from a glob that matched no files, after the prior
    rows had already been rolled back.

The staged-count invariant (staged variants must equal the record count
from the tabix index or the cyvcf2 pre-pass) is exercised by patching the
staging step to lose one row, and the staging-directory rules by a
failing ingest and a caller-supplied staging_dir.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from tests.conftest import bgzip_vcf, run_cli

TINY_VCF = Path(__file__).parent / "fixtures" / "tiny.vcf.gz"

HEADER = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=248956422>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
]


def _backend() -> str:
    env = os.environ.get("VCFCLICK_BACKEND")
    if env:
        return env
    from storage import backend

    return backend()


def _write_vcf(path: Path, records: list[tuple[int, str, str]]) -> Path:
    """Write chr1 records (pos, ref, alt) for three samples, bgzip + tabix."""
    plain = path.with_suffix("")
    assert plain.suffix == ".vcf", path
    lines = list(HEADER)
    for pos, ref, alt in records:
        lines.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t0/1\t1/1\t0/0")
    plain.write_text("\n".join(lines) + "\n")
    return bgzip_vcf(plain, path)


def _sparse_records() -> list[tuple[int, str, str]]:
    return [(1000 + i * 20_000, "A", "G") for i in range(150)]


def _tail_records() -> list[tuple[int, str, str]]:
    return [(1000 + i * 10, "A", "G") for i in range(20_000)]


def _scalar(home: Path, backend: str, db: str, sql: str) -> str:
    r = run_cli(home, backend, "db", "query", db, sql, "--format", "TSV")
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    return lines[-1].strip()


def _ingest(home: Path, backend: str, db: str, vcf: Path, *extra: str) -> None:
    run_cli(home, backend, "db", "create", db)
    run_cli(
        home,
        backend,
        "db",
        "ingest",
        db,
        str(vcf),
        "--cohort",
        "A",
        "--ingest-id",
        "t1",
        *extra,
    )


def _counts(home: Path, backend: str, db: str) -> dict[str, str]:
    return {
        "variants": _scalar(home, backend, db, "SELECT count(*) FROM variants"),
        "distinct": _scalar(
            home,
            backend,
            db,
            "SELECT count(*) FROM (SELECT DISTINCT chrom, pos, ref, alt FROM variants)",
        ),
        "genotypes": _scalar(home, backend, db, "SELECT count(*) FROM genotypes"),
        "n_variants": _scalar(
            home, backend, db, "SELECT sum(n_variants) FROM ingestions"
        ),
        "max_pos": _scalar(home, backend, db, "SELECT max(pos) FROM variants"),
    }


@pytest.mark.parametrize(
    "shape, records, expected_variants",
    [
        ("sparse", _sparse_records(), 150),
        ("tail", _tail_records(), 20_000),
        (
            "boundary_deletion",
            sorted(
                _tail_records()
                + [(99_998, "ACGTACGTA", "A"), (199_998, "ACGTACGTA", "A")]
            ),
            20_002,
        ),
    ],
    ids=["sparse_multi_mb", "contig_tail_window", "boundary_crossing_deletion"],
)
def test_parallel_matches_serial(
    tmp_path, vcfclick_home, shape, records, expected_variants
):
    backend = _backend()
    vcf = _write_vcf(tmp_path / f"{shape}.vcf.gz", records)

    _ingest(vcfclick_home, backend, "serial", vcf, "--serial")
    _ingest(vcfclick_home, backend, "par", vcf, "--workers", "4")

    serial = _counts(vcfclick_home, backend, "serial")
    parallel = _counts(vcfclick_home, backend, "par")
    assert serial["variants"] == str(expected_variants), serial
    assert serial["distinct"] == str(expected_variants), serial
    assert parallel == serial, {"serial": serial, "parallel": parallel}


def test_parallel_reingest_of_empty_vcf_records_zero_variants(tmp_path, vcfclick_home):
    backend = _backend()
    full = _write_vcf(tmp_path / "full.vcf.gz", _sparse_records())
    empty = _write_vcf(tmp_path / "empty.vcf.gz", [])

    _ingest(vcfclick_home, backend, "par", full)
    assert (
        _scalar(vcfclick_home, backend, "par", "SELECT count(*) FROM variants") == "150"
    )

    run_cli(
        vcfclick_home,
        backend,
        "db",
        "ingest",
        "par",
        str(empty),
        "--cohort",
        "A",
        "--ingest-id",
        "t1",
    )
    assert (
        _scalar(vcfclick_home, backend, "par", "SELECT count(*) FROM variants") == "0"
    )
    assert (
        _scalar(vcfclick_home, backend, "par", "SELECT count(*) FROM ingestions") == "1"
    )
    assert (
        _scalar(vcfclick_home, backend, "par", "SELECT n_variants FROM ingestions")
        == "0"
    )
    assert (
        _scalar(
            vcfclick_home,
            backend,
            "par",
            "SELECT count(DISTINCT sample_id) FROM samples",
        )
        == "3"
    )
    staging = vcfclick_home / "dbs" / "par" / "staging"
    assert not staging.exists() or not any(staging.iterdir()), list(staging.iterdir())


def test_failed_parallel_ingest_removes_its_staging_dir(tmp_path, vcfclick_home):
    backend = _backend()
    plain = tmp_path / "multi.vcf"
    lines = list(HEADER) + [
        "chr1\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\t0/0",
        "chr1\t2000\t.\tA\tG,T\t50\tPASS\t.\tGT\t0/1\t1/2\t0/0",
    ]
    plain.write_text("\n".join(lines) + "\n")
    vcf = bgzip_vcf(plain, tmp_path / "multi.vcf.gz")

    run_cli(vcfclick_home, backend, "db", "create", "par")
    r = run_cli(
        vcfclick_home,
        backend,
        "db",
        "ingest",
        "par",
        str(vcf),
        "--ingest-id",
        "t1",
        ok=False,
    )
    assert "Multi-allelic" in r.stderr + r.stdout
    staging = vcfclick_home / "dbs" / "par" / "staging"
    assert not staging.exists() or not any(staging.iterdir()), list(staging.iterdir())


_LOST_ROW_INGEST = """
import json
import sys
from pathlib import Path

import ingest.parallel as par
from ingest.parallel import ingest_parallel
from storage import apply_schema, db_path, get_session

db_path("par").mkdir(parents=True, exist_ok=True)
apply_schema()
ingest_parallel(sys.argv[1], cohort="A", ingest_id="t1", workers=2)

orig = par._stage_regions


def lossy(*a, **k):
    total, elapsed = orig(*a, **k)
    return total - 1, elapsed


par._stage_regions = lossy
try:
    ingest_parallel(sys.argv[2], cohort="A", ingest_id="t1", workers=2)
    print("RESULT " + json.dumps({"raised": None}))
except RuntimeError as exc:
    sess = get_session("par")
    raw = sess.query("SELECT count(*) FROM variants", "CSV").bytes().decode().strip()
    n = [ln for ln in raw.splitlines() if ln.strip()][-1]
    staging = db_path("par") / "staging"
    left = sorted(p.name for p in staging.iterdir()) if staging.exists() else []
    print("RESULT " + json.dumps({"raised": str(exc), "variants": n, "staging": left}))
"""


def test_staged_count_must_match_index_record_count(
    tmp_path, vcfclick_home, run_python
):
    """One staged row short of the index's record count fails the ingest
    before the commit, so the rows from the earlier ingest under the same
    id survive and no staging directory is left behind."""
    sparse = _write_vcf(tmp_path / "sparse.vcf.gz", _sparse_records())
    out = run_python(
        vcfclick_home,
        _LOST_ROW_INGEST,
        str(TINY_VCF),
        str(sparse),
        VCFCLICK_DB_NAME="par",
        VCFCLICK_BACKEND=_backend(),
    )
    line = [ln for ln in out.splitlines() if ln.startswith("RESULT ")][-1]
    result = json.loads(line.removeprefix("RESULT "))
    assert result["raised"], result
    assert "149" in result["raised"] and "150" in result["raised"], result["raised"]
    assert "--serial" in result["raised"], result["raised"]
    assert result["variants"] == "5", result
    assert result["staging"] == [], result


_USER_STAGING_INGEST = """
import json
import sys
from pathlib import Path

from ingest.parallel import ingest_parallel
from storage import apply_schema, db_path, get_session

db_path("par").mkdir(parents=True, exist_ok=True)
apply_schema()
user_dir = Path(sys.argv[2])
ingest_parallel(sys.argv[1], cohort="A", ingest_id="t1", workers=2, staging_dir=str(user_dir))

sess = get_session("par")
raw = sess.query("SELECT count(*) FROM variants", "CSV").bytes().decode().strip()
n = [ln for ln in raw.splitlines() if ln.strip()][-1]
print("RESULT " + json.dumps({
    "variants": n,
    "left": sorted(p.name for p in user_dir.iterdir()),
}))
"""


def test_user_staging_dir_is_neither_imported_nor_deleted(
    tmp_path, vcfclick_home, run_python
):
    """A caller-supplied staging_dir keeps its own files: a stale
    `variants_*.parquet` in it is not globbed into the import (it is not
    a Parquet file, so importing it would fail) and the directory itself
    survives the post-commit cleanup."""
    user_dir = tmp_path / "userdir"
    user_dir.mkdir()
    (user_dir / "precious.txt").write_text("keep me\n")
    (user_dir / "variants_stale_0000.parquet").write_text("not parquet\n")

    out = run_python(
        vcfclick_home,
        _USER_STAGING_INGEST,
        str(TINY_VCF),
        str(user_dir),
        VCFCLICK_DB_NAME="par",
        VCFCLICK_BACKEND=_backend(),
    )
    line = [ln for ln in out.splitlines() if ln.startswith("RESULT ")][-1]
    result = json.loads(line.removeprefix("RESULT "))
    assert result["variants"] == "5", result
    assert result["left"] == ["precious.txt", "variants_stale_0000.parquet"], result
