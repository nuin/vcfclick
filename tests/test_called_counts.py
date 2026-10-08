"""Called-genotype accounting: variants.n_called / an_called / ac_called
and the missing_genotypes table, across every ingest path, Parquet
dump/load and bundles, plus older databases and dumps.

Fixture: tests/fixtures/popgen.vcf (see conftest.popgen_vcf). The rows
used below, 11 samples A1..A4 B1..B4 C1 C2 U1:

  1:400   A1 ./.                         → n 10, an 20, ac 3
  1:500   B1 ./1 (partially missing)     → n 10, an 21, ac 3
  1:600   C1 haploid `1`                 → n 11, an 21, ac 2
  1:1100  B1 B2 B3 ./.                   → n  8, an 16, ac 2
  1:1300  U1 ./.                         → n 10, an 20, ac 6
  1:1400  phased, U1 1|1                 → n 11, an 22, ac 4
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from tests.conftest import run_cli

REPO = Path(__file__).resolve().parent.parent

EXPECTED_COUNTS = {
    400: [10, 20, 3],
    500: [10, 21, 3],
    600: [11, 21, 2],
    1100: [8, 16, 2],
    1300: [10, 20, 6],
    1400: [11, 22, 4],
}
EXPECTED_MISSING = [
    [400, "A1"],
    [1100, "B1"],
    [1100, "B2"],
    [1100, "B3"],
    [1300, "U1"],
]


def _vc(home: Path, backend: str, *args: str) -> subprocess.CompletedProcess:
    return run_cli(home, backend, *args)


def _rows(home: Path, backend: str, db: str, sql: str) -> list[list]:
    r = _vc(home, backend, "db", "query", db, sql, "--format", "JSONCompact")
    return json.loads(r.stdout)["data"]


def _columns(home, backend, db: str, table: str) -> list[str]:
    if backend == "duckdb":
        sql = (
            "SELECT column_name FROM information_schema.columns "
            f"WHERE table_name = '{table}' ORDER BY ordinal_position"
        )
    else:
        sql = (
            "SELECT name FROM system.columns WHERE database = currentDatabase() "
            f"AND table = '{table}' ORDER BY position"
        )
    return [r[0] for r in _rows(home, backend, db, sql)]


def _counts(home, backend, db="pg", ingest_id="kg") -> dict[int, list]:
    rows = _rows(
        home,
        backend,
        db,
        "SELECT pos, n_called, an_called, ac_called FROM variants "
        f"WHERE ingest_id = '{ingest_id}' AND chrom = '1' ORDER BY pos, alt",
    )
    return {r[0]: r[1:] for r in rows}


def _missing(home, backend, db="pg", ingest_id="kg") -> list[list]:
    return _rows(
        home,
        backend,
        db,
        "SELECT pos, sample_id FROM missing_genotypes "
        f"WHERE ingest_id = '{ingest_id}' ORDER BY pos, sample_id",
    )


@pytest.mark.parametrize("backend", ["duckdb", "chdb"])
@pytest.mark.parametrize("mode", [["--serial"], ["--workers", "2"]])
def test_ingest_records_called_counts_and_missing(
    vcfclick_home, popgen_vcf, backend, mode
):
    _vc(vcfclick_home, backend, "db", "create", "pg")
    _vc(
        vcfclick_home,
        backend,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "kg",
        *mode,
    )
    counts = _counts(vcfclick_home, backend)
    for pos, expected in EXPECTED_COUNTS.items():
        assert counts[pos] == expected, (pos, counts[pos])
    # Every fully called diploid site: 11 samples, 22 alleles.
    assert counts[100] == [11, 22, 8]
    assert _missing(vcfclick_home, backend) == EXPECTED_MISSING
    # `genotypes` keeps its meaning: no ./. rows, the partial ./1 is a het
    # (cyvcf2 gt_types), the haploid `1` is stored as gt=2.
    assert _rows(
        vcfclick_home,
        backend,
        "pg",
        "SELECT pos, sample_id, gt FROM genotypes WHERE pos IN (500, 600) "
        "AND sample_id IN ('B1', 'C1') ORDER BY pos",
    ) == [[500, "B1", 1], [600, "C1", 2]]


def test_no_record_missing_keeps_site_counts(vcfclick_home, popgen_vcf):
    _vc(vcfclick_home, "duckdb", "db", "create", "pg")
    _vc(
        vcfclick_home,
        "duckdb",
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "kg",
        "--serial",
        "--no-record-missing",
    )
    assert _missing(vcfclick_home, "duckdb") == []
    assert _counts(vcfclick_home, "duckdb")[400] == EXPECTED_COUNTS[400]


@pytest.mark.parametrize("backend", ["duckdb", "chdb"])
def test_dump_ingest_parquet_and_bundle_round_trip(
    vcfclick_home, popgen_vcf, tmp_path, backend
):
    home = vcfclick_home
    _vc(home, backend, "db", "create", "pg")
    _vc(
        home,
        backend,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "kg",
        "--serial",
    )
    _vc(home, backend, "db", "panel", "pg", str(REPO / "tests/fixtures/popgen.panel"))

    dump = tmp_path / "dump"
    _vc(home, backend, "db", "dump", "pg", "--out", str(dump))
    assert (dump / "missing_genotypes.parquet").exists()
    assert (dump / "populations.parquet").exists()

    _vc(home, backend, "db", "create", "copy")
    _vc(home, backend, "db", "ingest-parquet", "copy", str(dump), "--ingest-id", "kg")
    assert _counts(home, backend, "copy") == _counts(home, backend)
    assert _missing(home, backend, "copy") == EXPECTED_MISSING
    assert _rows(home, backend, "copy", "SELECT count(*) FROM populations") == [[10]]
    # Idempotent: re-loading the same dump under the same id replaces.
    _vc(home, backend, "db", "ingest-parquet", "copy", str(dump), "--ingest-id", "kg")
    assert _rows(home, backend, "copy", "SELECT count(*) FROM populations") == [[10]]
    assert _missing(home, backend, "copy") == EXPECTED_MISSING

    bundle = tmp_path / "pg.tar.gz"
    _vc(home, backend, "db", "push", "pg", str(bundle))
    _vc(home, backend, "db", "pull", "pulled", str(bundle))
    assert _counts(home, backend, "pulled") == _counts(home, backend)
    assert _missing(home, backend, "pulled") == EXPECTED_MISSING
    assert _rows(home, backend, "pulled", "SELECT count(*) FROM populations") == [[10]]


def _strip_to_old_dump(dump: Path) -> None:
    """Make `dump` look like one written before called-count tracking."""
    t = pq.read_table(dump / "variants.parquet")
    pq.write_table(
        t.drop_columns(["n_called", "an_called", "ac_called"]),
        dump / "variants.parquet",
    )
    for name in ("missing_genotypes", "populations"):
        (dump / f"{name}.parquet").unlink()


@pytest.mark.parametrize("backend", ["duckdb", "chdb"])
def test_older_dump_and_bundle_still_load(vcfclick_home, popgen_vcf, tmp_path, backend):
    home = vcfclick_home
    _vc(home, backend, "db", "create", "pg")
    _vc(
        home,
        backend,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "kg",
        "--serial",
    )
    dump = tmp_path / "dump"
    _vc(home, backend, "db", "dump", "pg", "--out", str(dump))
    _strip_to_old_dump(dump)

    _vc(home, backend, "db", "create", "old")
    _vc(home, backend, "db", "ingest-parquet", "old", str(dump), "--ingest-id", "kg")
    assert set(map(tuple, _counts(home, backend, "old").values())) == {
        (None, None, None)
    }
    assert _missing(home, backend, "old") == []

    import tarfile

    bundle = tmp_path / "old.tar.gz"
    with tarfile.open(bundle, "w:gz") as tar:
        for f in sorted(dump.iterdir()):
            tar.add(f, arcname=f.name)
    out = _vc(home, backend, "db", "pull", "oldpull", str(bundle)).stdout
    assert "no missing_genotypes.parquet in bundle" in out
    assert _counts(home, backend, "oldpull")[400] == [None, None, None]


@pytest.mark.parametrize("backend", ["duckdb", "chdb"])
def test_ingest_upgrades_an_older_database(vcfclick_home, popgen_vcf, backend):
    """A database created before these columns/tables existed gets them
    added in place on the next ingest; its old rows read as NULL."""
    home = vcfclick_home
    _vc(home, backend, "db", "create", "pg")
    _vc(
        home,
        backend,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "old",
        "--serial",
    )
    for sql in (
        "ALTER TABLE variants DROP COLUMN n_called",
        "ALTER TABLE variants DROP COLUMN an_called",
        "ALTER TABLE variants DROP COLUMN ac_called",
        "DROP TABLE missing_genotypes",
        "DROP TABLE populations",
    ):
        _vc(home, backend, "db", "query", "pg", sql)
    info = _vc(home, backend, "db", "info", "pg").stdout
    assert "missing_gt: (not present)" in info

    # Dump of the old database skips the absent tables instead of failing.
    r = _vc(home, backend, "db", "dump", "pg", "--out", str(home / "d"))
    assert not (home / "d" / "missing_genotypes.parquet").exists()

    r = _vc(
        home,
        backend,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "new",
        "--serial",
    )
    assert "upgraded schema" in r.stderr
    assert _counts(home, backend, ingest_id="new")[400] == EXPECTED_COUNTS[400]
    assert _counts(home, backend, ingest_id="old")[400] == [None, None, None]
    assert _missing(home, backend, ingest_id="new") == EXPECTED_MISSING
    # Upgraded and fresh databases agree on column order (DuckDB can only
    # append, so its DDL declares the new columns last too).
    _vc(home, backend, "db", "create", "fresh")
    for table in ("variants", "missing_genotypes", "populations"):
        assert _columns(home, backend, "pg", table) == _columns(
            home, backend, "fresh", table
        ), table
    # Upgrading is idempotent: nothing left to add on the next ingest.
    again = _vc(
        home,
        backend,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "new",
        "--serial",
    )
    assert "upgraded schema" not in again.stderr
    # Re-ingesting the old id now records its counts too.
    _vc(
        home,
        backend,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "old",
        "--serial",
    )
    assert _counts(home, backend, ingest_id="old")[400] == EXPECTED_COUNTS[400]


_SMALL_BATCH_INGEST = """
import json
import sys

import ingest.vcf_load as vl
from storage import get_session

if sys.argv[2] == "small":
    vl.MAX_BATCH_ROWS = 1  # flush after every site that has any row
vl.ingest(sys.argv[1], cohort="c", ingest_id="i1")
sess = get_session()


def n(table):
    raw = sess.query(f"SELECT count(*) FROM {table}", "CSV").bytes().decode()
    return [ln for ln in raw.splitlines() if ln.strip()][-1]


print("COUNTS " + json.dumps({t: n(t) for t in (
    "variants", "genotypes", "missing_genotypes"
)}))
"""


def test_row_count_flush_lands_the_same_rows(vcfclick_home, popgen_vcf, run_python):
    """Flushing staging batches on row count (not only every BATCH_SIZE
    sites) must not lose or duplicate rows."""
    counts = {}
    for mode in ("small", "default"):
        home = vcfclick_home / mode
        home.mkdir()
        (home / "dbs" / "rc").mkdir(parents=True)
        out = run_python(
            home,
            _SMALL_BATCH_INGEST,
            str(popgen_vcf),
            mode,
            VCFCLICK_DB_NAME="rc",
            VCFCLICK_BACKEND="duckdb",
        )
        line = [ln for ln in out.splitlines() if ln.startswith("COUNTS ")][-1]
        counts[mode] = json.loads(line.removeprefix("COUNTS "))
    assert counts["small"] == counts["default"]
    assert int(counts["default"]["missing_genotypes"]) > 0
