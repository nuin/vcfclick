"""Parallel-ingest regression tests.

Two silent failures the parallel path used to have, locked in here so
future refactors of the splitter or the dispatch loop have to keep them
fixed:

  * The tabix-derived splitter returns `[]` (empty list) for small VCFs
    whose linear index is too sparse to balance across workers. The
    caller used to check `if regions is None`, so an empty list silently
    became "zero workers, zero rows ingested." Both backends were
    affected. Fixed by switching to `if not regions:` so empty lists
    fall through to the cyvcf2 pre-pass splitter too.

  * `variant_density()` in `ingest/_tabix.py` used to set the FINAL
    16Kb linear-index bucket's byte_cost to 0, then filter out
    zero-cost buckets, which dropped the very last bucket of every
    contig. On chr21 of 1000G phase 3 that silently lost 134 variants
    at the distal tip. Both backends affected. Fixed by giving the
    trailing bucket a placeholder cost of 1 so it survives the filter.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
TINY_VCF = Path(__file__).parent / "fixtures" / "tiny.vcf.gz"


def _vc(home: Path, *args: str) -> subprocess.CompletedProcess:
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
    return r


def _scalar(home: Path, db: str, sql: str) -> str:
    return _vc(home, "db", "query", db, sql + " FORMAT TabSeparated").stdout.strip()


def test_parallel_ingest_falls_back_when_tabix_returns_empty(
    vcfclick_home, monkeypatch
):
    """Lock the empty-list fallback specifically. Without forcing the
    tabix splitter to return [], we can't reliably observe the
    regression — after the variant_density fix, the tabix path on the
    tiny fixture returns a non-empty list. To prove the `if not regions`
    branch is the load-bearing change, monkeypatch split_via_tbi to
    return [] and require the ingest to still succeed via the cyvcf2
    pre-pass fallback.

    Done in-process (not via the CLI subprocess) because monkeypatches
    do not cross the process boundary. `ingest_parallel()` lazy-imports
    `split_via_tbi` from `ingest._tabix` inside the function body, so
    the patch has to land on the source module.
    """
    monkeypatch.setenv("VCFCLICK_DB_NAME", "par")
    monkeypatch.setenv("VCFCLICK_BACKEND", os.environ.get("VCFCLICK_BACKEND", "chdb"))

    # Both backends enforce a single live session per (process, DB).
    # Clear any session cache state leaked from earlier tests in the
    # same pytest process so the chdb EmbeddedServer can re-init under
    # the tmp_path VCFCLICK_HOME this test was handed.
    import storage.db as sdb

    sdb._sessions.clear()

    from storage import apply_schema, db_path, get_session

    db_path("par").mkdir(parents=True, exist_ok=True)
    apply_schema()

    import ingest._tabix as tabix_mod
    from ingest.parallel import ingest_parallel

    monkeypatch.setattr(tabix_mod, "split_via_tbi", lambda *a, **k: [])

    # If the fallback is broken (regions == [] is not None-checked),
    # this call surfaces the failure either as an exception during the
    # bulk-import phase ("no files match the pattern variants_*.parquet")
    # or as zero rows landed when the workers loop got no work to do.
    # Either way the assertions below catch it.
    ingest_parallel(
        str(TINY_VCF),
        cohort="A",
        ingest_id="t1",
        workers=2,
    )

    # Query through the live in-process session — both backends enforce
    # one-writer-per-DB-file, so spawning a CLI subprocess to query
    # would collide with this process's session.
    sess = get_session("par")
    raw = sess.query("SELECT count(*) FROM variants", "CSV").bytes().decode().strip()
    last = [ln for ln in raw.splitlines() if ln.strip()]
    assert last and last[-1] == "5", f"variants count {last!r}"

    raw = sess.query("SELECT count(*) FROM genotypes", "CSV").bytes().decode().strip()
    last = [ln for ln in raw.splitlines() if ln.strip()]
    assert last and last[-1] == "10", f"genotypes count {last!r}"

    raw = (
        sess.query("SELECT count(DISTINCT sample_id) FROM samples", "CSV")
        .bytes()
        .decode()
        .strip()
    )
    last = [ln for ln in raw.splitlines() if ln.strip()]
    assert last and last[-1] == "3", f"samples count {last!r}"


def test_parallel_ingest_lands_same_rows_as_serial(vcfclick_home):
    """Smoke test: parallel default with 2 workers on the tiny fixture
    must land the same counts as the serial path. Covers the end-to-end
    worker → pool.map → bulk-import chain on both backends. Does NOT by
    itself lock the empty-list fallback — see the monkeypatched test
    above for that — because after the variant_density fix the tiny
    fixture's tabix splitter no longer returns []."""
    _vc(vcfclick_home, "db", "create", "par")
    _vc(
        vcfclick_home,
        "db",
        "ingest",
        "par",
        str(TINY_VCF),
        "--cohort",
        "A",
        "--ingest-id",
        "t1",
        "--workers",
        "2",
    )

    assert _scalar(vcfclick_home, "par", "SELECT count(*) FROM variants") == "5"
    assert (
        _scalar(vcfclick_home, "par", "SELECT count(DISTINCT sample_id) FROM samples")
        == "3"
    )
    assert _scalar(vcfclick_home, "par", "SELECT count(*) FROM genotypes") == "10"


def test_variant_density_includes_final_linear_bucket():
    """Direct unit test on the helper that used to drop the final
    16Kb linear-index bucket — the trailing-bucket regression that
    silently lost 134 variants on 1000G chr21.

    The bug-exposing shape: a run of empty (delta-zero) buckets,
    then ONE non-zero delta whose linear-index position maps to a
    different rolled-up position bucket than the trailing bucket.
    Under the old code the trailing bucket gets byte_cost=0 and is
    dropped by the zero-filter, so its position bucket vanishes.
    Under the fix the trailing bucket gets a placeholder cost of 1
    and its position bucket survives.
    """
    from ingest._tabix import variant_density

    # Seven identical entries (delta=0 from index 0 through 5), then a
    # jump at index 6, then the trailing entry at index 7.
    #
    #   linear_idx 0..5: delta = 0  → skipped under old code
    #   linear_idx 6:    delta = 10 → counted, pos_bucket = 6 * 16384 // 100_000 = 0
    #   linear_idx 7:    trailing   → pos_bucket = 7 * 16384 // 100_000 = 1
    #
    # Under the old code, only linear_idx 6 contributes (pos bucket 0).
    # Pos bucket 1 has no contributor and is absent from the density map.
    # Under the fix, linear_idx 7's placeholder cost of 1 makes pos
    # bucket 1 appear.
    offsets = [(10 << 16)] * 7 + [(20 << 16), (20 << 16)]
    assert len(offsets) == 9
    density = variant_density(offsets, position_bucket_size=100_000)

    # Sanity: pos bucket 0 was already present pre-fix (from linear_idx 6).
    assert 0 in density, f"pre-existing bucket vanished: {density}"
    # Load-bearing: pos bucket 1 only appears because the trailing
    # linear bucket gets a placeholder cost.
    assert (
        1 in density
    ), f"trailing-bucket fix regressed — pos bucket 1 missing: {density}"
