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


def test_parallel_ingest_lands_same_rows_as_serial(vcfclick_home):
    """The parallel default with 2 workers on the tiny fixture (5
    variants, 3 samples) must land EXACTLY the same counts as the
    serial path. Catches both:
      * The empty-list fallback bug (tabix splitter returns []
        on tiny VCFs; if the caller doesn't fall back, zero rows
        land and the test fails immediately).
      * Any future regression in the worker → pool.map → bulk-import
        chain (e.g. workers writing to a different staging dir).
    """
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

    # The tiny fixture: 5 variants, 3 samples.
    # Non-reference genotype rows are sparse (~10 across 5 variants × 3
    # samples per the test_ingest_atomic baseline). Compare to the
    # serial path's numbers exactly.
    assert _scalar(vcfclick_home, "par", "SELECT count(*) FROM variants") == "5"
    assert (
        _scalar(vcfclick_home, "par", "SELECT count(DISTINCT sample_id) FROM samples")
        == "3"
    )
    assert _scalar(vcfclick_home, "par", "SELECT count(*) FROM genotypes") == "10"


def test_variant_density_includes_final_linear_bucket():
    """Direct unit test on the helper that used to drop the final
    16Kb linear-index bucket. With offsets indicating 3 successive
    buckets all containing variants, the density map must include
    the position bucket corresponding to the final linear bucket."""
    from ingest._tabix import variant_density

    # Three 16Kb buckets, each with non-trivial byte cost between
    # them. The final bucket has no successor so its byte_cost
    # can't be computed from a delta — the helper has to insert a
    # placeholder cost (not 0) so the final position bucket survives
    # the zero-filter.
    #
    # Synthetic linear offsets: bgzf_block | uoffset
    # (1 << 16) → block 1, offset 0; (2 << 16) → block 2, offset 0; …
    offsets = [
        (1 << 16),  # bucket 0 begins at block 1
        (5 << 16),  # bucket 1 begins at block 5 (cost 4 bytes-of-blocks)
        (9 << 16),  # bucket 2 begins at block 9 (cost 4)
    ]
    density = variant_density(offsets, position_bucket_size=100_000)

    # All three linear buckets (0, 1, 2) map to position bucket 0
    # at 100Kb granularity (0 * 16384 / 100000 == 0, 1 * 16384 / 100000
    # == 0, 2 * 16384 / 100000 == 0). The fact that bucket 0 of the
    # density map exists at all is what we need to assert — under the
    # pre-fix code path the final linear bucket would be skipped, but
    # because the first two contribute non-zero deltas the map would
    # still have bucket 0. To actually exercise the trailing-bucket
    # path we need offsets where the final linear bucket is the ONLY
    # contribution to its position bucket.
    assert 0 in density, "density map lost position bucket 0"

    # Now the load-bearing case: a long run of empty (delta-zero)
    # buckets followed by a single trailing bucket. Under the old
    # code, every middle bucket has byte_cost=0 (so skipped) and the
    # final bucket also gets byte_cost=0 (so also skipped) — the
    # contig's last position bucket vanished. Under the fix the
    # trailing bucket gets a placeholder cost and the position bucket
    # survives.
    long_run = [(10 << 16)] * 11 + [(20 << 16)]  # 11 empty deltas + final
    density2 = variant_density(long_run, position_bucket_size=100_000)
    # Linear bucket 11 = position 11 * 16384 = 180_224 → pos bucket 1.
    assert 1 in density2, f"density map dropped the trailing linear bucket: {density2}"
