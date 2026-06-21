"""Parquet ingest tests.

The contract: `db dump` produces three Parquet files matching the
locked Arrow schemas in ingest/_arrow.py; `db ingest-parquet` reads
the same files back under a NEW (cohort, ingest_id) label and lands
the data byte-equivalent into the destination DB. That's the
"vcfclick is a citizen of the Arrow/Parquet stack" claim — tested
end-to-end here.

The trickier guarantees we lock in:
  * The source ingest_id and cohort columns are NOT honoured. The
    caller's --ingest-id and --cohort win, even if the source dump
    had different values. This is what makes round-tripping safe:
    moving a dump between cohorts/labels is what users will do.
  * Schema-mismatched input is rejected during Phase 1, before any
    chDB writes. A bad parquet file at re-ingest time must leave
    the prior ingest_id's data fully intact.
  * Empty dumps (no samples.parquet, no genotypes.parquet) are valid
    — that's the "cohort-summary AF table" use case where a tool
    produces variants without per-sample genotypes.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
TINY_VCF = Path(__file__).parent / "fixtures" / "tiny.vcf.gz"


def _vc(
    home: Path, *args: str, expect_failure: bool = False
) -> subprocess.CompletedProcess:
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
    if expect_failure:
        assert r.returncode != 0, f"expected failure but got rc=0:\n{r.stdout}"
    else:
        assert r.returncode == 0, (
            f"`vcfclick {' '.join(args)}` failed (rc={r.returncode}):\n"
            f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
        )
    return r


def _query_json(home: Path, db: str, sql: str) -> list[list]:
    r = _vc(home, "db", "query", db, sql, "--format", "JSONCompact")
    return json.loads(r.stdout)["data"]


def _ingest_tiny(home: Path, db: str, cohort: str, ingest_id: str) -> None:
    """Set up: create a db and load the tiny fixture VCF into it."""
    _vc(home, "db", "create", db)
    _vc(
        home,
        "db",
        "ingest",
        db,
        str(TINY_VCF),
        "--cohort",
        cohort,
        "--ingest-id",
        ingest_id,
        "--serial",
    )


# ─────────────────────── round-trip ───────────────────────


def test_dump_then_ingest_parquet_preserves_variant_and_sample_counts(vcfclick_home):
    """Foundational round-trip: VCF → dump → ingest-parquet. Variant
    and sample counts in the destination DB equal the source DB."""
    _ingest_tiny(vcfclick_home, "src", cohort="A", ingest_id="original")

    src_variants = _query_json(vcfclick_home, "src", "SELECT count() FROM variants")[0][
        0
    ]
    src_samples = _query_json(
        vcfclick_home, "src", "SELECT count(DISTINCT sample_id) FROM samples"
    )[0][0]

    dump_dir = vcfclick_home / "dump"
    _vc(vcfclick_home, "db", "dump", "src", "--out", str(dump_dir))
    # Sanity: the three core files actually got written.
    assert (dump_dir / "variants.parquet").exists()
    assert (dump_dir / "genotypes.parquet").exists()
    assert (dump_dir / "samples.parquet").exists()

    _vc(vcfclick_home, "db", "create", "dst")
    _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(dump_dir),
        "--cohort",
        "B",
        "--ingest-id",
        "imported",
    )

    dst_variants = _query_json(vcfclick_home, "dst", "SELECT count() FROM variants")[0][
        0
    ]
    dst_samples = _query_json(
        vcfclick_home, "dst", "SELECT count(DISTINCT sample_id) FROM samples"
    )[0][0]

    assert dst_variants == src_variants, (
        f"variant count drifted across round-trip: src={src_variants} dst={dst_variants}"
    )
    assert dst_samples == src_samples, (
        f"sample count drifted across round-trip: src={src_samples} dst={dst_samples}"
    )


def test_ingest_parquet_overrides_ingest_id_and_cohort(vcfclick_home):
    """The caller's --ingest-id and --cohort win. The source dump's
    label is NOT preserved — that's the contract that lets users move
    data between cohorts/runs without rewriting their parquet files."""
    _ingest_tiny(vcfclick_home, "src", cohort="ORIGINAL_COHORT", ingest_id="src_id")

    dump_dir = vcfclick_home / "dump"
    _vc(vcfclick_home, "db", "dump", "src", "--out", str(dump_dir))

    _vc(vcfclick_home, "db", "create", "dst")
    _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(dump_dir),
        "--cohort",
        "NEW_COHORT",
        "--ingest-id",
        "new_id",
    )

    # The destination DB must NOT contain the source labels anywhere.
    rows = _query_json(vcfclick_home, "dst", "SELECT DISTINCT ingest_id FROM variants")
    assert rows == [["new_id"]], f"variants ingest_id not overridden: {rows}"

    rows = _query_json(vcfclick_home, "dst", "SELECT DISTINCT ingest_id FROM genotypes")
    assert rows == [["new_id"]], f"genotypes ingest_id not overridden: {rows}"

    rows = _query_json(
        vcfclick_home,
        "dst",
        "SELECT DISTINCT ingest_id, cohort FROM samples",
    )
    assert rows == [["new_id", "NEW_COHORT"]], (
        f"samples ingest_id/cohort not overridden: {rows}"
    )

    # Belt-and-suspenders: explicitly check the source labels aren't there.
    rows = _query_json(
        vcfclick_home,
        "dst",
        "SELECT count() FROM variants WHERE ingest_id = 'src_id'",
    )
    assert rows == [[0]], "source ingest_id leaked into destination"
    rows = _query_json(
        vcfclick_home,
        "dst",
        "SELECT count() FROM samples WHERE cohort = 'ORIGINAL_COHORT'",
    )
    assert rows == [[0]], "source cohort leaked into destination"


def test_ingest_parquet_reingest_same_id_replaces_prior_data(vcfclick_home):
    """Re-ingesting under the same --ingest-id replaces the prior data
    under that id, matching the VCF-path rollback semantics. After
    two ingest-parquet calls under the same id, the destination row
    counts match a single ingest's, not double."""
    _ingest_tiny(vcfclick_home, "src", cohort="A", ingest_id="x")
    dump_dir = vcfclick_home / "dump"
    _vc(vcfclick_home, "db", "dump", "src", "--out", str(dump_dir))

    _vc(vcfclick_home, "db", "create", "dst")
    _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(dump_dir),
        "--cohort",
        "B",
        "--ingest-id",
        "stable",
    )
    after_first = _query_json(vcfclick_home, "dst", "SELECT count() FROM variants")[0][
        0
    ]

    _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(dump_dir),
        "--cohort",
        "B",
        "--ingest-id",
        "stable",
    )
    after_second = _query_json(vcfclick_home, "dst", "SELECT count() FROM variants")[0][
        0
    ]

    assert after_first == after_second, (
        f"re-ingest doubled rows instead of replacing: "
        f"first={after_first} second={after_second}"
    )


def test_ingest_parquet_writes_ingestions_catalog_row(vcfclick_home):
    """Every ingest must leave a row in the ingestions catalog so
    `db info` can show provenance. The vcf_path column for a Parquet
    ingest gets a `parquet://` prefix so the source type is obvious."""
    _ingest_tiny(vcfclick_home, "src", cohort="A", ingest_id="x")
    dump_dir = vcfclick_home / "dump"
    _vc(vcfclick_home, "db", "dump", "src", "--out", str(dump_dir))

    _vc(vcfclick_home, "db", "create", "dst")
    _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(dump_dir),
        "--cohort",
        "B",
        "--ingest-id",
        "imported",
    )

    rows = _query_json(
        vcfclick_home,
        "dst",
        "SELECT ingest_id, cohort, vcf_path FROM ingestions",
    )
    assert len(rows) == 1, f"expected one ingestions row, got {rows}"
    assert rows[0][0] == "imported"
    assert rows[0][1] == "B"
    assert rows[0][2].startswith("parquet://"), (
        f"expected parquet:// provenance prefix, got {rows[0][2]!r}"
    )


# ─────────────────────── partial-dump variants ───────────────────────


def test_ingest_parquet_works_without_samples_or_genotypes_file(
    vcfclick_home, tmp_path
):
    """A variants-only dump (cohort-summary AF table) is valid: an
    external tool that produces only allele-frequency variants without
    per-sample data should still land. No samples, no genotypes — just
    the variants table populated."""
    _ingest_tiny(vcfclick_home, "src", cohort="A", ingest_id="x")
    full_dump = vcfclick_home / "full_dump"
    _vc(vcfclick_home, "db", "dump", "src", "--out", str(full_dump))

    # Shave it down to just variants.parquet.
    variants_only = tmp_path / "variants_only"
    variants_only.mkdir()
    shutil.copy(full_dump / "variants.parquet", variants_only / "variants.parquet")

    _vc(vcfclick_home, "db", "create", "dst")
    _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(variants_only),
        "--cohort",
        "B",
        "--ingest-id",
        "var_only",
    )

    n_variants = _query_json(vcfclick_home, "dst", "SELECT count() FROM variants")[0][0]
    n_genotypes = _query_json(vcfclick_home, "dst", "SELECT count() FROM genotypes")[0][
        0
    ]
    n_samples = _query_json(vcfclick_home, "dst", "SELECT count() FROM samples")[0][0]

    assert n_variants > 0, "variants didn't land"
    assert n_genotypes == 0, f"genotypes table should be empty, got {n_genotypes}"
    assert n_samples == 0, f"samples table should be empty, got {n_samples}"


def test_ingest_parquet_derives_samples_from_genotypes_when_samples_file_absent(
    vcfclick_home, tmp_path
):
    """If the user supplies genotypes.parquet but not samples.parquet,
    the sample list is derived from `SELECT DISTINCT sample_id` against
    the genotypes file. This is the realistic case for external tools
    that produce genotype tables but don't bother with a samples table."""
    _ingest_tiny(vcfclick_home, "src", cohort="A", ingest_id="x")
    full_dump = vcfclick_home / "full_dump"
    _vc(vcfclick_home, "db", "dump", "src", "--out", str(full_dump))

    partial = tmp_path / "no_samples"
    partial.mkdir()
    shutil.copy(full_dump / "variants.parquet", partial / "variants.parquet")
    shutil.copy(full_dump / "genotypes.parquet", partial / "genotypes.parquet")
    # NOTE: deliberately no samples.parquet

    _vc(vcfclick_home, "db", "create", "dst")
    _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(partial),
        "--cohort",
        "DERIVED",
        "--ingest-id",
        "derived_samples",
    )

    src_samples = _query_json(
        vcfclick_home, "src", "SELECT count(DISTINCT sample_id) FROM samples"
    )[0][0]
    dst_samples = _query_json(
        vcfclick_home, "dst", "SELECT count(DISTINCT sample_id) FROM samples"
    )[0][0]
    assert dst_samples == src_samples, (
        f"derived sample count diverged from source: "
        f"src={src_samples} dst={dst_samples}"
    )
    # And the derived rows must carry the caller's cohort, not the source's.
    rows = _query_json(vcfclick_home, "dst", "SELECT DISTINCT cohort FROM samples")
    assert rows == [["DERIVED"]], f"derived samples got wrong cohort: {rows}"


# ─────────────────────── validation + safety ───────────────────────


def test_ingest_parquet_missing_variants_file_raises(vcfclick_home, tmp_path):
    """variants.parquet is the one required file. Without it the
    command must fail loudly rather than silently no-op."""
    empty = tmp_path / "empty_dump"
    empty.mkdir()
    _vc(vcfclick_home, "db", "create", "dst")
    r = _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(empty),
        "--cohort",
        "X",
        "--ingest-id",
        "y",
        expect_failure=True,
    )
    assert "variants.parquet" in (r.stderr + r.stdout), (
        f"expected error to mention variants.parquet, got:\nSTDOUT:{r.stdout}\nSTDERR:{r.stderr}"
    )


def test_ingest_parquet_rejects_schema_mismatch_before_touching_chdb(
    vcfclick_home, tmp_path
):
    """A parquet file with the wrong column set must be rejected in
    Phase 1, before any chDB write happens. Concretely: if the
    destination already has data under the same ingest_id from a prior
    good ingest, a bad re-ingest must leave that prior data intact."""
    # First, set up a destination DB with a known-good ingest under
    # ingest_id 'stable' so we have something to NOT-wipe.
    _ingest_tiny(vcfclick_home, "src", cohort="A", ingest_id="x")
    good_dump = vcfclick_home / "good_dump"
    _vc(vcfclick_home, "db", "dump", "src", "--out", str(good_dump))

    _vc(vcfclick_home, "db", "create", "dst")
    _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(good_dump),
        "--cohort",
        "B",
        "--ingest-id",
        "stable",
    )
    before = _query_json(vcfclick_home, "dst", "SELECT count() FROM variants")[0][0]
    assert before > 0, "set-up failed — no variants in dst"

    # Now construct a bad dump: a variants.parquet with the wrong
    # column set entirely. Should be rejected by _validate_parquet_schema.
    bad_dump = tmp_path / "bad_dump"
    bad_dump.mkdir()
    bogus_table = pa.table({"definitely_not_a_real_column": [1, 2, 3]})
    pq.write_table(bogus_table, bad_dump / "variants.parquet")

    r = _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(bad_dump),
        "--cohort",
        "B",
        "--ingest-id",
        "stable",  # SAME id as the prior good ingest — must not wipe it
        expect_failure=True,
    )
    assert (
        "schema" in (r.stderr + r.stdout).lower()
        or "missing columns" in (r.stderr + r.stdout).lower()
    ), f"expected schema-mismatch error, got:\nSTDOUT:{r.stdout}\nSTDERR:{r.stderr}"

    # The critical invariant: prior 'stable' data must still be there
    # because Phase 1 validation rejected the bad input BEFORE the
    # rollback_ingest in Phase 2 ran.
    after = _query_json(vcfclick_home, "dst", "SELECT count() FROM variants")[0][0]
    assert after == before, (
        f"bad re-ingest wiped prior data under same ingest_id: "
        f"before={before} after={after}"
    )


def test_ingest_parquet_validates_ingest_id_format(vcfclick_home, tmp_path):
    """ingest_id is interpolated into chDB SQL (in WHERE clauses for
    rollback, in INSERT SELECT for the override). validate_ingest_id
    must reject anything that isn't ASCII letters/digits/_/./-."""
    _ingest_tiny(vcfclick_home, "src", cohort="A", ingest_id="x")
    dump_dir = vcfclick_home / "dump"
    _vc(vcfclick_home, "db", "dump", "src", "--out", str(dump_dir))

    _vc(vcfclick_home, "db", "create", "dst")
    r = _vc(
        vcfclick_home,
        "db",
        "ingest-parquet",
        "dst",
        str(dump_dir),
        "--cohort",
        "B",
        "--ingest-id",
        "bad; DROP TABLE variants; --",
        expect_failure=True,
    )
    assert "ingest_id" in (r.stderr + r.stdout).lower(), (
        f"expected an ingest_id validation error, got:\nSTDOUT:{r.stdout}\nSTDERR:{r.stderr}"
    )
