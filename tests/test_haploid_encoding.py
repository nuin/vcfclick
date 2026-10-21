"""How haploid and partial calls are stored, pinned on both backends.

schema/02_genotypes.sql, schema/duckdb/02_genotypes.sql, docs/SCHEMA.md and
the MCP briefing all say the same thing: the encoder never writes a negative
`gt`; a haploid ALT call (a male on the non-PAR X, chrM) is stored as `gt = 2`
with `phased = 1`, a haploid REF call as `gt = 0` (only with
--keep-reference), a partial `./1` as `gt = 1`, and a full no-call is not
stored but recorded in `missing_genotypes`. Those comments once described a
`gt = -1` "hemizygous" code that nothing writes. This test is the evidence
the corrected wording rests on, run through the real ingest path.
"""

from __future__ import annotations

import json

import pytest

from tests.conftest import run_cli
from tests.fixtures._vcfgen import Site, write_vcf

BACKENDS = ("chdb", "duckdb")

SAMPLES = ["HAP_ALT", "HAP_REF", "HAP_MISSING", "PARTIAL", "HOM_ALT", "HET"]
# one chrX site, one GT per sample, in SAMPLES order
GTS = ("1", "0", ".", "./1", "1/1", "0/1")


def _rows(home, backend, sql: str) -> list:
    r = run_cli(home, backend, "db", "query", "t", sql, "--format", "JSONCompact")
    return json.loads(r.stdout)["data"]


@pytest.mark.parametrize("backend", BACKENDS)
def test_haploid_and_partial_calls_are_stored_as_documented(
    backend, vcfclick_home, tmp_path
):
    vcf = write_vcf(
        tmp_path / "hap.vcf", SAMPLES, [Site("chrX", 1000000, "A", "G", GTS)],
        contigs=("chrX",),
    )  # fmt: skip
    run_cli(vcfclick_home, backend, "db", "create", "t")
    run_cli(
        vcfclick_home, backend, "db", "ingest", "t", str(vcf),
        "--ingest-id", "i1", "--serial", "--keep-reference",
    )  # fmt: skip

    rows = _rows(
        vcfclick_home,
        backend,
        "SELECT sample_id, gt, phased FROM genotypes ORDER BY sample_id",
    )
    stored = {sample: (int(gt), int(phased)) for sample, gt, phased in rows}
    assert stored == {
        "HAP_ALT": (2, 1),  # haploid ALT: stored as 1/1, reads as phased
        "HAP_REF": (0, 1),  # haploid REF: a reference row (keep-reference)
        "PARTIAL": (1, 0),  # ./1: one ALT allele called
        "HOM_ALT": (2, 0),
        "HET": (1, 0),
        # HAP_MISSING (.) has no genotypes row
    }
    assert min(gt for gt, _ in stored.values()) >= 0, "the encoder never writes < 0"

    missing = _rows(vcfclick_home, backend, "SELECT sample_id FROM missing_genotypes")
    assert [m[0] for m in missing] == ["HAP_MISSING"]

    # called-allele counts are ploidy-aware: HAP_ALT 1 + HAP_REF 1 + PARTIAL 1
    # + HOM_ALT 2 + HET 2 = 7 alleles called; ALT alleles 1 + 0 + 1 + 2 + 1 = 5
    an_called, ac_called = _rows(
        vcfclick_home, backend, "SELECT an_called, ac_called FROM variants"
    )[0]
    assert (int(an_called), int(ac_called)) == (7, 5)


def test_the_schema_files_describe_that_encoding():
    """The DDL comments and the docs must not mention a negative code again."""
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent
    for rel in (
        "schema/02_genotypes.sql",
        "schema/duckdb/02_genotypes.sql",
        "docs/SCHEMA.md",
    ):
        text = (root / rel).read_text()
        assert "-1 = hemizygous" not in text, rel
        assert "`-1` = hemizygous" not in text, rel
        assert "haploid" in text.lower(), rel
        assert "never writes a negative" in text, rel
