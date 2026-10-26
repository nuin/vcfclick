"""chrX regions: the SQL predicate, the Python classifier and bcftools agree.

Pure Python plus an in-memory DuckDB (the predicate is dialect-neutral SQL);
no chDB. The bcftools comparison is skipped when bcftools is absent.
"""

from __future__ import annotations

import shutil
import subprocess

import pytest

from cli.xchrom import (
    AMBIGUOUS_BAND,
    GRCH37,
    GRCH37_CHRX_LENGTH,
    GRCH38,
    NPX_INTERVALS,
    UNKNOWN,
    X_NAMES,
    band_sql,
    is_hemizygous,
    npx_sql,
    region,
    resolve_build,
)
from tests.conftest import FIXTURES

BUILDS = (GRCH37, GRCH38, UNKNOWN)


def _edges() -> list[int]:
    out = {1, 10_000, 10_001, 60_000, 60_001}
    for table in (NPX_INTERVALS, AMBIGUOUS_BAND):
        for intervals in table.values():
            for lo, hi in intervals:
                out |= {lo - 1, lo, lo + 1, hi - 1, hi, hi + 1}
    return sorted(out)


@pytest.mark.parametrize("build", BUILDS)
def test_sql_predicate_equals_the_classifier(build):
    duckdb = pytest.importorskip("duckdb")
    con = duckdb.connect()
    positions = _edges()
    names = (*X_NAMES, "chrY", "chr1", "chrXY")
    rows = ", ".join(f"('{c}', {p})" for c in names for p in positions)
    got = con.execute(
        f"SELECT chrom, pos, {npx_sql(build)}, {band_sql(build)} "
        f"FROM (VALUES {rows}) t(chrom, pos)"
    ).fetchall()
    for chrom, pos, npx, band in got:
        want = region(chrom, pos, build)
        assert (npx, band) == (want == "npx", want == "ambiguous"), (chrom, pos)


def test_unknown_is_the_intersection_and_the_bands_its_complement():
    def points(intervals):
        return {p for p in _edges() if any(lo <= p <= hi for lo, hi in intervals)}

    n37, n38 = points(NPX_INTERVALS[GRCH37]), points(NPX_INTERVALS[GRCH38])
    # a position beyond the GRCh37 chrX length says nothing about GRCh37
    beyond = {p for p in _edges() if p > GRCH37_CHRX_LENGTH}
    assert points(NPX_INTERVALS[UNKNOWN]) == (n37 | beyond) & n38
    assert points(AMBIGUOUS_BAND[UNKNOWN]) == (n37 ^ n38) - beyond


def test_build_resolution_threshold():
    assert resolve_build("auto", GRCH37_CHRX_LENGTH) == (UNKNOWN, "unknown", False)
    assert resolve_build("auto", GRCH37_CHRX_LENGTH + 1) == (GRCH38, "data", False)
    assert resolve_build("auto", None) == (UNKNOWN, "unknown", False)
    assert resolve_build(GRCH37, GRCH37_CHRX_LENGTH + 1) == (GRCH37, "option", True)
    assert resolve_build(GRCH38, 1) == (GRCH38, "option", False)
    with pytest.raises(ValueError):
        resolve_build("CHM13", 1)


def test_hemizygous_only_for_a_male_on_npx():
    assert is_hemizygous("male", "chrX", 50_000_000, GRCH38)
    assert not is_hemizygous("female", "chrX", 50_000_000, GRCH38)
    assert not is_hemizygous("male", "chrX", 1_000_000, GRCH38)  # PAR1
    assert not is_hemizygous("male", "chr1", 50_000_000, GRCH38)
    assert not is_hemizygous(None, "X", 50_000_000, GRCH38)


def _bcftools_hemizygous(build: str) -> list[tuple[int, int]]:
    r = subprocess.run(
        ["bcftools", "+mendelian2", str(FIXTURES / "xlinked.vcf"),
         "-p", "1X:XSON,XDAD,XMUM", "--rules", f"{build}?"],
        capture_output=True, text=True,
    )  # fmt: skip
    out = []
    for line in (r.stdout + r.stderr).splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0] == "1X" and parts[1].startswith("chrX:"):
            lo, hi = parts[1].split(":")[1].split("-")
            out.append((int(lo), int(hi)))
    return out


@pytest.mark.parametrize(
    "build,par1_start,expected_diff",
    [
        # GRC puts PAR2 at 155,701,383; bcftools ends NPX one base earlier
        (GRCH38, 10_001, {(155_701_382, 155_701_382)}),
        # bcftools has no rule past GRCh37 PAR2 (155,260,561-155,270,560)
        (GRCH37, 60_001, {(155_260_561, 155_270_560)}),
    ],
)
def test_intervals_match_bcftools(build, par1_start, expected_diff):
    if not shutil.which("bcftools"):
        pytest.skip("bcftools not on PATH")
    theirs = [(lo, hi) for lo, hi in _bcftools_hemizygous(build) if lo >= par1_start]
    assert theirs, "no chrX haploid rule parsed from bcftools"
    ours = NPX_INTERVALS[build]
    # every bcftools NPX interval starts where ours does
    assert {lo for lo, _ in theirs} <= {lo for lo, _ in ours}
    covered = {(lo, hi) for lo, hi in theirs}
    extra = set()
    for lo, hi in ours:
        inside = [t for t in covered if t[0] == lo]
        if not inside:
            extra.add((lo, hi))
        elif inside[0][1] != hi:
            extra.add((inside[0][1] + 1, hi))
    assert extra == expected_diff
