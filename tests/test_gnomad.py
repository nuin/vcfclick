"""Tests for the gnomAD allele-frequency annotation.

tests/fixtures/gnomad_cftr.vcf.gz holds REAL gnomAD v4.1 genomes AF /
AF_grpmax for the CFTR golden alleles (source URL in the header). The
loader/lookup are validated against those published values.
"""

from __future__ import annotations

from pathlib import Path

from annotations.db import GnomadAF

FIX = Path(__file__).parent / "fixtures"
GNOMAD_VCF = FIX / "gnomad_cftr.vcf.gz"


def test_gnomad_loader_and_lookup(tmp_path, monkeypatch):
    """Loading the real gnomAD slice and looking a variant up returns its
    published AF and popmax."""
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    from annotations import gnomad_af
    from annotations.loaders.gnomad import load

    assert load(GNOMAD_VCF, replace=True) == 5

    g = gnomad_af("chr7", 117488888, "A", "C")
    assert g is not None
    assert abs(g.af - 0.548299) < 1e-4
    assert abs(g.af_grpmax - 0.827770) < 1e-4
    assert g.popmax == g.af_grpmax  # popmax prefers the group-max AF

    rare = gnomad_af("chr7", 117480471, "A", "G")
    assert abs(rare.popmax - 0.009072) < 1e-4  # genuinely rare


def test_gnomad_lookup_absent_is_none(tmp_path, monkeypatch):
    """A locus outside the loaded slice returns None — the caller treats
    that as rare, not as AF 0."""
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    from annotations import gnomad_af
    from annotations.loaders.gnomad import load

    load(GNOMAD_VCF, replace=True)
    assert gnomad_af("chr7", 1, "A", "C") is None
    assert gnomad_af("chr7", 117488888, "A", "T") is None  # wrong alt


def test_popmax_falls_back_to_overall_af():
    """popmax uses AF_grpmax when present, else the overall AF."""
    assert GnomadAF("chr1", 1, "A", "G", 0.02, 0.05).popmax == 0.05
    assert GnomadAF("chr1", 1, "A", "G", 0.02, None).popmax == 0.02
    assert GnomadAF("chr1", 1, "A", "G", None, None).popmax is None


def _seed_gnomad(rows):
    """Fill the (env-redirected) annotation store's gnomad_af table."""
    from annotations.db import get_connection

    conn = get_connection()
    conn.executemany("INSERT INTO gnomad_af VALUES (?, ?, ?, ?, ?, ?)", rows)
    conn.close()


def test_gnomad_af_many_answers_like_gnomad_af(tmp_path, monkeypatch):
    """One batched lookup returns exactly what one `gnomad_af` call per
    allele does: present alleles with their values (a NULL AF or popmax
    stays NULL), absent alleles (wrong alt, wrong position, another
    chromosome) missing, repeated alleles once, a position given as a string
    as JSON from chDB can hold it."""
    from annotations import gnomad_af, gnomad_af_many

    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    assert gnomad_af_many([("chr1", 10, "A", "C")]) == {}  # nothing loaded
    _seed_gnomad(
        [
            ("chr1", 10, "A", "C", 0.002, 0.004),
            ("chr1", 10, "A", "G", 0.3, None),  # same locus, other allele
            ("chr1", 20, "AT", "A", None, None),
            ("chrX", 5, "G", "T", 0.5, 0.6),
        ]
    )
    probes = [
        ("chr1", 10, "A", "C"),
        ("chr1", 10, "A", "G"),
        ("chr1", 10, "A", "T"),  # wrong alt
        ("chr1", 20, "AT", "A"),
        ("chr1", 21, "AT", "A"),  # wrong position
        ("chr2", 10, "A", "C"),  # other chromosome
        ("chrX", "5", "G", "T"),  # position as a string
        ("chr1", 10, "A", "C"),  # repeated
    ]
    got = gnomad_af_many(probes)
    want = {}
    for chrom, pos, ref, alt in probes:
        g = gnomad_af(chrom, int(pos), ref, alt)
        if g is not None:
            want[(chrom, int(pos), ref, alt)] = g
    assert got == want
    assert len(got) == 4
    assert got[("chr1", 10, "A", "G")].popmax == 0.3
    assert got[("chr1", 20, "AT", "A")].popmax is None
    assert gnomad_af_many([]) == {}


def test_gnomad_af_many_leaves_a_given_connection_open(tmp_path, monkeypatch):
    """A caller that holds the store open (the JSON trio document) passes
    its connection: it is used and not closed, and the registered batch does
    not outlive the call."""
    from annotations import gnomad_af_many
    from annotations.db import get_connection

    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    _seed_gnomad([("chr1", 10, "A", "C", 0.1, 0.2)])
    conn = get_connection()
    for _ in range(2):  # a second call must not trip over the first's batch
        got = gnomad_af_many([("chr1", 10, "A", "C")], conn)
        assert got[("chr1", 10, "A", "C")].popmax == 0.2
    assert conn.execute("SELECT count(*) FROM gnomad_af").fetchone()[0] == 1
    conn.close()


def test_point_lookups_miss_on_an_impossible_position(tmp_path, monkeypatch):
    """`gnomad_af` and `clinvar_lookup` cast the position to the column's
    UINTEGER so DuckDB can skip row groups. A position that is no UINTEGER
    (the MCP tools take any int) is a miss, not a conversion error."""
    from annotations import clinvar_lookup, gnomad_af

    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    _seed_gnomad([("chr1", 10, "A", "C", 0.1, 0.2)])
    assert gnomad_af("chr1", 10, "A", "C").popmax == 0.2
    for pos in (-1, 2**32, 2**40):
        assert gnomad_af("chr1", pos, "A", "C") is None
        assert clinvar_lookup("chr1", pos, "A", "C") is None
