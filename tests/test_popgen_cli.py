"""`vcfclick db popgen` end to end, on chDB and DuckDB.

Fixture: tests/fixtures/popgen.vcf + popgen.panel (YRI = A1..A4,
CEU = B1..B4, CHB = C1 C2; U1 is in the VCF but not the panel). With
default filters the retained sites are 1:100 200 300 1000 1200 1300 1400:

  dropped  X:100 (non-autosomal), 1:700 ×2 (split multi-allelic),
           1:800 (indel), 1:900 (LowQual), 1:500 (./1) and 1:600
           (haploid) — per-group counts not exact, 1:400 (YRI call
           rate 3/4) and 1:1100 (CEU call rate 1/4).

Every subcommand runs on both backends and must produce byte-identical
JSON. Values are checked three ways: (1) an independent reference
computed here directly from the VCF text with exact fractions, (2)
literals from scikit-allel 1.3.13, (3) literals from dadi 2.4.4 (for the
projected SFS / Tajima's D with missing data). (2) and (3) were run once
in a throwaway venv outside the project; they are not dependencies.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from tests.conftest import FIXTURES, bgzip_vcf, run_cli

BACKENDS = ("duckdb", "chdb")

POPS = {
    "YRI": ["A1", "A2", "A3", "A4"],
    "CEU": ["B1", "B2", "B3", "B4"],
    "CHB": ["C1", "C2"],
}
ALL = {"all": ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "U1"]}
RETAINED = [100, 200, 300, 1000, 1200, 1300, 1400]
RETAINED_ALL = [100, 200, 300, 400, 1000, 1200, 1300, 1400]


# ───────────────────────── independent reference ────────────────────────


def _vcf_rows() -> dict[int, dict[str, str]]:
    lines = (FIXTURES / "popgen.vcf").read_text().splitlines()
    header = next(ln for ln in lines if ln.startswith("#CHROM")).split("\t")
    out = {}
    for ln in lines:
        if ln.startswith("#"):
            continue
        f = ln.split("\t")
        if f[0] == "1" and f[3] in ("A", "C", "G", "T") and len(f[4]) == 1:
            out.setdefault(int(f[1]), dict(zip(header[9:], f[9:], strict=True)))
    return out


def _site(gts: dict[str, str], samples: list[str]) -> tuple[int, int, int, int]:
    """(n haplotypes, k alt, het individuals, called individuals)."""
    n = k = het = called = 0
    for s in samples:
        alleles = gts[s].replace("|", "/").split("/")
        if "." in alleles:
            continue
        called += 1
        n += len(alleles)
        k += alleles.count("1")
        het += len(set(alleles)) > 1
    return n, k, het, called


def _a(n: int) -> Fraction:
    return sum((Fraction(1, i) for i in range(1, n)), Fraction(0))


def reference(groups: dict[str, list[str]], positions: list[int]) -> dict[str, dict]:
    rows = _vcf_rows()
    out = {}
    for g, samples in groups.items():
        pi = theta_w = ho = he = Fraction(0)
        seg = 0
        for p in positions:
            n, k, het, called = _site(rows[p], samples)
            site_pi = Fraction(2 * k * (n - k), n * (n - 1))
            pi += site_pi
            if 0 < k < n:
                seg += 1
                theta_w += 1 / _a(n)
            ho += Fraction(het, called)
            he += site_pi
        L = len(positions)
        out[g] = {
            "sites": L,
            "segregating_sites": seg,
            "pi": float(pi),
            "theta_w": float(theta_w),
            "ho": float(ho / L),
            "he": float(he / L),
            "f": None if he == 0 else float(1 - ho / he),
        }
    return out


def reference_fst(a: list[str], b: list[str], positions: list[int]) -> float:
    rows = _vcf_rows()
    num = den = Fraction(0)
    for p in positions:
        n1, k1, *_ = _site(rows[p], a)
        n2, k2, *_ = _site(rows[p], b)
        p1, p2 = Fraction(k1, n1), Fraction(k2, n2)
        num += (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        den += p1 * (1 - p2) + p2 * (1 - p1)
    return float(num / den)


# ─────────────────────────────── harness ────────────────────────────────


def _vc(home: Path, backend: str, *args: str, ok: bool = True):
    return run_cli(home, backend, *args, ok=ok)


@pytest.fixture(scope="module")
def homes(tmp_path_factory) -> dict[str, Path]:
    """One populated database per backend, shared by the module."""
    root = tmp_path_factory.mktemp("popgen")
    vcf = bgzip_vcf(FIXTURES / "popgen.vcf", root / "popgen.vcf.gz")
    out = {}
    for b in BACKENDS:
        home = root / b
        home.mkdir()
        _vc(home, b, "db", "create", "pg")
        _vc(
            home,
            b,
            "db",
            "ingest",
            "pg",
            str(vcf),
            "--ingest-id",
            "kg",
            "--workers",
            "2",
        )
        _vc(home, b, "db", "panel", "pg", str(FIXTURES / "popgen.panel"))
        out[b] = home
    return out


def popgen(homes, *args: str) -> dict:
    """Run on both backends, assert identical output, return parsed JSON."""
    outs = {
        b: _vc(home, b, "db", "popgen", *args, "--format", "json").stdout
        for b, home in homes.items()
    }
    assert outs["duckdb"] == outs["chdb"], "backends disagree"
    return json.loads(outs["duckdb"])


def by_group(doc: dict) -> dict[str, dict]:
    return {r["group"]: r for r in doc["results"]}


# ─────────────────────────────── summary ────────────────────────────────


def test_summary_default_filters_and_report(homes):
    doc = popgen(homes, "summary", "pg")
    assert doc["grouping"] == "population"
    assert doc["groups"] == ["CEU", "CHB", "YRI"]
    assert doc["ancestral"] == "aa"  # INFO/AA present → default aa
    assert doc["missing_data_tracked"] is True
    assert doc["sites"] == {
        "in_scope": 16,
        "retained": 7,
        "dropped": {
            "non_autosomal": 1,
            "not_biallelic": 2,
            "variant_type": 1,
            "filter": 1,
            "inexact_group_counts": 2,
            "call_rate": 2,
            "maf": 0,
        },
        "polarised": 6,
        "unpolarised": 1,
        "samples": 11,
        "unlabelled_samples": 1,
    }
    assert any("inexact" in w or "exactly" in w for w in doc["warnings"])


def test_summary_matches_independent_reference(homes):
    got = by_group(popgen(homes, "summary", "pg"))
    want = reference(POPS, RETAINED)
    for g, w in want.items():
        for key, value in w.items():
            assert got[g][key] == pytest.approx(value, rel=1e-10, abs=1e-12), (g, key)
        assert got[g]["n_samples"] == len(POPS[g])
        assert got[g]["projection_n"] == 2 * len(POPS[g])
        assert got[g]["pi_per_site"] == pytest.approx(w["pi"] / 7)


def test_summary_matches_scikit_allel(homes):
    """scikit-allel 1.3.13 on the same seven sites (no missing data in any
    population there, so its fixed-n estimators coincide with ours)."""
    got = by_group(popgen(homes, "summary", "pg"))
    allel = {
        "CEU": (1.75, 1.54269972452, 0.5862256612704615, 0.25),
        "YRI": (1.92857142857, 1.92837465565, 0.0004627241059868953, 0.285714285714),
        "CHB": (1.0, 1.09090909091, None, 0.142857142857),
    }
    for g, (pi, tw, d, ho) in allel.items():
        assert got[g]["pi"] == pytest.approx(pi, rel=1e-10)
        assert got[g]["theta_w"] == pytest.approx(tw, rel=1e-10)
        assert got[g]["ho"] == pytest.approx(ho, rel=1e-10)
        if d is not None:  # allel returns nan below 3 segregating sites
            assert got[g]["tajima_d"] == pytest.approx(d, rel=1e-9)


def test_summary_whole_cohort_with_missing_data(homes):
    """--by all keeps 1:400 and 1:1300 (one missing call each among 11
    samples), so π uses per-site n and D the SFS projected to 20."""
    doc = popgen(homes, "summary", "pg", "--by", "all")
    assert doc["grouping"] == "all"
    assert doc["sites"]["retained"] == 8
    row = doc["results"][0]
    want = reference(ALL, RETAINED_ALL)["all"]
    for key, value in want.items():
        assert row[key] == pytest.approx(value, rel=1e-10), key
    assert row["projection_n"] == 20
    # scikit-allel mean_pairwise_difference (per-site n) and dadi 2.4.4
    # Spectrum.from_data_dict(..., projections=[20]).Tajima_D().
    assert row["pi"] == pytest.approx(2.40749601276, rel=1e-10)
    assert row["tajima_d"] == pytest.approx(1.34972003299, rel=1e-9)


def test_polarisation_modes(homes):
    def sites(*extra):
        return popgen(homes, "summary", "pg", *extra)["sites"]

    assert (sites()["polarised"], sites()["unpolarised"]) == (6, 1)
    high = sites("--ancestral", "aa-high")  # 1:200 (t), 1:1400 (a) are lower case
    assert (high["polarised"], high["unpolarised"]) == (4, 3)
    assert sites("--ancestral", "ref")["polarised"] == 7
    none = popgen(homes, "summary", "pg", "--ancestral", "none")
    assert none["sites"]["polarised"] == 0
    assert all(r["fay_wu_h"] is None for r in none["results"])
    # Folded statistics do not depend on polarisation.
    base = by_group(popgen(homes, "summary", "pg"))
    for r in none["results"]:
        assert r["tajima_d"] == base[r["group"]]["tajima_d"]


def test_site_filters(homes):
    def report(*extra):
        return popgen(homes, "summary", "pg", *extra)["sites"]

    assert report("--include-indels")["retained"] == 8  # + 1:800 AT>A
    assert report("--all-filters")["retained"] == 8  # + 1:900 LowQual
    # cohort MAF 0 at 1:1000 (all ref) and 1:1200 (all alt)
    maf = report("--maf", "0.1")
    assert maf["retained"] == 5 and maf["dropped"]["maf"] == 2
    low = report("--min-call-rate", "0.7")  # keeps 1:400 (YRI 3/4)
    assert low["retained"] == 8 and low["dropped"]["call_rate"] == 1
    region = popgen(homes, "summary", "pg", "--region", "chr1:1-500")
    assert region["sites"]["in_scope"] == 5  # bare "1" matched via chr1
    assert region["sites"]["retained"] == 3


def test_super_population_grouping(homes):
    sup = by_group(popgen(homes, "summary", "pg", "--by", "super_population"))
    pop = by_group(popgen(homes, "summary", "pg"))
    for s, p in (("AFR", "YRI"), ("EUR", "CEU"), ("EAS", "CHB")):
        assert {k: v for k, v in sup[s].items() if k != "group"} == {
            k: v for k, v in pop[p].items() if k != "group"
        }


# ───────────────────────────────── sfs ──────────────────────────────────


def test_sfs_populations(homes):
    got = by_group(popgen(homes, "sfs", "pg"))
    ceu = got["CEU"]
    assert ceu["projection_n"] == 8 and ceu["sites_used"] == 7
    # CEU ALT counts over the 7 sites: 1, 2, 5, 0, 8, 3, 0
    assert ceu["folded"] == [3, 1, 1, 2, 0]
    # derived counts (1:300 unpolarised; 1:200 and 1:1200 flipped):
    # 1, 8-2=6, -, 0, 8-8=0, 3, 0
    assert ceu["unfolded"] == [3, 1, 0, 1, 0, 0, 1, 0, 0]
    assert ceu["polarised_sites_used"] == 6


def test_sfs_projection_matches_dadi(homes):
    row = popgen(homes, "sfs", "pg", "--by", "all")["results"][0]
    # dadi 2.4.4 Spectrum.from_data_dict(dd, ["all"], [20], polarized=...)
    assert row["folded"] == pytest.approx(
        [2.0, 0.0, 0.025974025974, 1.311688311688, 0.727272727273,
         0.415584415584, 1.761904761905, 0.969696969697, 0.787878787879,
         0.0, 0.0],
        abs=1e-11,
    )  # fmt: skip
    assert row["unfolded"] == pytest.approx(
        [2.0, 0.0, 0.025974025974, 1.311688311688, 0.662337662338, 0.0,
         1.242424242424, 0.969696969697, 0.787878787879] + [0.0] * 12,
        abs=1e-11,
    )  # fmt: skip


def test_sfs_explicit_projection_drops_sites(homes):
    got = by_group(
        popgen(homes, "sfs", "pg", "--min-call-rate", "0.7", "--project", "8")
    )
    # 1:400 has 6 YRI haplotypes (< 8) → dropped from YRI's spectrum only
    assert got["YRI"]["sites_used"] == 7 and got["YRI"]["sites_dropped"] == 1
    assert got["CEU"]["sites_dropped"] == 0
    default = by_group(popgen(homes, "sfs", "pg", "--min-call-rate", "0.7"))
    assert default["YRI"]["projection_n"] == 6
    assert default["YRI"]["sites_used"] == 8


# ───────────────────────────────── fst ──────────────────────────────────


def test_fst_matches_reference_and_scikit_allel(homes):
    doc = popgen(homes, "fst", "pg")
    pairs = {(p["group1"], p["group2"]): p for p in doc["pairs"]}
    allel = {
        ("CEU", "CHB"): 0.45,
        ("CEU", "YRI"): 0.254972875226,
        ("CHB", "YRI"): -0.0649350649351,
    }
    for (a, b), value in allel.items():
        assert pairs[(a, b)]["fst"] == pytest.approx(value, rel=1e-10)
        assert pairs[(a, b)]["fst"] == pytest.approx(
            reference_fst(POPS[a], POPS[b], RETAINED), rel=1e-12
        )
        assert pairs[(a, b)]["sites"] == 7


def test_fst_windows(homes):
    doc = popgen(homes, "fst", "pg", "--window", "500")
    assert [(w["start"], w["end"], w["sites"]) for w in doc["windows"]] == [
        (1, 500, 3),
        (501, 1000, 1),
        (1001, 1400, 3),
    ]
    first = doc["windows"][0]
    assert first["fst_CEU_YRI"] == pytest.approx(
        reference_fst(POPS["CEU"], POPS["YRI"], [100, 200, 300])
    )
    assert doc["windows"][1]["fst_CEU_YRI"] is None  # monomorphic only


# ─────────────────────────────── windows ────────────────────────────────


def test_windows_tsv_and_json(homes):
    doc = popgen(homes, "windows", "pg", "--window", "700", "--step", "350", "--fst")
    rows = doc["results"]
    assert [(r["start"], r["end"], r["sites"]) for r in rows] == [
        (1, 700, 3),
        (351, 1050, 1),
        (701, 1400, 4),
    ]
    want = reference(POPS, [100, 200, 300])
    for g in POPS:
        assert rows[0][f"{g}_pi"] == pytest.approx(want[g]["pi"])
        assert rows[0][f"{g}_theta_w"] == pytest.approx(want[g]["theta_w"])
        assert rows[0][f"{g}_pi_per_bp"] == pytest.approx(want[g]["pi"] / 700)
    assert "fst_CEU_YRI" in rows[0]
    tsv = {
        b: _vc(home, b, "db", "popgen", "windows", "pg", "--window", "700").stdout
        for b, home in homes.items()
    }
    assert tsv["duckdb"] == tsv["chdb"]
    lines = tsv["duckdb"].splitlines()
    assert lines[0].startswith("chrom\tstart\tend\tsites\tCEU_S")
    assert len(lines) == 3  # header + windows 1-700, 701-1400


def test_table_output_runs(homes):
    for b, home in homes.items():
        out = _vc(home, b, "db", "popgen", "summary", "pg").stdout
        assert "missing_data_tracked: true" in out and "tajima_d" in out


# ─────────────────────────────── refusals ───────────────────────────────


def test_refusals(homes):
    home = homes["duckdb"]
    r = _vc(
        home,
        "duckdb",
        "db",
        "popgen",
        "summary",
        "pg",
        "--include-sex-chroms",
        ok=False,
    )
    assert "not supported" in r.stderr
    r = _vc(
        home,
        "duckdb",
        "db",
        "popgen",
        "summary",
        "pg",
        "--region",
        "chrX:1-500",
        ok=False,
    )
    assert "autosomes only" in r.stderr
    r = _vc(home, "duckdb", "db", "popgen", "fst", "pg", "--by", "all", ok=False)
    assert "at least two groups" in r.stderr
    r = _vc(
        home, "duckdb", "db", "popgen", "summary", "pg", "--ingest-id", "x", ok=False
    )
    assert "no ingestion 'x'" in r.stderr


# ───────────────────────── older databases degrade ──────────────────────


@pytest.mark.parametrize("backend", BACKENDS)
def test_older_database_warns_and_treats_missing_as_reference(
    vcfclick_home, popgen_vcf, backend
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
    for sql in (
        "ALTER TABLE variants DROP COLUMN n_called",
        "ALTER TABLE variants DROP COLUMN an_called",
        "ALTER TABLE variants DROP COLUMN ac_called",
        "DROP TABLE missing_genotypes",
        "DROP TABLE populations",
    ):
        _vc(home, backend, "db", "query", "pg", sql)
    r = _vc(home, backend, "db", "popgen", "summary", "pg", "--format", "json")
    doc = json.loads(r.stdout)
    assert doc["missing_data_tracked"] is False
    assert r.stderr.count("predates called-genotype tracking") == 1  # warned once
    assert doc["grouping"] == "all"  # no panel table → one cohort group
    # Without tracking, ./. reads as 0/0: 1:400, 1:500 and 1:1100 are kept
    # and every site has 22 haplotypes. 1:600 too (haploid stored as gt=2).
    assert doc["sites"]["dropped"]["inexact_group_counts"] == 0
    assert doc["results"][0]["projection_n"] == 22


def test_no_record_missing_drops_sites_with_missing_calls(vcfclick_home, popgen_vcf):
    home, b = vcfclick_home, "duckdb"
    _vc(home, b, "db", "create", "pg")
    _vc(
        home,
        b,
        "db",
        "ingest",
        "pg",
        str(popgen_vcf),
        "--ingest-id",
        "kg",
        "--serial",
        "--no-record-missing",
    )
    r = _vc(home, b, "db", "popgen", "summary", "pg", "--format", "json")
    doc = json.loads(r.stdout)
    assert doc["missing_data_tracked"] is False
    assert "--no-record-missing" in r.stderr
    # 1:400, 1:1100 and 1:1300 have ./. calls → counts not exact → dropped
    assert doc["sites"]["dropped"]["inexact_group_counts"] == 5
    assert doc["sites"]["retained"] == 6
