"""The X-linked model, the X guard and sex / build handling, through the CLI.

Fixture: tests/fixtures/xlinked.vcf (make_xlinked_vcf.py). Every expectation
is checked two ways: against the planted cases (numbers derived by hand from
the pattern table, see the generator) and against a pure-Python reading of
the VCF (`reference`), which applies the same gates to the genotypes as
ingest stores them. The reference uses `cli.xlinked.classify`, whose table
tests/test_xlinked_oracle.py checks against allele transmission.

One ingest per database per backend (module scope); never a chDB session in
the pytest process.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from unittest import mock

import pytest

from cli.xchrom import region
from cli.xlinked import BY_NAME, classify
from tests.conftest import FIXTURES, bgzip_vcf, run_cli

BACKENDS = ("chdb", "duckdb")
VCF = FIXTURES / "xlinked.vcf"
PROBANDS = ("XSON", "XDAU", "YSON", "ZKID")
PARENTS = {
    "XSON": ("XDAD", "XMUM"),
    "XDAU": ("XDAD", "XMUM"),
    "YSON": ("YDAD", "YMUM"),
    "ZKID": ("ZDAD", "ZMUM"),
}
SEX = {"XSON": "male", "XDAU": "female", "YSON": "male", "ZKID": None}
GENES = (("XCH1", 50_000_000, 50_999_999), ("XCH2", 51_000_000, 51_999_999),
         ("XCH3", 52_000_000, 52_999_999))  # fmt: skip
GRCH37_LEN = 155_270_560


# ───────────────────────────── the reference ─────────────────────────────


@dataclass(frozen=True)
class Call:
    gt: int | None
    gq: int | None
    dp: int | None
    ad: tuple[int | None, int | None]


def stored_gt(gt: str, keep_reference: bool) -> int | None:
    """The `genotypes.gt` ingest stores for a GT string, None for no row."""
    alleles = gt.replace("|", "/").split("/")
    if all(a == "." for a in alleles):
        return None
    if len(alleles) == 1:
        value = 2 if alleles[0] == "1" else 0
    else:
        value = sum(a == "1" for a in alleles)
    if value == 0 and not keep_reference:
        return None
    return value


def parse(path: Path, keep_reference: bool) -> list[dict]:
    rows = []
    samples: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith("##"):
            continue
        cols = line.split("\t")
        if line.startswith("#"):
            samples = cols[9:]
            continue
        info = dict(kv.split("=") for kv in cols[7].split(";"))
        keys = cols[8].split(":")
        calls = {}
        for s, cell in zip(samples, cols[9:], strict=True):
            f = dict(zip(keys, cell.split(":"), strict=True))
            gt = stored_gt(f["GT"], keep_reference)
            if gt is None:
                continue
            ad = f.get("AD", ".")
            ref, alt = (None, None) if ad == "." else map(int, ad.split(","))
            calls[s] = Call(gt, int(f["GQ"]), int(f["DP"]), (ref, alt))
        rows.append(
            {
                "chrom": cols[0],
                "pos": int(cols[1]),
                "name": cols[2],
                "af": float(info["AF"]),
                "calls": calls,
            }
        )
    return rows


@dataclass(frozen=True)
class Opts:
    min_gq: int = 20
    min_dp: int = 10
    max_af: float = 0.01
    min_ab: float = 0.25
    max_ab: float = 0.75
    min_male_x_dp: int | None = None
    max_af_x: float = 0.001
    min_hemi_ab: float = 0.75
    max_parent_alt_reads: int = 1

    def args(self) -> list[str]:
        out = []
        for name, default in (
            ("min_male_x_dp", None),
            ("max_af_x", 0.001),
            ("max_parent_alt_reads", 1),
        ):
            value = getattr(self, name)
            if value != default:
                out += [f"--{name.replace('_', '-')}", str(value)]
        return out


def _frac(c: Call) -> float | None:
    ref, alt = c.ad
    if ref is None or alt is None or ref + alt == 0:
        return None
    return alt / (ref + alt)


def _quality(c: Call | None, o: Opts, dp: int) -> bool:
    return c is None or (
        (c.gq is None or c.gq >= o.min_gq) and (c.dp is None or c.dp >= dp)
    )


def _het_ok(c: Call | None, o: Opts) -> bool:
    if c is None or c.gt != 1:
        return True
    f = _frac(c)
    return f is None or o.min_ab <= f <= o.max_ab


def reference(rows, proband, sex, build, o: Opts, has_ref=True):
    """(candidates {(pos, pattern)}, excluded counters) for one proband."""
    father, mother = PARENTS[proband]
    male_dp = o.min_male_x_dp if o.min_male_x_dp is not None else o.min_dp // 2
    cap_af = min(o.max_af, o.max_af_x) + 1e-7
    found, excluded = set(), dict.fromkeys(
        ("male_het_calls", "father_het_calls", "ambiguous_band",
         "mendel_consistent_unpatterned", "mendel_error_other"), 0,
    )  # fmt: skip
    for r in rows:
        reg = region(r["chrom"], r["pos"], build)
        p = r["calls"].get(proband)
        if reg == "ambiguous" and p is not None and p.gt > 0:
            excluded["ambiguous_band"] += 1
        if reg != "npx" or p is None or sex is None:
            continue
        f, m = r["calls"].get(father), r["calls"].get(mother)
        fg, mg = (None if c is None else c.gt for c in (f, m))
        cls = classify(sex, p.gt, fg, mg)
        key = {
            "male_het": "male_het_calls",
            "father_het": "father_het_calls",
            "consistent_unpatterned": "mendel_consistent_unpatterned",
            "mendel_error_other": "mendel_error_other",
        }.get(cls or "")
        if key:
            excluded[key] += 1
        if cls not in BY_NAME:
            continue
        pat = BY_NAME[cls]
        if pat.mode in ("denovo", "female_het") and not has_ref:
            continue
        roles = {"p": (p, male_dp if sex == "male" else o.min_dp), "m": (m, o.min_dp)}
        if sex == "female":
            roles["f"] = (f, male_dp)
        ok = all(_quality(c, o, dp) for c, dp in roles.values())
        ok = ok and _het_ok(p, o) and _het_ok(m, o)
        if sex == "male" and p.gt == 2 and _frac(p) is not None:
            ok = ok and _frac(p) >= o.min_hemi_ab
        for role in pat.cap:
            c = {"f": f, "m": m}[role]
            ok = ok and (c.ad[1] is None or c.ad[1] <= o.max_parent_alt_reads)
        ok = ok and r["af"] <= cap_af
        if ok:
            found.add((r["pos"], cls))
    return found, excluded


# ───────────────────────────── the databases ─────────────────────────────


@dataclass
class World:
    home: Path
    backend: str
    ann: Path


@pytest.fixture(scope="module", params=BACKENDS)
def backend(request):
    if request.param == "chdb" and os.environ.get("VCFCLICK_BACKEND") == "duckdb":
        pytest.skip("chDB variant skipped under VCFCLICK_BACKEND=duckdb")
    return request.param


def _le37(path: Path) -> Path:
    """The fixture without the records beyond the GRCh37 chrX length."""
    lines = [
        ln
        for ln in VCF.read_text().splitlines()
        if ln.startswith("#") or int(ln.split("\t")[1]) <= GRCH37_LEN
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture(scope="module")
def world(backend, tmp_path_factory):
    from annotations.db import get_connection

    base = tmp_path_factory.mktemp(f"xlinked-{backend}")
    home = base / "home"
    home.mkdir()
    ann = base / "ann.duckdb"
    with mock.patch.dict(os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(ann)}):
        conn = get_connection()  # in-process DuckDB, never a chDB session
        for sym, start, end in GENES:
            conn.execute(
                "INSERT INTO refseq_genes VALUES "
                "(?, 'chrX', ?, ?, '+', 'id', 'protein_coding')",
                [sym, start, end],
            )
        conn.close()
    w = World(home, backend, ann)
    le37 = _le37(base / "xlinked_le37.vcf")
    for name, vcf, keep in (("k", VCF, True), ("s", VCF, False), ("le37", le37, True)):
        vc(w, "db", "create", name)
        args = ["db", "ingest", name, str(vcf), "--ingest-id", "i1", "--serial"]
        vc(w, *(args + ["--keep-reference"] if keep else args))
        vc(w, "db", "ped", name, str(FIXTURES / "xlinked.ped"))
    return w


def vc(w: World, *args: str, ok: bool = True):
    with mock.patch.dict(os.environ, {"VCFCLICK_ANNOTATIONS_DB": str(w.ann)}):
        return run_cli(w.home, w.backend, *args, ok=ok)


def trios(w: World, db: str, *extra: str) -> dict[str, dict]:
    """{proband: trio document} for every proband of the pedigree."""
    out = vc(w, "db", "trio", db, "--format", "json", *extra).stdout
    return {d["trio"]["proband"]: d for d in json.loads(out)["trios"]}


def xl(doc: dict) -> dict:
    return doc["models"]["xlinked"]


def hits(doc: dict) -> set[tuple[int, str]]:
    return {(c["pos"], c["pattern"]) for c in xl(doc)["candidates"]}


def modes(doc: dict) -> dict[str, int]:
    return {m: v["count"] for m, v in xl(doc)["modes"].items() if v["count"]}


ROWS = {True: parse(VCF, True), False: parse(VCF, False)}


# ───────────────────────────── planted cases ─────────────────────────────


def test_planted_counts_under_defaults(world):
    docs = trios(world, "k")
    son, dau = docs["XSON"], docs["XDAU"]
    assert xl(son)["count"] == 29 and modes(son) == {"recessive": 18, "denovo": 11}
    assert xl(dau)["count"] == 23
    assert modes(dau) == {"recessive": 3, "female_het": 14, "denovo": 6, "comphet": 1}
    assert [g["gene"] for g in xl(dau)["genes"]] == ["XCH1"]
    patterns = {c["pattern"] for c in xl(dau)["candidates"] if c["mode"] == "denovo"}
    assert patterns == {"x_denovo_daughter", "x_denovo_daughter_mat"}
    ex_son, ex_dau = xl(son)["context"]["excluded"], xl(dau)["context"]["excluded"]
    assert ex_son["male_het_calls"] == 5  # A6, G1, G2, L3, L6
    assert ex_dau["father_het_calls"] == 2  # L1, L2
    assert ex_dau["mendel_error_other"] == 10
    for doc in (son, dau):
        ctx = xl(doc)["context"]
        assert (ctx["build"], ctx["build_source"]) == ("GRCh38", "data")
        assert ctx["proband_sex_source"] == "pedigree" and ctx["warnings"] == []
        assert xl(doc)["blocked"] is False and xl(doc)["blocked_reason"] is None


def test_planted_positives_and_negatives(world):
    son = {c["pos"]: c for c in xl(trios(world, "k")["XSON"])["candidates"]}
    a4 = son[10_300_001]  # father hemizygous ALT: reported, flagged
    assert (a4["pattern"], a4["father_flag"]) == ("xlr_son", "hemizygous_alt")
    assert a4["hemizygous"] == ["proband", "father"] and a4["parent_affected"] is False
    assert son[13_100_001]["father_flag"] == "anomalous_het"  # L2
    assert son[10_200_001]["father_gt"] == -1  # A3: father no-call, not used
    assert son[40_200_001]["pattern"] == "x_denovo_son"  # D3: father carries
    for pos in (10_500_001, 10_600_001, 10_700_001, 10_900_001, 12_000_001,
                13_300_001, 13_400_001, 13_500_001, 1_100_001, 2_781_479,
                155_701_384, 2_750_000):  # fmt: skip
        assert pos not in son, pos
    for pos in (2_781_480, 155_701_381, 154_950_000, 10_800_001, 11_000_001):
        assert son[pos]["pattern"] == "xlr_son", pos


@pytest.mark.parametrize(
    "extra,son,dau",
    [
        ((), 29, 23),
        (("--min-male-x-dp", "8"), 28, 23),  # A9 (son DP 6) drops
        (("--max-af-x", "0.01"), 30, 24),  # + R1, + S1
        (("--max-parent-alt-reads", "10"), 30, 23),  # + L5
    ],
)
def test_option_sets_and_the_reference_agree(world, extra, son, dau):
    o = Opts()
    for flag, value in zip(extra[::2], extra[1::2], strict=True):
        name = flag.lstrip("-").replace("-", "_")
        o = Opts(**{**o.__dict__, name: type(getattr(Opts(), name) or 0)(value)})
    for db, keep in (("k", True), ("s", False)):
        docs = trios(world, db, *extra)
        if keep:
            assert (xl(docs["XSON"])["count"], xl(docs["XDAU"])["count"]) == (son, dau)
        for p in ("XSON", "XDAU", "YSON"):
            want, excluded = reference(ROWS[keep], p, SEX[p], "GRCh38", o, keep)
            assert hits(docs[p]) == want, (db, p)
            assert xl(docs[p])["context"]["excluded"] == excluded, (db, p)


def test_haploid_and_diploid_coded_males_give_the_same_sites(world):
    docs = trios(world, "k")
    assert {p for p, _ in hits(docs["XSON"])} == {p for p, _ in hits(docs["YSON"])}
    assert xl(docs["YSON"])["context"]["excluded"]["male_het_calls"] == 5


def test_sparse_ingest_blocks_the_modes_that_need_a_parent_0_0(world):
    docs = trios(world, "s")
    son, dau = xl(docs["XSON"]), xl(docs["XDAU"])
    assert son["count"] == 18 and son["blocked"] is False
    assert son["modes"]["recessive"]["blocked"] is False
    assert son["modes"]["denovo"]["blocked_reason"] == "needs_keep_reference"
    assert dau["count"] == 3
    for mode in ("dominant", "female_het", "denovo", "comphet"):
        assert dau["modes"][mode]["blocked_reason"] == "needs_keep_reference", mode
    text = vc(world, "db", "trio", "s", "--proband", "XSON").stdout
    assert "  xlinked        18  (recessive 18; denovo needs --keep-reference)" in text


def test_unresolved_sex_blocks_xlinked_and_keeps_the_classic_models(world):
    doc = trios(world, "k")["ZKID"]
    x = xl(doc)
    assert (x["blocked"], x["blocked_reason"], x["count"]) == (True, "needs_sex", 0)
    assert x["context"]["proband_sex"] is None and x["candidates"] == []
    assert {"denovo", "recessive", "dominant", "comphet"} <= set(doc["models"])
    text = vc(world, "db", "trio", "k", "--proband", "ZKID").stdout
    assert "  xlinked         -  (needs proband sex" in text


def test_text_lines(world):
    text = vc(world, "db", "trio", "k", "--proband", "XSON").stdout
    assert "  xlinked        29  (recessive 18, denovo 11)" in text
    text = vc(world, "db", "trio", "k", "--proband", "XDAU").stdout
    assert "  xlinked        23  (recessive 3, female_het 14, denovo 6)" in text
    r = vc(world, "db", "trio", "k", "--proband", "XSON", "--category", "xlinked")
    assert (
        "xlinked candidates: 29  (proband=XSON male [pedigree], father=XDAD male "
        "[pedigree], mother=XMUM female [pedigree]; GRCh38 [data])"
    ) in r.stdout
    assert "chrX:10300001 C>T  proband_gt=2 father_gt=2 mother_gt=1" in r.stdout
    assert "mode=recessive pattern=xlr_son father=hemizygous_alt" in r.stdout
    assert "excluded: male_het_calls=5" in r.stdout
    r = vc(world, "db", "trio", "k", "--proband", "XDAU", "--category", "xlinked")
    assert "X compound-het genes: 1" in r.stdout and "  XCH1" in r.stdout


# ───────────────────────────── the X guard ─────────────────────────────


def _legacy(doc: dict) -> dict:
    return {m: doc["models"][m]["count"] for m in ("denovo", "recessive", "dominant")}


def test_x_guard_and_no_x_guard(world):
    on, off = trios(world, "k", "--x-guard"), trios(world, "k")
    assert _legacy(on["XSON"]) == {"denovo": 6, "recessive": 3, "dominant": 0}
    assert _legacy(off["XSON"]) == {"denovo": 7, "recessive": 4, "dominant": 4}
    assert _legacy(on["XDAU"]) == {"denovo": 6, "recessive": 0, "dominant": 15}
    assert _legacy(off["XDAU"]) == {"denovo": 6, "recessive": 1, "dominant": 16}
    assert on["XSON"]["x_guard"] is True and off["XSON"]["x_guard"] is False
    excluded = {m: on["XSON"]["models"][m]["x_guard_excluded"] for m in on["XSON"]["models"]
                if m != "xlinked"}  # fmt: skip
    assert excluded == {"denovo": 1, "recessive": 1, "dominant": 4, "comphet": 1}
    assert off["XSON"]["models"]["denovo"]["x_guard_excluded"] is None
    # comphet: the male-het pair over XCH3 is gone with the guard
    assert on["XSON"]["models"]["comphet"]["count"] == 0
    assert [g["gene"] for g in off["XSON"]["models"]["comphet"]["genes"]] == ["XCH3"]
    # a daughter's rows with a hemizygous father are untouched
    dom_on = {c["pos"] for c in on["XDAU"]["models"]["dominant"]["candidates"]}
    dom_off = {c["pos"] for c in off["XDAU"]["models"]["dominant"]["candidates"]}
    assert dom_off - dom_on == {13_000_001}  # L1, the het father
    # the X model itself does not depend on the guard
    assert hits(on["XSON"]) == hits(off["XSON"])


def test_af_exactly_at_the_cap_is_kept_on_both_backends(world):
    """INFO/AF 0.001 is a Float32 just above 0.001: the gate's epsilon keeps it
    on chDB as on DuckDB."""
    doc = trios(world, "k", "--max-af", "0.001", "--x-guard")["XDAU"]
    assert doc["models"]["dominant"]["count"] == 14  # S1 (0.005) out, L1 guarded


# ───────────────────────────── build ─────────────────────────────


def test_build_auto_unknown_and_explicit(world):
    auto = trios(world, "le37")
    ctx = xl(auto["XSON"])["context"]
    assert (ctx["build"], ctx["build_source"]) == ("unknown", "unknown")
    assert [w["code"] for w in ctx["warnings"]] == ["build_unknown"]
    assert ctx["excluded"]["ambiguous_band"] == 3  # P2, P6, P7
    assert xl(auto["XSON"])["count"] == 27
    assert 154_950_000 not in {p for p, _ in hits(auto["XSON"])}
    b38 = trios(world, "le37", "--build", "GRCh38")
    b37 = trios(world, "le37", "--build", "GRCh37")
    assert xl(b38["XSON"])["count"] == 28  # P7 joins
    assert xl(b37["XSON"])["count"] == 29  # P2, P6 join; P7 leaves
    assert xl(b37["XSON"])["context"]["build_source"] == "option"
    for docs, build in ((auto, "unknown"), (b38, "GRCh38"), (b37, "GRCh37")):
        assert xl(docs["XDAU"])["count"] == 23
        rows = parse(_le37_path(world), True)
        for p in ("XSON", "XDAU"):
            want, excluded = reference(rows, p, SEX[p], build, Opts())
            assert hits(docs[p]) == want, (build, p)
            assert xl(docs[p])["context"]["excluded"] == excluded, (build, p)
    conflict = xl(trios(world, "k", "--build", "GRCh37")["XSON"])
    assert (conflict["blocked"], conflict["blocked_reason"]) == (True, "build_conflict")


def _le37_path(world: World) -> Path:
    return _le37(world.home.parent / "xlinked_le37_ref.vcf")


# ───────────────────────────── sex ─────────────────────────────


def test_sex_paths(world, tmp_path):
    def load(ped: str) -> dict:
        vc(world, "db", "ped", "k", str(FIXTURES / ped))
        return trios(world, "k")

    try:
        docs = load("xlinked_nosex.ped")
        for p, sex in (("XSON", "male"), ("XDAU", "female")):
            ctx = xl(docs[p])["context"]
            assert (ctx["proband_sex"], ctx["proband_sex_source"]) == (sex, "inferred")
            assert ctx["father_sex_source"] == "inferred"
        assert xl(docs["XSON"])["count"] == 29 and xl(docs["XDAU"])["count"] == 23

        docs = load("xlinked_conflict.ped")
        ctx = xl(docs["XSON"])["context"]
        assert (ctx["proband_sex"], ctx["proband_sex_source"]) == ("female", "pedigree")
        assert [w["code"] for w in ctx["warnings"]] == ["sex_conflict"]

        docs = load("xlinked_badfather.ped")
        for p in ("XSON", "XDAU"):
            x = xl(docs[p])
            assert x["blocked_reason"] == "parent_sex_conflict" and x["count"] == 0
            assert docs[p]["models"]["recessive"]["count"] > 0  # legacy still runs
        text = vc(world, "db", "trio", "k", "--proband", "XSON").stdout
        assert "(parent sex conflict: XDAD is the father, declared female)" in text

        docs = load("xlinked_affected.ped")
        assert modes(docs["XDAU"]) == {
            "recessive": 3, "dominant": 14, "denovo": 6, "comphet": 1,
        }  # fmt: skip
        paternal = [c for c in xl(docs["XDAU"])["candidates"]
                    if c["pattern"] == "xhet_paternal"]  # fmt: skip
        assert paternal and all(c["parent_affected"] for c in paternal)
    finally:
        vc(world, "db", "ped", "k", str(FIXTURES / "xlinked.ped"))


def test_db_ped_flags_never_make_dominant_and_sex_comes_from_chrx(world):
    """The desktop apps' path: `db ped` flags leave parents' phenotype unknown,
    so a het daughter stays female_het; her sex is inferred until --sex."""
    vc(world, "db", "create", "f")
    vc(world, "db", "ingest", "f", str(VCF), "--ingest-id", "i1", "--serial",
       "--keep-reference")  # fmt: skip
    vc(world, "db", "ped", "f", "--proband", "XDAU", "--father", "XDAD",
       "--mother", "XMUM")  # fmt: skip
    doc = json.loads(
        vc(world, "db", "trio", "f", "--proband", "XDAU", "--format", "json").stdout
    )
    x = xl(doc)
    assert x["context"]["proband_sex_source"] == "inferred"
    assert (
        x["modes"]["dominant"]["count"] == 0 and x["modes"]["female_het"]["count"] == 14
    )
    vc(world, "db", "ped", "f", "--proband", "XDAU", "--father", "XDAD",
       "--mother", "XMUM", "--sex", "XDAU=female")  # fmt: skip
    doc = json.loads(
        vc(world, "db", "trio", "f", "--proband", "XDAU", "--format", "json").stdout
    )
    assert xl(doc)["context"]["proband_sex_source"] == "pedigree"
    assert xl(doc)["count"] == 23


def test_db_qc_x_nonpar_reports_the_sex_db_trio_uses(world):
    plain = json.loads(vc(world, "db", "qc", "k", "--format", "json").stdout)
    r = vc(world, "db", "qc", "k", "--x-nonpar", "--format", "json")
    npx = {s["sample_id"]: s for s in json.loads(r.stdout)}
    assert "non-PAR X (GRCh38 [data])" in r.stderr
    assert [list(s) for s in plain] == [list(npx[s["sample_id"]]) for s in plain]
    for sample, sex in (("XDAD", "male"), ("XSON", "male"), ("XMUM", "female"),
                        ("XDAU", "female"), ("ZKID", "unknown")):  # fmt: skip
        assert npx[sample]["inferred_sex"] == sex, sample
    # the PAR hets of the diploid males are left out
    by = {s["sample_id"]: s for s in plain}
    assert npx["XDAD"]["chrx_het_frac"] < by["XDAD"]["chrx_het_frac"]
    r = vc(world, "db", "qc", "k", "--x-nonpar", "--build", "GRCh37", ok=False)
    assert "not GRCh37" in r.stderr


# ───────────────────────────── bcftools (optional) ─────────────────────────────


def _merr(vcf: Path, pfm: str) -> dict[int, int]:
    r = subprocess.run(
        ["bcftools", "+mendelian2", str(vcf), "-p", pfm, "--rules", "GRCh38",
         "-m", "a"],
        capture_output=True, text=True, check=True,
    )  # fmt: skip
    out = {}
    for line in r.stdout.splitlines():
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        info = dict(kv.split("=") for kv in cols[7].split(";") if "=" in kv)
        out[int(cols[1])] = int(info["MERR"])
    return out


def test_candidates_agree_with_bcftools_mendelian2(world):
    if not shutil.which("bcftools"):
        pytest.skip("bcftools not on PATH")
    docs = trios(world, "k")
    for proband, pfm in (("XSON", "1X:XSON,XDAD,XMUM"), ("XDAU", "2X:XDAU,XDAD,XMUM")):
        merr = _merr(VCF, pfm)
        for c in xl(docs[proband])["candidates"]:
            want = 1 if c["mode"] == "denovo" else 0
            assert merr[c["pos"]] == want, (proband, c["pos"], c["pattern"])


def test_inferred_sexes_agree_with_bcftools_guess_ploidy(world, tmp_path):
    """On the samples with chrX evidence. guess-ploidy calls a sample with no
    X data male (ZMUM, ZKID), which is why inference needs 20 calls."""
    if not shutil.which("bcftools"):
        pytest.skip("bcftools not on PATH")
    gz = bgzip_vcf(VCF, tmp_path / "xlinked.vcf.gz")
    r = subprocess.run(
        ["bcftools", "+guess-ploidy", "-g", "hg38", "-t", "GT", str(gz)],
        capture_output=True, text=True, check=True,
    )  # fmt: skip
    theirs = dict(line.split("\t") for line in r.stdout.strip().splitlines())
    out = vc(world, "db", "qc", "k", "--x-nonpar", "--format", "json").stdout
    ours = {s["sample_id"]: s["inferred_sex"] for s in json.loads(out)}
    for sample in ("XDAD", "XMUM", "XSON", "XDAU", "YDAD", "YMUM", "YSON"):
        assert {"M": "male", "F": "female"}[theirs[sample]] == ours[sample], sample
