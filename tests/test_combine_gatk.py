"""GATK CombineVariants parity options for `combine` (all opt-in; default
behaviour unchanged). Fixture shapes mirror the real shadow-comparison sites
without the clinical data.
"""

from __future__ import annotations

from pathlib import Path

from cyvcf2 import VCF

from ingest.combine import combine_vcfs

_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1>\n"
    '##FILTER=<ID=PASS,Description="x">\n'
    '##FILTER=<ID=AFB,Description="x">\n'
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="d">\n'
    '##INFO=<ID=AO,Number=A,Type=Integer,Description="a">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def _wv(path: Path, rows, contig="chr1") -> str:
    """rows: (pos, ref, alt, qual, filter, info, gt)."""
    hdr = _HDR.replace("ID=chr1", f"ID={contig}")
    body = "".join(
        f"{contig}\t{p}\t.\t{r}\t{a}\t{q}\t{flt}\t{info}\tGT\t{gt}\n"
        for (p, r, a, q, flt, info, gt) in rows
    )
    path.write_text(hdr + body)
    return str(path)


def _sites(path: Path) -> dict:
    return {rec.POS: rec for rec in VCF(str(path))}


def test_pass_only_counts_only_pass_calls(tmp_path):
    # chr3:10091049 shape: one PASS caller + one filtered caller at the same allele.
    a = _wv(tmp_path / "a.vcf", [(100, "C", "T", ".", "PASS", ".", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(100, "C", "T", ".", "AFB", ".", "0/1")])

    # Default: both inputs count -> kept at min_callsets 2.
    d = tmp_path / "def.vcf"
    combine_vcfs([a, b], d, min_callsets=2)
    assert 100 in _sites(d)

    # --pass-only: only 1 PASS input -> dropped at min_callsets 2.
    p = tmp_path / "po.vcf"
    combine_vcfs([a, b], p, min_callsets=2, pass_only=True)
    assert 100 not in _sites(p)


def test_pass_only_names_filtered_inputs_filterin(tmp_path):
    a = _wv(tmp_path / "a.vcf", [(100, "C", "T", ".", "PASS", ".", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(100, "C", "T", ".", "AFB", ".", "0/1")])
    out = tmp_path / "po1.vcf"
    combine_vcfs(
        [a, b], out, min_callsets=1, pass_only=True, names=["gatk3", "octopus"]
    )
    assert _sites(out)[100].INFO.get("set") == "gatk3-filterInoctopus"


def test_filter_dot_treated_as_pass(tmp_path):
    # §6: an unfiltered record written as FILTER=. counts as PASS.
    a = _wv(tmp_path / "a.vcf", [(100, "C", "T", ".", ".", ".", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(100, "C", "T", ".", ".", ".", "0/1")])
    out = tmp_path / "o.vcf"
    combine_vcfs([a, b], out, min_callsets=2, pass_only=True)
    assert 100 in _sites(out)  # both '.' count as PASS


def _alts_at(path: Path, pos: int) -> list[str]:
    return sorted(rec.ALT[0] for rec in VCF(str(path)) if rec.POS == pos and rec.ALT)


def test_count_by_site_keeps_all_alleles_at_shared_position(tmp_path):
    # chr1:241663902 shape: two callers, each with a DIFFERENT allele at one
    # position. Allele-count drops both at min_callsets 2; site-count keeps both.
    a = _wv(tmp_path / "a.vcf", [(100, "T", "A", ".", "PASS", ".", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(100, "T", "C", ".", "PASS", ".", "0/1")])

    d = tmp_path / "allele.vcf"
    combine_vcfs([a, b], d, min_callsets=2)  # default = allele
    assert _alts_at(d, 100) == []  # each allele seen by 1 input -> dropped

    s = tmp_path / "site.vcf"
    combine_vcfs([a, b], s, min_callsets=2, count_by="site", names=["fb", "gatk"])
    assert _alts_at(s, 100) == ["A", "C"]  # position has 2 inputs -> both kept
    # each allele labelled with the position's set= (Intersection: both present)
    for rec in VCF(str(s)):
        if rec.POS == 100:
            assert rec.INFO.get("set") == "Intersection"


def test_carry_info_takes_priority_input_qual_filter_info(tmp_path):
    # Priority input A (higher) supplies QUAL/FILTER/INFO for the shared allele.
    a = _wv(tmp_path / "a.vcf", [(100, "C", "T", "55", "PASS", "DP=30;AO=5", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(100, "C", "T", "10", "AFB", "DP=8", "0/1")])

    out = tmp_path / "carry.vcf"
    combine_vcfs([a, b], out, carry_info=True, names=["fb", "gatk"])
    rec = _sites(out)[100]
    assert rec.QUAL == 55
    assert rec.FILTER is None  # PASS carried
    assert rec.INFO.get("DP") == 30 and rec.INFO.get("AO") == 5
    assert rec.INFO.get("set") == "Intersection"

    # Default (no carry): QUAL/FILTER are '.', INFO holds only set=.
    d = tmp_path / "nocarry.vcf"
    combine_vcfs([a, b], d)
    rd = _sites(d)[100]
    assert rd.QUAL is None
    assert rd.INFO.get("DP") is None


def test_reference_normalizes_padded_indels_to_one_site(tmp_path):
    # chr9:135773000 shape: two spellings of the same homopolymer insertion
    # left-align to one key -> merge -> Intersection.
    import pytest

    pytest.importorskip("pyfaidx")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nCAAAAAT\n")
    a = _wv(tmp_path / "a.vcf", [(1, "C", "CA", ".", "PASS", ".", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(3, "A", "AA", ".", "PASS", ".", "0/1")])
    out = tmp_path / "o.vcf"
    combine_vcfs([a, b], out, reference=str(ref), min_callsets=2)
    recs = list(VCF(str(out)))
    assert len(recs) == 1  # both representations merged
    assert recs[0].INFO.get("set") == "Intersection"
    assert (recs[0].POS, recs[0].REF, recs[0].ALT[0]) == (1, "C", "CA")


def test_reference_splits_multiallelic_instead_of_refusing(tmp_path):
    import pytest

    pytest.importorskip("pyfaidx")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    a = _wv(tmp_path / "a.vcf", [(1, "A", "C,G", ".", "PASS", ".", "1/2")])
    b = _wv(tmp_path / "b.vcf", [(1, "A", "C", ".", "PASS", ".", "0/1")])
    out = tmp_path / "o.vcf"
    combine_vcfs([a, b], out, reference=str(ref))
    alts = sorted(rec.ALT[0] for rec in VCF(str(out)))
    assert alts == ["C", "G"]  # split: A>C (both), A>G (a only)


def test_cli_wires_all_four_options(tmp_path):
    from click.testing import CliRunner

    from cli.combine import combine

    a = _wv(tmp_path / "a.vcf", [(100, "C", "T", "9", "PASS", "DP=30", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(100, "C", "T", "1", "AFB", "DP=8", "0/1")])
    out = tmp_path / "o.vcf"
    res = CliRunner().invoke(
        combine,
        [
            str(a),
            str(b),
            "-o",
            str(out),
            "--pass-only",
            "--count-by",
            "site",
            "--carry-info",
            "--name",
            "gatk3",
            "--name",
            "octopus",
            "--min-callsets",
            "1",
        ],
    )
    assert res.exit_code == 0, res.output
    rec = _sites(out)[100]
    assert rec.INFO.get("set") == "gatk3-filterInoctopus"  # pass-only naming
    assert rec.QUAL == 9 and rec.INFO.get("DP") == 30  # carry from priority


# --- codex-review fixes ---------------------------------------------------

import pytest  # noqa: E402

from ingest.combine import CombineError  # noqa: E402


def _raw_col(path: Path, pos: int, col: int) -> str:
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        c = line.split("\t")
        if int(c[1]) == pos:
            return c[col]
    return ""


def test_carry_info_filter_dot_stays_verbatim(tmp_path):
    # #7: carried FILTER must be verbatim '.', not rewritten to PASS.
    a = _wv(tmp_path / "a.vcf", [(100, "C", "T", "9", ".", "DP=30", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(100, "C", "T", "1", ".", "DP=8", "0/1")])
    out = tmp_path / "o.vcf"
    combine_vcfs([a, b], out, carry_info=True)
    assert _raw_col(out, 100, 6) == "."  # FILTER column verbatim


def test_carry_info_strips_existing_set(tmp_path):
    # #2: an input that already has set= must not yield a duplicate INFO key.
    a = _wv(tmp_path / "a.vcf", [(100, "C", "T", "9", "PASS", "set=old;DP=30", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(100, "C", "T", "1", "PASS", "DP=8", "0/1")])
    out = tmp_path / "o.vcf"
    combine_vcfs([a, b], out, carry_info=True)
    assert _raw_col(out, 100, 7).count("set=") == 1


def test_reference_handles_contig_start_indel(tmp_path):
    # #3: POS=1 insertion must not crash left-align.
    pytest.importorskip("pyfaidx")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    a = _wv(tmp_path / "a.vcf", [(1, "A", "AA", ".", "PASS", ".", "0/1")])
    b = _wv(tmp_path / "b.vcf", [(1, "A", "AA", ".", "PASS", ".", "0/1")])
    out = tmp_path / "o.vcf"
    combine_vcfs([a, b], out, reference=str(ref))  # must not raise
    assert 1 in _sites(out)


def test_reference_carry_info_projects_per_allele(tmp_path):
    # #1: per-alt INFO (Number=A) must be projected to each split allele.
    pytest.importorskip("pyfaidx")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    a = _wv(tmp_path / "a.vcf", [(1, "A", "C,G", ".", "PASS", "AO=3,7", "1/2")])
    b = _wv(tmp_path / "b.vcf", [(1, "A", "C", ".", "PASS", "AO=5", "0/1")])
    out = tmp_path / "o.vcf"
    combine_vcfs([a, b], out, reference=str(ref), carry_info=True)
    recs = {rec.ALT[0]: rec for rec in VCF(str(out))}
    assert recs["C"].INFO.get("AO") == 3
    assert recs["G"].INFO.get("AO") == 7


_HDR_AD = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
    '##FILTER=<ID=PASS,Description="x">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def test_reference_splits_multiallelic_ad_per_allele(tmp_path):
    # #4: AD (Number=R) must be subset to [ref, alt_j] per split allele.
    pytest.importorskip("pyfaidx")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    a = tmp_path / "a.vcf"
    a.write_text(_HDR_AD + "chr1\t1\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t1/2:10,3,7\n")
    b = tmp_path / "b.vcf"
    b.write_text(_HDR_AD + "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT:AD\t0/1:12,5\n")
    out = tmp_path / "o.vcf"
    combine_vcfs([str(a), str(b)], out, reference=str(ref))
    recs = {rec.ALT[0]: rec for rec in VCF(str(out))}
    assert list(recs["C"].format("AD")[0]) == [10, 3]
    assert list(recs["G"].format("AD")[0]) == [10, 7]


def test_reference_keeps_reference_record_for_site_consensus(tmp_path):
    # #5: an ALT='.' (monomorphic) record must still count at the position.
    pytest.importorskip("pyfaidx")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    a = _wv(tmp_path / "a.vcf", [(1, "A", ".", ".", "PASS", ".", "0/0")])
    b = _wv(tmp_path / "b.vcf", [(1, "A", "C", ".", "PASS", ".", "0/1")])
    out = tmp_path / "o.vcf"
    combine_vcfs([a, b], out, reference=str(ref), count_by="site", min_callsets=2)
    assert "C" in _alts_at(out, 1)


def test_carry_info_rejects_incompatible_header(tmp_path):
    # #6: same INFO ID with a different Type across inputs must error, not
    # silently misparse.
    ha = _HDR
    hb = _HDR.replace("ID=DP,Number=1,Type=Integer", "ID=DP,Number=1,Type=String")
    a = tmp_path / "a.vcf"
    a.write_text(ha + "chr1\t100\t.\tC\tT\t.\tPASS\tDP=30\tGT\t0/1\n")
    b = tmp_path / "b.vcf"
    b.write_text(hb + "chr1\t200\t.\tC\tT\t.\tPASS\tDP=9\tGT\t0/1\n")
    out = tmp_path / "o.vcf"
    with pytest.raises(CombineError, match="incompatible"):
        combine_vcfs([str(a), str(b)], out, carry_info=True)
