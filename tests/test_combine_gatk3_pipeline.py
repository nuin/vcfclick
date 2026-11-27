"""`vcfclick combine --gatk3` against frozen GATK 3.8-1 CombineVariants output.

The four synthetic single-sample VCFs in tests/fixtures/gatk3_pipeline are the
raw output shape of FreeBayes, GATK4 HaplotypeCaller, GATK3 HaplotypeCaller
and Octopus for one sample: multi-allelic records, padded indels, non-PASS
FILTER values. GATK 3.8-1 (Java 8) merged them with the clinical pipeline
command; the output is frozen in the same directory (GOLDEN_COMMANDS.txt
holds the commands). Every header line that carries data and every record
must match. GATK needs Java 8 and is not run in CI.

Sites in the fixtures (chr1 unless noted) and what each pins:
  100  all four callers, GT agree, INFO from the first of an AC tie
  150  GT disagree; INFO from the first of a max-AC tie (gatk, octopus)
  200  one caller only: dropped by minimumN 2
  250  exactly two callers
  350, 400, 450, 455, 460  filtered records: counted out of minimumN
  470, 480  FILTER PASS and . merge
  500, 550, 600, 765  multi-allelic records: allele union, PL/AD dropped
  650  first caller ./. wins over a later called genotype
  760, 765  same-start deletions of different length: REF extension
  799/804, 899/901  the same indel spelled two ways: two records
  820, 850, 880, chr2:200  SNP/indel/MNP at one position: split by type
  950  no AC anywhere: GATK3 drops set=
  1000  ID union; 1050-1070 QUAL; 1100 AC recomputed from a 0/0 genotype
  1150  phased genotypes; 1200 GL converted to PL; 1020-1045 GQ rounding
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from cli.main import cli
from ingest.combine import CombineError, combine_vcfs
from ingest.combine_gatk3 import (
    format_qual,
    format_vcf_double,
    java_hashmap_order,
    java_string_hash,
)

FIX = Path(__file__).parent / "fixtures" / "gatk3_pipeline"
NAMES = ["freebayes", "gatk", "gatk3", "octopus"]
INPUTS = [
    FIX / "S_freebayes.final.vcf",
    FIX / "S_GATK.vcf",
    FIX / "S_GATK3.vcf",
    FIX / "S_octopus.vcf",
]
REF = FIX / "ref.fasta"

DATA_HEADERS = ("##fileformat", "##INFO", "##FORMAT", "##FILTER", "##contig", "#CHROM")


def data_lines(path: Path) -> tuple[list[str], list[str]]:
    """(header lines that carry data, records) of a VCF."""
    header, records = [], []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            if line.startswith(DATA_HEADERS):
                header.append(line)
        elif line:
            records.append(line)
    return header, records


def merge(tmp_path, **kw) -> Path:
    out = tmp_path / "merged.vcf"
    combine_vcfs(INPUTS, out, names=NAMES, gatk3=True, reference=REF, **kw)
    return out


def at(records: list[str], pos: int, chrom: str = "chr1") -> list[list[str]]:
    return [r.split("\t") for r in records if r.startswith(f"{chrom}\t{pos}\t")]


def assert_same_as_golden(out: Path, golden: str) -> None:
    got_h, got_r = data_lines(out)
    want_h, want_r = data_lines(FIX / golden)
    assert got_h == want_h
    assert len(got_r) == len(want_r)
    diffs = [(g, w) for g, w in zip(got_r, want_r) if g != w]
    assert not diffs, diffs[:3]


def test_pipeline_command_matches_gatk3(tmp_path):
    """--genotypemergeoption UNSORTED --mergeInfoWithMaxAC --minimumN 2."""
    out = merge(
        tmp_path,
        genotype_merge="unsorted",
        merge_info_with_max_ac=True,
        min_callsets=2,
    )
    assert_same_as_golden(out, "golden_pipeline.gatk3.vcf")


def test_default_info_merge_minimum_one_matches_gatk3(tmp_path):
    out = merge(tmp_path, genotype_merge="unsorted", min_callsets=1)
    assert_same_as_golden(out, "golden_unsorted_minn1.gatk3.vcf")


def test_prioritize_with_priority_and_filtered_uncalled_matches_gatk3(tmp_path):
    out = merge(
        tmp_path,
        genotype_merge="prioritize",
        priority=["octopus", "gatk3", "gatk", "freebayes"],
        merge_info_with_max_ac=True,
        min_callsets=2,
        filtered_are_uncalled=True,
    )
    assert_same_as_golden(out, "golden_prioritize_uncalled.gatk3.vcf")


def test_prioritize_in_input_order_equals_unsorted(tmp_path):
    """GATK3 UNSORTED keeps the first input's record for a sample, and so does
    PRIORITIZE with the priority list equal to the input order."""
    a = tmp_path / "a.vcf"
    b = tmp_path / "b.vcf"
    kw = dict(names=NAMES, gatk3=True, reference=REF, merge_info_with_max_ac=True)
    combine_vcfs(INPUTS, a, genotype_merge="unsorted", min_callsets=2, **kw)
    combine_vcfs(
        INPUTS, b, genotype_merge="prioritize", priority=NAMES, min_callsets=2, **kw
    )
    assert data_lines(a) == data_lines(b)


@pytest.fixture(scope="module")
def pipeline_records(tmp_path_factory):
    out = tmp_path_factory.mktemp("gatk3") / "m.vcf"
    combine_vcfs(
        INPUTS,
        out,
        names=NAMES,
        gatk3=True,
        reference=REF,
        genotype_merge="unsorted",
        merge_info_with_max_ac=True,
        min_callsets=2,
    )
    return data_lines(out)[1]


def test_genotype_disagreement_keeps_first_input(pipeline_records):
    # freebayes 0/1, gatk 1/1, gatk3 0/1, octopus 1/1 at chr1:150
    (rec,) = at(pipeline_records, 150)
    assert rec[9].startswith("0/1:")
    # INFO from gatk (AC=2 tie with octopus, first wins), AC recomputed to 1
    assert "MQ=60;QD=20" in rec[7] and "AC=1;AF=0.500" in rec[7]
    assert rec[7].endswith("set=Intersection")


def test_first_input_no_call_wins(pipeline_records):
    (rec,) = at(pipeline_records, 650)
    assert rec[9].startswith("./.:")
    assert "AC=0;AF=0.00;AN=0" in rec[7]


def test_min_callsets_counts_records_minus_filtered(pipeline_records):
    pos = {int(r.split("\t")[1]) for r in pipeline_records if r.startswith("chr1\t")}
    assert 200 not in pos  # one caller
    assert 400 not in pos  # one PASS record and one filtered
    assert 450 not in pos  # one PASS record and two filtered
    assert 460 not in pos  # four filtered records
    assert 455 in pos  # two PASS and two filtered
    assert 250 in pos  # exactly two callers
    (rec,) = at(pipeline_records, 455)
    assert rec[6] == "PASS"
    assert rec[7].endswith("set=freebayes-filterIngatk-gatk3-filterInoctopus")
    # one input with a SNP and an indel at 880 reaches two records
    assert len(at(pipeline_records, 880)) == 2


def test_indel_spelled_differently_stays_two_records(pipeline_records):
    (left,) = at(pipeline_records, 799)
    (right,) = at(pipeline_records, 804)
    assert left[7].endswith("set=gatk-gatk3")
    assert right[7].endswith("set=freebayes-octopus")


def test_same_start_deletions_extend_reference(pipeline_records):
    (rec,) = at(pipeline_records, 760)
    assert rec[3] == "CGT" and rec[4] == "CT,C"
    # allele lists differ from the merged one, so AD and PL are dropped
    assert rec[8] == "GT:AO:DP:GQ:RO"


def test_snp_and_indel_at_one_position_split_by_type(pipeline_records):
    snp, ins = at(pipeline_records, 820)
    assert snp[4] == "C" and ins[4] == "ACC"
    assert snp[7].endswith("set=freebayes") and ins[7].endswith("set=gatk3")


def test_no_ac_drops_set(pipeline_records):
    (rec,) = at(pipeline_records, 950)
    assert rec[7] == "AC=1;AF=0.500;AN=2"


def test_qual_id_and_gq_rules(pipeline_records):
    assert at(pipeline_records, 1000)[0][2] == "rs11,rs22"
    assert at(pipeline_records, 1050)[0][5] == "60.50"  # first non-missing QUAL
    assert at(pipeline_records, 1070)[0][5] == "0"
    gq = {p: at(pipeline_records, p)[0][9] for p in (1020, 1030, 1040, 1045)}
    assert gq[1020].split(":")[3] == "99"  # 150.3 capped; DP missing so absent
    assert gq[1030].split(":")[3] == "99"  # integer 100 capped
    assert ":13:" in gq[1040] and ":12:" in gq[1045]  # 12.5 rounds up, 12.4 down


def test_gl_becomes_pl(pipeline_records):
    (rec,) = at(pipeline_records, 1200)
    assert rec[9].split(":")[-2] == "2,0,24"  # GL -0.3,-0.1,-2.5


def test_filter_union_of_a_filtered_type_group(pipeline_records):
    snp, ins = at(pipeline_records, 1010)
    assert snp[6] == "PASS" and ins[6] == "LowQual;q10"
    assert ins[7].endswith("set=FilteredInAll")


def test_header_merge_follows_java_hashmap_order():
    # GATK3 3.8-1 merges headers in HashMap order of the rod names; the first
    # header in that order wins a description conflict.
    order = java_hashmap_order([java_string_hash(n) for n in NAMES])
    assert [NAMES[i] for i in order] == ["gatk", "octopus", "freebayes", "gatk3"]
    order = java_hashmap_order([java_string_hash(n) for n in ["x1", "x2", "x3", "x4"]])
    assert order == [0, 1, 2, 3]


def test_java_number_formats():
    assert format_vcf_double(0.0) == "0.00"
    assert format_vcf_double(0.5) == "0.500"
    assert format_vcf_double(2 / 3) == "0.667"
    assert format_vcf_double(1.0) == "1.00"
    assert format_vcf_double(0.0005) == "5.000e-04"
    assert format_qual(-6.05) == "60.50"
    assert format_qual(-9.9) == "99"
    assert format_qual(-0.0004) == "0"


def test_gatk3_rejects_normalising_options(tmp_path):
    with pytest.raises(CombineError, match="--pass-only"):
        combine_vcfs(INPUTS, tmp_path / "x.vcf", gatk3=True, pass_only=True)
    with pytest.raises(CombineError, match="--priority"):
        combine_vcfs(
            INPUTS,
            tmp_path / "x.vcf",
            gatk3=True,
            genotype_merge="unsorted",
            priority=NAMES,
        )


def test_incomplete_header_is_an_error(tmp_path):
    base = (
        '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )
    for n in "ab":
        (tmp_path / f"{n}.vcf").write_text(
            base + "chr1\t100\t.\tT\tC\t10\tLowQ\t.\tGT\t0/1\n"
        )
    with pytest.raises(CombineError, match="LowQ"):
        combine_vcfs(
            [tmp_path / "a.vcf", tmp_path / "b.vcf"],
            tmp_path / "o.vcf",
            gatk3=True,
            min_callsets=1,
        )


def test_cli_gatk3_flag(tmp_path):
    out = tmp_path / "cli.vcf"
    args = ["combine", *map(str, INPUTS), "-o", str(out), "--gatk3"]
    for n in NAMES:
        args += ["--name", n]
    args += [
        "--reference",
        str(REF),
        "--genotype-merge",
        "unsorted",
        "--merge-info-with-max-ac",
        "--min-callsets",
        "2",
    ]
    res = CliRunner().invoke(cli, args)
    assert res.exit_code == 0, res.output
    assert_same_as_golden(out, "golden_pipeline.gatk3.vcf")


# --- the normalising engine: --genotype-merge unsorted, --merge-info-with-max-ac

_MINI_HDR = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="a">\n'
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="d">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
)


def _mini(tmp_path, name, rows):
    p = tmp_path / f"{name}.vcf"
    p.write_text(
        _MINI_HDR
        + "".join(
            f"chr1\t{pos}\t.\tC\tT\t{q}\tPASS\t{info}\tGT\t{gt}\n"
            for pos, q, info, gt in rows
        )
    )
    return p


def _first_record(path):
    line = next(ln for ln in path.read_text().splitlines() if not ln.startswith("#"))
    return line.split("\t")


def test_unsorted_keeps_a_no_call_from_the_first_input(tmp_path):
    a = _mini(tmp_path, "a", [(100, 10, "AC=1;DP=10", "./.")])
    b = _mini(tmp_path, "b", [(100, 20, "AC=3;DP=20", "0/1")])
    prio = tmp_path / "prio.vcf"
    combine_vcfs([a, b], prio)
    assert _first_record(prio)[9] == "0/1"
    uns = tmp_path / "uns.vcf"
    combine_vcfs([a, b], uns, genotype_merge="unsorted")
    assert _first_record(uns)[9] == "./."


def test_max_ac_takes_info_from_the_highest_ac_input(tmp_path):
    a = _mini(tmp_path, "a", [(100, 10, "AC=1;DP=10", "0/1")])
    b = _mini(tmp_path, "b", [(100, 20, "AC=3;DP=20", "1/1")])
    first = tmp_path / "first.vcf"
    combine_vcfs([a, b], first, carry_info=True)
    rec = _first_record(first)
    assert rec[5] == "10" and "DP=10" in rec[7]
    mx = tmp_path / "max.vcf"
    combine_vcfs([a, b], mx, merge_info_with_max_ac=True)
    rec = _first_record(mx)
    assert rec[5] == "20" and "DP=20" in rec[7] and "AC=3" in rec[7]


def test_max_ac_tie_goes_to_the_first_input(tmp_path):
    a = _mini(tmp_path, "a", [(100, 10, "AC=2;DP=10", "1/1")])
    b = _mini(tmp_path, "b", [(100, 20, "AC=2;DP=20", "1/1")])
    out = tmp_path / "tie.vcf"
    combine_vcfs([a, b], out, merge_info_with_max_ac=True)
    assert "DP=10" in _first_record(out)[7]


def test_gatk3_gz_output_is_bgzf_and_indexed(tmp_path):
    import shutil

    if not (shutil.which("bgzip") and shutil.which("tabix")):
        pytest.skip("needs htslib")
    out = tmp_path / "merged.vcf.gz"
    combine_vcfs(
        INPUTS,
        out,
        names=NAMES,
        gatk3=True,
        reference=REF,
        genotype_merge="unsorted",
        merge_info_with_max_ac=True,
        min_callsets=2,
    )
    assert (tmp_path / "merged.vcf.gz.tbi").exists()
    import gzip

    text = tmp_path / "plain.vcf"
    text.write_bytes(gzip.open(out).read())
    assert data_lines(text) == data_lines(FIX / "golden_pipeline.gatk3.vcf")


def test_header_description_keeps_escaped_quotes(tmp_path):
    hdr = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth, \\"quoted\\"">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )
    row = "chr1\t100\t.\tT\tC\t10\tPASS\tDP=3\tGT\t0/1\n"
    for n in "ab":
        (tmp_path / f"{n}.vcf").write_text(hdr + row)
    out = tmp_path / "o.vcf"
    combine_vcfs(
        [tmp_path / "a.vcf", tmp_path / "b.vcf"], out, gatk3=True, min_callsets=1
    )
    assert 'Description="Depth, \\"quoted\\""' in out.read_text()
