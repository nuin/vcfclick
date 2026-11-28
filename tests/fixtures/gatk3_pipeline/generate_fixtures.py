"""Write the four synthetic single-sample caller VCFs for the GATK3 golden test.

Run from this directory: python3 generate_fixtures.py

The reference (ref.fasta) is random sequence with three planted repeats.
Every record below targets one CombineVariants behaviour, named in the
comment on its row. REF alleles are read from ref.fasta, so the VCFs stay
consistent with it.
"""

from pathlib import Path

HERE = Path(__file__).parent


def load_ref() -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    for line in (HERE / "ref.fasta").read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


REF = load_ref()


def b(chrom: str, pos: int, n: int = 1) -> str:
    return REF[chrom][pos - 1 : pos - 1 + n]


def alt(chrom: str, pos: int, k: int = 0) -> str:
    """The k-th base (ACGT order) that differs from the reference at pos."""
    return [x for x in "ACGT" if x != b(chrom, pos)][k]


# --- headers ---------------------------------------------------------------
CONTIGS = "##contig=<ID=chr1,length=1400>\n##contig=<ID=chr2,length=400>\n"
COMMON_INFO = (
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">\n'
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">\n'
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">\n'
)
HEADERS = {
    "freebayes": (
        "##fileformat=VCFv4.2\n##source=freebayes v1.3.6\n"
        + COMMON_INFO
        + '##INFO=<ID=AB,Number=1,Type=Float,Description="Allele balance">\n'
        '##INFO=<ID=TYPE,Number=A,Type=String,Description="Variant type">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
        '##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference count">\n'
        '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate count">\n'
        '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods">\n'
    ),
    "gatk": (
        "##fileformat=VCFv4.2\n"
        '##FILTER=<ID=LowQual,Description="Low quality">\n'
        + COMMON_INFO
        + '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">\n'
        '##INFO=<ID=QD,Number=1,Type=Float,Description="Qual by depth">\n'
        '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">\n'
        "##source=HaplotypeCaller\n" + CONTIGS
    ),
    "gatk3": (
        "##fileformat=VCFv4.2\n"
        '##FILTER=<ID=LowQual,Description="Low quality">\n'
        + COMMON_INFO
        + '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">\n'
        '##INFO=<ID=FS,Number=1,Type=Float,Description="Strand bias">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">\n'
        + CONTIGS
    ),
    "octopus": (
        "##fileformat=VCFv4.2\n"
        '##FILTER=<ID=q10,Description="Quality below 10">\n'
        '##FILTER=<ID=LowQual,Description="Low quality">\n'
        + COMMON_INFO
        + '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">\n'
        '##INFO=<ID=MP,Number=1,Type=Float,Description="Model posterior">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n' + CONTIGS
    ),
}
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"

# per-caller FORMAT strings
FMT = {
    "freebayes": "GT:GQ:DP:AD:RO:AO:GL",
    "gatk": "GT:AD:DP:GQ:PL",
    "gatk3": "GT:AD:DP:GQ:PL",
    "octopus": "GT:GQ:DP:PS",
}

rows: dict[str, list[tuple]] = {k: [] for k in HEADERS}


def add(
    caller,
    chrom,
    pos,
    ref,
    alts,
    gt,
    info,
    qual="50",
    flt="PASS",
    rid=".",
    ad=None,
    pl=None,
    dp=30,
    gq=60,
    gl=None,
    ps=".",
):
    """One record. alts is a comma string of ALT alleles; ad/pl default to a
    plausible shape for the allele count."""
    n_alt = len(alts.split(",")) if alts != "." else 0
    nal = n_alt + 1
    ad = ad or ",".join(["15"] + ["15"] * n_alt)
    npl = nal * (nal + 1) // 2
    pl = pl or ",".join(str(i * 10) for i in range(npl))
    if caller == "freebayes":
        gl = gl or ",".join(f"-{i + 1}.5" for i in range(npl))
        ao = ",".join(["15"] * max(n_alt, 1))
        cell = f"{gt}:{gq}:{dp}:{ad}:15:{ao}:{gl}"
    elif caller in ("gatk", "gatk3"):
        cell = f"{gt}:{ad}:{dp}:{gq}:{pl}"
    else:
        cell = f"{gt}:{gq}:{dp}:{ps}"
    rows[caller].append(
        (chrom, pos, rid, ref, alts, qual, flt, info, FMT[caller], cell)
    )


def snp(caller, pos, gt, ac, k=0, chrom="chr1", **kw):
    ref = b(chrom, pos)
    info = (
        f"AC={ac};AF={ac / 2:g};AN=2;DP={kw.get('dp', 30)}"
        if ac is not None
        else f"DP={kw.get('dp', 30)}"
    )
    extra = kw.pop("extra", "")
    if extra:
        info += ";" + extra
    kw.pop("dp", None)
    fdp = kw.pop("fdp", 30)
    add(caller, chrom, pos, ref, alt(chrom, pos, k), gt, info, dp=fdp, **kw)


# --- 1. site in all four callers, GTs agree, AC ties: INFO from first
for c, ex in (
    ("freebayes", "AB=0.5"),
    ("gatk", "MQ=60;QD=12.5"),
    ("gatk3", "MQ=59;FS=1.2"),
    ("octopus", "MP=30.2"),
):
    snp(
        c,
        100,
        "0/1",
        1,
        extra=ex,
        dp={"freebayes": 20, "gatk": 31, "gatk3": 29, "octopus": 25}[c],
        rid=".",
    )

# --- 2. GT disagreement + max-AC tie between gatk and octopus (both AC=2)
snp("freebayes", 150, "0/1", 1, extra="AB=0.4", qual="77.7")
snp("gatk", 150, "1/1", 2, extra="MQ=60;QD=20", qual="900.5")
snp("gatk3", 150, "0/1", 1, extra="MQ=58", qual="120")
snp("octopus", 150, "1/1", 2, extra="MP=40", qual="500")

# --- 3. site in one caller only (minimumN drops it)
snp("gatk3", 200, "0/1", 1, extra="MQ=60")

# --- 4. site in exactly two callers
snp("freebayes", 250, "0/1", 1, extra="AB=0.52")
snp("octopus", 250, "0/1", 1, extra="MP=22")

# --- 5. site in three callers
snp("gatk", 300, "0/1", 1, extra="MQ=60")
snp("gatk3", 300, "0/1", 1, extra="MQ=60;FS=0.5")
snp("octopus", 300, "1/1", 2, extra="MP=31")

# --- 6. one caller filtered, two PASS (merged FILTER PASS, filterIngatk in set)
snp("freebayes", 350, "0/1", 1)
snp("gatk", 350, "0/1", 1, flt="LowQual", extra="MQ=30")
snp("gatk3", 350, "0/1", 1, extra="MQ=60")

# --- 7. one PASS + one filtered: counted 1 of 2, dropped by minimumN 2
snp("freebayes", 400, "0/1", 1)
snp("gatk3", 400, "0/1", 1, flt="LowQual")

# --- 8. one PASS + two filtered: dropped
snp("freebayes", 450, "0/1", 1)
snp("gatk", 450, "0/1", 1, flt="LowQual")
snp("octopus", 450, "0/1", 1, flt="q10")

# --- 9. two PASS + two filtered: kept, FILTER PASS, filterIn labels for two
snp("freebayes", 455, "0/1", 1)
snp("gatk", 455, "0/1", 1, flt="LowQual", extra="MQ=30")
snp("gatk3", 455, "0/1", 1)
snp("octopus", 455, "0/1", 1, flt="q10")

# --- 10. all four filtered: dropped
for c, f in (
    ("freebayes", "LowQual"),
    ("gatk", "LowQual"),
    ("gatk3", "LowQual"),
    ("octopus", "q10"),
):
    snp(c, 460, "0/1", 1, flt=f)

# --- 11. FILTER '.' vs PASS: merged PASS; and both '.' : '.'
snp("freebayes", 470, "0/1", 1, flt=".")
snp("gatk", 470, "0/1", 1)
snp("freebayes", 480, "0/1", 1, flt=".")
snp("gatk", 480, "0/1", 1, flt=".")

# --- 12. multi-allelic: fb biallelic, gatk two alts, octopus other alt
snp("freebayes", 500, "0/1", 1, extra="AB=0.5")
add(
    "gatk",
    "chr1",
    500,
    b("chr1", 500),
    f"{alt('chr1', 500, 0)},{alt('chr1', 500, 1)}",
    "1/2",
    "AC=1,1;AF=0.5,0.5;AN=2;DP=40;MQ=60",
)
snp("octopus", 500, "1/1", 2, k=1, extra="MP=10")

# --- 13. same two alts, opposite order: allele lists disagree -> PL/AD stripped
add(
    "gatk",
    "chr1",
    550,
    b("chr1", 550),
    f"{alt('chr1', 550, 0)},{alt('chr1', 550, 1)}",
    "1/2",
    "AC=1,1;AF=0.5,0.5;AN=2;DP=40",
)
add(
    "gatk3",
    "chr1",
    550,
    b("chr1", 550),
    f"{alt('chr1', 550, 1)},{alt('chr1', 550, 0)}",
    "1/2",
    "AC=1,1;AF=0.5,0.5;AN=2;DP=38",
)

# --- 14. multi-allelic in first caller, higher-AC biallelic in another
add(
    "freebayes",
    "chr1",
    600,
    b("chr1", 600),
    f"{alt('chr1', 600, 0)},{alt('chr1', 600, 1)}",
    "1/2",
    "AC=1,1;AF=0.5,0.5;AN=2;DP=33;AB=0.5",
)
snp("octopus", 600, "1/1", 2, k=0, extra="MP=44")

# --- 15. first caller no-call, second called
add(
    "freebayes",
    "chr1",
    650,
    b("chr1", 650),
    alt("chr1", 650),
    "./.",
    "AC=0;AF=0;AN=0;DP=5",
    gq=0,
    ad="0,0",
    gl="-1,-1,-1",
)
snp("gatk", 650, "0/1", 1)

# --- 16. same-start deletions of different length (REF extension)
add(
    "freebayes",
    "chr1",
    760,
    b("chr1", 760, 2),
    b("chr1", 760),
    "0/1",
    "AC=1;AF=0.5;AN=2;DP=30;TYPE=del",
)
add(
    "gatk",
    "chr1",
    760,
    b("chr1", 760, 3),
    b("chr1", 760),
    "0/1",
    "AC=1;AF=0.5;AN=2;DP=30;MQ=60",
)

# --- 17. same indel spelled two ways: left-aligned (gatk, gatk3) vs right (freebayes, octopus)
for c in ("gatk", "gatk3"):
    add(
        c,
        "chr1",
        799,
        b("chr1", 799, 2),
        b("chr1", 799),
        "0/1",
        "AC=1;AF=0.5;AN=2;DP=30;MQ=60",
    )
for c in ("freebayes", "octopus"):
    add(
        c,
        "chr1",
        804,
        b("chr1", 804, 2),
        b("chr1", 804),
        "0/1",
        "AC=1;AF=0.5;AN=2;DP=30",
    )
# the same with a dinucleotide repeat
for c in ("gatk", "gatk3"):
    add(
        c,
        "chr1",
        899,
        b("chr1", 899, 3),
        b("chr1", 899),
        "0/1",
        "AC=1;AF=0.5;AN=2;DP=22",
    )
for c in ("freebayes", "octopus"):
    add(
        c, "chr1", 901, b("chr1", 901, 3), b("chr1", 901), "1/1", "AC=2;AF=1;AN=2;DP=22"
    )

# --- 18. SNP and insertion at one position, different callers: two records, both kept
snp("freebayes", 820, "0/1", 1)
add(
    "gatk3",
    "chr1",
    820,
    b("chr1", 820),
    b("chr1", 820) + "CC",
    "0/1",
    "AC=1;AF=0.5;AN=2;DP=30",
)

# --- 19. MNP vs SNP at one position
add(
    "freebayes",
    "chr1",
    850,
    b("chr1", 850, 2),
    alt("chr1", 850, 0) + alt("chr1", 851, 0),
    "0/1",
    "AC=1;AF=0.5;AN=2;DP=30;TYPE=mnp",
)
snp("gatk", 850, "0/1", 1)
snp("octopus", 850, "0/1", 1)

# --- 20. no AC on any record: GATK3 drops the set= label with max-AC
snp("gatk3", 950, "0/1", None)
snp("octopus", 950, "0/1", None)

# --- 21. IDs
snp("freebayes", 1000, "0/1", 1, rid="rs11")
snp("gatk", 1000, "0/1", 1, rid="rs22")
snp("gatk3", 1000, "0/1", 1, rid="rs11")
snp("octopus", 1000, "0/1", 1, rid=".")

# --- 22. QUAL: first non-missing wins; formatting
snp("freebayes", 1050, "0/1", 1, qual=".")
snp("gatk", 1050, "0/1", 1, qual="60.5")
snp("gatk3", 1050, "0/1", 1, qual="1234.567")
snp("freebayes", 1060, "0/1", 1, qual="99")
snp("gatk", 1060, "0/1", 1, qual="12.3456")
snp("freebayes", 1070, "0/1", 1, qual="0.004")
snp("gatk", 1070, "0/1", 1, qual="3")

# --- 23. hom-ref genotype in first caller; AC recomputed from genotypes
snp("freebayes", 1100, "0/0", 1)
snp("gatk", 1100, "0/1", 1)

# --- 24. phased genotypes
snp("freebayes", 1150, "0|1", 1)
snp("gatk", 1150, "1|0", 1)

# --- 25. FreeBayes GL to PL; GATK PL passthrough
snp("freebayes", 1200, "0/1", 1, gl="-0.3,-0.1,-2.5")
snp("gatk", 1200, "0/1", 1, pl="30,0,90")

# --- 26. max-AC picks a later caller
snp("freebayes", 1250, "0/1", 1)
snp("gatk3", 1250, "1/1", 2, extra="MQ=60")
snp("octopus", 1250, "1/1", 2, extra="MP=12")

# --- 27. one filtered of four, set label order
snp("freebayes", 1300, "0/1", 1)
snp("gatk", 1300, "0/1", 1)
snp("gatk3", 1300, "0/1", 1)
snp("octopus", 1300, "0/1", 1, flt="q10")

# --- 28. flag INFO key and ordering of INFO keys
snp("gatk", 1350, "0/1", 1, extra="DB;MQ=60;QD=3")
snp("gatk3", 1350, "0/1", 1, extra="MQ=60")

# --- 29. second contig, ordering
snp("freebayes", 100, "0/1", 1, chrom="chr2")
snp("gatk3", 100, "0/1", 1, chrom="chr2")
snp("octopus", 50, "0/1", 1, chrom="chr2")
snp("gatk", 50, "0/1", 1, chrom="chr2")
snp("freebayes", 300, "0/1", 1, chrom="chr2")

# --- 30. one input with a SNP and an indel at one position: two records count
# toward minimumN (GATK counts records, not inputs), so both are kept
snp("gatk", 880, "0/1", 1)
add(
    "gatk",
    "chr1",
    880,
    b("chr1", 880, 2),
    b("chr1", 880),
    "0/1",
    "AC=1;AF=0.5;AN=2;DP=30;MQ=60",
)

# --- 31. a type group that is filtered in every record: FILTER union, set label
snp("freebayes", 1010, "0/1", 1)
snp("gatk", 1010, "0/1", 1)
add(
    "gatk3",
    "chr1",
    1010,
    b("chr1", 1010),
    b("chr1", 1010) + "CT",
    "0/1",
    "AC=1;AF=0.5;AN=2;DP=30",
    flt="LowQual",
)
add(
    "octopus",
    "chr1",
    1010,
    b("chr1", 1010),
    b("chr1", 1010) + "CT",
    "0/1",
    "AC=1;AF=0.5;AN=2;DP=30",
    flt="q10",
)

# --- 32. GQ rounding and the 99 cap, DP missing
snp("freebayes", 1020, "0/1", 1, gq="150.3", fdp=".")
snp("gatk", 1020, "0/1", 1)
snp("gatk", 1030, "0/1", 1, gq=100)
snp("gatk3", 1030, "0/1", 1, gq=100)
snp("freebayes", 1040, "0/1", 1, gq="12.5")
snp("gatk", 1040, "0/1", 1)
snp("freebayes", 1045, "0/1", 1, gq="12.4")
snp("gatk", 1045, "0/1", 1)

# --- 33. reference-only records (no ALT)
add("gatk3", "chr1", 1380, b("chr1", 1380), ".", "0/0", "DP=30")
snp("octopus", 1380, "0/1", 1)
add("gatk3", "chr1", 1385, b("chr1", 1385), ".", "0/0", "DP=30")
add("octopus", "chr1", 1385, b("chr1", 1385), ".", "0/0", "DP=28")

# --- 34. mixed SNP + insertion record vs a SNP and an insertion from others
add(
    "freebayes",
    "chr2",
    200,
    b("chr2", 200),
    f"{alt('chr2', 200, 0)},{b('chr2', 200)}CC",
    "1/2",
    "AC=1,1;AF=0.5,0.5;AN=2;DP=30",
)
snp("gatk", 200, "0/1", 1, chrom="chr2")
add(
    "gatk3",
    "chr2",
    200,
    b("chr2", 200),
    b("chr2", 200) + "CC",
    "0/1",
    "AC=1;AF=0.5;AN=2;DP=30",
)


# --- 35. SNP inside the REF span of the deletions at 760
snp("gatk", 761, "0/1", 1)
snp("gatk3", 761, "1/1", 2)

# --- 36. multi-allelic with a longer REF elsewhere: alleles extend and reorder
add(
    "freebayes",
    "chr1",
    765,
    b("chr1", 765, 2),
    f"{b('chr1', 765)},{b('chr1', 765, 2)}AT,{b('chr1', 765)}G",
    "1/2",
    "AC=1,1,0;AF=0.5,0.5,0;AN=2;DP=30",
)
add(
    "gatk",
    "chr1",
    765,
    b("chr1", 765, 3),
    b("chr1", 765),
    "0/1",
    "AC=1;AF=0.5;AN=2;DP=30;MQ=60",
)


def main() -> None:
    names = {
        "freebayes": "S_freebayes.final.vcf",
        "gatk": "S_GATK.vcf",
        "gatk3": "S_GATK3.vcf",
        "octopus": "S_octopus.vcf",
    }
    for caller, fname in names.items():
        recs = sorted(rows[caller], key=lambda r: (r[0], r[1]))
        # callers emit coordinate order; stable on ties
        body = "".join("\t".join(str(x) for x in r) + "\n" for r in recs)
        (HERE / fname).write_text(HEADERS[caller] + COLUMNS + body)
        print(fname, len(recs))


if __name__ == "__main__":
    main()
