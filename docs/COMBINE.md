# Combining call sets (`vcfclick combine`)

`vcfclick combine` merges multiple VCF **call sets** into one, recording
where each variant came from. It is a native reimplementation of GATK3's
`CombineVariants` — a tool GATK4 dropped and never replaced, and which
`bcftools` does not cover.

## `combine` vs `merge`

vcfclick has two ways to put several VCFs together; they solve different
problems.

| | `vcfclick merge` | `vcfclick combine` |
|---|---|---|
| Inputs | per-sample VCFs, **disjoint** samples | call sets that may **share** samples |
| Use case | assemble a trio / cohort joint VCF | two callers, or pre/post-filter, of the same cohort |
| Same sample in two inputs | rejected (must be disjoint) | resolved by **priority** |
| Engine | wraps `bcftools merge` | native (cyvcf2 read + write) |
| Provenance | — | `set=` INFO field per record |

Reach for `merge` to build one multi-sample VCF from separate samples;
reach for `combine` to reconcile multiple *callings* of the same samples.

## Usage

```bash
# Two callers over the same cohort; keep all sites, annotate origin.
vcfclick combine gatk.vcf.gz deepvariant.vcf.gz -o all.vcf.gz

# Consensus: keep only sites called by at least 2 of 3 callers.
vcfclick combine gatk.vcf.gz dv.vcf.gz strelka.vcf.gz \
    -o consensus.vcf.gz --min-callsets 2

# Name the sets explicitly (default names come from the filenames).
vcfclick combine a.vcf.gz b.vcf.gz -o out.vcf.gz --name gatk --name dv
```

Inputs must be on the same reference and **decomposed** (one ALT per
record) — the same precondition as ingest. Decompose first with
`bcftools norm -m -` if needed. (`--gatk3` and `--reference` lift this; see
below.)

A plain `.vcf` output is fully native (no external tools). A `.vcf.gz`
output is BGZF-compressed and tabix-indexed via htslib's `bgzip` /
`tabix`, so it is region-queryable and usable by parallel ingest like
every other vcfclick VCF — install htslib (`brew install htslib` /
`conda install -c bioconda htslib`) or write a plain `.vcf`.

## What it does

**Union of sites.** Every `(chrom, pos, ref, alt)` across all inputs
appears once in the output.

**`set=` provenance.** Each output record's INFO carries `set=`:

- `set=Intersection` — present in every input;
- `set=<name>` — present in one input only;
- `set=<a>-<b>` — present in those inputs (priority order), but not all.

**PRIORITY genotype resolution.** Input order is priority, highest
first. When a sample appears in more than one input at a site, its
genotype is taken from the **first input with a non-missing call**. A
`./.` in a higher-priority input does not win — the next input fills it.

**Consensus filter.** `--min-callsets N` keeps only sites present in at
least N inputs, for "agreed by ≥N callers" workflows.

## GATK CombineVariants parity options

Four opt-in options bring `combine` closer to GATK3 `CombineVariants`.
**All default to off — current behaviour is unchanged unless you ask.**

- **`--pass-only`** — count only PASS calls toward `--min-callsets`. A
  filtered input that has the site is still named `filterIn<name>` in
  `set=` (GATK's convention), e.g. `freebayes-gatk-gatk3-filterInoctopus`.
  A `FILTER` of `.` counts as PASS.
- **`--count-by site|allele`** (default `allele`) — `allele` counts inputs
  that call the exact `(POS, REF, ALT)`; `site` counts inputs with any
  record at the position and keeps every allele there (GATK's `--minimumN`
  semantics), each labelled with the position's `set=`.
- **`--reference REF.fa`** — split multi-allelic records and left-align +
  trim against the reference internally (the `bcftools norm -m - -f`
  equivalent), instead of refusing multi-allelic input. This lets raw
  caller VCFs — multi-allelic, padded indels, even missing `##contig`
  headers — combine with **no separate `bcftools` step**; contig order is
  taken from the reference `.fai`. Needs `pyfaidx` (`vcfclick[benchmark]`).
- **`--atomize`** (requires `--reference`) — after left-alignment, split
  complex alleles into primitives (SNPs plus at most one indel), so a
  caller that packs substitutions into one record converges with one that
  emits them separately. Note the division of labour: differently
  *anchored* spellings of the same indel converge via reference
  left-alignment (`--reference`); `--atomize` handles the packed-substitution
  case that left-alignment alone cannot.
- **`--carry-info`** — carry `QUAL`, `FILTER`, and `INFO` (with their
  `##INFO`/`##FILTER` header lines) from the highest-priority input that
  called each allele, rather than writing `QUAL`/`FILTER` as `.` and
  `INFO` as only `set=`.

```bash
# Combine raw caller VCFs with GATK-like semantics, no bcftools step:
vcfclick combine freebayes.vcf gatk.vcf gatk3.vcf octopus.vcf \
    -o merged.vcf --reference hg19.fasta \
    --min-callsets 2 --pass-only --count-by site --carry-info \
    --name freebayes --name gatk --name gatk3 --name octopus
```

## Reproducing GATK 3.8-1 CombineVariants (`--gatk3`)

`--gatk3` runs a separate engine that follows GATK 3.8-1's
`CombineVariants.map` and `GATKVariantContextUtils.simpleMerge` step by
step. Raw caller output goes in (multi-allelic records, padded indels,
non-PASS FILTER, VCFv4.3 headers) and the output matches GATK3 record for
record and header line for header line.

The clinical pipeline command

```bash
java -jar GenomeAnalysisTK.jar -T CombineVariants -R ref.fasta \
    --variant:freebayes S_freebayes.final.vcf --variant:gatk S_GATK.vcf \
    --variant:gatk3 S_GATK3.vcf --variant:octopus S_octopus.vcf \
    -o S_merged.vcf --genotypemergeoption UNSORTED --mergeInfoWithMaxAC --minimumN 2
```

is, in vcfclick:

```bash
vcfclick combine --gatk3 \
    S_freebayes.final.vcf S_GATK.vcf S_GATK3.vcf S_octopus.vcf \
    -o S_merged.vcf --reference ref.fasta \
    --name freebayes --name gatk --name gatk3 --name octopus \
    --genotype-merge unsorted --merge-info-with-max-ac --min-callsets 2
```

Input order is GATK's `--variant` order. `--reference` only supplies the
contig order and the `##contig` lines (GATK writes them from the reference
dictionary and ignores the contig lines of the inputs); it needs the `.fai`
and no pyfaidx. Without it the contig order comes from the input headers.
`.vcf.gz` output works as in the normalising engine. The other GATK options
map as `-priority a,b` to `--priority a,b` (with `--genotype-merge
prioritize`), `-filteredAreUncalled` to `--filtered-are-uncalled`, and
`--minimumN` to `--min-callsets`. `--gatk3` refuses `--pass-only`,
`--count-by site`, `--carry-info` and `--atomize`: GATK does none of them.

### What GATK 3.8-1 does (read from the source, then confirmed on its output)

1. **`--genotypemergeoption UNSORTED`.** The genotype of a sample is the one
   in the first record, in `--variant` order, that has the sample. A `./.`
   there wins over a later called genotype. UNSORTED skips the sort that
   PRIORITIZE does, so PRIORITIZE with the priority list equal to the input
   order gives identical output (`test_prioritize_in_input_order_equals_unsorted`).
   The normalising engine's default skips a `./.` and takes the next call;
   `--genotype-merge unsorted` switches it to GATK's rule.
2. **`--mergeInfoWithMaxAC`.** INFO is the INFO of the record with the highest
   INFO `AC` (the largest value of a multi-allelic list; the first record wins
   a tie; a record without `AC` never wins). `AC`, `AF` and `AN` are then
   recomputed from the merged genotypes, and INFO keys sort alphabetically.
   When no record has an `AC`, GATK3 writes no `set=` and no other INFO key,
   only the recomputed `AC`/`AF`/`AN` (Octopus and DeepVariant write no `AC`,
   so this happens). `--merge-info-with-max-ac` also works without `--gatk3`:
   it makes `--carry-info` take QUAL, FILTER and INFO from the input with the
   highest `AC` for the allele.
3. **`--minimumN 2`.** GATK counts *records* that start at the position,
   minus records whose FILTER is neither `PASS` nor `.`. It does not count
   inputs and does not look at alleles. Filtered records do not count:
   3.8-1 subtracts them (`vcs.size() - numFilteredRecords < minimumN`), and
   `-filteredAreUncalled` only decides whether they take part in the merge.
   One input with a SNP and an indel at one position therefore satisfies
   `--minimumN 2`. `--count-by site` in the normalising engine counts inputs
   with a record at the position, so it differs on that case, and on filtered
   records unless `--pass-only` is also set.
4. **No normalisation.** Records merge only when they start at the same
   position. Two spellings of one indel at different start positions stay two
   records (`--reference` in the normalising engine would join them).
   Multi-allelic records stay multi-allelic. The normalising engine without
   `--reference` also does not left-align, but it refuses multi-allelic input
   and writes one record per exact allele, where GATK3 puts the alleles of one
   position into one record (A>C and A>G from two inputs become `A C,G`).
   Records at one position are split
   by type (SNP, MNP, INDEL, SYMBOLIC, MIXED) with GATK's allele-subset rule;
   each type group is one output record, in that type order. `*` is typed by
   length like any allele. When REFs differ in length at one position the
   longest REF is kept and shorter records' ALTs are extended with its extra
   bases.
5. **Record fields.**
   - QUAL: the first non-missing QUAL, printed with two decimals and the
     trailing `.00` trimmed (`60.50`, `99`).
   - FILTER: `PASS` when any merged record is not filtered, else the sorted
     union; `.` when no record had FILTER applied.
   - ID: the comma-joined union of non-`.` IDs.
   - `set=`: `Intersection` (no filtered record and every input present),
     `FilteredInAll`, `ReferenceInAll`, else the names in record order,
     `filterIn<name>` for a filtered record.
   - FORMAT: `GT` then the other keys alphabetically; `GQ` is rounded and
     capped at 99; `GL` becomes `PL`; trailing missing values are trimmed;
     `AD` and `PL` are dropped on every sample when any record's allele list
     differs from the merged one.
6. **Header.** Lines merge in the iteration order of a Java `HashMap` keyed
   by input name (so `gatk` wins before `freebayes` for the names above),
   first line wins a description, an `Integer`/`Float` clash keeps the first
   line, a `Number` clash becomes `.`, standard INFO/FORMAT lines whose type
   or number is wrong are replaced by the standard line (a `Float` `GQ`
   becomes `Integer`), and `AC`/`AF`/`AN`/`set` lines are added. Lines are
   sorted as text. GATK3 stops on a FILTER, INFO or FORMAT key that no input
   header defines; so does `--gatk3`.

### Validation

Synthetic golden sets (`tests/fixtures/gatk3_pipeline`, GATK 3.8-1 under Java
8, docker `broadinstitute/gatk3:3.8-1`): four single-sample VCFs, one
reference, three GATK3 runs. Every data-carrying header line (`##fileformat`,
`##INFO`, `##FORMAT`, `##FILTER`, `##contig`, `#CHROM`) and every record match
in all three (48, 53 and 47 records; 0 differ).
`tests/test_combine_gatk3_pipeline.py` asserts it.

Public data, GIAB HG001 (NA12878) Garvan HiSeq exome, library NIST7035 lane 1,
hg19 `chr20`, reads sliced from the public BAM. Callers: GATK 4.6.0.0
HaplotypeCaller, GATK 3.8-1 HaplotypeCaller, FreeBayes 1.3.6, Octopus (the
`dancooke/octopus` image) and DeepVariant 1.6.1 (WES model). Octopus writes
`VCFv4.3`, which GATK3 refuses, so its header line was changed to `VCFv4.2`
for the GATK3 run only (`--gatk3` reads both).

| Region | Callers merged | Records GATK3 / vcfclick | Differing |
|---|---|---|---|
| chr20:3.6-3.9 Mb (18,568 reads) | freebayes, gatk4, gatk3, octopus, `--minimumN 2` | 76 / 76 | 0 |
| chr20:3.6-3.9 Mb | freebayes, gatk4, gatk3, deepvariant, `--minimumN 2` | 77 / 77 | 0 |
| chr20:3.6-3.9 Mb | freebayes, gatk4, deepvariant, `--minimumN 2` (2 of 3) | 75 / 75 | 0 |
| chr20 three windows, 1.2 Mb (82,988 reads) | freebayes, gatk4, gatk3, octopus, `--minimumN 2` | 350 / 350 | 0 |
| same | freebayes, gatk4, gatk3, deepvariant, `--minimumN 2` | 353 / 353 | 0 |
| same | freebayes, gatk4, deepvariant, `--minimumN 2` (2 of 3) | 346 / 346 | 0 |

The merged real sets include multi-allelic records (3 per set), a `*`
spanning-deletion allele (Octopus sets), 13 to 50 records with a `filterIn` label, DeepVariant
`RefCall` records, and 6 positions with two records of different type. Headers
match line for line (91 and 79 data-carrying lines, same order).

Rule comparison on the DeepVariant sets (GATK3 output). The 2-of-3 rule
(freebayes, gatk4, deepvariant) loses records the four-caller rule keeps
because gatk4 and gatk3 are one engine and count as two callers: 2 of 77 sites
at 3.6-3.9 Mb and 7 of 353 in the 1.2 Mb set, all labelled `gatk-gatk3` (6) or
`gatk-gatk3-filterIndeepvariant` (1). The 2-of-3 rule keeps nothing the
four-caller rule drops. Replacing Octopus by DeepVariant changes 1 record
(3.6-3.9 Mb) and 9 (3 gone, 6 new, 1.2 Mb), all sites called by freebayes and
one of the two.

To repeat the public-data check (Docker, `samtools`, about 20 MB of
download; the BAM is read over HTTPS, only the slice is fetched):

```bash
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome
B=$G/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam
# hg19 chr20 from UCSC, named chr20; BAM header cut to chr20
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz
gunzip -c chr20.fa.gz | sed 's/^>.*/>chr20/' > ref.fa && samtools faidx ref.fa
# ref.dict: one @SQ line (SN, LN, M5 = md5 of the sequence)
{ printf '@HD\tVN:1.5\tSO:coordinate\n@SQ\tSN:chr20\tLN:63025520\n'
  samtools view -H $B | grep '^@RG'
  samtools view $B chr20:3600000-3900000; } | samtools view -b - | samtools sort -o region.bam -
samtools index region.bam
freebayes -f ref.fa -r chr20:3600000-3900000 region.bam > freebayes.vcf
gatk HaplotypeCaller -R ref.fa -I region.bam -L chr20:3600000-3900000 -O gatk4.vcf
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ref.fa -I region.bam \
    -L chr20:3600000-3900000 -o gatk3.vcf
octopus -R ref.fa -I region.bam -T chr20:3600000-3900000 -o octopus.vcf
sed -i 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' octopus.vcf   # GATK3 only
run_deepvariant --model_type=WES --ref=ref.fa --reads=region.bam \
    --regions=chr20:3600000-3900000 --output_vcf=dv.vcf.gz
```

Then run the GATK3 command above and the `--gatk3` command on the same
inputs and compare the records.

### Remaining differences from GATK3

- GATK3 writes `##GATKCommandLine` and `##reference`; `--gatk3` writes
  neither and adds `##vcfclickCombine=gatk3-compatible`. Other header text
  (`##source`, `##ALT`, ...) is kept the way GATK3 keeps it.
- A malformed number (`DP=.`, `AC=.` in INFO, a non-integer genotype `DP`)
  stops GATK3 with a `NumberFormatException`; `--gatk3` stops with a
  `CombineError` naming the record.
- GATK3 rejects `VCFv4.3` inputs and INFO/FORMAT lines with `Source=` or
  `Version=`; `--gatk3` reads them. Structured header lines are written back
  with only `ID`, `Number`, `Type` and `Description`.
- Alleles are upper-cased on read, as htsjdk does.
- Not covered by a golden run: multi-sample inputs beyond a two-file check
  (samples sort alphabetically, a sample absent from every record at a site
  is `./.`), `--filteredRecordsMergeType` other than the default,
  `--multipleAllelesMergeType MIX_TYPES`, `-minimalVCF`,
  `-excludeNonVariants`, and UNIQUIFY. They are not implemented.
- The normalising engine (without `--gatk3`) still differs from GATK3 by
  design: it normalises and splits, refuses an INFO header type clash that
  GATK3 resolves (Octopus `MQ` Integer against GATK `MQ` Float), and writes a
  different header. On the DeepVariant 1.2 Mb set the documented legacy
  command (`--reference --pass-only --count-by site --carry-info`) shares
  329 of 353 records with GATK3 by position and alleles, and none of the
  shared records is identical in INFO.

## Output and limits

- Output is a fresh VCF with the union of input samples (priority
  order), `GT` per sample, the `set=` INFO field, and the `GQ`/`DP`/`AD`
  FORMAT fields.
- **Quality travels with the genotype.** When a sample's genotype comes
  from a given input, that same input's `GQ`/`DP`/`AD` come with it — so
  the depth and allele balance describe the call that was actually kept,
  and a combined VCF feeds straight into the trio quality gates (`db
  trio --min-gq/--min-dp/--min-ab`). Output FORMAT lists only the
  passthrough fields that some input actually carries (just `GT` if none
  do).
- A sample absent from every input at a site is `./.` with `.` for every
  FORMAT field, never silently `0/0` or a borrowed quality value —
  `combine` does not assert hom-reference where no input has a record.
- The site union is held in memory, which suits the typical "a few call
  sets" use; whole-genome many-caller merges of very large cohorts are
  not the target.

## Comparison with GATK3 CombineVariants

`combine` is a reimplementation of GATK3's `CombineVariants`. The two
produce the same site union, the same `set=` provenance, and the same
PRIORITIZE genotype resolution. The equivalent invocations:

```bash
# GATK3 (removed in GATK4) — needs Java 8 and a reference
java -jar GenomeAnalysisTK.jar -T CombineVariants -R ref.fasta \
    -V:first first.vcf.gz -V:second second.vcf.gz \
    -o combined.vcf \
    -genotypeMergeOptions PRIORITIZE -priority first,second

# vcfclick — no Java, no reference; input order is the priority list
vcfclick combine first.vcf second.vcf -o combined.vcf \
    --name first --name second
```

On two call sets that share sample `S2` (where `S2` is `0/1` in `first`
and `1/1` in `second` at chr1:202), both tools produce the same core —
`set=` and genotypes, with PRIORITIZE keeping `first`'s call for `S2`:

| pos | `set=` | S1 | S2 | S3 |
|---|---|---|---|---|
| 101 | `first` | 0/1 | 0/1 | ./. |
| 202 | `Intersection` | 1/1 | **0/1** | 0/1 |
| 303 | `Intersection` | 0/1 | 1/1 | 0/0 |
| 404 | `second` | ./. | 0/0 | 1/1 |

The `-V:name` binding GATK3 uses for the `set=` label is `--name` in
vcfclick; `-priority a,b` is just the input order. Without `--gatk3`, differences are representational: GATK3 also recomputes
`AC/AF/AN` INFO, whereas vcfclick emits `set=` plus the `GQ/DP/AD`
passthrough. `--gatk3` removes them (see above).

This equivalence is pinned by a test against **real** GATK 3.8-1 output:
`tests/test_combine_gatk3.py` runs `combine` on the checked-in call sets
and asserts it reproduces the frozen GATK3 golden file
(`tests/fixtures/gatk3/combine_priority.gatk3.vcf`, whose header records
the exact generating command). GATK3 needs Java 8 and is not run in CI;
the golden output is captured once and asserted against.
