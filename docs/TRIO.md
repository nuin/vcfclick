# Trio / family analysis

vcfclick can report candidate variants under Mendelian inheritance
models for a trio (affected child + two parents), with the genotype
quality gates and population-frequency rarity filter that rare-disease
analysis uses in practice.

This is **candidate filtering, not variant calling**. It works
downstream of your caller, over an already-genotyped cohort. Quality
gating is only as strong as the FORMAT fields your VCF carries
(`gq`, `dp`, `ad`); where they are absent the gates pass through.

## The pipeline

```bash
# 1. Trios usually arrive as three per-sample VCFs. Combine them into
#    one joint multi-sample VCF (decomposed, ingest-ready).
vcfclick merge proband.vcf.gz father.vcf.gz mother.vcf.gz -o trio.vcf.gz

# 2. Ingest as ONE cohort. Use --keep-reference so a parent's
#    confident 0/0 calls are stored (see "Why --keep-reference" below).
vcfclick db create fam1
vcfclick db ingest fam1 trio.vcf.gz \
    --cohort trio --ingest-id fam1 --keep-reference

# 3. Load the pedigree (standard PED/FAM).
vcfclick db ped fam1 fam1.ped

# 4. Report candidates.
vcfclick db trio fam1 --proband PROBAND
vcfclick db trio fam1 --proband PROBAND --category denovo
vcfclick db trio fam1 --proband PROBAND --category recessive,comphet
vcfclick db trio fam1 --proband PROBAND --preset rare-recessive
```

If your trio was already joint-called into one multi-sample VCF, skip
step 1 and ingest it directly.

## The PED file

Standard six-column PED/FAM, whitespace-delimited:

```
# family  individual  paternal  maternal  sex  phenotype
FAM1  CHILD   FATHER  MOTHER  1  2
FAM1  FATHER  0       0       1  1
FAM1  MOTHER  0       0       2  1
```

- paternal/maternal `0` = no parent declared (a founder, or the unknown side
  of a duo)
- sex: 1 = male, 2 = female
- phenotype: 1 = unaffected, 2 = affected, `0` or `-9` = unknown (stored as
  NULL, which is *unknown*, not unaffected)

`db ped` validates that every individual and named parent is a real
sample in the cohort, and that the pedigree has no cycle, no one who is
their own parent and no parent in a different family, before it changes
anything. Loading a PED file replaces the prior pedigree; the family column
lets one file (or several `--merge` loads) hold more than one family, and a
family can span generations. See [Families, generations and
probands](#families-generations-and-probands).

## Inheritance models

`db trio` reports, per model, candidate variants where the proband and
parent genotypes fit the pattern. Genotype encoding: `gt` 0 = hom-ref,
1 = het, 2 = hom-alt.

| Model | Pattern | Needs `--keep-reference`? |
|---|---|---|
| `denovo` | proband carries; both parents provably hom-ref | yes |
| `recessive` | proband hom-alt; both parents heterozygous carriers | no |
| `dominant` | proband het; exactly one parent carries, other hom-ref | yes |
| `comphet` | two rare proband hets in the same gene, one from each parent | yes |

`--category all` (default) prints per-model counts; a specific
`--category` prints the variants.

### Compound heterozygous

`--category comphet` finds genes where the proband carries **two** rare
heterozygous variants in *trans* — one inherited from the father, one
from the mother — so both copies of the gene are hit (a recessive
mechanism that per-variant filters miss). Each candidate variant is a
rare proband het with a clear parent-of-origin (one parent carries, the
other is provably hom-ref, hence `--keep-reference`); a gene is reported
only when it has both a paternal and a maternal one. Results are grouped
**per gene**, not per variant.

Because genes live in the annotation store (separate from the cohort),
this needs gene coordinates loaded first:

```bash
vcfclick annotations load          # GENCODE gene coordinates
vcfclick db trio fam1 --proband CHILD --category comphet
```

Parent-of-origin here is genotype-based (no read-backed phasing), so it
is honest candidate filtering: it cannot distinguish true *trans* from
two variants that happen to be on the same parental copy. The annotation
store path can be overridden with `VCFCLICK_ANNOTATIONS_DB`.

## Quality gates

Defaults follow common rare-disease practice and are all tunable:

| Flag | Default | Meaning |
|---|---|---|
| `--min-gq` | 20 | genotype quality (≈99% accurate) |
| `--min-dp` | 10 | read depth |
| `--min-ab` / `--max-ab` | 0.25 / 0.75 | het allele balance `ad_alt/(ad_ref+ad_alt)`, expected ≈0.5 |
| `--max-af` | 0.01 | keep variants with the VCF's own `info_AF` at or below this (a site exactly at the threshold is kept; see [below](#cohort-derived-infoaf)) |
| `--gnomad-max-af` | off | additionally drop candidates whose **gnomAD popmax** AF exceeds this; needs `vcfclick annotations load-gnomad`. Variants absent from the loaded gnomAD slice are kept as rare. |
| `--strict-missing` | off | a missing GQ, DP or (for a het call) AD on the proband or a parent counts as a **fail**, not a pass, for each gate that is asked for (`--min-gq` / `--min-dp` above 0, an allele-balance band narrower than 0..1, the parent-allele screen); a field nobody asked about is not held against a call |
| `--pass-only` | off | keep only variants whose VCF `FILTER` is `PASS` or unset (on a VCF that labels no variant otherwise it removes nothing, and says so) |
| `--parent-alt-screen` / `--max-parent-alt-frac` | off / 0.05 | de novo only: a parent's stored 0/0 is trusted only while its alt reads stay under the fraction (two or more alt reads and a higher fraction exclude the site); shows the parents' read counts |
| `--require-bound-gates` | off | an error, instead of a warning, when a requested gate cannot bind in this data (see [Presets](PRESETS.md#does-the-gate-actually-bind)) |

`--max-af` filters on the cohort VCF's `info_AF` (often missing on small
cohorts); `--gnomad-max-af` adds a true population-frequency filter from
gnomAD, the defensible way to require a variant be rare in the general
population. Load gnomAD AFs for your regions of interest first:

```bash
# pull a region from the public gnomAD bucket and load it
tabix https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr7.vcf.bgz \
    chr7:117480000-117670000 | bgzip > cftr.gnomad.vcf.gz
vcfclick annotations load-gnomad cftr.gnomad.vcf.gz
vcfclick db trio fam1 --proband CHILD --gnomad-max-af 0.001
```

Gates are lenient on NULL: a VCF lacking `gq`/`dp`/`ad` is not silently
zeroed out — the gate simply doesn't filter on the missing field. (That is
the CLI's rule. The generic `gq >= 20 AND dp >= 10` an SQL query would write
drops those rows instead; to reproduce a `db trio` count by hand, use
`(gq IS NULL OR gq >= 20) AND (dp IS NULL OR dp >= 10)`.)

Whether a gate does anything depends on the data: with no GQ, DP or allele
depths the tolerant gates remove nothing. A run that uses a [filter
preset](PRESETS.md) (`--preset NAME`) says, per gate, whether the field it
reads exists in this data, and the opt-in `--strict-missing` and
`--require-bound-gates` turn that into a refusal. The four built-in presets are
starting points, not validated standards.

### The `--max-af` boundary

`variants.info_AF` is a 32-bit float on chDB and a `REAL` on DuckDB. A site
whose AF is exactly `--max-af` (say `0.001`) is stored as
`0.0010000000474…`, which is greater than the double `0.001`, so the plain
comparison `info_AF <= 0.001` dropped it on chDB and kept it on DuckDB. The
gate now compares with `--max-af + 1e-7`: far below any real frequency
step, above float32 rounding, so a site **at** the threshold is kept on both
backends. This is the only change to the SQL of the four classic models; it
can add a candidate sitting exactly on the threshold, never remove one.

### Cohort-derived `INFO/AF`

`--max-af` assumes `INFO/AF` is a *population* frequency. A VCF joint-called
from just the family (GATK-style) instead carries the callset's own
frequency, and a private het in a trio is `1/6 = 0.167`. The default
`--max-af 0.01` would then silently remove every private candidate, which
is most of them.

`db trio` checks for this from the called-genotype counts stored at ingest
(`ac_called / an_called`): if `INFO/AF` equals the cohort's own frequency at
90% or more of at least 20 sites, it is treated as cohort-derived. "Equals"
is relative, to 1% of the cohort's frequency, plus 0.0005 for a frequency
of 0.01 or more (a writer that prints three decimals turns 1/6 into
0.167). It has to be relative: in a cohort of a hundred or more samples
the cohort's own frequency of a rare site is 0.005 or less, and a real
population frequency such as 0.0001 must not be mistaken for it. Then:

- without `--max-af` on the command line (a `--max-af` that a preset supplies
  does not count as typed), the AF gate is **skipped** and a
  `note:` on stderr says so (and points at `--gnomad-max-af`, which is a real
  population frequency);
- with `--max-af` given explicitly, the gate is **applied as asked** and the
  note says it may remove real candidates.

The outcome is reported as `gates.af_gate` in the JSON (see below):
`applied`, `cohort_derived_suspected` (looks cohort-derived, applied because
`--max-af` was explicit), `skipped_cohort_derived`, or `unavailable` (the
probe could not run: fewer than 20 sites carry `INFO/AF` and called counts,
or the database was ingested before those columns existed; the gate runs as it
always did). The web Trio tab does the same: leave *max AF* blank for the
default and the auto-skip, or enter a value to always apply it.

## Why `--keep-reference` (and the de-novo caveat)

vcfclick's `genotypes` table is **sparse** — only non-reference calls
are stored, so a sample absent at a site is `0/0` (reference) **or**
`./.` (no-call), indistinguishable. De novo means "the child carries a
variant neither parent has," which requires the parents to be
*confidently* hom-reference. Absence alone can't prove that: a no-call
parent is not evidence.

`--keep-reference` additionally stores confident `0/0` calls (gt=0) at
variant sites, while still dropping no-calls. Then de-novo's join
(`f.gt = 0 AND m.gt = 0`) matches only when both parents have a stored
reference call — a no-call parent has no row, so that site is correctly
**excluded** from de novo rather than falsely reported.

Without `--keep-reference`, `db trio` still runs but de-novo and
dominant return nothing (no stored parent reference rows) and print a
note to re-ingest with the flag. Recessive works either way.

Rigorous de-novo *calling* (PL-based Bayesian posteriors using genotype
likelihoods) is future work; this is honest candidate filtering.

## Validation against the GIAB benchmark trio

The inheritance models are tested against **real** genotypes from the
Genome in a Bottle Ashkenazi trio — HG002 (son), HG003 (father), HG004
(mother) — taken from the NIST v4.2.1 GRCh38 small-variant benchmark, not
hand-built fixtures. In a CFTR slice the benchmark genotypes reproduce
every model: two recessive sites (proband hom-alt, both parents het), a
paternal-origin and maternal-origin het that together make CFTR a
compound-het candidate gene, and no de novo. The golden fixture
(`tests/fixtures/giab/cftr_trio.vcf.gz`) records its GIAB source URLs in
the header.

De novo is validated harder, against an **independent** truth: GIAB's
high-confidence BED. A one-time scan of a 5 Mb chr20 region found HG002
carries 50 variants absent from both parents (naive de novos), but only
**4** have both parents confidently hom-ref per the BED; the other 46
have a no-call parent. The checked-in test exercises a representative
7-site subset — vcfclick's `--keep-reference` de novo recovers exactly
the 4 confident sites and excludes the 3 sampled no-call sites, and
`bcftools +mendelian2` independently flags the same 4. See
[Validation](VALIDATION.md).

## Families, generations and probands

A pedigree can hold several families, and a family can span generations:
a child's parents may have parents of their own in the pedigree. Each
individual is one row of `pedigree` (`family_id`, `father_id`, `mother_id`,
`sex`, `affected`); `0` means no parent declared. Every individual and
every named parent must be a sample ingested under the same `ingest_id`.
Each proband is analysed as its own trio against its own parents, so a de
novo call in one sibling and not the other is the expected picture; shared
candidates have to be read across the reports.

### Declaring a family without a PED file

```bash
# a trio, or siblings who share the same parents
vcfclick db ped fam1 --proband SIB1 --proband SIB2 --father DAD --mother MUM

# grandparents, an unaffected sibling and an affected half sibling, in a
# named family
vcfclick db ped fam1 --family famA --proband C1 --father F --mother M \
    --parents F=PGF,PGM --parents C2=F,M --parents H1=F,M2 \
    --unaffected C2 --affected H1 --sex C1=male
```

| Flag | Meaning |
|---|---|
| `--proband S` | an affected child of `--father` and `--mother` (repeatable) |
| `--father F`, `--mother M` | the parents of the probands; `0` for a parent that is not known (a *duo*), never both |
| `--family FID` | the family every listed sample joins (default `fam1`) |
| `--parents CHILD=F,M` | the parents of any other member, for grandparents and half siblings (repeatable; `0` for an unknown side) |
| `--affected S`, `--unaffected S`, `--unknown S` | set the phenotype of any listed sample (repeatable) |
| `--sex ID=male\|female\|unknown` | set the sex of any listed sample (repeatable); `--proband-sex` does the same for a single proband |
| `--merge` | change only the people this call names and keep everybody else as stored (see below) |

Probands are stored as affected and **everyone else as unknown**, parents
included: a parent's status is an assertion (it decides reduced
penetrance), so only a PED file or an explicit `--affected`/`--unaffected`
makes one. A father is stored male and a mother female unless `--sex` says
otherwise (reported as a warning). A PED file carries each individual's sex
and phenotype itself. (Under `--merge`, these are what a person who is not
stored yet gets; someone already stored keeps what is stored unless a flag
says otherwise.)

Everything is checked before the stored pedigree is touched: a sample that is
not in the cohort, an individual who is their own parent, a cycle of
ancestors, a parent in a different family, or a sample moving to another
family under `--merge` is refused, and nothing is changed. The error names
the family it is in. Under `--merge` only a problem the call would *introduce*
refuses it: a cycle that an earlier vcfclick let into another family (or
into this one, somewhere the call does not touch) is reported as a `warning:`
and left as it is, so one broken family does not make the rest unreadable
to edit; `db families` lists what is wrong.

### Replacing and merging

A PED file is a whole pedigree: loading one replaces what is stored for the
`ingest_id`. The flags describe **a few people** of one family, so `db ped`
treats them more carefully:

- `--merge` with flags rewrites **only the people the flags name**, and only
  what the flags say about them. Everybody else, in the same family or
  another, stays exactly as stored. A person the flags name without saying
  more (a father given only by `--father`) keeps their own parents, sex and
  phenotype, so `--family famA --merge --proband C1 --father F --mother M`
  on a family of 13 rewrites C1, F and M as they were and changes nothing.
  What the flags do say wins: the proband's parents, `--parents`,
  `--affected`/`--unaffected`/`--unknown`, `--sex`. Those last flags may also
  name a person who is stored but not otherwise in the call
  (`--unaffected K2`), to change just that. A sample that is stored in one
  family cannot join another this way; a person cannot be removed with flags
  either (use a PED file).
- `--merge` with a PED file replaces the families the file lists **as
  wholes** and keeps every other family. A member of such a family that the
  file does not list is removed, and `db ped` says so
  (`removed from famB: K2, XB`).
- Without `--merge`, a flag call that would delete individuals it does not
  name **fails** and says which, and what to add (`--family FID --merge`
  when the samples it names already live in one family). The desktop apps
  show nothing from a call that succeeds, so a quiet replace would lose
  people without anyone noticing. There is no exception for a small family:
  an unaffected sibling the call does not name is as protected as a
  grandparent, and so is the stored proband when a call names another one.
  To change the people a call names, add `--merge`; to replace everything on
  purpose, load a PED file. A call that names every stored individual
  replaces them and removes nobody.
- When a flag call does replace the stored phenotype of an individual that was
  marked affected or unaffected with a different one it was not asked for
  (probands default to affected, others to unknown), a `warning:` on stderr
  names them.

### Probands and `db families`

A **proband** is an analysable child, that is a child with both parents
declared and in the data (a trio), who is **not marked unaffected** and who
**either has no children of their own or is affected**. Two consequences: in
a three-generation family the parents in the middle are trios but not
probands unless they are affected, and a child marked unaffected (a healthy
sibling) is not a proband. One rule is used everywhere: `db families`,
`db trio` without `--proband`, and the proband list of the web UI.

```bash
vcfclick db families fam1                  # text
vcfclick db families fam1 --format json    # what other programs read
```

`db families` is read-only. It lists each family's members by generation
(founders are generation 0), the probands, every trio with its siblings and
why a child is not a proband, and any problem it finds in the pedigree: an
error (a cycle, a parent in another family, an individual as their own
parent), a warning (a father recorded as female, an individual that is not a
sample in the data) or an info (a duo). It reads tolerantly, so a pedigree
that an older vcfclick let in is described, not refused, and the command exits
0 either way. JSON shape:

```json
{"ingest_id": "i1",
 "families": [{"family_id": "famA",
               "members": [{"sample_id": "C1", "family_id": "famA", "father": "F",
                            "mother": "M", "sex": "male", "affected": "affected",
                            "genotyped": true, "generation": 2}],
               "probands": ["C1"],
               "trios": [{"child": "C1", "father": "F", "mother": "M",
                          "kind": "trio", "is_proband": true,
                          "siblings": [{"sample": "C2", "kind": "full"}]}]}],
 "problems": [{"family_id": "famA", "level": "info", "code": "duo",
               "message": "..."}]}
```

`sex`, `affected`, `father`, `mother` and `generation` are `null` when unknown
(a member on or below a cycle has no generation). `problems[].code` is one of
`cycle`, `self_parent`, `same_parent`, `cross_family_parent`, `parent_sex`,
`not_genotyped`, `parent_without_row` or `duo`.

### Running `db trio` over a family

```bash
vcfclick db trio fam1                      # the probands of every family
vcfclick db trio fam1 --family famB        # the probands of one family
vcfclick db trio fam1 --all-trios          # every child with both parents declared
vcfclick db trio fam1 --proband SIB2       # one
```

Reports are ordered by family and then sample id; with more than one family
the text output names the family above each trio. A note on stderr lists the
children that are not probands (and why) and the duos that were skipped, and
points at `--all-trios`. `--family` and `--all-trios` cannot be combined with
`--proband`. If a family has trios but no proband (every child is marked
unaffected, say) the error says so and names `--all-trios`.

Before probands were defined, `db trio` without `--proband` ran every child
with both parents declared. A PED file with an unaffected sibling, or with a
middle generation, now gets fewer reports by default, and `--all-trios`
restores the old set. A one-generation family whose children are not marked
unaffected gets exactly the output it always did.

### Siblings

Siblings share a declared parent: **full** siblings share both,
**half** siblings one (`half_paternal` or `half_maternal`). Each trio of a
`db trio --format json` document carries its siblings (see below), so a
consumer reading shared candidates across reports knows who is related how.
A half sibling is evidence about a variant only if it is a child of the
parent the variant came from; a half sibling who lacks the allele says nothing
about an allele from the other parent.

### Duos

A child with one parent declared (`--father 0`, or a PED row with a `0`
parent) is a **duo**. It is stored, shown by `db families` (`kind: "duo"`, an
`info` problem), and not offered as a proband: duo analysis is a separate
feature. `db trio` skips duos without `--proband` and says so; naming one
still gives the "missing a parent" error.

### Older databases

A database that predates the pedigree table says `no pedigree loaded` (from
`db families`, `db trio` and the web UI) instead of failing in SQL.

## Natural language

The MCP server teaches an LLM the pedigree table, the inheritance
models, and the de-novo sparse caveat, so you can ask trio questions in
English and get the SQL back to inspect. See [MCP](MCP.md).

## Machine-readable output

`vcfclick db trio NAME --proband ID --format json` reports every model in one
document: per model the candidate count, whether it is blocked by a missing
`--keep-reference` ingest, and the candidates with their trio genotypes,
population AF and, when those annotations are loaded, gene, gnomAD popmax and
ClinVar significance. Compound hets come grouped by gene with the paternal and
maternal hits. `--limit` caps candidates per model (default 1000). The desktop
apps use this. Without `--proband` the output is
`{"trios": [...], "families": [...]}`: one such document per proband (see
[Probands](#probands-and-db-families)), ordered by family and sample id, and
the ids of the families they belong to.

A family can be declared without a PED file:
`vcfclick db ped NAME --proband HG002 --father HG003 --mother HG004`
(`--proband` repeats for siblings; see [Declaring a
family](#declaring-a-family-without-a-ped-file) for the rest of the flags).

### The document

```json
{
  "trio": {"ingest_id": "fam1", "proband": "CHILD", "father": "FATHER", "mother": "MOTHER",
           "family_id": "FAM1", "kind": "trio",
           "siblings": [{"sample": "SIB", "kind": "full"}]},
  "keep_reference": true,
  "gates": {"min_gq": 20, "min_dp": 10, "max_af": 0.01, "min_ab": 0.25,
            "max_ab": 0.75, "gnomad_max_af": null, "af_gate": "applied"},
  "annotations": {"genes": true, "clinvar": false, "gnomad": false},
  "models": {
    "denovo": {"count": 1, "blocked": false, "blocked_reason": null,
               "truncated": false, "candidates": [ ... ]},
    "recessive": { ... }, "dominant": { ... },
    "comphet": {"count": 0, "blocked": false, "blocked_reason": null,
                "needs_genes": false, "truncated": false, "genes": []}
  }
}
```

- `trio.family_id` is the proband's family in the pedigree (`null` if the
  proband has no row of its own), `trio.kind` is `trio` or `duo`, and
  `trio.siblings` lists everyone who shares a declared parent, each as
  `{"sample", "kind"}` with `kind` one of `full`, `half_paternal` or
  `half_maternal` (an empty list for an only child). They come after the four
  original keys, in that order.
- `models` lists the models of `--category` (a name, `all`, or a comma list
  such as `recessive,comphet`), in the order `denovo`, `recessive`, `dominant`,
  `comphet`.
- `blocked` is true when the model cannot be decided because the database has
  no stored hom-reference calls (see *Why `--keep-reference`*). `blocked_reason`
  says why a model has no result, and is `null` when nothing is wrong. It comes
  from one closed vocabulary:

  | `blocked_reason` | Meaning | Reported by |
  |---|---|---|
  | `needs_keep_reference` | no stored hom-ref rows to prove a parent is 0/0 (re-ingest with `--keep-reference`) | `denovo`, `dominant`, `comphet` |
  | `needs_genes` | gene annotations are not loaded (`vcfclick annotations load`) | `comphet` |

  The vocabulary also reserves `needs_phase`, `needs_parents`, `needs_sex`,
  `needs_ad`, `needs_parent_ad`, `no_autosomes`, `no_x_data`,
  `parent_sex_conflict` and `build_conflict` for models that are not in this
  release; a reader should map an unknown value to a generic "blocked" message
  rather than fail. When both apply, `needs_keep_reference` wins (fix the first
  requirement first), and `comphet` keeps `blocked: false` for a missing gene
  annotation, which it reports through `needs_genes` and `blocked_reason` alone.
- `gates.af_gate` is the outcome of the `--max-af` gate: `applied`,
  `cohort_derived_suspected`, `skipped_cohort_derived` or `unavailable` (see
  [Cohort-derived `INFO/AF`](#cohort-derived-infoaf)).
- `gates` also reports `strict_missing` (then `strict_exempt`, the gates
  `--no-require-bound-gates` did not apply it to, only when there are any),
  `pass_only`, `parent_alt_screen` and `max_parent_alt_frac` (after `af_gate`),
  but only for a gate that is on; the seven original keys carry the values that
  took effect after a preset and the flags were merged.
- With a preset (`--preset`) the document ends with a `preset` object (what ran:
  name, source, revision, fingerprint, the filters requested and the value each
  filter took with where it came from, what a flag overrode, what was skipped)
  and an `effective` object (per gate, how much of the data carries its field,
  and the warnings). `--strict-missing`, `--require-bound-gates`,
  `--parent-alt-screen` and `--pass-only` without a preset add `effective` alone
  (for `--parent-alt-screen` or `--pass-only` alone its warnings are about that
  gate only, and the `bound` of every other gate, `gnomad.bound` and
  `population_filter` are `null`: not examined, which is not the same as
  bound). With none
  of them neither key is present. See [Presets](PRESETS.md#what-the-output-carries).
- With a preset or `--parent-alt-screen`, each `denovo` candidate also carries
  `parent_alt_reads`: `{"father": {"ref", "alt"}, "mother": {...}}`, the allele
  depths of the parents' stored 0/0 rows (`null` counts where the call has no AD).
- A comphet gene is `{"gene", "paternal": [...], "maternal": [...]}`. A reader
  must treat `paternal` and `maternal` as **optional**: a model that cannot tell
  parent of origin leaves them out, and may add `tier`, `k_carried` and
  `k_absent`.

### Compatibility policy

The JSON is a contract with scripts and with the desktop apps, which decode it
with strict typed decoders and show a blank panel if a field they require goes
missing or changes type. So:

- **Additive only.** A key already in the output is never removed, renamed or
  retyped, and the relative order of existing keys does not change. New keys
  are `snake_case`, lower case, with digits only inside a word (never leading a
  segment after an underscore).
- **Old keys are byte-stable.** `tests/test_trio_goldens.py` captures the CLI
  output of the four classic models on `related.vcf`, `trio.vcf.gz` and the
  GIAB trios, text and JSON, on both backends; after removing the keys the test
  whitelists as new (`blocked_reason`, `gates.af_gate`, `trio.family_id`,
  `trio.kind`, `trio.siblings`, the top-level `families` of a no-`--proband`
  list, and the `preset` and `effective` objects of a run that uses a preset)
  it must match byte for byte. A change that needs an old key to move has
  to say so in its commit and regenerate those goldens deliberately. The
  complete document of one proband, new keys included, is pinned separately
  (`tests/goldens/trio_full/`), and the `db families` document too
  (`tests/goldens/families/`).
- **Readers ignore keys they do not know,** and treat every field a model might
  one day leave out as optional.
- Every model object has an integer `count`, a boolean `blocked` and a boolean
  `truncated`; every candidate has a `genes` list of strings and integer
  genotypes; `annotations` is `{string: bool}`. `tests/test_trio_contract.py`
  asserts this on real output and is reusable by later models.

Two reader rules are reserved now, so consumers can be written once, although
no model emits either today: a parent genotype in a candidate may be `null`
(that parent was not sequenced), and `-1` as a genotype means *a sequenced
parent with no stored row at the site*. The stored `genotypes.gt` is never
negative (see [SCHEMA.md](SCHEMA.md#identity--genotype-encoding)), so a `-1`
here is a display state, not a genotype: **an exporter must not do arithmetic
on it or write it into a VCF**; render it as `./.`.

### Web UI and MCP

The Trio tab of `vcfclick web` offers the same four models (its list is
generated from the CLI's, so the two cannot drift), `/api/trio` returns
`af_gate` and a `note` for the keep-reference and cohort-AF cases, and the MCP
briefing gives an LLM the same model names, the NULL-tolerant gate wording,
the cohort-derived-AF warning and how families, `affected = NULL` and
probands work. `/api/meta` lists the probands by the same rule as `db trio`
and a `families` list (per family its probands and its trios with siblings).
See [WEB](WEB.md) and [MCP](MCP.md).

