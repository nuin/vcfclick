# Filter presets

A **preset** is a named set of `db trio` filter values: quality gates, a rarity
ceiling, which inheritance models to run. Instead of retyping
`--min-gq 30 --min-dp 20 --max-af 0.001 ...` you write
`db trio fam1 --proband C1 --preset strict-denovo`, and you can save your own,
hand them to a colleague, and cite the exact version in a methods section.

A preset holds filter values and nothing else: no sample, family, path or
database. The same preset applies to any trio.

> Research preview, not clinical. The four presets that ship with vcfclick are
> **the author's starting points, not validated clinical standards**. Every one
> is labelled `starting point, not validated` wherever it is shown, and the
> numbers need review by a clinical geneticist before anyone relies on them.

## Quick start

```bash
vcfclick preset list                              # the built-ins and your own
vcfclick preset show strict-denovo                # values, what it covers, what it does not

vcfclick db trio fam1 --proband C1 --preset broad-first-pass
vcfclick db trio fam1 --proband C1 --preset rare-recessive --max-af 0.005   # the flag wins

vcfclick preset save id-panel --from-preset rare-recessive --gnomad-max-af 0.002 \
    --description "ID/ASD recessive first pass"
vcfclick preset export id-panel -o id-panel.vcfclick-preset.json            # share it
vcfclick preset import id-panel.vcfclick-preset.json                        # a colleague
```

## The built-in presets

| Preset | Models | GQ / DP | Allele balance | `--max-af` | `--gnomad-max-af` | Other |
|---|---|---|---|---|---|---|
| `broad-first-pass` | all | 15 / 8 | 0.20-0.80 | 0.01 | 0.01 | warn-only |
| `strict-denovo` | denovo | 30 / 20 | 0.30-0.70 | 0.001 | 0.0001 | `strict_missing`, `require_bound`, `pass_only`, `parent_alt_screen` |
| `rare-recessive` | recessive, comphet | 20 / 10 | 0.25-0.75 | 0.01 | 0.01 | warn-only |
| `rare-dominant` | dominant, denovo | 25 / 12 | 0.25-0.75 | 0.001 | 0.0001 | warn-only |

`preset show NAME` prints, for each, two short paragraphs: **covers** (what the
numbers select, in words) and **does not cover** (what they cannot see). Before
there is a consequence filter, none of them excludes non-coding variants; none
looks for chrX in males, UPD or mosaic variants; `rare-dominant` keeps a variant
inherited from either parent because parent phenotype is not considered;
`rare-recessive` takes a single 1% ceiling, which is the usual carrier-frequency
limit but can still hide alleles near 1% popmax (some *CFTR* and *GJB2* alleles;
check against the gnomAD release in use). All four say that a variant absent from
the loaded gnomAD slice is **kept** as rare, and that the `INFO/AF` gate steps
aside when `INFO/AF` is the cohort's own frequency; `strict-denovo` also says that
a homozygous-alt proband call is not allele-balance gated and a missing allele
depth does not fail it, and that a parent 0/0 with fewer than two alt reads is
trusted whatever the fraction.

Built-ins are read-only. A built-in only changes when a release adds a model or
a filter that belongs in it; the change bumps its `revision`, adds a line to its
history with the new **fingerprint**, and is in the CHANGELOG. Anything else
gets a new name. The revision and fingerprint are in every result.

## Using a preset

```bash
vcfclick db trio fam1 --proband C1 --preset NAME
vcfclick db trio fam1 --proband C1 --preset NAME@a1b2c3d4e5f6   # fail if the content changed
vcfclick db trio fam1 --proband C1 --preset ./lab/id-panel.json  # a file (command line only)
```

**Precedence is: a flag you type > the preset > the default.** Flags and preset
can be given in any order. A flag equal to the preset's value still counts as
yours, and the result records it (`from: "flag"`).

- `--preset-ignore KEY` switches one preset key off for one run (repeatable):
  `--preset id-panel --preset-ignore gnomad_max_af`. The key is recorded in
  `preset.overridden` with `"flag": "off"`.
- `--preset NAME@FINGERPRINT` (or `--preset-fingerprint`, which also works with a
  path) makes the run fail if the preset's content is not the one you pinned. For
  a built-in whose revision has moved on, the message says which revision had that
  fingerprint.
- A key in the preset that this vcfclick does not know (a preset written by a
  newer version, stored with `preset import --allow-unknown-keys`) is an error;
  `--preset-skip-unavailable` applies the rest and says, in the output, that
  those filters were **not** applied.
- A key the command does not take (a trio-only gate on another command) is
  skipped, and because skipping it would widen the search the output carries a
  `WARNING`.
- The `--max-af` that a preset supplies is **not** counted as typed. If the VCF's
  `INFO/AF` turns out to be the cohort's own frequency (GATK-style family VCFs),
  the gate is skipped exactly as it is for the default (see
  [TRIO.md](TRIO.md#cohort-derived-infoaf)), and the run says so; typing
  `--max-af` forces the gate.
- `--category` takes a comma list: `--category recessive,comphet`. `all` and a
  single name behave as they always have. With a list, the text output has one
  summary line per model.

## Does the gate actually bind?

`db trio`'s quality gates are tolerant of missing data: a call with no GQ passes
`--min-gq`. So a preset called "strict" can silently be a no-op on a VCF that
carries no GQ, DP or allele depths, or when `INFO/AF` is the cohort's own
frequency and gnomAD is not loaded. A run with a preset therefore prints an
**effective-gate report**: for each gate, whether the field it reads exists in
this data. A gate **binds** when its field is present on at least 90% of the rows
it reads (GQ, DP and AD on the **candidate rows** of each selected model, below;
`INFO/AF` on the variants; gnomAD when it can judge the cohort, below).

"The rows it reads" is neither all of a sample's rows nor all of one genotype class.
A caller can fill AD on its inherited calls and leave it empty on the new ones, and
then AD is on 95% of a parent's 0/0 rows and on none of the 0/0 rows a de novo call
is made on; a caller can leave GQ empty on the proband's de novo hets and fill it on
his 200 inherited hom-alts, and GQ is on 95% of his non-reference rows. So each
GQ / DP / AD gate is judged on the **candidate rows** of every model the run
selects: the proband's and the parents' stored rows at the sites where that model's
genotype pattern holds (for `denovo`, the sites where the proband carries and both
parents are 0/0), counted over the same join `db trio` runs but before any quality
or rarity gate (one grouped count per model). Of those rows the gate reads these genotype classes, and the weakest
model decides:

| Model | Proband rows | Parent rows |
|---|---|---|
| `denovo` | 0/1 and 1/1 | 0/0 (both) |
| `recessive` | 1/1 | 0/1 (both) |
| `dominant`, `comphet` | 0/1 | 0/0, 0/1 or 1/1 |

`min_gq` and `min_dp` read all of those rows. `min_ab` reads the het rows (the
proband's; with `strict_missing` a parent's het call must carry depths too, so the
parents' 0/1 rows of `recessive`, `dominant` and `comphet`). `parent_alt_screen` reads
the parents' 0/0 rows, in `denovo` only. Only a model with no candidate rows of those
classes (no candidate at all, or no stored 0/0 rows because the database was ingested
without `--keep-reference`) is judged on the sample's stored rows of those classes, and
a role with none of those on its whole sample. The warning says which rows fell short
(`only 0% of father/mother 0/0 rows at denovo candidate sites carry AD`), and a
refusal's remedy names the field (`GQ is missing where the gates read it`). The
`gq` / `dp` / `ad` shares in `effective` stay the whole sample's, a fact about the
data and not the verdict; `bound` is the verdict of the gates that read the field
(`ad.bound` needs both `min_ab` and `parent_alt_screen`). The
`FILTER` gate (`pass_only`) binds when some variant in the database carries a
`FILTER` label other than `PASS` or unset: the column is null for both `PASS` and
`.`, so on a VCF that labels nothing the gate removes nothing, which the run says
(`pass_only_no_labels`) without refusing, because a VCF that was hard-filtered
before it was ingested is all `PASS` and the gate is right to remove nothing.

The text puts the `preset:` and `gates:` lines under the `trio:` line, then one
`WARNING:` line per problem, on stdout:

```
trio: proband=C1 father=F mother=M
preset: broad-first-pass (builtin, rev 1, cd3726cc3853, starting point)
  gates: category=all [preset]  min_gq=15 [preset]  ... max_af=0.01 [preset; skipped: INFO/AF is cohort-derived] ...
  WARNING: min_gq=15 cannot be enforced: only 31% of proband/father/mother rows carry GQ
  WARNING: NO population frequency filter applied: INFO/AF is this cohort's own frequency, so the --max-af gate was skipped, and no gnomAD data is loaded
```

| Warning code | When |
|---|---|
| `gate_unbound` | a requested gate's field is present on under 90% of its rows |
| `gate_partial` | the gate binds but not on every row (a note on stderr, `info`) |
| `thin_coverage` | the proband's median DP is under 1.5 x `min_dp` |
| `pass_only_no_labels` | `--pass-only` (or the preset's) is on and no variant carries a `FILTER` label other than `PASS` or unset, so it removes nothing (a warning, never a refusal: a VCF hard-filtered before ingest is all `PASS`) |
| `gnomad_not_loaded` | `--gnomad-max-af` was asked for and no gnomAD data is loaded |
| `gnomad_no_coverage` | gnomAD is loaded but cannot judge the cohort: it has records on chromosomes that hold under 90% of the cohort's variants (another chromosome naming, a slice elsewhere), or a record for under 10% of a sample of the cohort's own variants (a one-row test load, a slice of one gene, another genome build) |
| `gnomad_no_record` | after the run: the gnomAD filter found a record for none of the candidates, so all are kept as rare (novel variants, which is what de novo candidates are; the text says how much of the cohort gnomAD does reach) |
| `no_population_frequency_filter` | neither an `INFO/AF` gate nor gnomAD narrows by frequency |
| `model_blocked` | a model the run selected cannot run on this database: it has no stored hom-reference calls (`denovo`, `dominant` and `comphet` need `--keep-reference`) or no gene table is loaded (`comphet`). The model's count is 0 because of the database, not the genome. Named per model, in the text and in `effective.warnings`; not reported in a run limited to the gates it typed |
| `skipped_widening` | a preset key was skipped and skipping widens the search |
| `preset_max_af_not_applied` | the preset's `max_af` was skipped as cohort-derived (`info`) |

The same information is in the JSON `effective` object (below).

### Making a gate insist

Warnings do not stop a run, so a broad preset still works on data without GQ or
DP. Two opt-in gates turn the warnings into refusals:

- `--require-bound-gates` (`require_bound`): an error, with the remedy, when a
  requested gate cannot bind, when a model the run selected cannot run on this
  database (`model_blocked`), or when no population-frequency filter is in
  effect.
- `--strict-missing` (`strict_missing`): a missing GQ, DP or (for a het call) AD
  on the proband or a parent **counts as a fail**, not a pass, for each of those
  gates that is asked for: `min_gq` or `min_dp` above 0, an allele-balance band
  narrower than 0..1, the parent-allele screen on (a parent 0/0 without allele
  depths cannot be screened). A field nobody asked about is not held against a
  call, so `--strict-missing --min-gq 0 --min-dp 0 --min-ab 0 --max-ab 1` is the
  same run as without it. A requested gate that cannot bind is then an error too
  (a field missing everywhere would fail everything).

`strict-denovo` turns both on, so on a family VCF with no gnomAD loaded, or on a
database ingested without `--keep-reference` (where its one model cannot prove a
parent is 0/0 and would print `denovo candidates: 0`), it **refuses to run** and
says why. `--no-require-bound-gates` overrides both
refusals for one run; the warnings stay in the output. For that run,
`strict_missing` is **not applied to the gates that cannot bind**: a call
without GQ passes the `min_gq` gate rather than failing it, so the run reports
what the other gates leave and not "no candidates" for a reason that has nothing
to do with the genome. The warnings say so ("strict_missing is not applied to it
for this run"), and `gates.strict_exempt` lists the gates in the JSON:

```
$ vcfclick db trio fam1 --proband C1 --preset strict-denovo
Error: preset 'strict-denovo' needs every gate to bind, and these cannot in this data:
  - gnomad_max_af=0.0001 cannot be enforced: no gnomAD data is loaded (...)
  - NO population frequency filter applied: INFO/AF is this cohort's own frequency, ...
Remedy: load gnomAD for this cohort's chromosomes and genome build with `vcfclick annotations load-gnomad` and set --gnomad-max-af, or use broad-first-pass.
To run anyway, pass --no-require-bound-gates: the warnings stay in the output, and a call whose field is missing passes the gate instead of failing it.
```

gnomAD **binds** when it can judge the cohort, which is tested before the run on
the cohort's own variants, not on whether the store has a row:

- it has records on the chromosomes that hold at least 90% of the cohort's
  variants (a slice on other chromosomes, or under another naming, `1` against
  `chr1`, judges nothing), **and**
- the exact allele (chromosome, position, ref, alt) is in it for at least 10% of
  a sample of up to 5,000 of the cohort's variants (the same lookup the
  `--gnomad-max-af` filter does, so a one-row test load, a slice of one gene or a
  genome-build mismatch finds almost none). A real load has a record for most of a
  cohort's variants; the floor is low because a cohort's private variants are
  legitimately absent.

A gnomAD that does not bind applies no frequency filter that can be relied on,
and `strict-denovo` refuses it as it refuses an unloaded one (the JSON carries
`gnomad.cohort_variants_sampled` / `cohort_variants_with_record`, and the message
says the share). A slice of one gene on an exome therefore reads as not covering,
and `--no-require-bound-gates` runs it anyway.

What is left for after the run is a fact about the *candidates*: a de novo
candidate is novel by nature, so one that gnomAD has no record for is the
expected case and not a reason to refuse. The candidates the filter kept for want
of a record are counted (`gnomAD: no record for 6 of 8 candidates (kept as
rare); it has a record for 62.7% of a sample of the cohort's variants`, and
`gnomad.candidates_checked` / `candidates_without_record` in `effective`), and
when it is all of them the text ends with a `WARNING:` line and `effective`
carries `gnomad_no_record`. That is a warning and not a refusal because gnomAD
does judge the cohort, which the probe has already established.

`--parent-alt-screen` and `--pass-only` typed on their own are told about too:
on a VCF whose parents have no allele depths the screen screens nothing (every
parent 0/0 is trusted), and on a VCF whose `FILTER` column carries no label the
FILTER gate removes nothing, so the run prints that `WARNING:` line and carries an
`effective` object whose only warning is the one for the gate that was typed
(`gate_unbound` for `parent_alt_screen`, `pass_only_no_labels` for `pass_only`).
The other gates were not asked for, so they are not reported, and the object does
not claim anything about them either: gnomAD is not probed in this run, so `gnomad.bound`, `info_af.bound`,
`gq.bound`, `dp.bound` and `population_filter` are `null` (not examined; not the
same as bound or unbound), and `ad.bound` is the screen's own verdict. The shares
of rows, the medians and `gnomad.loaded` are facts about the data and stay.

### The other new gates

All default off, and with none of them (and no preset) the output and the SQL
are exactly what they were.

| Flag | Effect |
|---|---|
| `--pass-only` | keep only variants whose VCF `FILTER` is `PASS` or unset. It binds only where some variant carries another label: `FILTER` is stored as null for both `PASS` and `.`, so on a VCF that labels nothing it removes nothing, and the run says so (`pass_only_no_labels`) |
| `--parent-alt-screen` | de novo only: a parent's stored 0/0 counts as "does not carry it" only when its alt reads stay under `--max-parent-alt-frac` (default 0.05; a parent 0/0 with two or more alt reads and a higher fraction is excluded). The most common false de novo is a parent called 0/0 with real alt support. |

When a preset or `--parent-alt-screen` is on, each de novo candidate also shows
the parents' read counts (`parent_alt_reads` in the JSON; `parents' alt reads:
father 4/40, mother 0/38` in the text), so you see at a glance that a "de novo"
has 4 alt reads of 40 in the father.

## Making your own

```bash
vcfclick preset save id-panel --from-preset rare-recessive --max-af 0.003   # copy and change
vcfclick preset save mine --min-gq 25 --min-dp 12 --category recessive,comphet
vcfclick preset set id-panel max_af=0.002 strict_missing=true
vcfclick preset unset id-panel gnomad_max_af
vcfclick preset edit id-panel               # $EDITOR, validated before it is saved
vcfclick preset rename id-panel id-panel-v2
vcfclick preset rm id-panel-v2 --yes
```

`preset save` takes exactly the filter flags `db trio` takes (they are generated
from one registry), and saves only the ones you type. There is no inheritance:
`--from-preset` copies at save time. `--from-result RESULT.json` rebuilds the
preset from a saved `db trio --format json` result that used a preset (the
values that came from the preset or a flag), and says so when a built-in it used
has changed revision since. `--from-json FILE|-` reads `{"filters": {...}}`.

Presets live in `$VCFCLICK_HOME/presets/NAME.json`. Names are lowercase letters,
digits, `.`, `_` and `-` (up to 63), and may not end in `.json`, be a Windows
device name or a built-in's name. Writes are atomic. A file an export wrote
(`NAME.vcfclick-preset.json`) dropped into that directory is the same preset as
`NAME.json`: `save`, `import` and `rename` refuse to shadow it without `--force`,
and `--force` replaces that file rather than adding a second one beside it.

The file name is what a stored preset is called. A colleague's `q1.json` copied in
as `lab-q.json` still says `"name": "q1"` inside: it is listed with a
`name_mismatch` warning and used as `lab-q`, and `set`, `edit`, `rename` and
`export` all write `lab-q`, so an export imports back as `lab-q` and `edit` is not
asked to rename it. (A file with no `name` at all gets the same treatment.)

### The file

```json
{
  "schema": "vcfclick.preset",
  "schema_version": 1,
  "name": "id-panel",
  "title": "ID/ASD recessive first pass",
  "description": "...",
  "applies_to": ["trio"],
  "assay": "exome",
  "filters": {"category": ["recessive", "comphet"], "min_gq": 20, "max_af": 0.003},
  "meta": {"created_with": "0.15.0", "updated": "2026-10-21T09:00:00Z",
           "imported_from": null, "original_fingerprint": null,
           "approved_by": null, "approved_on": null, "locked": false}
}
```

- `filters` keys are the registered ones: `category`, `min_gq` (0-99), `min_dp`
  (0-10000), `max_af`, `min_ab`, `max_ab`, `gnomad_max_af` (0-1), `strict_missing`,
  `pass_only`, `parent_alt_screen`, `max_parent_alt_frac`, `require_bound`. A key
  not in that list is an error with a did-you-mean (and the vcfclick version
  recorded in `meta.created_with`, if different); `x_...` keys are kept and
  ignored. A `null`, a boolean where a number goes, `NaN`, an infinity (or a
  number too large for a float, anywhere in the file, `x_` keys included) or an
  out-of-range number is an error. `min_ab` above `max_ab` is checked on the resolved values, so
  a flag can repair it.
- **Never in a preset** (an error with the reason): which sample or family is
  analysed, output control (`format`, `limit`), paths and environment, safety
  checks (`parentage`, the X guard) and attestations. This holds for any key of
  the file, not only those under `filters`: `proband`, `family`, `db`, `path`,
  `sample` and the like at the top level, in `meta` or under `filters` are errors
  too, whatever the case or separator (`Patient`, `sample-id`, `sampleId`); other
  unknown top-level keys are kept for forward compatibility. A preset can only make a
  run stricter or equal, so nothing in it can switch off a check.
- `applies_to` (`trio`, `segregate`, `case`) and `assay` (`genome`, `exome`,
  `panel`) are hints for pickers and for the median-depth warning. A `title`,
  `description`, `applies_to` or `assay` that is `null` is the same as leaving it
  out; one of the wrong type is an error, and `preset list` shows the file with
  that error instead of failing.
- `schema_version` changes only when the envelope changes; adding a filter key
  never bumps it. A newer file is rejected with "made by a newer vcfclick; upgrade".
  Unknown envelope keys are preserved when a preset is edited.
- `meta.locked`, `approved_by` and `approved_on` are **advisory**: they are shown
  and echoed in results, and `set`/`unset`/`edit`/`rename`/`rm` refuse a locked
  preset without `--force`, but anyone with write access can edit the file. They
  are not signatures.

### Fingerprint and revision

The **fingerprint** is the first 12 hex digits of a SHA-256 over the
`schema_version` and the `filters` in canonical form (keys sorted, model names
sorted, numbers as written by Python's `repr`); `meta`, the title and the
description are left out. The same filters always give the same fingerprint,
however the file was written, and a changed value always changes it. Every result
produced with a preset carries its name, source, revision (built-ins), fingerprint
and the resolved values, which is what makes a result reproducible. A built-in's
revision and fingerprint history is in `preset show`.

## Sharing

```bash
vcfclick preset export id-panel -o id-panel.vcfclick-preset.json
vcfclick preset export a b --bundle -o lab-presets.json      # several in one file
vcfclick preset import id-panel.vcfclick-preset.json [--as NAME] [--force]
```

`export -o DIR` writes `DIR/NAME.vcfclick-preset.json`; in a preset directory
that file is the preset `NAME` (a preset cannot itself be named `...vcfclick-preset`,
and if `NAME.json` sits beside it, `NAME.json` is the one used and the other is
listed with a `duplicate_name` warning).

Export leaves out patient-specific keys unless `--with-patient-data`. Import is
all or nothing: every file is parsed and validated, every name checked (a
built-in name and an existing name are refused without `--as` / `--force`)
before anything is written. It records only the source file's name in
`meta.imported_from`, and the fingerprint of what was imported in
`meta.original_fingerprint`.

A lab can keep validated presets in a git folder and put it on
`VCFCLICK_PRESET_PATH` (`:`-separated; `;` on Windows). Those presets are listed as
`shared` and are read-only; a preset of yours with the same name wins, with a
note. Pin them with `--preset NAME@FINGERPRINT`. `preset path` prints where
presets are looked up, `preset validate FILE` checks a file (exit 1 on error).

A preset file is data. It is parsed and validated, never executed; its values
reach the database only as typed numbers and names; no error message quotes the
file's content. A preset is always addressed by a name: a value with a path
separator or a leading dot never reaches the file system for `set`, `unset`,
`edit`, `rename`, `rm`, `export` or `save --from-preset`, so `preset rm ../x`
cannot leave the presets directory. `--preset ""` (an unset shell variable) is an
error, not "no preset".

## What the output carries

With a preset, `db trio --format json` adds two top-level objects after
`models`; with none, the document is byte-identical to what it was.

- `preset`: `name`, `source` (`builtin`, `user`, `shared`, `file`), `revision`,
  `status` (`starting_point` for built-ins), `schema_version`, `fingerprint`,
  `pinned`, `meta`, `path` (a file name only), `requested` (the preset's
  filters), `resolved` (per filter: `value`, `from` = `preset` / `flag` /
  `default`, and `bound` for a gate), `overridden` (preset keys a flag or
  `--preset-ignore` changed), `skipped`, `models_not_run`, `covers`,
  `does_not_cover`, `notes`.
- `effective`: `rows`, `gq`, `dp`, `ad` (the share of rows carrying each, per
  sample, and `bound`), `median_dp`, `info_af` (`state` is `gates.af_gate`),
  `filter_column` (`present`, and `non_pass_variants`: the variants whose `FILTER`
  is a label other than `PASS` or unset), `gnomad` (`loaded`, `bound`, `candidates_checked`,
  `candidates_without_record`, `cohort_variants_sampled`,
  `cohort_variants_with_record`), `population_filter` and `warnings` (`{code,
  level, key, message}`). `gnomad.bound` is given only for a filter that was
  asked for and probed: with no `--gnomad-max-af` (or the key switched off with
  `--preset-ignore`) a loaded store was never looked at and `bound` is `null`, and
  `candidates_checked` / `candidates_without_record` are `null` while no gnomAD is
  loaded (over an empty store every candidate would be "without a record"). In a
  run limited to some gates (`--parent-alt-screen` typed alone) every `bound` of a
  gate outside the scope, and `population_filter`, is `null`: it was not
  examined.

`gates` also reports `strict_missing` (and, after it, `strict_exempt`: the gates
a `--no-require-bound-gates` run did not apply it to), `pass_only`,
`parent_alt_screen` and `max_parent_alt_frac`, but only when they are on. See
[TRIO.md](TRIO.md#machine-readable-output).

## Not in this release

Gene lists and regions as a preset scope, the `--funnel` / `--near-miss`
diagnostics and `preset diff`, the web and MCP preset endpoints, and the desktop
picker come in later steps; none of them is needed to use or share a preset
today. A duo (a proband with one parent declared) is not analysed by `db trio`
yet, so there is no duo warning either; it comes with duo analysis. gnomAD
homozygote counts are not stored, so there is no `gnomad_max_hom`.
