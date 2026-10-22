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
check against the gnomAD release in use).

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
it reads (GQ, DP and AD on the stored rows of the proband and both parents;
`INFO/AF` on the variants; gnomAD when any is loaded).

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
| `gnomad_not_loaded` | `--gnomad-max-af` was asked for and no gnomAD data is loaded |
| `gnomad_no_coverage` | gnomAD is loaded but has records on chromosomes that hold under 90% of the cohort's variants (another chromosome naming, a slice elsewhere): it judges nothing |
| `gnomad_no_record` | after the run: the gnomAD filter found a record for none of the candidates (a slice that does not reach them, another genome build); all are kept as rare |
| `no_population_frequency_filter` | neither an `INFO/AF` gate nor gnomAD narrows by frequency |
| `skipped_widening` | a preset key was skipped and skipping widens the search |
| `duo_denovo` | one parent: the de novo answer is `possible_denovo`, weaker evidence |
| `preset_max_af_not_applied` | the preset's `max_af` was skipped as cohort-derived (`info`) |

The same information is in the JSON `effective` object (below).

### Making a gate insist

Warnings do not stop a run, so a broad preset still works on data without GQ or
DP. Two opt-in gates turn the warnings into refusals:

- `--require-bound-gates` (`require_bound`): an error, with the remedy, when a
  requested gate cannot bind or when no population-frequency filter is in effect.
- `--strict-missing` (`strict_missing`): a missing GQ, DP or (for a het call) AD
  on the proband or a parent **counts as a fail**, not a pass, for each of those
  gates that is asked for: `min_gq` or `min_dp` above 0, an allele-balance band
  narrower than 0..1, the parent-allele screen on (a parent 0/0 without allele
  depths cannot be screened). A field nobody asked about is not held against a
  call, so `--strict-missing --min-gq 0 --min-dp 0 --min-ab 0 --max-ab 1` is the
  same run as without it. A requested gate that cannot bind is then an error too
  (a field missing everywhere would fail everything).

`strict-denovo` turns both on, so on a family VCF with no gnomAD loaded it
**refuses to run** and says why. `--no-require-bound-gates` overrides both
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

gnomAD **binds** when it is loaded and has records on the chromosomes that hold
at least 90% of the cohort's variants. A slice on other chromosomes, or under
another naming (`1` against `chr1`), applies no frequency filter at all, and
`strict-denovo` refuses it as it refuses an unloaded one. Whether a slice that
covers the right chromosomes has a record for any candidate (a build mismatch, a
slice of one gene) cannot be known before the query: the run says so afterwards,
with a `WARNING:` line at the end of the text (and `gnomad_no_record` plus
`gnomad.candidates_checked` / `candidates_without_record` in `effective`), and
the candidates the filter kept for want of a record are counted
(`gnomAD: no record for 6 of 8 candidates (kept as rare)`). This one is a
warning, not a refusal, because the run is already done.

### The other new gates

All default off, and with none of them (and no preset) the output and the SQL
are exactly what they were.

| Flag | Effect |
|---|---|
| `--pass-only` | keep only variants whose VCF `FILTER` is `PASS` or unset |
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
device name or a built-in's name. Writes are atomic.

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
  `sample` and the like at the top level or in `meta` are errors too (other
  unknown top-level keys are kept for forward compatibility). A preset can only make a
  run stricter or equal, so nothing in it can switch off a check.
- `applies_to` (`trio`, `segregate`, `case`) and `assay` (`genome`, `exome`,
  `panel`) are hints for pickers and for the median-depth warning.
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
  `filter_column`, `gnomad` (`loaded`, `bound`, `candidates_checked`,
  `candidates_without_record`), `population_filter` and `warnings` (`{code,
  level, key, message}`).

`gates` also reports `strict_missing` (and, after it, `strict_exempt`: the gates
a `--no-require-bound-gates` run did not apply it to), `pass_only`,
`parent_alt_screen` and `max_parent_alt_frac`, but only when they are on. See
[TRIO.md](TRIO.md#machine-readable-output).

## Not in this release

Gene lists and regions as a preset scope, the `--funnel` / `--near-miss`
diagnostics and `preset diff`, the web and MCP preset endpoints, and the desktop
picker come in later steps; none of them is needed to use or share a preset
today. gnomAD homozygote counts are not stored, so there is no `gnomad_max_hom`.
