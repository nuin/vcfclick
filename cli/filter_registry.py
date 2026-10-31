"""The one place a filter is declared.

Every value `db trio` (and, later, `db segregate` and the case command) takes
as a gate or a scope is a `FilterKey` here: its key in a preset file, its
Click parameter and flag (three different words more often than not: file
`require_bound`, parameter `require_bound`, flag `--require-bound-gates`),
type and range, default, the step that owns it and the commands it applies
to. The Click options of a command are generated from this table
(`cli.preset_options.gate_options`), so a typed flag and a preset value are
converted and checked by the same rules, and `tests/test_filter_registry_drift.py`
fails when a Click option of a filter-bearing command is neither registered
here nor listed in `NEVER_PRESET`, when a Click default differs from the
registry's, or when a registered key is a parameter of no command. A later
step registers its own keys in its own change; that test is what makes it.

There are no reserved keys: a key this binary does not know is an error with a
did-you-mean, never a silent no-op.

Pure: no click, no SQL, no file system. `resolve_filters_from_values` is the
precedence rule (flag over preset over default) for a caller that has no
Click context (the web route, the MCP tools, `preset save`).
"""

from __future__ import annotations

import difflib
import math
import re
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from typing import Any

# The commands a filter can apply to. `case` is the multi-sample case command
# (`db case`), not an alias for every multi-sample command.
COMMANDS: tuple[str, ...] = ("trio", "segregate", "case")

# The default for `db trio --max-af` (and the web API's max_af).
DEFAULT_MAX_AF = 0.01

VALUE = "value"  # a filter the command applies, one Click parameter
META = "meta"  # a descriptive companion (a gene-list label): no parameter

# Value types whose value is a list of tokens (gene symbols or ids, regions).
LIST_TYPES: tuple[str, ...] = ("genes", "regions")
# Caps on a list in a preset (and on what a command line or a file may give).
MAX_GENES = 20_000
MAX_REGIONS = 5_000
GENE_TOKEN_RE = re.compile(r"^[A-Za-z0-9._@:-]{1,40}$")
# 1-based inclusive `chr7:117100000-117400000` (or `7:...`).
REGION_RE = re.compile(r"^(?:chr)?([0-9]+|X|Y|M|MT):([0-9]+)-([0-9]+)$")
MAX_TEXT = 200  # a META text value (a gene-list label)


class _Derived:
    """A default that is computed from other values, not a constant. A preset
    must not freeze it; the echo shows `derived`."""

    def __repr__(self) -> str:
        return "DERIVED"


DERIVED = _Derived()


def _category_choices() -> tuple[str, ...]:
    # Late: cli.trio_analysis reads DEFAULTS from this module.
    from cli.trio_analysis import CATEGORY_CHOICES

    return CATEGORY_CHOICES


@dataclass(frozen=True)
class FilterKey:
    file_key: str  # the key in a preset file's `filters`
    param: str | None  # the Click parameter name (None for kind=meta)
    flag: str | None  # the declaration, `--min-gq` or `--strict-missing/--no-...`
    type: str  # int | float | bool | category
    default: Any  # a constant, or DERIVED
    owner: str  # the plan step that registered it
    commands: tuple[str, ...]
    kind: str = VALUE
    lo: float | None = None  # inclusive bounds of a preset value (int/float)
    hi: float | None = None
    presettable: bool = True
    patient_specific: bool = False  # stripped on export (HPO terms, ...)
    # Skipping this key would make the search wider, so a command that cannot
    # apply it says so loudly.
    widening: bool = True
    help: str = ""
    choices: Callable[[], tuple[str, ...]] | None = None  # category / text tokens
    # A list key (`genes`, `regions`): the flag repeats, and a list typed on
    # the command line REPLACES the preset's list.
    multiple: bool = False
    # A command-line-only companion that reads the list from a file
    # (`--gene-file` for `genes`), as (flag, Click parameter). Its tokens are
    # merged into the key before resolution: typed together with the flag, the
    # two are unioned. The parameter is in NEVER_PRESET (a preset stores the
    # list itself, never a path).
    file_option: tuple[str, str] | None = None
    # Left out of a preset's fingerprint: a purely descriptive key (a gene
    # list's label) that changes neither which rows a run keeps nor how its
    # values are read. `regions_build` is NOT exempt: it decides how regions
    # are checked, so a pin must notice a change to it.
    fingerprint_exempt: bool = False

    @property
    def is_bool(self) -> bool:
        return self.type == "bool"

    @property
    def is_list(self) -> bool:
        return self.type in LIST_TYPES


REGISTRY: dict[str, FilterKey] = {}


def register(key: FilterKey) -> FilterKey:
    """Add a key after the existing ones. Registration order is the order of
    the generated Click options and of the echo."""
    if key.file_key in REGISTRY:
        raise ValueError(f"filter key {key.file_key!r} is already registered")
    if key.kind == VALUE and not (key.param and key.flag):
        raise ValueError(f"{key.file_key!r}: a value key needs a param and a flag")
    if key.param and any(k.param == key.param for k in REGISTRY.values()):
        raise ValueError(f"Click parameter {key.param!r} is already registered")
    unknown = set(key.commands) - set(COMMANDS)
    if unknown:
        raise ValueError(f"{key.file_key!r}: unknown commands {sorted(unknown)}")
    REGISTRY[key.file_key] = key
    return key


# What a preset may never carry, by file key or Click parameter name, and why.
# `preset save` refuses these with the reason; the drift test accepts a Click
# option that is listed here instead of registered.
_SELECTION = "selects who is analysed; a preset never changes which trio"
_OUTPUT = "output control, not a filter"
_ENVIRONMENT = "an environment setting or a path; a preset file carries neither"
_SAFETY = "a safety check; a shared preset must not switch it off"
_ATTEST = "an attestation about a patient; a preset can never assert one"
_EVIDENCE = "would weaken the evidence; a preset may only make a run stricter"
_PRESET = "an option of the preset mechanism itself; presets do not nest"

NEVER_PRESET: dict[str, str] = {
    **dict.fromkeys(
        (
            "name",
            "proband",
            "family",
            "father",
            "mother",
            "all_trios",
            "ingest_id",
            "mz_pair",
            "affected",
            "unaffected",
            "compare",
            "group",
            "case",
            "write_case",
            "since",
            "max_samples",
        ),
        _SELECTION,
    ),
    **dict.fromkeys(
        (
            "format",
            "out_format",
            "limit",
            "funnel",
            "near_miss",
            "view",
            "class",
            "carrier_set",
            "explain",
        ),
        _OUTPUT,
    ),
    **dict.fromkeys(("build", "fasta", "hpo_file", "curated_genes"), _ENVIRONMENT),
    **dict.fromkeys(("parentage", "x_guard", "allow_unresolved_genes"), _SAFETY),
    **dict.fromkeys(
        ("gene_file", "exclude_gene_file", "region_file"),
        "a file read on the command line; a preset stores the list itself "
        "(genes, exclude_genes, regions), never a path",
    ),
    **dict.fromkeys(
        (
            "preset",
            "preset_fingerprint",
            "preset_ignore",
            "preset_skip_unavailable",
        ),
        _PRESET,
    ),
    **dict.fromkeys(("attest_affected", "attest_affected_from_pedigree"), _ATTEST),
    "calling": _EVIDENCE + " (it declares how the cohort was called)",
}


# ───────────────────────────── the keys of Pr-a ─────────────────────────────

_ALL = COMMANDS
_AB = "allele balance"

register(
    FilterKey(
        "category",
        "category",
        "--category",
        "category",
        "all",
        "Pr-a",
        ("trio",),
        choices=_category_choices,
        help="Inheritance model: a name, a comma list (recessive,comphet) or "
        "'all'. 'all' prints per-model candidate counts.",
    )
)
register(
    FilterKey(
        "min_gq",
        "min_gq",
        "--min-gq",
        "int",
        20,
        "Pr-a",
        _ALL,
        lo=0,
        hi=99,
        help="Minimum genotype quality (GQ) of the proband and both parents. "
        "A call with no GQ passes unless --strict-missing.",
    )
)
register(
    FilterKey(
        "min_dp",
        "min_dp",
        "--min-dp",
        "int",
        10,
        "Pr-a",
        _ALL,
        lo=0,
        hi=10000,
        help="Minimum read depth (DP) of the proband and both parents. "
        "A call with no DP passes unless --strict-missing.",
    )
)
register(
    FilterKey(
        "max_af",
        "max_af",
        "--max-af",
        "float",
        DEFAULT_MAX_AF,
        "Pr-a",
        _ALL,
        lo=0.0,
        hi=1.0,
        help="Keep variants with population info_AF <= this (rarity filter). "
        "If INFO/AF turns out to be this cohort's own frequency (GATK-style "
        "family VCFs), the filter is skipped with a note unless this option is "
        "given explicitly.",
    )
)
register(
    FilterKey(
        "min_ab",
        "min_ab",
        "--min-ab",
        "float",
        0.25,
        "Pr-a",
        _ALL,
        lo=0.0,
        hi=1.0,
        help=f"Lowest {_AB} (alt reads / all reads) a proband het call may have.",
    )
)
register(
    FilterKey(
        "max_ab",
        "max_ab",
        "--max-ab",
        "float",
        0.75,
        "Pr-a",
        _ALL,
        lo=0.0,
        hi=1.0,
        help=f"Highest {_AB} a proband het call may have.",
    )
)
register(
    FilterKey(
        "gnomad_max_af",
        "gnomad_max_af",
        "--gnomad-max-af",
        "float",
        None,
        "Pr-a",
        _ALL,
        lo=0.0,
        hi=1.0,
        help="Additionally drop candidates whose gnomAD popmax AF exceeds this "
        "(needs `vcfclick annotations load-gnomad`). Variants absent from the "
        "loaded gnomAD slice are kept as rare.",
    )
)
register(
    FilterKey(
        "strict_missing",
        "strict_missing",
        "--strict-missing/--no-strict-missing",
        "bool",
        False,
        "Pr-a",
        _ALL,
        help="A missing GQ, DP or (for a het call) AD on the proband or a "
        "parent counts as a fail, not a pass, for each of those gates that is "
        "asked for. A gate whose field is absent from more than 10% of the "
        "rows the selected models read at their candidate sites cannot be "
        "enforced and is an error (--no-require-bound-gates runs anyway, with "
        "that gate tolerant for the run).",
    )
)
register(
    FilterKey(
        "pass_only",
        "pass_only",
        "--pass-only/--no-pass-only",
        "bool",
        False,
        "Pr-a",
        _ALL,
        help="Keep only variants whose VCF FILTER is PASS or unset.",
    )
)
register(
    FilterKey(
        "parent_alt_screen",
        "parent_alt_screen",
        "--parent-alt-screen/--no-parent-alt-screen",
        "bool",
        False,
        "Pr-a",
        ("trio",),
        help="De novo only: a parent's stored 0/0 is trusted only when its alt "
        "reads stay under --max-parent-alt-frac (a parent 0/0 with >= 2 alt "
        "reads and a high alt fraction is not proof of absence). Shows the "
        "parents' read counts on de novo candidates.",
    )
)
register(
    FilterKey(
        "max_parent_alt_frac",
        "max_parent_alt_frac",
        "--max-parent-alt-frac",
        "float",
        0.05,
        "Pr-a",
        _ALL,
        lo=0.0,
        hi=1.0,
        help="With --parent-alt-screen: the parent alt-read fraction at and "
        "above which a stored 0/0 is not trusted.",
    )
)
register(
    FilterKey(
        "require_bound",
        "require_bound",
        "--require-bound-gates/--no-require-bound-gates",
        "bool",
        False,
        "Pr-a",
        _ALL,
        help="Error, instead of a warning, when a requested gate cannot bind "
        "in this data (the field is missing on too many rows), a selected model "
        "cannot run on this database, or no population-frequency filter is in "
        "effect. --strict-missing implies it for the GQ, DP and AD gates; "
        "--no-require-bound-gates overrides both for this run.",
    )
)

# ───────────────────────────── the keys of X ─────────────────────────────

register(
    FilterKey(
        "min_male_x_dp",
        "min_male_x_dp",
        "--min-male-x-dp",
        "int",
        DERIVED,
        "X",
        ("trio",),
        lo=0,
        hi=10000,
        help="xlinked: minimum DP of a hemizygous call (a male proband, and the "
        "father of a daughter) on non-PAR chrX. Default: half of --min-dp, "
        "rounded down (one X copy gives half the reads).",
    )
)
register(
    FilterKey(
        "max_af_x",
        "max_af_x",
        "--max-af-x",
        "float",
        0.001,
        "X",
        ("trio",),
        lo=0.0,
        hi=1.0,
        help="xlinked: rarity cap on non-PAR chrX; the model uses the lower of "
        "this and --max-af. A hemizygous male is affected at the allele "
        "frequency q, not q squared, so X needs a stricter ceiling.",
    )
)
register(
    FilterKey(
        "min_hemi_ab",
        "min_hemi_ab",
        "--min-hemi-ab",
        "float",
        0.75,
        "X",
        ("trio",),
        lo=0.0,
        hi=1.0,
        help="xlinked: lowest alt fraction (alt reads / all reads) of a male "
        "proband's hemizygous ALT call. A call without allele depths passes.",
    )
)
register(
    FilterKey(
        "max_parent_alt_reads",
        "max_parent_alt_reads",
        "--max-parent-alt-reads",
        "int",
        1,
        "X",
        ("trio",),
        lo=0,
        hi=10000,
        help="xlinked de novo: most alt reads a parent's stored 0/0 may carry "
        "(a missing allele depth passes).",
    )
)

# ───────────────────────────── the keys of Pr-b ─────────────────────────────
# The scope: which genes and regions a run looks at (cli/scope.py).

_BUILDS = ("GRCh37", "GRCh38")

register(
    FilterKey(
        "genes",
        "genes",
        "--genes",
        "genes",
        (),
        "Pr-b",
        _ALL,
        multiple=True,
        file_option=("--gene-file", "gene_file"),
        help="Only variants inside these genes: symbols (any case), Ensembl gene "
        "ids (ENSG...) and, after `vcfclick annotations load-hgnc`, previous "
        "symbols, HGNC ids (HGNC:1234) and Entrez ids. Comma list, repeatable. "
        "A gene is its whole span (introns and UTRs); promoters need "
        "--gene-flank. Needs `vcfclick annotations load`. --gene-file FILE reads "
        "one per line (# comments); typed with --genes, the two are joined.",
    )
)
register(
    FilterKey(
        "genes_label",
        "genes_label",
        "--genes-label",
        "text",
        None,
        "Pr-b",
        _ALL,
        kind=META,
        widening=False,
        fingerprint_exempt=True,
        help="A label for the gene list, such as the panel and its version "
        "(`PanelApp ID v4.2`). Descriptive only.",
    )
)
register(
    FilterKey(
        "exclude_genes",
        "exclude_genes",
        "--exclude-genes",
        "genes",
        (),
        "Pr-b",
        _ALL,
        multiple=True,
        file_option=("--exclude-gene-file", "exclude_gene_file"),
        help="Leave out variants inside these genes (same tokens as --genes). "
        "--exclude-gene-file FILE reads one per line.",
    )
)
register(
    FilterKey(
        "gene_flank",
        "gene_flank",
        "--gene-flank",
        "int",
        0,
        "Pr-b",
        _ALL,
        lo=0,
        hi=1_000_000,
        widening=False,
        help="Widen every gene of --genes and --exclude-genes by this many base "
        "pairs on each side (promoters, upstream elements).",
    )
)
register(
    FilterKey(
        "regions",
        "region",
        "--region",
        "regions",
        (),
        "Pr-b",
        _ALL,
        multiple=True,
        file_option=("--region-file", "region_file"),
        help="Only variants inside this region, chrN:FROM-TO, 1-based and "
        "inclusive (repeatable). With --genes, a variant in either is kept. "
        "--region-file BED reads 0-based half-open BED intervals.",
    )
)
register(
    FilterKey(
        "exclude_regions",
        "exclude_region",
        "--exclude-region",
        "regions",
        (),
        "Pr-b",
        _ALL,
        multiple=True,
        help="Leave out variants inside this region, chrN:FROM-TO (repeatable).",
    )
)
register(
    FilterKey(
        "regions_build",
        "regions_build",
        "--regions-build",
        "text",
        None,
        "Pr-b",
        _ALL,
        kind=META,
        widening=False,
        choices=lambda: _BUILDS,
        help="The assembly the regions are written in (GRCh37 or GRCh38), "
        "checked against the database's build when it is known.",
    )
)

DEFAULTS: dict[str, Any] = {
    k.file_key: k.default for k in REGISTRY.values() if k.default is not DERIVED
}


# ──────────────────────────────── lookups ────────────────────────────────


def keys_for(command: str, *, kind: str | None = VALUE) -> list[FilterKey]:
    """The keys that apply to `command`, in registration order."""
    return [
        k
        for k in REGISTRY.values()
        if command in k.commands and (kind is None or k.kind == kind)
    ]


def by_param(command: str) -> dict[str, FilterKey]:
    """Click parameter name -> key, for one command."""
    return {k.param: k for k in keys_for(command) if k.param}


def suggest(name: str) -> list[str]:
    """Registered keys close to a misspelt one (hyphens count as underscores)."""
    wanted = name.replace("-", "_").lower()
    pool = [*REGISTRY, *NEVER_PRESET]
    if wanted in REGISTRY and wanted != name:
        return [wanted]
    return difflib.get_close_matches(wanted, pool, n=3, cutoff=0.6)


# ───────────────────────────── values and checks ─────────────────────────────


def category_tokens(value: Any) -> list[str] | None:
    """A category value as a list of unique tokens, or None if it is neither a
    string (maybe comma-separated) nor a list of strings."""
    if isinstance(value, str):
        tokens: list[Any] = [t.strip() for t in value.split(",")]
    elif isinstance(value, (list, tuple)):
        tokens = list(value)
    else:
        return None
    if not tokens or not all(isinstance(t, str) for t in tokens):
        return None
    return list(dict.fromkeys(tokens))


def check_value(key: FilterKey, value: Any) -> tuple[Any, str | None]:
    """Validate and normalise one preset value: (value, None) or (None, why).

    `bool` is refused for numbers (Click's int type accepts True), NaN and
    infinities are refused (a JSON `1e999` parses to inf), and category tokens
    must be the registry's. Ranges are those of a preset value; a typed flag
    keeps the looser rule it always had."""
    if value is None:
        return None, "null is not allowed; omit the key to inherit the default"
    if key.type == "bool":
        if isinstance(value, bool):
            return value, None
        return None, "must be true or false"
    if key.type == "category":
        tokens = category_tokens(value)
        if tokens is None:
            return None, "must be a model name or a list of model names"
        allowed = (key.choices or _category_choices)()
        bad = [t for t in tokens if t not in allowed]
        if bad:
            return None, (
                f"{bad[0]!r} is not one of {', '.join(allowed)}"
                + _near(bad[0], allowed)
            )
        return tokens, None
    if key.type in LIST_TYPES:
        return check_list(key, value)
    if key.type == "text":
        if not isinstance(value, str) or not value.strip():
            return None, "must be text"
        if len(value) > MAX_TEXT:
            return None, f"must be at most {MAX_TEXT} characters"
        allowed = key.choices() if key.choices else None
        if allowed and value not in allowed:
            return None, f"must be one of {', '.join(allowed)}"
        return value, None
    if key.type == "int":
        if isinstance(value, bool) or not isinstance(value, int):
            return None, "must be an integer"
        number: float = value  # compared as it is: an int of any size is exact
    elif key.type == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return None, "must be a number"
        try:
            number = float(value)  # an int of 400 digits does not fit a float
        except OverflowError:
            return None, "must be a finite number"
        if not math.isfinite(number):
            return None, "must be a finite number"
        value = number
    else:
        return None, f"unsupported type {key.type!r}"
    below = key.lo is not None and number < key.lo
    above = key.hi is not None and number > key.hi
    if below or above:
        return None, f"must be between {key.lo:g} and {key.hi:g}"
    return value, None


def list_tokens(key: FilterKey, value: Any) -> list[str] | None:
    """A list key's value as its tokens: a string or a list of strings, each
    split on commas (genes only: a region never holds one), stripped, empty
    ones dropped, duplicates removed in order (a gene without regard to case,
    the first spelling kept). None for a value that is no list of text."""
    if isinstance(value, str):
        items: list[Any] = [value]
    elif isinstance(value, (list, tuple)):
        items = list(value)
    else:
        return None
    if not all(isinstance(t, str) for t in items):
        return None
    tokens: list[str] = []
    for item in items:
        parts = item.split(",") if key.type == "genes" else [item]
        tokens += [p.strip() for p in parts if p.strip()]
    seen: dict[str, str] = {}
    for t in tokens:
        seen.setdefault(t.upper() if key.type == "genes" else t, t)
    return list(seen.values())


def parse_region(text: str) -> tuple[str, int, int] | None:
    """`chr7:100-200` as (contig as written, 100, 200), or None when it is not
    a region (bad shape, FROM above TO, a zero start)."""
    m = REGION_RE.match(text)
    if m is None:
        return None
    lo, hi = int(m.group(2)), int(m.group(3))
    if lo < 1 or lo > hi:
        return None
    return text.split(":", 1)[0], lo, hi


def check_list(key: FilterKey, value: Any) -> tuple[Any, str | None]:
    tokens = list_tokens(key, value)
    if tokens is None:
        return None, "must be a list of text values"
    if key.type == "genes":
        cap, what = MAX_GENES, "gene"
        bad = [t for t in tokens if not GENE_TOKEN_RE.match(t)]
        why = "is not a gene symbol or id"
    else:
        cap, what = MAX_REGIONS, "region"
        bad = [t for t in tokens if parse_region(t) is None]
        why = "is not a region chrN:FROM-TO with FROM <= TO"
    if len(tokens) > cap:
        return None, f"at most {cap:,} {what}s"
    if bad:
        return None, f"{bad[0][:40]!r} {why}"
    return tokens, None


def _near(token: str, allowed: Iterable[str]) -> str:
    close = difflib.get_close_matches(token, list(allowed), n=1)
    return f" (did you mean {close[0]!r}?)" if close else ""


def normalise_value(key: FilterKey, value: Any) -> Any:
    """A value as the echo and the fingerprint show it: category as a list,
    numbers as they are. (Not a validator.)"""
    if key.type == "category":
        tokens = category_tokens(value)
        return tokens if tokens is not None else value
    if key.type in LIST_TYPES:
        tokens = list_tokens(key, value)
        return tokens if tokens is not None else value
    return value


# ───────────────────────────── precedence ─────────────────────────────

FROM_PRESET, FROM_FLAG, FROM_DEFAULT, FROM_DERIVED = (
    "preset",
    "flag",
    "default",
    "derived",
)
OFF = "off"  # `overridden[].flag` for a key switched off by --preset-ignore


@dataclass(frozen=True)
class Resolved:
    value: Any
    source: str


def resolve_filters_from_values(
    command: str,
    preset_filters: Mapping[str, Any],
    explicit: Mapping[str, Any],
    ignore: Iterable[str] = (),
) -> tuple[dict[str, Resolved], list[dict]]:
    """Explicit value > preset > default, for every key `command` takes.

    `explicit` holds only what the caller sent (a typed flag, a request
    parameter); `ignore` the keys to switch off for this run. Returns
    {file_key: Resolved} in registry order, and `overridden`: one
    {"key", "preset", "flag"} per key the preset set whose final source is a
    flag (its value) or "off". A preset key this command does not take is not
    resolved here (the caller reports it as skipped)."""
    ignored = set(ignore)
    resolved: dict[str, Resolved] = {}
    overridden: list[dict] = []
    for key in keys_for(command):
        fk = key.file_key
        has_preset = fk in preset_filters
        if fk in explicit:
            value = normalise_value(key, explicit[fk])
            resolved[fk] = Resolved(value, FROM_FLAG)
            if has_preset:
                overridden.append(
                    {
                        "key": fk,
                        "preset": normalise_value(key, preset_filters[fk]),
                        "flag": value,
                    }
                )
        elif has_preset and fk in ignored:
            resolved[fk] = _default(key)
            overridden.append(
                {
                    "key": fk,
                    "preset": normalise_value(key, preset_filters[fk]),
                    "flag": OFF,
                }
            )
        elif has_preset:
            resolved[fk] = Resolved(
                normalise_value(key, preset_filters[fk]), FROM_PRESET
            )
        else:
            resolved[fk] = _default(key)
    return resolved, overridden


def _default(key: FilterKey) -> Resolved:
    if key.default is DERIVED:
        return Resolved(None, FROM_DERIVED)
    return Resolved(normalise_value(key, key.default), FROM_DEFAULT)


def cross_key_problems(values: Mapping[str, Any]) -> list[str]:
    """Rules that span keys, checked on RESOLVED values so that a flag can
    repair a preset: the allele-balance band must not be empty."""
    out = []
    lo, hi = values.get("min_ab"), values.get("max_ab")
    if lo is not None and hi is not None and lo > hi:
        out.append(f"min_ab ({lo:g}) is above max_ab ({hi:g})")
    return out
