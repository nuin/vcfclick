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
    choices: Callable[[], tuple[str, ...]] | None = None  # category tokens

    @property
    def is_bool(self) -> bool:
        return self.type == "bool"


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
    **dict.fromkeys(("parentage", "x_guard"), _SAFETY),
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
        "parent counts as a fail, not a pass. A gate whose field is absent "
        "from more than 10% of the rows cannot be enforced and is an error "
        "(--no-require-bound-gates runs anyway).",
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
        "in this data (the field is missing on too many rows) or no "
        "population-frequency filter is in effect. --strict-missing implies it "
        "for the GQ, DP and AD gates; --no-require-bound-gates overrides both "
        "for this run.",
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
    if key.type == "int":
        if isinstance(value, bool) or not isinstance(value, int):
            return None, "must be an integer"
        number: float = value
    elif key.type == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return None, "must be a number"
        number = float(value)
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


def _near(token: str, allowed: Iterable[str]) -> str:
    close = difflib.get_close_matches(token, list(allowed), n=1)
    return f" (did you mean {close[0]!r}?)" if close else ""


def normalise_value(key: FilterKey, value: Any) -> Any:
    """A value as the echo and the fingerprint show it: category as a list,
    numbers as they are. (Not a validator.)"""
    if key.type == "category":
        tokens = category_tokens(value)
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
