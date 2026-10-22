"""What a preset file is: its envelope, validation, versioning, fingerprint.

A preset is a small JSON document of filter values, with no sample, family,
path or database identity (a key of that kind anywhere in the envelope, at the
top level or in `meta`, is an error, not a thing to preserve):

    {"schema": "vcfclick.preset", "schema_version": 1, "name": "id-panel",
     "title": "...", "description": "...", "applies_to": ["trio"],
     "assay": "exome",
     "filters": {"min_gq": 20, "max_af": 0.005, "category": ["recessive"]},
     "meta": {"created_with": "0.15.0", "updated": "2026-10-21T09:00:00Z",
              "imported_from": null, "original_fingerprint": null,
              "approved_by": null, "approved_on": null, "locked": false}}

The file is data: parsed, validated, never executed, and used in SQL only as
typed values. Nothing here touches the file system or Click; the store is
`cli.presets`. Error messages never quote file content beyond a key name.

Versioning. `schema_version` changes only when the envelope changes; adding a
filter key never bumps it. A reader accepts versions up to SCHEMA_VERSION
(`migrate` is the chain, empty at 1) and rejects a newer file by name. Unknown
envelope keys are preserved on load and written back on save. `meta` is not
part of the fingerprint.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

from cli import filter_registry as fr

SCHEMA = "vcfclick.preset"
BUNDLE_SCHEMA = "vcfclick.preset-bundle"
SCHEMA_VERSION = 1  # the newest envelope this vcfclick reads
MAX_BYTES = 2 * 1024 * 1024
ASSAYS = ("genome", "exome", "panel")
EXTENSION_PREFIX = "x_"  # a free namespace: kept in the file, never applied
# `NAME.vcfclick-preset.json` is what `preset export -o DIR` writes: the name of
# the preset in such a file is NAME, so a name may not end in the stem itself.
EXPORT_STEM = ".vcfclick-preset"

NAME_RE = re.compile(r"^[a-z0-9][a-z0-9._-]{0,62}$")
_WINDOWS_DEVICES = frozenset(
    {"con", "prn", "aux", "nul"}
    | {f"com{i}" for i in range(1, 10)}
    | {f"lpt{i}" for i in range(1, 10)}
)

ERROR = "error"
WARN = "warning"

# The order the envelope is written in; keys not listed follow in file order.
_ENVELOPE_ORDER = (
    "schema",
    "schema_version",
    "name",
    "title",
    "description",
    "applies_to",
    "assay",
    "filters",
    "meta",
)
# Envelope or meta keys that would carry who or where: never kept, whatever the
# file says (NEVER_PRESET's keys are refused too, bar the envelope's own `name`).
IDENTITY_KEYS = frozenset(
    {
        "db",
        "database",
        "path",
        "sample",
        "samples",
        "sample_id",
        "proband_id",
        "family_id",
        "patient",
        "patient_id",
        "mrn",
        "vcf",
        "ped",
    }
)
# The same keys with case, spaces, hyphens, underscores and dots taken out, so
# `Patient`, `sample-id`, `sampleId` and `SAMPLE_ID` are caught as the keys they
# spell. (NEVER_PRESET is matched the same way.)
_KEY_NOISE = re.compile(r"[\s_.\-]+")


def flat_key(key: str) -> str:
    return _KEY_NOISE.sub("", key).casefold()


_IDENTITY_FLAT = frozenset(flat_key(k) for k in IDENTITY_KEYS)
_NEVER_FLAT = {flat_key(k): k for k in fr.NEVER_PRESET}
_META_STRINGS = (
    "created_with",
    "updated",
    "imported_from",
    "original_fingerprint",
    "approved_by",
    "approved_on",
)


class PresetError(Exception):
    """A problem with a preset that the user can act on. The message is safe
    to show (it quotes no file content)."""


@dataclass(frozen=True)
class Problem:
    level: str  # error | warning
    code: str
    message: str
    key: str | None = None

    def to_json(self) -> dict:
        return {
            "level": self.level,
            "code": self.code,
            "message": self.message,
            "key": self.key,
        }


# ─────────────────────────────── names ───────────────────────────────


def name_problem(name: object) -> str | None:
    """Why `name` cannot be a preset name, or None. Names are lowercase and
    never look like a path: a value with a separator or ending .json is read
    as a file on the command line, so a preset may not have one."""
    if not isinstance(name, str) or not NAME_RE.match(name):
        return (
            "a preset name is 1-63 characters: lowercase letters, digits, "
            "'.', '_' or '-', starting with a letter or digit"
        )
    if name.endswith(".json"):
        return "a preset name may not end in .json (it would be read as a file)"
    if name.endswith(EXPORT_STEM):
        return (
            f"a preset name may not end in {EXPORT_STEM} (that is the suffix of an "
            "exported file, whose name is what comes before it)"
        )
    if name.split(".")[0] in _WINDOWS_DEVICES:
        return "a preset name may not be a Windows device name"
    return None


def identity_problem(key: str, *, envelope: bool = True) -> str | None:
    """Why a top-level (`envelope`), `meta` or filter key may not be in a preset
    file, or None: it names a sample, family, path or database (a preset carries
    none), or is one of the keys a preset may never carry. Matched without regard
    to case, spaces, hyphens, underscores and dots (`Patient`, `sample-id`). The
    envelope's own keys (`name`, `filters`, ...) are not in question."""
    if envelope and key in _ENVELOPE_ORDER:
        return None
    flat = flat_key(key)
    never = key if key in fr.NEVER_PRESET else _NEVER_FLAT.get(flat)
    if never is not None:
        return f"can never be part of a preset: {fr.NEVER_PRESET[never]}"
    if flat in _IDENTITY_FLAT:
        return (
            "can never be part of a preset: it would carry a sample, family, "
            "path or database identity"
        )
    return None


def filter_identity_problem(key: str) -> str | None:
    """`identity_problem` for a key under `filters`: an extension (`x_`) key and
    a registered filter are what they are; any other key that spells an identity
    is refused (an unknown key that does not is `unknown_key`, a warning in a
    stored file)."""
    if key.startswith(EXTENSION_PREFIX) or key in fr.REGISTRY:
        return None
    return identity_problem(key, envelope=False)


# ───────────────────────────── parsing ─────────────────────────────


def _reject_constant(token: str) -> Any:
    raise ValueError(f"{token} is not allowed in a preset")


def _reject_nonfinite(text: str) -> float:
    """`1e999` parses to infinity, which `parse_constant` never sees."""
    value = float(text)
    if not math.isfinite(value):
        raise ValueError("a number is out of range")
    return value


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict:
    out: dict = {}
    for key, value in pairs:
        if key in out:
            raise ValueError(f"duplicate key {key[:40]!r}")
        out[key] = value
    return out


def parse_json(data: bytes | str) -> Any:
    """Parse preset JSON: no NaN or Infinity (or a number too large for a
    float, which `1e999` is), no duplicate keys, at most MAX_BYTES. Raises
    PresetError (the message names the position, not the content)."""
    raw = data.encode() if isinstance(data, str) else data
    if len(raw) > MAX_BYTES:
        raise PresetError(f"file is larger than {MAX_BYTES // (1024 * 1024)} MiB")
    try:
        return json.loads(
            raw.decode("utf-8"),
            parse_constant=_reject_constant,
            parse_float=_reject_nonfinite,
            object_pairs_hook=_no_duplicates,
        )
    except UnicodeDecodeError:
        raise PresetError("not valid UTF-8 JSON") from None
    except json.JSONDecodeError as e:
        raise PresetError(
            f"not valid JSON: {e.msg} (line {e.lineno}, column {e.colno})"
        ) from None
    except (ValueError, RecursionError) as e:
        raise PresetError(f"not a valid preset: {e}") from None


# ───────────────────────────── versioning ─────────────────────────────

_MIGRATIONS: dict[int, Callable[[dict], dict]] = {}  # from_version -> to the next


def migrate(doc: dict, from_version: int) -> dict:
    """Bring a document of an older envelope up to SCHEMA_VERSION. Empty at
    version 1; a step is registered here by the change that bumps the version."""
    version = from_version
    while version < SCHEMA_VERSION:
        doc = _MIGRATIONS[version](doc)
        version += 1
    return doc


# ───────────────────────────── validation ─────────────────────────────


@dataclass(frozen=True)
class Validated:
    doc: dict  # the normalised envelope (unknown keys kept)
    problems: tuple[Problem, ...]
    unknown_keys: tuple[str, ...]  # filter keys this binary does not know

    @property
    def errors(self) -> list[Problem]:
        return [p for p in self.problems if p.level == ERROR]

    @property
    def ok(self) -> bool:
        return not self.errors


def vcfclick_version() -> str:
    try:
        from importlib.metadata import version

        return version("vcfclick")
    except Exception:  # noqa: BLE001 - not installed as a distribution
        return "unknown"


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _author_note(doc: dict) -> str:
    meta = doc.get("meta")
    made_by = meta.get("created_with") if isinstance(meta, dict) else None
    if isinstance(made_by, str) and made_by and made_by != vcfclick_version():
        return f" (written by vcfclick {made_by[:20]}; upgrade)"
    return ""


def validate_document(doc: object, *, allow_unknown: bool = False) -> Validated:
    """Check a parsed preset and return it normalised, with every problem.
    Never raises: a document the checks themselves cannot digest is reported
    as not a readable preset, so a bad file in the presets directory can be
    listed (`preset list` shows it with its problem) instead of crashing every
    command that reads the directory.

    A filter key this binary does not know is an error with a did-you-mean
    (`allow_unknown` downgrades it to a warning, for a file that was stored
    on purpose with `import --allow-unknown-keys` and for reading a stored
    preset, where such a key is "unavailable" at apply time instead). A key in
    `NEVER_PRESET`, a wrong type, an out-of-range number and a null are errors
    in every mode. So is any key, at the top level or in `meta`, that would
    carry a sample, family, path or database identity, and a number anywhere
    in the file that is not finite."""
    try:
        return _validate_document(doc, allow_unknown)
    except Exception:  # noqa: BLE001 - whatever the file held, it is a problem
        return Validated(
            {},
            (Problem(ERROR, "not_a_preset", "not a readable preset"),),
            (),
        )


def _has_nonfinite(value: Any) -> bool:
    """Does a NaN or infinity sit anywhere in this parsed value?"""
    stack = [value]
    while stack:
        v = stack.pop()
        if isinstance(v, float) and not math.isfinite(v):
            return True
        if isinstance(v, dict):
            stack.extend(v.values())
        elif isinstance(v, (list, tuple)):
            stack.extend(v)
    return False


def _validate_document(doc: object, allow_unknown: bool) -> Validated:
    problems: list[Problem] = []

    def bad(code: str, message: str, key: str | None = None, level: str = ERROR):
        problems.append(Problem(level, code, message, key))

    if not isinstance(doc, dict):
        bad("not_a_preset", "not a vcfclick preset (the top level is not an object)")
        return Validated({}, tuple(problems), ())
    if doc.get("schema") != SCHEMA:
        bad("not_a_preset", f"not a vcfclick preset (schema must be {SCHEMA!r})")
        return Validated({}, tuple(problems), ())
    version = doc.get("schema_version")
    if isinstance(version, bool) or not isinstance(version, int) or version < 1:
        bad("unsupported_schema", "not a readable preset: schema_version is missing")
        return Validated({}, tuple(problems), ())
    if version > SCHEMA_VERSION:
        bad(
            "newer_schema",
            f"made by a newer vcfclick (schema {version}); upgrade to read it",
        )
        return Validated({}, tuple(problems), ())
    doc = migrate(doc, version)

    out: dict = {"schema": SCHEMA, "schema_version": SCHEMA_VERSION}
    _check_text_fields(doc, out, bad)
    out["filters"], unknown = _check_filters(doc, allow_unknown, bad)
    out["meta"] = _check_meta(doc.get("meta"), bad)
    for key, value in doc.items():  # unknown envelope keys survive a rewrite
        if key in _ENVELOPE_ORDER:
            continue  # the checks above own these: a wrong one is dropped, not kept
        why = identity_problem(key)
        if why:
            bad("never_presettable", f"{key[:40]} {why}", key[:40])
            continue
        out.setdefault(key, value)
    ordered = {k: out[k] for k in _ENVELOPE_ORDER if k in out}
    ordered.update({k: v for k, v in out.items() if k not in ordered})
    for key, value in ordered.items():  # x_ values, unknown and meta keys too
        if _has_nonfinite(value):
            bad("bad_value", f"{key[:40]} holds a number that is not finite", key[:40])

    if not any(p.level == ERROR for p in problems):
        for message in fr.cross_key_problems(ordered["filters"]):
            bad("inconsistent", message)
    return Validated(ordered, tuple(problems), tuple(unknown))


def _check_text_fields(doc: dict, out: dict, bad: Callable) -> None:
    name = doc.get("name")
    if name is not None:
        why = name_problem(name)
        if why:
            bad("bad_name", why, "name")
        else:
            out["name"] = name
    for key in ("title", "description"):
        value = doc.get(key)
        if value is None:
            continue
        if not isinstance(value, str):
            bad("bad_type", f"{key} must be text", key)
        else:
            out[key] = value
    applies = doc.get("applies_to")
    if applies is not None:
        if not isinstance(applies, list) or not all(
            isinstance(a, str) for a in applies
        ):
            bad("bad_type", "applies_to must be a list of command names", "applies_to")
        else:
            out["applies_to"] = list(applies)
            for a in applies:
                if a not in fr.COMMANDS:
                    bad(
                        "unknown_command",
                        f"applies_to names {a[:20]!r}, which is not one of "
                        f"{', '.join(fr.COMMANDS)}",
                        "applies_to",
                        WARN,
                    )
    assay = doc.get("assay")
    if assay is not None:
        if assay not in ASSAYS:
            bad("bad_value", f"assay must be one of {', '.join(ASSAYS)}", "assay")
        else:
            out["assay"] = assay


def _check_filters(
    doc: dict, allow_unknown: bool, bad: Callable
) -> tuple[dict, list[str]]:
    raw = doc.get("filters")
    if not isinstance(raw, dict):
        bad("bad_type", "filters must be an object", "filters")
        return {}, []
    known: dict[str, Any] = {}
    other: dict[str, Any] = {}
    unknown: list[str] = []
    for key, value in raw.items():
        shown = key[:40]
        if key.startswith(EXTENSION_PREFIX):
            other[key] = value
        elif key in fr.REGISTRY:
            spec = fr.REGISTRY[key]
            if not spec.presettable:
                bad("never_presettable", f"{shown} cannot be part of a preset", key)
                continue
            checked, why = fr.check_value(spec, value)
            if why:
                bad("bad_value", f"{shown}: {why}", key)
            else:
                known[key] = checked
        elif key in fr.NEVER_PRESET:
            bad(
                "never_presettable",
                f"{shown} can never be part of a preset: {fr.NEVER_PRESET[key]}",
                key,
            )
        elif who := filter_identity_problem(key):
            bad("never_presettable", f"{shown} {who}", key)
        else:
            near = fr.suggest(key)
            hint = f" (did you mean {', '.join(near)}?)" if near else ""
            bad(
                "unknown_key",
                f"unknown filter key {shown!r}{hint}{_author_note(doc)}",
                key,
                WARN if allow_unknown else ERROR,
            )
            unknown.append(key)
            other[key] = value
    ordered = {k: known[k] for k in fr.REGISTRY if k in known}
    ordered.update(other)
    return ordered, unknown


def _check_meta(meta: object, bad: Callable) -> dict:
    out: dict = {}
    if meta is None:
        return out
    if not isinstance(meta, dict):
        bad("bad_type", "meta must be an object", "meta")
        return out
    for key, value in meta.items():
        if key in _META_STRINGS:
            if value is not None and not isinstance(value, str):
                bad("bad_type", f"meta.{key} must be text or null", f"meta.{key}")
                continue
        elif key == "locked":
            if not isinstance(value, bool):
                bad("bad_type", "meta.locked must be true or false", "meta.locked")
                continue
        else:
            why = identity_problem(key, envelope=False)
            if why:
                bad(
                    "never_presettable",
                    f"meta.{key[:40]} {why}",
                    f"meta.{key[:40]}",
                )
                continue
        out[key] = value
    return out


# ───────────────────────────── fingerprint ─────────────────────────────


def _canonical_value(key: str, value: Any) -> Any:
    spec = fr.REGISTRY.get(key)
    if spec is not None and spec.type == "category":
        tokens = fr.category_tokens(value)
        return sorted(tokens) if tokens is not None else value
    if spec is not None and spec.type == "float" and isinstance(value, int):
        return float(value)
    if isinstance(value, float) and not math.isfinite(value):
        raise PresetError("a filter value is not finite")
    return value


def fingerprint(doc: dict) -> str:
    """First 12 hex of the sha256 of the canonical JSON of {schema_version,
    filters}: keys sorted, category tokens sorted, floats by repr, extension
    (x_) keys and `meta` left out. The same filters give the same fingerprint
    whatever the key order or the descriptive text."""
    filters = {
        k: _canonical_value(k, v)
        for k, v in doc.get("filters", {}).items()
        if not k.startswith(EXTENSION_PREFIX)
    }
    payload = {
        "schema_version": doc.get("schema_version", SCHEMA_VERSION),
        "filters": filters,
    }
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()[:12]


def named(doc: dict, name: str) -> dict:
    """`doc` with its `name` set, in the envelope's place (after the schema)."""
    out: dict = {}
    placed = False
    for key, value in doc.items():
        if key == "name":
            continue
        if not placed and key not in ("schema", "schema_version"):
            out["name"] = name
            placed = True
        out[key] = value
    if not placed:
        out["name"] = name
    return out


def to_text(doc: dict) -> str:
    """The file text: indented, UTF-8, one trailing newline. A NaN or infinity
    is never written (it would be a file this program then refuses to read)."""
    try:
        return json.dumps(doc, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    except ValueError:
        raise PresetError("the preset holds a number that is not finite") from None
