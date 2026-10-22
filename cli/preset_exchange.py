"""Moving presets between people and machines: export, import, from a result.

A preset is shared as a file (`NAME.vcfclick-preset.json`) or a bundle of
several (`{"schema": "vcfclick.preset-bundle", "presets": [...]}`). Export
strips the patient-specific keys (HPO terms, ...) unless asked; import is
all-or-nothing: every file is parsed and validated, every name checked
(a built-in name is refused, a collision refused without --force), and only
then is anything written, with a rollback if a write fails. An imported
preset records the source file's BASENAME and the fingerprint of what was
imported, never a path.
"""

from __future__ import annotations

import contextlib
from dataclasses import dataclass
from pathlib import Path

from cli import filter_registry as fr
from cli.builtin_presets import BUILTINS
from cli.preset_schema import (
    BUNDLE_SCHEMA,
    EXPORT_STEM,
    SCHEMA_VERSION,
    PresetError,
    filter_identity_problem,
    fingerprint,
    identity_problem,
    name_problem,
    now_utc,
    parse_json,
    to_text,
    validate_document,
    vcfclick_version,
)
from cli.presets import (
    Preset,
    find_preset,
    user_file,
    user_path,
    write_text_atomic,
)

EXPORT_SUFFIX = EXPORT_STEM + ".json"  # NAME.vcfclick-preset.json
BUNDLE_NAME = "presets.vcfclick-bundle.json"


# ───────────────────────────────── export ─────────────────────────────────


def export_document(preset: Preset, *, with_patient_data: bool = False) -> dict:
    """The preset as a file to hand to a colleague: without the patient-specific
    keys (unless asked), without where it was imported from, with its
    approval metadata."""
    # a loaded preset has no identity keys, at the top level or under filters
    # (validation refuses them); drop any anyway, so nothing can leave through
    # an export
    doc = {k: v for k, v in preset.doc.items() if identity_problem(k) is None}
    doc["filters"] = {
        k: v
        for k, v in preset.filters.items()
        if filter_identity_problem(k) is None
        and (
            with_patient_data
            or not (k in fr.REGISTRY and fr.REGISTRY[k].patient_specific)
        )
    }
    meta = preset.meta
    doc["meta"] = {
        "created_with": meta.get("created_with") or vcfclick_version(),
        "updated": meta.get("updated") or now_utc(),
        "imported_from": None,
        "original_fingerprint": None,
        "approved_by": meta.get("approved_by"),
        "approved_on": meta.get("approved_on"),
        "locked": bool(meta.get("locked", False)),
    }
    return doc


def bundle_document(docs: list[dict]) -> dict:
    return {"schema": BUNDLE_SCHEMA, "schema_version": SCHEMA_VERSION, "presets": docs}


def export_target(
    names: list[str], out: str | None, bundle: bool
) -> tuple[bool, Path | None]:
    """Where an export goes: (as one bundle?, the file or directory; None is
    stdout). A directory gets one file per preset, or the bundle with --bundle."""
    target = Path(out).expanduser() if out else None
    if target is not None and target.is_dir():
        return (True, target / BUNDLE_NAME) if bundle else (False, target)
    if len(names) > 1 and not bundle:
        raise PresetError(
            "several presets need --bundle (one file) or -o DIRECTORY (one file each)"
        )
    return bundle, target


def load_for_export(name: str) -> Preset:
    if name_problem(name):
        raise PresetError(f"{name[:40]!r} is not a preset name")
    preset, _ = find_preset(name)
    if preset is None:
        raise PresetError(f"no preset named {name!r}")
    if preset.has_errors:
        raise PresetError(
            f"preset {name!r} cannot be exported: {preset.problems[0].message}"
        )
    return preset


# ───────────────────────────────── import ─────────────────────────────────


@dataclass(frozen=True)
class Incoming:
    name: str
    doc: dict  # normalised envelope, meta set for storage
    source: str  # basename it came from


def _envelopes(path: Path) -> list[dict]:
    try:
        data = path.read_bytes()
    except OSError as e:
        raise PresetError(f"{path.name}: cannot read the file ({e.strerror})") from None
    try:
        doc = parse_json(data)
    except PresetError as e:
        raise PresetError(f"{path.name}: {e}") from None
    if isinstance(doc, dict) and doc.get("schema") == BUNDLE_SCHEMA:
        items = doc.get("presets")
        if not isinstance(items, list) or not items:
            raise PresetError(f"{path.name}: the bundle holds no presets")
        return items
    return [doc]


def read_incoming(
    paths: list[str], *, as_name: str | None, allow_unknown: bool
) -> list[Incoming]:
    """Parse and validate every file and bundle; name each preset. Raises
    PresetError on the first problem, naming the file."""
    pending: list[tuple[dict, str]] = []
    for raw in paths:
        path = Path(raw).expanduser()
        for env in _envelopes(path):
            pending.append((env, path.name))
    if as_name is not None and len(pending) != 1:
        raise PresetError("--as names one preset; this import holds several")
    out: list[Incoming] = []
    for env, source in pending:
        checked = validate_document(env, allow_unknown=allow_unknown)
        if not checked.ok:
            raise PresetError(f"{source}: {checked.errors[0].message}")
        doc = checked.doc
        name = as_name or doc.get("name")
        if name is None:
            stem = source.removesuffix(".json").removesuffix(".vcfclick-preset")
            name = stem
        why = name_problem(name)
        if why:
            raise PresetError(f"{source}: {why}")
        if name in BUILTINS:
            raise PresetError(
                f"{source}: {name!r} is a built-in preset's name; import it under "
                "another name with --as"
            )
        meta = dict(doc.get("meta", {}))
        meta.update(
            imported_from=source,
            original_fingerprint=fingerprint(doc),
            updated=now_utc(),
        )
        meta.setdefault("created_with", None)
        out.append(Incoming(name, {**doc, "name": name, "meta": meta}, source))
    names = [i.name for i in out]
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        raise PresetError(f"the import names {dupes[0]!r} more than once")
    return out


def import_presets(incoming: list[Incoming], *, force: bool) -> list[Path]:
    """Write every incoming preset as a user preset, or none: a name that is
    taken is refused without `force` before anything is written, and a failure
    while writing puts back what the import changed."""
    taken = [i.name for i in incoming if user_file(i.name) is not None]
    if taken and not force:
        raise PresetError(
            f"preset {taken[0]!r} already exists; pass --force to replace it, or "
            "import it under another name with --as"
        )
    previous: dict[Path, str | None] = {}
    written: list[Path] = []
    try:
        for item in incoming:
            # the file that is there, which may be an exported NAME.vcfclick-preset.json
            path = user_file(item.name) or user_path(item.name)
            previous[path] = path.read_text("utf-8") if path.exists() else None
            write_text_atomic(path, to_text(item.doc))
            written.append(path)
    except BaseException:
        for path in written:
            old = previous[path]
            with contextlib.suppress(OSError):
                if old is None:
                    path.unlink()
                else:
                    write_text_atomic(path, old)
        raise
    return written


# ─────────────────────────────── from a result ───────────────────────────────


def filters_from_result(text: str | bytes) -> tuple[dict, dict]:
    """Rebuild a preset's filters from a saved `db trio --format json` result:
    the values that took effect from the preset or from a flag (a default is
    not frozen into a preset). Returns (filters, the result's `preset` header:
    name, source, revision, fingerprint). A result of several probands uses the
    first document that has a `preset` block."""
    doc = parse_json(text)
    docs = doc.get("trios") if isinstance(doc, dict) and "trios" in doc else [doc]
    block = next(
        (d["preset"] for d in docs if isinstance(d, dict) and "preset" in d), None
    )
    if not isinstance(block, dict) or not isinstance(block.get("resolved"), dict):
        raise PresetError(
            "this result has no preset block (it was produced without --preset); "
            "for a flag-only run use `vcfclick preset save NAME` with the same "
            "flags"
        )
    filters = {
        key: entry["value"]
        for key, entry in block["resolved"].items()
        if isinstance(entry, dict) and entry.get("from") in ("preset", "flag")
    }
    return filters, {
        k: block.get(k) for k in ("name", "source", "revision", "fingerprint")
    }


def revision_change(header: dict) -> str | None:
    """`rare-recessive rev 1 a1b2c3 -> now rev 2 d4e5f6` when the built-in a
    result used has changed since; None otherwise."""
    name = header.get("name")
    if header.get("source") != "builtin" or name not in BUILTINS:
        return None
    now = BUILTINS[name]
    current = now.history[-1]["fingerprint"]
    if current == header.get("fingerprint"):
        return None
    return (
        f"{name} rev {header.get('revision')} {header.get('fingerprint')} -> "
        f"now rev {now.revision} {current}"
    )
