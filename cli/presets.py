"""The preset store: where presets live, how one is found and applied.

Files, not tables. A preset is a JSON file (`cli.preset_schema` says what is in
one) in one of these places, looked up in this order for a bare name:

  * the user's directory, `$VCFCLICK_HOME/presets/<name>.json` (read-write);
  * every directory on `VCFCLICK_PRESET_PATH` (os.pathsep-separated; a lab's git
    folder; read-only here), where `<name>.vcfclick-preset.json` (what
    `preset export -o DIR` writes) is the same preset as `<name>.json`;
  * the built-ins compiled into the CLI (`cli.builtin_presets`; read-only).

A built-in name belongs to the built-in: a user or shared file with that name
is listed with an error and never used. A name in both the user directory and
a shared one resolves to the user's with a note. An explicit file path
(`--preset ./x.json`) is `source=file`, skips the lookup and exists on the
command line only: the name-only resolver the web and MCP front ends use
(`load_named`) rejects anything that is not a name before it touches the file
system, and no error message quotes a file's content.

Every lookup, rewrite and delete takes a NAME that has passed `name_problem`,
whatever the command: a value with a separator, a leading dot or an absolute
path never reaches the file system, so `preset rm ../x` cannot leave the
presets directory.

Writes are atomic (a temp file in the same directory, fsync, os.replace); there
is no lock and the last writer wins. A corrupt file is reported by `list` and
never deleted.

The pure rules are in `cli.preset_schema` and `cli.filter_registry`; the
Click side is `cli.preset_options` (applying a preset) and `cli.preset_cmds`
(the `preset` command group).
"""

from __future__ import annotations

import contextlib
import difflib
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from cli import filter_registry as fr
from cli.builtin_presets import BUILTINS, BuiltinPreset
from cli.preset_schema import (
    ERROR,
    EXPORT_STEM,
    EXTENSION_PREFIX,
    MAX_BYTES,
    SCHEMA,
    SCHEMA_VERSION,
    WARN,
    PresetError,
    Problem,
    fingerprint,
    name_problem,
    named,
    now_utc,
    parse_json,
    to_text,
    validate_document,
    vcfclick_version,
)

BUILTIN, USER, SHARED, FILE = "builtin", "user", "shared", "file"
FINGERPRINT_LEN = 12


# ───────────────────────────── a located preset ─────────────────────────────


@dataclass(frozen=True)
class Preset:
    """A preset found somewhere. `problems` may hold errors: such a preset is
    listed, never applied."""

    name: str
    source: str  # builtin | user | shared | file
    doc: dict  # the normalised envelope ({} when the file was unreadable)
    path: Path | None = None
    problems: tuple[Problem, ...] = ()
    unknown_keys: tuple[str, ...] = ()
    builtin: BuiltinPreset | None = field(default=None, compare=False)

    @property
    def filters(self) -> dict:
        return self.doc.get("filters", {})

    @property
    def meta(self) -> dict:
        return self.doc.get("meta", {})

    @property
    def has_errors(self) -> bool:
        return any(p.level == ERROR for p in self.problems)

    @property
    def available(self) -> bool:
        """Can this binary apply every key in it?"""
        return not self.has_errors and not self.unknown_keys

    @property
    def fingerprint(self) -> str | None:
        return None if self.has_errors or not self.doc else fingerprint(self.doc)

    # The text fields are read defensively: a document that failed validation is
    # still listed (with its problem), so whatever a file held must not crash
    # `preset list` or `show`.
    @property
    def title(self) -> str:
        value = self.doc.get("title")
        return value if isinstance(value, str) and value else self.name

    @property
    def description(self) -> str:
        value = self.doc.get("description")
        return value if isinstance(value, str) else ""

    @property
    def applies_to(self) -> list[str]:
        value = self.doc.get("applies_to")
        if not isinstance(value, list):
            return []
        return [a for a in value if isinstance(a, str)]

    @property
    def assay(self) -> str | None:
        value = self.doc.get("assay")
        return value if isinstance(value, str) else None

    @property
    def locked(self) -> bool:
        return bool(self.meta.get("locked", False))

    @property
    def revision(self) -> int | None:
        return self.builtin.revision if self.builtin else None

    @property
    def status(self) -> str | None:
        return self.builtin.status if self.builtin else None

    @property
    def categories(self) -> list[str]:
        return fr.category_tokens(self.filters.get("category", "all")) or ["all"]


def builtin_preset(name: str) -> Preset:
    b = BUILTINS[name]
    checked = validate_document(b.document())
    return Preset(name, BUILTIN, checked.doc, None, checked.problems, (), b)


# ───────────────────────────── where they live ─────────────────────────────


def user_dir() -> Path:
    from storage.db import _vcfclick_home

    return _vcfclick_home() / "presets"


def shared_dirs() -> list[Path]:
    raw = os.environ.get("VCFCLICK_PRESET_PATH", "")
    return [Path(p).expanduser() for p in raw.split(os.pathsep) if p.strip()]


def preset_name_of(path: Path) -> str:
    """The name a preset file stands for: its name without `.json`, and
    without `.vcfclick-preset` (the suffix `preset export -o DIR` adds)."""
    return path.stem.removesuffix(EXPORT_STEM)


def candidate_files(directory: Path, name: str) -> list[Path]:
    """The files in `directory` that are the preset `name`, the plain one
    first. `name` must already have passed `name_problem`."""
    return [directory / f"{name}.json", directory / f"{name}{EXPORT_STEM}.json"]


def read_preset_file(path: Path, source: str) -> Preset:
    """Read and validate one file. Never raises for a bad file: the problems
    are on the result. (Unknown filter keys are warnings here: the preset is
    stored, and `plan_apply` decides what to do about them.)"""
    name = preset_name_of(path)
    problems: list[Problem] = []

    def fail(code: str, message: str) -> Preset:
        problems.append(Problem(ERROR, code, message))
        return Preset(name, source, {}, path, tuple(problems))

    try:
        if path.stat().st_size > MAX_BYTES:
            return fail("too_large", f"file is larger than {MAX_BYTES >> 20} MiB")
        raw = path.read_bytes()
    except OSError as e:
        return fail("unreadable", f"cannot read the file ({e.strerror or 'error'})")
    try:
        doc = parse_json(raw)
    except PresetError as e:
        return fail("invalid_json", str(e))
    checked = validate_document(doc, allow_unknown=True)
    problems.extend(checked.problems)
    doc = checked.doc
    if source in (USER, SHARED):
        why = name_problem(name)
        if why:
            problems.append(Problem(ERROR, "bad_name", f"file name: {why}"))
        elif name in BUILTINS:
            problems.append(
                Problem(
                    ERROR,
                    "builtin_name",
                    f"{name!r} is a built-in preset's name; the built-in is used "
                    "and this file is ignored",
                )
            )
        declared = checked.doc.get("name")
        if declared and declared != name:
            problems.append(
                Problem(
                    WARN,
                    "name_mismatch",
                    f"the file says its name is {declared!r}; it is used as {name!r}",
                    "name",
                )
            )
        if doc and not why:
            # the file name decides what a stored preset is called, so every
            # writer (set, edit, rename, export) must write that name, not the
            # one a copied file still declares
            doc = named(doc, name)
    return Preset(name, source, doc, path, tuple(problems), checked.unknown_keys)


def _scan(directory: Path, source: str) -> list[Preset]:
    try:
        files = sorted(
            p
            for p in directory.glob("*.json")
            if p.is_file() and not p.name.startswith(".")
        )
    except OSError:
        return []
    out: list[Preset] = []
    taken: dict[str, str] = {}
    for p in files:
        found = read_preset_file(p, source)
        if found.name in taken:  # NAME.json and NAME.vcfclick-preset.json
            found = _with_problem(
                found,
                Problem(
                    WARN,
                    "duplicate_name",
                    f"the same preset name as {taken[found.name]}, which is used",
                ),
            )
        else:
            taken[found.name] = p.name
        out.append(found)
    return out


def list_presets(source: str = "all") -> list[Preset]:
    """Every preset, built-ins first, then the user's, then each shared
    directory's. A shared preset shadowed by a user one of the same name
    carries a warning."""
    out: list[Preset] = []
    if source in ("all", BUILTIN):
        out += [builtin_preset(n) for n in BUILTINS]
    user = _scan(user_dir(), USER)
    if source in ("all", USER):
        out += user
    if source in ("all", SHARED):
        seen = {p.name for p in user}
        shadowed = Problem(
            WARN,
            "shadowed",
            "shadowed by a preset of the same name earlier on the lookup path",
        )
        for directory in shared_dirs():
            for p in _scan(directory, SHARED):
                out.append(_with_problem(p, shadowed) if p.name in seen else p)
                seen.add(p.name)
    return out


def _with_problem(p: Preset, problem: Problem) -> Preset:
    return Preset(
        p.name, p.source, p.doc, p.path, (*p.problems, problem), p.unknown_keys,
        p.builtin,
    )  # fmt: skip


def all_names() -> list[str]:
    return [p.name for p in list_presets()]


def _first_file(directory: Path, name: str) -> Path | None:
    for candidate in candidate_files(directory, name):
        if candidate.is_file():
            return candidate
    return None


def find_preset(name: str) -> tuple[Preset | None, list[str]]:
    """The preset a bare name means, and notes to show: built-in, else the
    user's, else the first shared directory that has it. A value that is not a
    valid name finds nothing, without touching the file system."""
    if name_problem(name):
        return None, []
    if name in BUILTINS:
        return builtin_preset(name), []
    notes: list[str] = []
    mine = _first_file(user_dir(), name)
    user = read_preset_file(mine, USER) if mine is not None else None
    shared = None
    for directory in shared_dirs():
        candidate = _first_file(directory, name)
        if candidate is not None:
            shared = read_preset_file(candidate, SHARED)
            break
    if user is not None and shared is not None:
        notes.append(
            f"preset {name!r} exists in your presets and in a shared directory; "
            "using yours"
        )
    return (user or shared), notes


# ───────────────────────────── finding one to use ─────────────────────────────


@dataclass(frozen=True)
class Located:
    preset: Preset
    pinned: bool  # a fingerprint pin was given and matched
    notes: tuple[str, ...] = ()


def looks_like_path(value: str) -> bool:
    return "/" in value or "\\" in value or value.endswith(".json")


def _usable(preset: Preset, label: str) -> None:
    errors = [p for p in preset.problems if p.level == ERROR]
    if errors:
        raise PresetError(f"{label} cannot be used: {errors[0].message}")


def _check_pin(preset: Preset, pin: str) -> None:
    if len(pin) != FINGERPRINT_LEN or any(c not in "0123456789abcdef" for c in pin):
        raise PresetError(
            f"{pin[:20]!r} is not a fingerprint ({FINGERPRINT_LEN} hex digits)"
        )
    actual = preset.fingerprint
    if actual == pin:
        return
    message = (
        f"preset {preset.name!r} has fingerprint {actual}, not the pinned {pin}: "
        "its filters are not the ones that were pinned"
    )
    if preset.builtin is not None:
        for entry in preset.builtin.history:
            if entry["fingerprint"] == pin:
                message += (
                    f". Built-in {preset.name!r} revision {entry['revision']} "
                    f"had that fingerprint; the current revision is "
                    f"{preset.builtin.revision}"
                )
    raise PresetError(message)


def _not_found(name: str) -> PresetError:
    if name_problem(name):
        return PresetError(f"{name[:40]!r} is not a preset name")
    near = difflib.get_close_matches(name, all_names(), n=3, cutoff=0.6)
    hint = f" Did you mean {', '.join(near)}?" if near else ""
    return PresetError(
        f"no preset named {name!r}.{hint} `vcfclick preset list` shows the ones "
        "there are."
    )


def load_named(value: str, *, pin: str | None = None) -> Located:
    """The name-only resolver (web, MCP, and the CLI's name form): `NAME` or
    `NAME@FINGERPRINT`. Anything that is not a name matching the name pattern
    is refused BEFORE the file system is touched, so a path (`/etc/passwd`,
    `../x`, `a/b`, `x.json`) can never be read this way."""
    name, at, tail = value.partition("@")
    if at:
        if pin is not None and pin != tail:
            raise PresetError("two different fingerprints were given")
        pin = tail
    if name_problem(name):
        raise PresetError(f"{name[:40]!r} is not a preset name")
    preset, notes = find_preset(name)
    if preset is None:
        raise _not_found(name)
    _usable(preset, f"preset {name!r}")
    if pin is not None:
        _check_pin(preset, pin)
    return Located(preset, pin is not None, tuple(notes))


def load_file(path: str, *, pin: str | None = None) -> Located:
    """A preset from an explicit file (CLI only)."""
    p = Path(path).expanduser()
    if not p.is_file():
        raise PresetError(f"no preset file {p.name!r}")
    preset = read_preset_file(p, FILE)
    _usable(preset, f"preset file {p.name!r}")
    if pin is not None:
        _check_pin(preset, pin)
    return Located(preset, pin is not None)


def locate(value: str, *, pin: str | None = None, allow_path: bool = True) -> Located:
    """`--preset VALUE` on the command line: a file when it looks like a path
    (and `allow_path`), otherwise a name, optionally `NAME@FINGERPRINT`."""
    if looks_like_path(value):
        if not allow_path:
            raise PresetError(f"{value[:40]!r} is not a preset name")
        return load_file(value, pin=pin)
    return load_named(value, pin=pin)


# ───────────────────────────── applying one ─────────────────────────────


@dataclass(frozen=True)
class ApplyPlan:
    """What a command does with a preset's keys."""

    applicable: dict  # file key -> value, every key this command takes
    ignored: tuple[str, ...]  # of those, switched off by --preset-ignore
    not_applicable: tuple[str, ...]  # known keys this command does not take
    unavailable: tuple[str, ...]  # keys this binary does not know

    @property
    def defaults(self) -> dict:
        """The values to feed Click as defaults."""
        return {k: v for k, v in self.applicable.items() if k not in self.ignored}


def check_ignore_keys(keys: tuple[str, ...]) -> None:
    for key in keys:
        if key not in fr.REGISTRY:
            near = fr.suggest(key)
            hint = f" Did you mean {', '.join(near)}?" if near else ""
            raise PresetError(f"--preset-ignore {key[:40]!r}: not a filter key.{hint}")


def plan_apply(
    preset: Preset,
    command: str,
    *,
    ignore: tuple[str, ...] = (),
    skip_unavailable: bool = False,
) -> ApplyPlan:
    """Split a preset's keys for `command`. A key this binary does not know is
    an error unless `skip_unavailable`; a known key the command does not take
    is skipped (and the caller says so, loudly when it would widen the
    search). Extension (x_) keys and descriptive (meta) keys are not applied."""
    check_ignore_keys(ignore)
    if preset.unknown_keys and not skip_unavailable:
        keys = ", ".join(repr(k[:40]) for k in preset.unknown_keys)
        raise PresetError(
            f"preset {preset.name!r} has keys this vcfclick does not know: "
            f"{keys}. Upgrade vcfclick, or pass --preset-skip-unavailable to "
            "apply the rest (those filters are then NOT applied)"
        )
    applicable: dict[str, Any] = {}
    not_applicable: list[str] = []
    for key, value in preset.filters.items():
        if key.startswith(EXTENSION_PREFIX) or key in preset.unknown_keys:
            continue
        spec = fr.REGISTRY[key]
        if spec.kind != fr.VALUE:
            continue
        if command not in spec.commands:
            not_applicable.append(key)
        else:
            applicable[key] = value
    ignored = tuple(k for k in dict.fromkeys(ignore) if k in applicable)
    return ApplyPlan(
        applicable, ignored, tuple(not_applicable), tuple(preset.unknown_keys)
    )


# ───────────────────────────── writing ─────────────────────────────


def write_text_atomic(path: Path, text: str) -> None:
    """Write `text` to `path` so that a reader sees the old file or the new
    one, never half of one: a temp file in the same directory, fsync, then
    os.replace."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".tmp-", suffix=".json")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        with contextlib.suppress(OSError):
            os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


def user_path(name: str) -> Path:
    """Where the user's preset `name` is written. Refuses anything that is not
    a preset name, so no caller can aim a write or a delete outside the
    presets directory."""
    why = name_problem(name)
    if why:
        raise PresetError(f"{str(name)[:40]!r} is not a preset name: {why}")
    return user_dir() / f"{name}.json"


def user_file(name: str) -> Path | None:
    """The file the user's preset `name` is in, whichever spelling it has
    (`NAME.json`, or the `NAME.vcfclick-preset.json` an export writes), or None.
    This, not `user_path(name).exists()`, is what "does that preset exist" asks:
    a write that only looked for `NAME.json` would put a second file beside an
    exported one and shadow it."""
    why = name_problem(name)
    if why:
        raise PresetError(f"{str(name)[:40]!r} is not a preset name: {why}")
    return _first_file(user_dir(), name)


def new_document(
    name: str,
    filters: dict,
    *,
    title: str | None = None,
    description: str | None = None,
    applies_to: list[str] | None = None,
    assay: str | None = None,
    extra: dict | None = None,
) -> dict:
    """A fresh envelope, validated (PresetError on the first problem)."""
    doc: dict[str, Any] = {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "name": name,
    }
    for key, value in (
        ("title", title),
        ("description", description),
        ("applies_to", applies_to),
        ("assay", assay),
    ):
        if value is not None:
            doc[key] = value
    doc["filters"] = dict(filters)
    doc["meta"] = {"created_with": vcfclick_version(), "updated": now_utc()}
    doc.update(extra or {})
    return validated_or_raise(doc)


def validated_or_raise(doc: dict, *, allow_unknown: bool = False) -> dict:
    checked = validate_document(doc, allow_unknown=allow_unknown)
    if not checked.ok:
        raise PresetError(checked.errors[0].message)
    return checked.doc


def save_user_preset(doc: dict, *, force: bool = False) -> Path:
    """Write `doc` as the user's preset of its name. Refuses a built-in name
    and, without `force`, an existing one."""
    name = doc.get("name")
    why = name_problem(name)
    if why:
        raise PresetError(why)
    if name in BUILTINS:
        raise PresetError(
            f"{name!r} is a built-in preset; choose another name (copy it with "
            f"`vcfclick preset save NEWNAME --from-preset {name}`)"
        )
    existing = user_file(name)
    if existing is not None and not force:
        raise PresetError(f"preset {name!r} already exists; pass --force to replace it")
    path = existing or user_path(name)  # replace the file that is there, not a twin
    write_text_atomic(path, to_text(doc))
    return path


def require_editable(preset: Preset, *, force: bool = False) -> None:
    """Built-in and shared presets are read-only; a locked one needs --force."""
    if preset.source == BUILTIN:
        raise PresetError(
            f"{preset.name!r} is a built-in preset and cannot be changed; copy "
            f"it with `vcfclick preset save NEWNAME --from-preset {preset.name}`"
        )
    if preset.source != USER:
        raise PresetError(
            f"{preset.name!r} is a shared preset (read-only here); copy it with "
            f"`vcfclick preset save NEWNAME --from-preset {preset.name}`"
        )
    if preset.locked and not force:
        raise PresetError(
            f"preset {preset.name!r} is locked (meta.locked); pass --force to change it"
        )


def get_user_preset(name: str) -> Preset:
    """The preset `name` for a change: a built-in or shared one is returned
    so `require_editable` can explain; a missing one is an error."""
    preset, _ = find_preset(name)
    if preset is None:
        raise _not_found(name)
    return preset


def touch(doc: dict) -> dict:
    """`doc` with meta.updated set to now and created_with kept."""
    meta = dict(doc.get("meta", {}))
    meta["updated"] = now_utc()
    meta.setdefault("created_with", vcfclick_version())
    return {**doc, "meta": meta}
