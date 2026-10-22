"""`vcfclick preset ...` — named, versioned, shareable filter sets.

    vcfclick preset list
    vcfclick preset show strict-denovo
    vcfclick preset save id-panel --from-preset rare-recessive --max-af 0.003
    vcfclick preset set id-panel max_af=0.002
    vcfclick preset export id-panel -o id-panel.vcfclick-preset.json
    vcfclick preset import id-panel.vcfclick-preset.json

None of these opens a database. The model is in `cli.presets` (the store) and
`cli.preset_schema` (the file); using a preset is `db trio --preset NAME`
(`cli.preset_options`). docs/PRESETS.md is the guide.
"""

from __future__ import annotations

import json
import sys
from collections.abc import Callable
from contextlib import contextmanager
from pathlib import Path

import click
from click.core import ParameterSource

from cli import filter_registry as fr
from cli.builtin_presets import BUILTINS
from cli.filter_options import gate_options
from cli.main import cli
from cli.preset_exchange import (
    EXPORT_SUFFIX,
    bundle_document,
    export_document,
    export_target,
    filters_from_result,
    import_presets,
    load_for_export,
    read_incoming,
    revision_change,
)
from cli.preset_schema import (
    ERROR,
    EXTENSION_PREFIX,
    PresetError,
    Problem,
    fingerprint,
    name_problem,
    parse_json,
    to_text,
    validate_document,
)
from cli.presets import (
    FILE,
    Preset,
    find_preset,
    get_user_preset,
    list_presets,
    looks_like_path,
    new_document,
    read_preset_file,
    require_editable,
    save_user_preset,
    shared_dirs,
    touch,
    user_dir,
    user_path,
    validated_or_raise,
    write_text_atomic,
)

_SOURCES = ("builtin", "user", "shared", "all")


@contextmanager
def _errors():
    """A PresetError is the user's to fix: report it like any CLI error."""
    try:
        yield
    except PresetError as e:
        raise click.ClickException(str(e)) from None


@cli.group(name="preset")
def preset() -> None:
    """Filter presets: named sets of `db trio` gates you can share.

    A preset holds filter values only (no samples, paths or database). Use one
    with `db trio --preset NAME`; a flag typed next to it wins. `preset list`
    shows the built-in starting points and your own."""


# ─────────────────────────────── list / show ───────────────────────────────


def _summary(p: Preset) -> dict:
    keys = [k for k in p.filters if not k.startswith(EXTENSION_PREFIX)]
    return {
        "name": p.name,
        "title": p.title,
        "description": p.description,
        "source": p.source,
        "revision": p.revision,
        "status": p.status,
        "fingerprint": p.fingerprint,
        "applies_to": p.applies_to,
        "assay": p.assay,
        "categories": p.categories,
        "keys": keys,
        "patient_specific_keys": [
            k for k in keys if k in fr.REGISTRY and fr.REGISTRY[k].patient_specific
        ],
        "available": p.available,
        "locked": p.locked,
        "problems": [x.to_json() for x in p.problems],
    }


@preset.command(
    name="list", short_help="List the built-in, your own and shared presets."
)
@click.option(
    "--source",
    type=click.Choice(_SOURCES),
    default="all",
    show_default=True,
    help="Which presets: the built-ins, your own, the shared directories "
    "(VCFCLICK_PRESET_PATH), or all.",
)
@click.option(
    "--command",
    "command",
    type=click.Choice(fr.COMMANDS),
    default=None,
    help="Only presets meant for this command (a hint the author set).",
)
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
)
def list_cmd(source: str, command: str | None, out_format: str) -> None:
    """List presets. Always exits 0: a file that cannot be read is listed with
    its problem, never hidden."""
    items = [
        p
        for p in list_presets(source)
        if command is None or not p.applies_to or command in p.applies_to
    ]
    if out_format == "json":
        click.echo(
            json.dumps(
                {"schema_version": 1, "presets": [_summary(p) for p in items]},
                indent=2,
            )
        )
        return
    if not items:
        click.echo("no presets")
        return
    width = max(len(p.name) for p in items)
    for p in items:
        mark = "" if p.available else "  [unavailable]"
        status = "starting point, not validated" if p.status else p.source
        click.echo(
            f"{p.name:<{width}}  {p.source:<8}  {p.fingerprint or '-' * 12}  "
            f"{status}{mark}"
        )
        if p.title and p.title != p.name:
            click.echo(f"{'':<{width}}  {p.title}")
        for problem in p.problems:
            click.echo(f"{'':<{width}}  {problem.level}: {problem.message}")


def _locate_for_show(target: str) -> Preset:
    if looks_like_path(target):
        path = Path(target).expanduser()
        if not path.is_file():
            raise PresetError(f"no preset file {path.name!r}")
        return read_preset_file(path, FILE)
    if name_problem(target):
        raise PresetError(f"{target[:40]!r} is not a preset name")
    found, _ = find_preset(target)
    if found is None:
        raise PresetError(f"no preset named {target!r}; see `vcfclick preset list`")
    return found


def _flag_of(key: str) -> str | None:
    spec = fr.REGISTRY.get(key)
    return spec.flag if spec else None


def _value_text(value) -> str:
    return json.dumps(value) if not isinstance(value, str) else value


@preset.command(short_help="Show a preset, its limits and history.")
@click.argument("target", metavar="NAME|PATH")
@click.option(
    "--command",
    "command",
    type=click.Choice(fr.COMMANDS),
    default=None,
    help="Mark the keys this command does not take.",
)
@click.option("--resolve", is_flag=True, help="Show the flag each key maps to.")
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
)
def show(target: str, command: str | None, resolve: bool, out_format: str) -> None:
    """Show one preset: its values, what it covers and does not, its revision
    history and fingerprint."""
    with _errors():
        p = _locate_for_show(target)
    unused = [
        k
        for k in p.filters
        if command and k in fr.REGISTRY and command not in fr.REGISTRY[k].commands
    ]
    if out_format == "json":
        doc = {
            **_summary(p),
            "filters": p.filters,
            "covers": p.builtin.covers if p.builtin else None,
            "does_not_cover": p.builtin.does_not_cover if p.builtin else None,
            "excludes": dict(p.builtin.excludes) if p.builtin else {},
            "history": [dict(h) for h in p.builtin.history] if p.builtin else [],
            "meta": p.meta,
        }
        if resolve:
            doc["flags"] = {k: _flag_of(k) for k in p.filters if _flag_of(k)}
        if command:
            doc["not_applicable"] = unused
        click.echo(json.dumps(doc, indent=2))
        return
    head = f"{p.name}  ({p.source}"
    if p.revision is not None:
        head += f", revision {p.revision}"
    head += f", fingerprint {p.fingerprint or 'n/a'}"
    if p.status:
        head += ", starting point, not validated"
    click.echo(head + ")")
    if p.title != p.name:
        click.echo(p.title)
    if p.description:
        click.echo(p.description)
    click.echo("\nfilters:")
    for key, value in p.filters.items():
        line = f"  {key:<22} {_value_text(value)}"
        if resolve and _flag_of(key):
            line += f"    {_flag_of(key)}"
        if key in unused:
            line += f"    (not used by `{command}`)"
        click.echo(line)
    if p.builtin:
        click.echo(f"\ncovers: {p.builtin.covers}")
        click.echo(f"does not cover: {p.builtin.does_not_cover}")
        if p.builtin.excludes:
            click.echo("models not run:")
            for model, why in p.builtin.excludes.items():
                click.echo(f"  {model}: {why}")
        click.echo("history:")
        for h in p.builtin.history:
            click.echo(f"  revision {h['revision']}  {h['fingerprint']}  {h['note']}")
    meta = p.meta
    click.echo(
        f"\napproved by: {meta.get('approved_by') or '-'}  on: "
        f"{meta.get('approved_on') or '-'}  locked: {'yes' if p.locked else 'no'}"
    )
    for problem in p.problems:
        click.echo(f"{problem.level}: {problem.message}")


# ───────────────────────────────── save ─────────────────────────────────


def _typed_filters(ctx: click.Context, params: dict) -> dict:
    """The filters the user typed on `preset save`, by file key."""
    out = {}
    for key in fr.REGISTRY.values():
        if key.param and ctx.get_parameter_source(key.param) == (
            ParameterSource.COMMANDLINE
        ):
            out[key.file_key] = fr.normalise_value(key, params[key.param])
    return out


def _read_text(source: str) -> str:
    if source == "-":
        return sys.stdin.read()
    try:
        return Path(source).expanduser().read_text("utf-8")
    except (OSError, UnicodeDecodeError) as e:
        raise PresetError(f"cannot read {Path(source).name!r} ({e})") from None


_JSON_KEYS = ("title", "description", "applies_to", "assay")
_ENVELOPE_OWN = {"schema", "schema_version", "name", "filters", "meta", *_JSON_KEYS}


@preset.command()
@click.argument("name")
@click.option(
    "--from-preset",
    "from_preset",
    metavar="BASE",
    help="Start from this preset's filters (copied; there is no inheritance).",
)
@click.option(
    "--from-result",
    "from_result",
    metavar="RESULT.json",
    help="Rebuild the filters from a saved `db trio --format json` result that "
    "used a preset: the values that came from the preset or a flag.",
)
@click.option(
    "--from-json",
    "from_json",
    metavar="FILE|-",
    help='Read {"filters": {...}, "title": ..., "description": ...} from a '
    "file or stdin.",
)
@click.option("--title", default=None, help="A short title.")
@click.option("--description", default=None, help="What the preset is for.")
@click.option("--force", is_flag=True, help="Replace a preset of that name.")
@gate_options(None, defaults=False)
@click.pass_context
def save(
    ctx: click.Context,
    name: str,
    from_preset: str | None,
    from_result: str | None,
    from_json: str | None,
    title: str | None,
    description: str | None,
    force: bool,
    **flags,
) -> None:
    """Save a preset under NAME in your presets directory.

    The filter flags are the ones `db trio` takes: only those you type are
    saved (a value you leave out is inherited from the defaults when the
    preset is used). Built-in names are taken; so is an existing name unless
    --force."""
    with _errors():
        chosen = [s for s in (from_preset, from_result, from_json) if s]
        if len(chosen) > 1:
            raise click.UsageError(
                "choose one of --from-preset, --from-result and --from-json"
            )
        why = name_problem(name)
        if why:
            raise PresetError(why)
        base: dict = {}
        header = {"title": title, "description": description}
        extra: dict = {}
        if from_preset:
            src, _ = find_preset(from_preset)
            if src is None or src.has_errors:
                raise PresetError(f"no usable preset named {from_preset!r}")
            base = dict(src.filters)
            header["applies_to"] = src.applies_to or None
            header["assay"] = src.assay
        elif from_result:
            base, was = filters_from_result(_read_text(from_result))
            change = revision_change(was)
            if change:
                click.echo(f"note: {change}")
        elif from_json:
            doc = parse_json(_read_text(from_json))
            if not isinstance(doc, dict) or not isinstance(doc.get("filters"), dict):
                raise PresetError('--from-json needs an object with a "filters" object')
            base = doc["filters"]
            for k in _JSON_KEYS:
                if header.get(k) is None and doc.get(k) is not None:
                    header[k] = doc[k]
            extra = {k: v for k, v in doc.items() if k not in _ENVELOPE_OWN}
        filters = {**base, **_typed_filters(ctx, flags)}
        doc = new_document(
            name,
            filters,
            title=header.get("title"),
            description=header.get("description"),
            applies_to=header.get("applies_to"),
            assay=header.get("assay"),
            extra=extra,
        )
        path = save_user_preset(doc, force=force)
    click.echo(f"saved preset {name!r} (fingerprint {fingerprint(doc)}) to {path}")


# ─────────────────────────── set / unset / edit ───────────────────────────


def _parse_assignment(text: str) -> tuple[str, object]:
    key, eq, raw = text.partition("=")
    key, raw = key.strip(), raw.strip()
    if not eq or not key:
        raise PresetError(f"{text[:40]!r}: expected KEY=VALUE")
    spec = _known_key(key)
    if spec.type == "bool":
        words = {
            "true": True,
            "yes": True,
            "on": True,
            "false": False,
            "no": False,
            "off": False,
        }
        if raw.lower() not in words:
            raise PresetError(f"{key}: must be true or false")
        value: object = words[raw.lower()]
    elif spec.type == "category":
        value = json.loads(raw) if raw.startswith("[") else raw
    else:
        try:
            value = json.loads(raw)
        except ValueError:
            raise PresetError(f"{key}: {raw[:20]!r} is not a number") from None
    checked, bad = fr.check_value(spec, value)
    if bad:
        raise PresetError(f"{key}: {bad}")
    return key, checked


def _known_key(key: str) -> fr.FilterKey:
    if key in fr.REGISTRY:
        return fr.REGISTRY[key]
    if key in fr.NEVER_PRESET:
        raise PresetError(
            f"{key} can never be part of a preset: {fr.NEVER_PRESET[key]}"
        )
    near = fr.suggest(key)
    hint = f" Did you mean {', '.join(near)}?" if near else ""
    raise PresetError(f"unknown filter key {key[:40]!r}.{hint}")


def _file_of(p: Preset) -> Path:
    """The file a user preset lives in (`NAME.json`, or the exported
    `NAME.vcfclick-preset.json` someone dropped in the directory)."""
    if p.path is None:
        raise PresetError(f"preset {p.name!r} has no file to change")
    return p.path


def _change(name: str, force: bool, edit: Callable[[dict], dict]) -> Preset:
    """Apply `edit` to the filters of a user preset and write it back."""
    p = get_user_preset(name)
    require_editable(p, force=force)
    if p.has_errors:
        raise PresetError(
            f"preset {name!r} cannot be changed because it is not valid "
            f"({p.problems[0].message}); `vcfclick preset edit {name}` can fix it"
        )
    doc = touch({**p.doc, "filters": edit(dict(p.filters))})
    checked = validated_or_raise(doc, allow_unknown=True)
    write_text_atomic(_file_of(p), to_text(checked))
    return p


@preset.command(name="set", short_help="Set values in a preset of yours.")
@click.argument("name")
@click.argument("assignments", nargs=-1, required=True, metavar="KEY=VALUE...")
@click.option("--force", is_flag=True, help="Change a locked preset.")
def set_cmd(name: str, assignments: tuple[str, ...], force: bool) -> None:
    """Set filter values in a preset of yours: `preset set id-panel max_af=0.003`."""
    with _errors():
        pairs = [_parse_assignment(a) for a in assignments]
        _change(name, force, lambda f: {**f, **dict(pairs)})
    click.echo(f"set {', '.join(k for k, _ in pairs)} in {name!r}")


@preset.command(short_help="Remove keys from a preset of yours.")
@click.argument("name")
@click.argument("keys", nargs=-1, required=True)
@click.option("--force", is_flag=True, help="Change a locked preset.")
def unset(name: str, keys: tuple[str, ...], force: bool) -> None:
    """Remove filter keys from a preset of yours (they fall back to defaults)."""
    with _errors():
        for key in keys:
            if key not in fr.REGISTRY and not key.startswith(EXTENSION_PREFIX):
                _known_key(key)
        _change(name, force, lambda f: {k: v for k, v in f.items() if k not in keys})
    click.echo(f"unset {', '.join(keys)} in {name!r}")


@preset.command(short_help="Edit a preset of yours in $EDITOR.")
@click.argument("name")
@click.option("--force", is_flag=True, help="Edit a locked preset.")
def edit(name: str, force: bool) -> None:
    """Open a preset of yours in $EDITOR; it is validated before it is saved."""
    with _errors():
        p = get_user_preset(name)
        require_editable(p, force=force)
        path = _file_of(p)
        try:
            text = path.read_text("utf-8")
        except OSError as e:
            raise PresetError(f"cannot read the preset ({e.strerror})") from None
        try:
            edited = click.edit(text, extension=".json", require_save=True)
        except click.ClickException:
            click.echo(f"no editor could be started; edit {path} yourself")
            return
        if edited is None:
            click.echo("no changes")
            return
        checked = validate_document(parse_json(edited), allow_unknown=True)
        if not checked.ok:
            shown = "; ".join(e.message for e in checked.errors[:3])
            raise PresetError(f"not saved: {shown}")
        declared = checked.doc.get("name", name)
        if declared != name:
            raise PresetError(
                "not saved: the name cannot change here; use `vcfclick preset rename`"
            )
        doc = touch({**checked.doc, "name": name})
        write_text_atomic(path, to_text(doc))
    click.echo(f"saved {name!r} (fingerprint {fingerprint(doc)})")


@preset.command()
@click.argument("old")
@click.argument("new")
@click.option("--force", is_flag=True, help="Rename a locked preset.")
def rename(old: str, new: str, force: bool) -> None:
    """Rename a preset of yours."""
    with _errors():
        p = get_user_preset(old)
        require_editable(p, force=force)
        why = name_problem(new)
        if why:
            raise PresetError(why)
        if new in BUILTINS:
            raise PresetError(f"{new!r} is a built-in preset's name; choose another")
        if user_path(new).exists():
            raise PresetError(f"a preset named {new!r} already exists")
        doc = validated_or_raise(touch({**p.doc, "name": new}), allow_unknown=True)
        write_text_atomic(user_path(new), to_text(doc))
        _file_of(p).unlink()
    click.echo(f"renamed {old!r} to {new!r}")


@preset.command(name="rm", short_help="Delete a preset of yours.")
@click.argument("name")
@click.option("--yes", is_flag=True, help="Do not ask.")
@click.option("--force", is_flag=True, help="Delete a locked preset.")
def rm_cmd(name: str, yes: bool, force: bool) -> None:
    """Delete a preset of yours. `rm` and `delete` are the same command."""
    with _errors():
        p = get_user_preset(name)
        require_editable(p, force=force)
        if not yes:
            click.confirm(f"Delete preset {name!r}?", abort=True)
        _file_of(p).unlink()
    click.echo(f"deleted {name!r}")


preset.add_command(rm_cmd, name="delete")


# ───────────────────────────── export / import ─────────────────────────────


@preset.command()
@click.argument("names", nargs=-1, required=True)
@click.option(
    "-o",
    "--output",
    "out",
    metavar="FILE|DIR",
    help="Where to write: a file, or a directory (one file per preset). "
    "Without it the preset is printed.",
)
@click.option("--bundle", is_flag=True, help="Put several presets in one file.")
@click.option(
    "--with-patient-data",
    is_flag=True,
    help="Keep patient-specific keys (such as HPO terms); they are left out "
    "of an export by default.",
)
def export(
    names: tuple[str, ...], out: str | None, bundle: bool, with_patient_data: bool
) -> None:
    """Write presets to a file to hand to a colleague."""
    with _errors():
        bundle, target = export_target(list(names), out, bundle)
        docs = [
            export_document(load_for_export(n), with_patient_data=with_patient_data)
            for n in names
        ]
        if bundle:
            text, path = to_text(bundle_document(docs)), target
        elif target is not None and target.is_dir():
            for n, doc in zip(names, docs, strict=True):
                dest = target / f"{n}{EXPORT_SUFFIX}"
                write_text_atomic(dest, to_text(doc))
                click.echo(f"wrote {dest}")
            return
        else:
            text, path = to_text(docs[0]), target
        if path is None:
            click.echo(text, nl=False)
        else:
            write_text_atomic(path, text)
            click.echo(f"wrote {path}")


@preset.command(name="import")
@click.argument("files", nargs=-1, required=True)
@click.option("--as", "as_name", metavar="NAME", help="Store it under this name.")
@click.option("--force", is_flag=True, help="Replace a preset of the same name.")
@click.option(
    "--allow-unknown-keys",
    is_flag=True,
    help="Store a preset with filter keys this vcfclick does not know. Using it "
    "is then an error unless --preset-skip-unavailable.",
)
def import_cmd(
    files: tuple[str, ...], as_name: str | None, force: bool, allow_unknown_keys: bool
) -> None:
    """Import presets from files or bundles. All or nothing."""
    with _errors():
        incoming = read_incoming(
            list(files), as_name=as_name, allow_unknown=allow_unknown_keys
        )
        import_presets(incoming, force=force)
    for item in incoming:
        click.echo(f"imported {item.name!r} (fingerprint {fingerprint(item.doc)})")


# ───────────────────────────── validate / path ─────────────────────────────


@preset.command()
@click.argument("file", metavar="FILE")
@click.option(
    "--command",
    "command",
    type=click.Choice(fr.COMMANDS),
    default=None,
    help="Also report the keys this command does not take.",
)
@click.option(
    "--format",
    "out_format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
)
def validate(file: str, command: str | None, out_format: str) -> None:
    """Check a preset file. Exits 1 when it has an error."""
    problems: list[Problem] = []
    doc: dict = {}
    path = Path(file).expanduser()
    try:
        raw = path.read_bytes()
        checked = validate_document(parse_json(raw))
        problems = list(checked.problems)
        doc = checked.doc
    except OSError as e:
        problems.append(
            Problem(ERROR, "unreadable", f"cannot read the file ({e.strerror})")
        )
    except PresetError as e:
        problems.append(Problem(ERROR, "invalid", str(e)))
    if command and doc:
        for key in doc.get("filters", {}):
            spec = fr.REGISTRY.get(key)
            if spec and command not in spec.commands:
                problems.append(
                    Problem(
                        "warning",
                        "not_applicable",
                        f"{key} is not a filter of `{command}`",
                        key,
                    )
                )
    valid = not any(p.level == ERROR for p in problems)
    if out_format == "json":
        click.echo(
            json.dumps(
                {
                    "valid": valid,
                    "fingerprint": fingerprint(doc) if valid and doc else None,
                    "problems": [p.to_json() for p in problems],
                },
                indent=2,
            )
        )
    else:
        if valid:
            click.echo(f"ok: fingerprint {fingerprint(doc)}")
        for p in problems:
            click.echo(f"{p.level}: {p.message}")
    if not valid:
        sys.exit(1)


@preset.command(name="path", short_help="Print where presets are looked up.")
def path_cmd() -> None:
    """Print where presets are looked up: your directory, then each shared one
    (VCFCLICK_PRESET_PATH)."""
    click.echo(str(user_dir()))
    for directory in shared_dirs():
        click.echo(str(directory))
