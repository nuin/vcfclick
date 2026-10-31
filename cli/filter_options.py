"""Click options generated from the filter registry.

    @gate_options("trio")      # every filter `db trio` takes, from the registry

A command's filter flags are not written by hand: they are generated from
`cli/filter_registry.py`, so a flag typed on the command line and a value in a
preset are converted and checked by the same rules, and
`tests/test_filter_registry_drift.py` fails when the two drift apart.

`CategoryList` is the `--category` type (a model name, `all`, or a comma list),
and `is_explicit` the one test of "did the user type this?": a value that came
from a preset (Click's DEFAULT_MAP) or from the default does not count, so a
preset's max_af cannot defeat the cohort-derived-AF skip. See
`cli/preset_options.py` for how a preset is fed to Click.
"""

from __future__ import annotations

from typing import Any

import click
from click.core import ParameterSource

from cli import filter_registry as fr

# Parameter names a subclass of PresetCommand (the case command's case file)
# can mark explicit.
EXPLICIT_META = "vcfclick.explicit"


class CategoryList(click.Choice):
    """`--category`: one model, `all`, or a comma list (`recessive,comphet`).

    It is a click.Choice for the metavar, the help text and the wording of an
    invalid value (checked token by token), and converts to a tuple of names.
    A preset's list arrives here too, as a list."""

    def convert(self, value: Any, param, ctx) -> tuple[str, ...]:
        tokens = fr.category_tokens(value)
        if tokens is None:
            tokens = [""]  # an empty or odd value: Choice says what is allowed
        choice = super()
        return tuple([choice.convert(token, param, ctx) for token in tokens])


def is_explicit(ctx: click.Context, name: str) -> bool:
    """Did the user type this parameter (or set it in the environment)? A
    value that came from a preset's default map, or from the default, did not."""
    if name in ctx.meta.get(EXPLICIT_META, ()):
        return True
    return ctx.get_parameter_source(name) in (
        ParameterSource.COMMANDLINE,
        ParameterSource.ENVIRONMENT,
    )


def _option(key: fr.FilterKey, *, defaults: bool) -> click.Option:
    kwargs: dict[str, Any] = {"help": key.help or None}
    if key.type == "category":
        from cli.trio_analysis import CATEGORY_CHOICES

        kwargs["type"] = CategoryList(CATEGORY_CHOICES)
    elif key.type in ("int", "float"):
        kwargs["type"] = int if key.type == "int" else float
    elif key.type == "text" and key.choices:
        kwargs["type"] = click.Choice(key.choices())
    if key.multiple:
        kwargs["multiple"] = True
        kwargs["metavar"] = "SYMBOL[,SYMBOL...]" if key.type == "genes" else "REGION"
    if defaults:
        kwargs["default"] = None if key.default is fr.DERIVED else key.default
        kwargs["show_default"] = not key.is_bool and not key.multiple
    else:
        kwargs["default"] = None
    return click.option(key.flag, key.param, **kwargs)


def _file_option(key: fr.FilterKey) -> click.Option:
    flag, param = key.file_option
    what = "BED" if key.type == "regions" else "FILE"
    return click.option(
        flag,
        param,
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        metavar=what,
        help=f"Read {key.flag} values from a file (see {key.flag}).",
    )


def _keys(command: str | None) -> list[fr.FilterKey]:
    if command:
        return fr.keys_for(command)
    # `preset save`: every filter, and the descriptive keys a flag can set
    return [k for k in fr.REGISTRY.values() if k.kind == fr.VALUE or k.flag]


def gate_options(command: str | None, *, defaults: bool = True):
    """The Click options of every filter `command` takes, generated from the
    registry (`None`: every registered filter, and the descriptive keys such as
    `--genes-label`, for `preset save`), each list key followed by its file
    companion (`--gene-file`). `defaults=False` is for `preset save`, which
    must tell a typed flag from an absent one and so gives every option a
    default of None."""
    keys = _keys(command)

    def decorate(fn):
        for key in reversed(keys):
            if key.file_option:
                fn = _file_option(key)(fn)
            fn = _option(key, defaults=defaults)(fn)
        return fn

    return decorate


def file_option_tokens(key: fr.FilterKey, path: str) -> list[str]:
    """The tokens a list key's file companion holds: a gene file is one token
    per line (`#` starts a comment; commas and blanks also separate), a region
    file is BED (0-based half-open, turned into 1-based inclusive regions;
    `track`, `browser` and `#` lines are skipped). Raises click.BadParameter."""
    from pathlib import Path

    flag = key.file_option[0] if key.file_option else key.flag
    try:
        text = Path(path).read_text("utf-8")
    except (OSError, UnicodeDecodeError) as e:
        raise click.BadParameter(f"cannot read {Path(path).name!r} ({e})") from None
    out: list[str] = []
    for n, line in enumerate(text.splitlines(), 1):
        body = line.split("#", 1)[0].strip()
        if not body or body.startswith(("track", "browser")):
            continue
        if key.type == "genes":
            out += [t for t in body.replace(",", " ").split() if t]
            continue
        fields = body.split()
        try:
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        except (IndexError, ValueError):
            raise click.BadParameter(
                f"{Path(path).name} line {n} is not a BED interval (chrom start end)",
                param_hint=flag,
            ) from None
        if start < 0 or end <= start:
            raise click.BadParameter(
                f"{Path(path).name} line {n}: start must be >= 0 and below end",
                param_hint=flag,
            )
        out.append(f"{chrom}:{start + 1}-{end}")
    return out


def typed_list_values(ctx: click.Context, command: str | None) -> dict[str, list]:
    """For each list key of `command` whose flag or file companion was typed:
    the tokens typed, the flag's and the file's together (a file and the flag
    typed together are unioned). By file key."""
    out: dict[str, list] = {}
    for key in _keys(command):
        if not key.is_list:
            continue
        values: list = []
        typed = False
        if key.param in ctx.params and is_explicit(ctx, key.param):
            values += list(ctx.params[key.param] or ())
            typed = True
        if key.file_option:
            fparam = key.file_option[1]
            path = ctx.params.get(fparam)
            if path is not None and is_explicit(ctx, fparam):
                values += file_option_tokens(key, path)
                typed = True
        if typed:
            out[key.file_key] = fr.normalise_value(key, values)
    return out
