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
    if defaults:
        kwargs["default"] = None if key.default is fr.DERIVED else key.default
        kwargs["show_default"] = not key.is_bool
    else:
        kwargs["default"] = None
    return click.option(key.flag, key.param, **kwargs)


def gate_options(command: str | None, *, defaults: bool = True):
    """The Click options of every filter `command` takes, generated from the
    registry (`None`: every registered filter, for `preset save`).
    `defaults=False` is for `preset save`, which must tell a typed flag from an
    absent one and so gives every option a default of None."""
    keys = (
        fr.keys_for(command)
        if command
        else [k for k in fr.REGISTRY.values() if k.kind == fr.VALUE]
    )

    def decorate(fn):
        for key in reversed(keys):
            fn = _option(key, defaults=defaults)(fn)
        return fn

    return decorate
