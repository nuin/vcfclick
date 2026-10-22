"""Click glue for filters and presets: one decorator pair, one precedence rule.

    @db.command(name="trio", cls=PresetCommand)
    @gate_options("trio")      # every filter of the command, from the registry
    @preset_option()           # --preset, --preset-fingerprint, --preset-ignore, ...
    @command_options(TrioOptions)
    def db_trio(options): ...

Precedence is: a flag typed on the command line (or in the environment) >
the preset > the default. Click already has the machinery, so a preset is fed
to it as `ctx.default_map`: `PresetCommand.parse_args` pre-scans the four
preset options, locates the preset, and sets the map BEFORE the real parse,
which means option order on the command line does not matter and
`--preset-ignore` is always seen first. The real parse still declares the
options (so `--help` and the apps' capability probes see them), converts a
preset value with the same ParamType as a typed flag, and tells the sources
apart: COMMANDLINE / ENVIRONMENT is a flag, DEFAULT_MAP a preset, DEFAULT a
default.

A preset-supplied value is never "explicit". `is_explicit(ctx, name)` is the
one test (replacing the `!= ParameterSource.DEFAULT` that would count
DEFAULT_MAP as typed): a preset's max_af must not defeat the cohort-derived-AF
skip that a default max_af is subject to, and only a typed `--max-af` forces
the gate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import click

from cli import filter_registry as fr
from cli.filter_options import (
    EXPLICIT_META,
    CategoryList,
    gate_options,
    is_explicit,
)
from cli.preset_schema import PresetError
from cli.presets import ApplyPlan, Located, locate, plan_apply

# gate_options, CategoryList and is_explicit live in cli/filter_options.py (they
# need no preset store); this module is where a command takes them from.
__all__ = [
    "EXPLICIT_META",
    "CategoryList",
    "PresetCommand",
    "PresetState",
    "gate_options",
    "is_explicit",
    "preset_option",
    "preset_state",
    "resolve_filters",
    "scan_preset_args",
]

PRESET_META = "vcfclick.preset"

_VALUE_OPTIONS = ("--preset", "--preset-fingerprint", "--preset-ignore")
_FLAG_OPTIONS = ("--preset-skip-unavailable",)


@dataclass(frozen=True)
class PresetScan:
    preset: str | None = None
    fingerprint: str | None = None
    ignore: tuple[str, ...] = ()
    skip_unavailable: bool = False

    @property
    def any(self) -> bool:
        return bool(
            self.preset or self.fingerprint or self.ignore or self.skip_unavailable
        )


def scan_preset_args(args: list[str]) -> PresetScan:
    """Find the preset options in a raw argument list without parsing the
    rest: the last --preset and --preset-fingerprint win, --preset-ignore
    repeats. `--` ends the scan."""
    found: dict[str, Any] = {"ignore": []}
    i = 0
    while i < len(args):
        arg = args[i]
        if arg == "--":
            break
        name, eq, value = arg.partition("=") if arg.startswith("--") else (arg, "", "")
        if name in _VALUE_OPTIONS:
            if not eq:
                i += 1
                if i >= len(args):
                    break  # the real parse reports the missing value
                value = args[i]
            if name == "--preset":
                found["preset"] = value
            elif name == "--preset-fingerprint":
                found["fingerprint"] = value
            else:
                found["ignore"].append(value)
        elif name in _FLAG_OPTIONS:
            found["skip_unavailable"] = True
        i += 1
    return PresetScan(
        found.get("preset"),
        found.get("fingerprint"),
        tuple(found["ignore"]),
        found.get("skip_unavailable", False),
    )


@dataclass(frozen=True)
class PresetState:
    """The preset a run uses, located and split for its command."""

    command: str
    located: Located
    plan: ApplyPlan
    ignore: tuple[str, ...]
    skip_unavailable: bool = False
    notes: tuple[str, ...] = field(default=())


class PresetCommand(click.Command):
    """A command that takes `--preset`. Its name is the command id the
    registry knows (`trio`, `segregate`, `case`)."""

    def default_layer(self, ctx: click.Context, args: list[str]) -> dict:
        """Defaults that sit between the built-in defaults and a preset (the
        case file of the case command). Keyed by Click parameter name. A
        subclass may also add parameter names to `ctx.meta[EXPLICIT_META]`."""
        return {}

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        defaults = dict(self.default_layer(ctx, args))
        helping = any(a in ctx.help_option_names for a in args)
        scan = PresetScan() if helping else scan_preset_args(args)
        if scan.any:
            if scan.preset is None:
                raise click.UsageError(
                    "--preset-fingerprint, --preset-ignore and "
                    "--preset-skip-unavailable need --preset",
                    ctx,
                )
            try:
                state = build_state(self.name or "", scan)
            except PresetError as e:
                raise click.ClickException(str(e)) from None
            ctx.meta[PRESET_META] = state
            for key, value in state.plan.defaults.items():
                defaults[fr.REGISTRY[key].param] = value
        if defaults:
            ctx.default_map = {**(ctx.default_map or {}), **defaults}
        return super().parse_args(ctx, args)


def build_state(command: str, scan: PresetScan) -> PresetState:
    """Locate the preset named by `scan` and split its keys for `command`.
    Raises PresetError."""
    if scan.preset is None:
        raise PresetError("no preset was named")
    located = locate(scan.preset, pin=scan.fingerprint)
    plan = plan_apply(
        located.preset,
        command,
        ignore=scan.ignore,
        skip_unavailable=scan.skip_unavailable,
    )
    return PresetState(
        command,
        located,
        plan,
        scan.ignore,
        scan.skip_unavailable,
        located.notes,
    )


def preset_state(ctx: click.Context | None = None) -> PresetState | None:
    """The preset this run uses, or None."""
    ctx = ctx or click.get_current_context()
    return ctx.meta.get(PRESET_META)


def resolve_filters(
    ctx: click.Context, command: str, state: PresetState | None
) -> tuple[dict[str, fr.Resolved], list[dict]]:
    """Walk the registry: for each key of `command`, its final value and where
    it came from (flag, preset, default), plus the preset keys a flag or
    `--preset-ignore` overrode. The same function the web and MCP resolvers use
    (`fr.resolve_filters_from_values`), fed from this Click context."""
    explicit = {
        key.file_key: ctx.params[key.param]
        for key in fr.keys_for(command)
        if is_explicit(ctx, key.param)
    }
    filters = state.plan.applicable if state else {}
    ignore = state.plan.ignored if state else ()
    resolved, overridden = fr.resolve_filters_from_values(
        command, filters, explicit, ignore
    )
    if state:
        problems = fr.cross_key_problems({k: r.value for k, r in resolved.items()})
        if problems:
            raise click.UsageError(f"{'; '.join(problems)} (after applying the preset)")
    return resolved, overridden


# ───────────────────────────── the decorators ─────────────────────────────


def preset_option():
    """`--preset NAME|NAME@FINGERPRINT|PATH` and its three companions."""

    def decorate(fn):
        options = [
            click.option(
                "--preset",
                "preset",
                default=None,
                metavar="NAME|PATH",
                help="Apply a filter preset (`vcfclick preset list`): its values "
                "are the defaults, and a flag typed here wins. NAME@FINGERPRINT "
                "fails the run if the preset's content changed; a value with a "
                "path separator or ending .json is a file.",
            ),
            click.option(
                "--preset-fingerprint",
                default=None,
                metavar="HEX12",
                help="Fail unless the preset's fingerprint is this (also for a "
                "preset given as a path).",
            ),
            click.option(
                "--preset-ignore",
                "preset_ignore",
                multiple=True,
                metavar="KEY",
                help="Switch one preset key off for this run (repeatable).",
            ),
            click.option(
                "--preset-skip-unavailable",
                is_flag=True,
                help="Apply the preset keys this vcfclick knows and skip the "
                "rest (they are NOT applied, and the output says so); "
                "the default is an error.",
            ),
        ]
        for option in reversed(options):
            fn = option(fn)
        return fn

    return decorate
