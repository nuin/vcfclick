"""The filter registry and the Click options cannot drift apart.

`cli/filter_registry.py` is the one place a filter is declared; the options of
`db trio` are generated from it. This file fails when a later step adds a
Click option to a filter-bearing command without registering it (or listing it
as never presettable, with the reason), when a Click default differs from the
registry's, or when a registered key is a parameter of no command. In-process;
nothing here opens a database.
"""

from __future__ import annotations

import click
import pytest

import cli.main  # noqa: F401  (registers the commands)
from cli import filter_registry as fr
from cli.main import db


def _command(command_id: str) -> click.Command | None:
    return db.commands.get(command_id)


def _implemented() -> list[str]:
    return [c for c in fr.COMMANDS if _command(c) is not None]


def unaccounted_options(cmd: click.Command, command_id: str) -> list[str]:
    """Click options of `cmd` that are neither a registered filter of the
    command nor in NEVER_PRESET."""
    known = set(fr.by_param(command_id)) | set(fr.NEVER_PRESET)
    return [
        p.name
        for p in cmd.params
        if isinstance(p, click.Option) and p.name != "help" and p.name not in known
    ]


def test_trio_is_a_filter_command_and_there_is_something_to_check():
    assert "trio" in _implemented()


@pytest.mark.parametrize("command_id", _implemented())
def test_every_option_is_registered_or_never_presettable(command_id):
    assert unaccounted_options(_command(command_id), command_id) == []


def test_the_check_is_not_vacuous():
    """An option added to the command without a registry entry is caught."""
    cmd = click.Command(
        "trio", params=[*_command("trio").params, click.Option(["--min-qual"])]
    )
    assert unaccounted_options(cmd, "trio") == ["min_qual"]


@pytest.mark.parametrize("command_id", _implemented())
def test_click_defaults_are_the_registry_defaults(command_id):
    cmd = _command(command_id)
    options = {p.name: p for p in cmd.params if isinstance(p, click.Option)}
    for param, key in fr.by_param(command_id).items():
        if key.default is fr.DERIVED:
            continue
        assert options[param].default == key.default, key.file_key


def test_every_registered_value_key_is_a_parameter_of_a_command():
    implemented = _implemented()
    for key in fr.REGISTRY.values():
        if key.kind != fr.VALUE:
            continue
        homes = [c for c in key.commands if c in implemented]
        assert homes, f"{key.file_key} is registered for no command that exists"
        for c in homes:
            names = {p.name for p in _command(c).params}
            assert key.param in names, (key.file_key, c)


def test_the_flag_is_declared_by_the_option():
    options = {p.name: p for p in _command("trio").params}
    for key in fr.keys_for("trio"):
        declared = set(options[key.param].opts) | set(options[key.param].secondary_opts)
        assert set(key.flag.split("/")) == declared, key.file_key
