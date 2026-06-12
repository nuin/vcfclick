from __future__ import annotations

import sys

from click.testing import CliRunner

from cli.main import cli


def test_tui_command_is_registered():
    runner = CliRunner()
    result = runner.invoke(cli, ["tui", "--help"])
    assert result.exit_code == 0
    assert "Launch the optional terminal UI" in result.output
    assert "--db" in result.output


def test_tui_command_reports_missing_textual(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "textual.app":
            raise ModuleNotFoundError("No module named 'textual'", name="textual")
        return real_import(name, *args, **kwargs)

    sys.modules.pop("tui.app", None)
    monkeypatch.setattr(builtins, "__import__", fake_import)

    runner = CliRunner()
    result = runner.invoke(cli, ["tui"])

    assert result.exit_code != 0
    assert "Install the TUI extra" in result.output


def test_tui_command_does_not_mask_missing_internal_module(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "tui.app":
            raise ModuleNotFoundError("No module named 'tui.app'", name="tui.app")
        return real_import(name, *args, **kwargs)

    sys.modules.pop("tui.app", None)
    monkeypatch.setattr(builtins, "__import__", fake_import)

    runner = CliRunner()
    result = runner.invoke(cli, ["tui"])

    assert isinstance(result.exception, ModuleNotFoundError)
    assert result.exception.name == "tui.app"
    assert "Install the TUI extra" not in result.output
