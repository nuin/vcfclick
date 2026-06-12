from __future__ import annotations

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
        if name == "tui.app":
            raise ModuleNotFoundError("No module named 'textual'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    runner = CliRunner()
    result = runner.invoke(cli, ["tui"])

    assert result.exit_code != 0
    assert "Install the TUI extra" in result.output
