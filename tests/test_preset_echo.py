"""The preset echo as text: the `preset:` and `gates:` lines and the header.

Pure functions of `cli/preset_echo.py`; the words a clinician reads under the
`trio:` line. In-process, no database.
"""

from __future__ import annotations

from cli import filter_registry as fr
from cli.preset_echo import gates_line, header_lines, preset_line
from cli.warning_codes import warning


def _resolved(**sources):
    out = {}
    for key in fr.keys_for("trio"):
        out[key.file_key] = fr.Resolved(
            None if key.default is fr.DERIVED else fr.normalise_value(key, key.default),
            "default",
        )
    for key, (value, source) in sources.items():
        out[key] = fr.Resolved(value, source)
    return out


def test_the_preset_line():
    base = {"name": "rare-recessive", "source": "builtin", "revision": 1,
            "fingerprint": "a1b2c3d4e5f6", "status": "starting_point", "pinned": False}  # fmt: skip
    assert preset_line(base) == (
        "preset: rare-recessive (builtin, rev 1, a1b2c3d4e5f6, starting point)"
    )
    user = {**base, "source": "user", "revision": None, "status": None, "pinned": True}
    assert preset_line(user) == "preset: rare-recessive (user, a1b2c3d4e5f6, pinned)"


def test_the_gates_line_shows_every_value_with_its_source():
    resolved = _resolved(
        min_gq=(30, "preset"), max_af=(0.005, "flag"), strict_missing=(True, "preset"),
        gnomad_max_af=(None, "default"),
    )  # fmt: skip
    over = [
        {"key": "max_af", "preset": 0.01, "flag": 0.005},
        {"key": "gnomad_max_af", "preset": 0.0001, "flag": "off"},
    ]
    line = gates_line(resolved, over, "applied")
    assert line.startswith("  gates: ")
    assert "min_gq=30 [preset]" in line
    assert "min_dp=10 [default]" in line
    assert "max_af=0.005 [flag; preset 0.01]" in line
    assert "gnomad_max_af=none [off; preset 0.0001]" in line
    assert "strict_missing=on [preset]" in line
    assert "pass_only" not in line and "max_parent_alt_frac" not in line  # off defaults
    skipped = gates_line(_resolved(), [], "skipped_cohort_derived")
    assert "max_af=0.01 [default; skipped: INFO/AF is cohort-derived]" in skipped


def test_the_header_prints_warnings_not_notes():
    block = {"name": "p", "source": "user", "revision": None, "fingerprint": "f" * 12,
             "status": None, "pinned": False, "overridden": []}  # fmt: skip
    lines = header_lines(
        block, _resolved(), "applied",
        [warning("gate_unbound", "min_gq", "a"), warning("gate_partial", "min_dp", "b")],
    )  # fmt: skip
    assert lines[0].startswith("preset: p (user")
    assert lines[2:] == ["  WARNING: a"]


def test_the_gates_line_shows_the_parent_fraction_with_the_screen_and_survives_a_command_without_it():
    on = _resolved(parent_alt_screen=(True, "flag"))
    assert "max_parent_alt_frac=0.05 [default]" in gates_line(on, [], "applied")
    # a command that has no parent_alt_screen key (it is trio-only) still prints
    del on["parent_alt_screen"]
    assert "max_parent_alt_frac" not in gates_line(on, [], "applied")
    on["max_parent_alt_frac"] = fr.Resolved(0.1, "flag")
    assert "max_parent_alt_frac=0.1 [flag]" in gates_line(on, [], "applied")
