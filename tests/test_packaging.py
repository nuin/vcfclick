"""Checks that a build from the sdist ships every package the wheel does.

`uv build` builds the sdist first and the wheel from that sdist, so a
package missing from the sdist include list is missing from the
published wheel too. Bioconda builds from the sdist as well.
"""

from __future__ import annotations

import shutil
import subprocess
import tarfile
import tomllib
import zipfile
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parent.parent

PACKAGE_DIRS = {
    "annotations",
    "benchmark",
    "cli",
    "export",
    "ingest",
    "popgen",
    "schema",
    "storage",
    "tui",
    "vcfclick_mcp",
    "vcfclick_web",
}


def _hatch_targets() -> dict:
    data = tomllib.loads((REPO / "pyproject.toml").read_text())
    return data["tool"]["hatch"]["build"]["targets"]


def test_wheel_packages_are_the_top_level_packages():
    """Every top-level package directory is in the wheel package list."""
    assert set(_hatch_targets()["wheel"]["packages"]) == PACKAGE_DIRS


def test_sdist_includes_every_wheel_package():
    """Each wheel package appears in the sdist include list, root-anchored."""
    targets = _hatch_targets()
    include = set(targets["sdist"]["include"])
    missing = [pkg for pkg in targets["wheel"]["packages"] if f"/{pkg}" not in include]
    assert missing == []


def test_sdist_include_patterns_are_anchored():
    """A bare `schema` pattern matches at any depth (a nested worktree copy
    under .claude/worktrees/ would ship); `/schema` matches the root only."""
    include = _hatch_targets()["sdist"]["include"]
    unanchored = [pattern for pattern in include if not pattern.startswith("/")]
    assert unanchored == []


def test_gitignore_excludes_agent_working_dirs():
    lines = {line.strip() for line in (REPO / ".gitignore").read_text().splitlines()}
    assert ".claude/" in lines
    assert ".cairn" in lines or ".cairn/" in lines


def _export_tracked_tree(dest: Path) -> None:
    """Copy the tracked files (working-tree content) into dest."""
    listing = subprocess.run(
        ["git", "ls-files", "-z"], cwd=REPO, check=True, capture_output=True
    ).stdout
    for rel in listing.decode().split("\0"):
        if not rel:
            continue
        src = REPO / rel
        if not src.is_file():
            continue
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, target)


@pytest.mark.skipif(
    shutil.which("uv") is None or shutil.which("git") is None,
    reason="needs uv and git on PATH",
)
def test_build_from_export_ships_web_and_mcp_and_no_nested_copies(tmp_path):
    """Build sdist + wheel from an exported tree with a planted nested
    worktree. The wheel must carry vcfclick_web and vcfclick_mcp; the
    sdist must not carry the nested copy."""
    export = tmp_path / "export"
    _export_tracked_tree(export)
    planted = export / ".claude" / "worktrees" / "agent-x"
    (planted / "schema").mkdir(parents=True)
    (planted / "vcfclick_mcp").mkdir(parents=True)
    (planted / "schema" / "x.sql").write_text("select 1\n")
    (planted / "vcfclick_mcp" / "p.py").write_text("x = 1\n")

    subprocess.run(["uv", "build", "-q"], cwd=export, check=True, timeout=600)

    (sdist,) = (export / "dist").glob("*.tar.gz")
    with tarfile.open(sdist) as tar:
        sdist_names = tar.getnames()
    (wheel,) = (export / "dist").glob("*.whl")
    with zipfile.ZipFile(wheel) as zf:
        wheel_names = zf.namelist()

    assert "vcfclick_web/app.py" in wheel_names
    assert "vcfclick_mcp/server.py" in wheel_names
    assert "cli/main.py" in wheel_names
    assert any(name.endswith("/vcfclick_web/app.py") for name in sdist_names)
    assert [name for name in sdist_names if ".claude" in name] == []
