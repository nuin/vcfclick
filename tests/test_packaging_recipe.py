"""Checks for the bioconda recipe smoke commands."""

from __future__ import annotations

import re
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_NAME_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9_-]{0,62}$")


def test_bioconda_smoke_db_name_is_valid():
    """The recipe's smoke DB name must pass storage.db_path validation."""
    recipe = REPO / "packaging" / "bioconda" / "meta.yaml"
    text = recipe.read_text()
    match = re.search(r"vcfclick db create (?P<name>\S+)", text)

    assert match is not None
    assert DB_NAME_RE.fullmatch(match.group("name"))
