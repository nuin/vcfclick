"""The legacy trio SQL, pinned as text.

The four classic models (denovo, recessive, dominant, comphet) must keep
producing exactly the SQL they produced before the family program: the
refactors that route them through shared helpers (`sample_join`, the
registry) are only allowed to change how the text is assembled, never the
text. The AF epsilon is the one deliberate change, made in its own commit
(it edits the `info_AF <=` line of tests/goldens/trio_sql/*.sql and nothing
else).

Regenerate (and say why in the commit):

    VCFCLICK_UPDATE_GOLDENS=1 .venv/bin/python -m pytest \\
        tests/test_trio_sql_text.py -p no:cacheprovider
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

GOLDEN_DIR = Path(__file__).parent / "goldens" / "trio_sql"
UPDATE = os.environ.get("VCFCLICK_UPDATE_GOLDENS") == "1"


def _statements() -> dict[str, str]:
    from cli.trio_analysis import Gates, Trio, _comphet_sql, trio_sql

    trio = Trio("i1", "C1", "F", "M")
    gates = Gates(20, 10, 0.01, 0.25, 0.75)
    out = {
        cat: trio_sql(cat, trio, gates, count_only=False)
        for cat in ("denovo", "recessive", "dominant")
    }
    out["denovo_count"] = trio_sql("denovo", trio, gates, count_only=True)
    out["comphet_candidates"] = _comphet_sql(trio, gates)
    return out


def _render(statements: dict[str, str]) -> str:
    return "".join(f"-- {name}\n{sql}\n\n" for name, sql in statements.items())


@pytest.mark.parametrize("backend", ["chdb", "duckdb"])
def test_legacy_trio_sql_text_is_unchanged(backend, monkeypatch):
    monkeypatch.setenv("VCFCLICK_BACKEND", backend)
    got = _render(_statements())
    path = GOLDEN_DIR / f"{backend}.sql"
    if UPDATE:
        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        path.write_text(got)
        return
    assert got == path.read_text()
