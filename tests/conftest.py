"""Shared pytest fixtures.

Every test gets its own VCFCLICK_HOME under tmp_path so nothing touches
the user's real ~/.vcfclick/dbs/. The fixture VCF is a tiny committed
file under tests/fixtures/ — 5 variants, 3 samples, bgzip+tabix indexed.
"""

from __future__ import annotations

from pathlib import Path

import pytest


FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def vcfclick_home(tmp_path, monkeypatch) -> Path:
    """Isolated VCFCLICK_HOME for one test.

    Tests invoke the CLI via subprocess (see test_cli._vc), so they pick
    this env var up at process startup — no module reload needed.
    """
    home = tmp_path / "vcfclick"
    home.mkdir()
    monkeypatch.setenv("VCFCLICK_HOME", str(home))
    monkeypatch.delenv("VCFCLICK_DB_NAME", raising=False)
    return home


@pytest.fixture
def tiny_vcf() -> Path:
    """Path to the committed 5-variant / 3-sample fixture VCF."""
    p = FIXTURES / "tiny.vcf.gz"
    assert p.exists(), f"missing fixture: {p}"
    assert (FIXTURES / "tiny.vcf.gz.tbi").exists(), "missing tabix index"
    return p
