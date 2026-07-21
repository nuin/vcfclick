"""Spec-derived hap.py / GA4GH intermediate-VCF conformance checks.

Parses the authoritative vocabulary carried by
tests/fixtures/benchmark/happy_intermediate_header.vcf (the ``##vcfclick_spec_*``
lines, compiled from hap.py @master source — SPEC-DERIVED, not from a real
hap.py run) and asserts that ``benchmark/constants.py`` string values are all
present in / consistent with that vocabulary, casing included.

A constant that disagrees with the spec is a real finding: the corresponding
test FAILS on purpose (constants.py is left unedited).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from benchmark.constants import (
    BD_FN,
    BD_FP,
    BD_N,
    BD_TP,
    BK_AM,
    BK_GM,
    BK_LM,
    BK_NONE,
    BLT_HAPLOID,
    BLT_HET,
    BLT_HETALT,
    BLT_HOMALT,
    BLT_HOMREF,
    BLT_NOCALL,
    VT_COMPLEX,
    VT_INDEL,
    VT_NOCALL,
    VT_SNP,
)

_FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "benchmark"
    / "happy_intermediate_header.vcf"
)

_PREFIX = "##vcfclick_spec_"


def _spec_vocab() -> dict[str, set[str]]:
    """Parse ``##vcfclick_spec_<FIELD>=v1,v2,...`` lines into value sets."""
    vocab: dict[str, set[str]] = {}
    for line in _FIXTURE.read_text().splitlines():
        if not line.startswith(_PREFIX):
            continue
        field, _, values = line[len(_PREFIX) :].partition("=")
        vocab[field] = {v for v in values.split(",") if v}
    return vocab


VOCAB = _spec_vocab()


def test_fixture_carries_all_expected_fields() -> None:
    assert set(VOCAB) == {"BD", "BK", "BVT", "BLT", "BI"}


def test_header_reproduces_happy_format_lines() -> None:
    """The verbatim hap.py FORMAT/INFO definitions must be present."""
    text = _FIXTURE.read_text()
    for needle in (
        '##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision for call (TP/FP/FN/N)">',
        '##FORMAT=<ID=BK,Number=1,Type=String,Description="Sub-type for decision (match/mismatch type)">',
        '##FORMAT=<ID=BI,Number=1,Type=String,Description="Additional comparison information">',
        '##FORMAT=<ID=BVT,Number=1,Type=String,Description="High-level variant type (SNP|INDEL).">',
        "##INFO=<ID=BS,Number=.,Type=Integer,",
    ):
        assert needle in text, needle


@pytest.mark.parametrize(
    "const_name, value",
    [
        ("BD_TP", BD_TP),
        ("BD_FP", BD_FP),
        ("BD_FN", BD_FN),
        ("BD_N", BD_N),
    ],
)
def test_bd_values_in_spec(const_name: str, value: str) -> None:
    assert value in VOCAB["BD"], f"{const_name}={value!r} not in spec BD {VOCAB['BD']}"


@pytest.mark.parametrize(
    "const_name, value",
    [
        ("BK_GM", BK_GM),
        ("BK_AM", BK_AM),
        ("BK_LM", BK_LM),
        ("BK_NONE", BK_NONE),
    ],
)
def test_bk_values_in_spec(const_name: str, value: str) -> None:
    assert value in VOCAB["BK"], f"{const_name}={value!r} not in spec BK {VOCAB['BK']}"


@pytest.mark.parametrize(
    "const_name, value",
    [
        ("VT_SNP", VT_SNP),
        ("VT_INDEL", VT_INDEL),
        ("VT_NOCALL", VT_NOCALL),
        ("VT_COMPLEX", VT_COMPLEX),
    ],
)
def test_bvt_values_in_spec(const_name: str, value: str) -> None:
    assert value in VOCAB["BVT"], (
        f"{const_name}={value!r} not in spec BVT {sorted(VOCAB['BVT'])}"
    )


@pytest.mark.parametrize(
    "const_name, value",
    [
        ("BLT_HET", BLT_HET),
        ("BLT_HOMALT", BLT_HOMALT),
        ("BLT_HETALT", BLT_HETALT),
        ("BLT_HOMREF", BLT_HOMREF),
        ("BLT_NOCALL", BLT_NOCALL),
        ("BLT_HAPLOID", BLT_HAPLOID),
    ],
)
def test_blt_values_in_spec(const_name: str, value: str) -> None:
    assert value in VOCAB["BLT"], (
        f"{const_name}={value!r} not in spec BLT {sorted(VOCAB['BLT'])}"
    )
