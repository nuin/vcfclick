from __future__ import annotations

import pytest

from tui.services import LocusInputError, ParsedLocus, parse_locus_input


def test_parse_range_with_commas():
    locus = parse_locus_input("chr17:43,044,295-43,125,483")
    assert locus == ParsedLocus(
        kind="range",
        label="chr17:43044295-43125483",
        chrom="chr17",
        start_pos=43044295,
        end_pos=43125483,
        gene_symbol=None,
    )


def test_parse_gene_symbol():
    locus = parse_locus_input("BRCA1")
    assert locus == ParsedLocus(
        kind="gene",
        label="BRCA1",
        chrom=None,
        start_pos=None,
        end_pos=None,
        gene_symbol="BRCA1",
    )


@pytest.mark.parametrize(
    "raw",
    ["", "chr17", "chr17:431-430", "chr17:start-end", "BRCA1 BRCA2"],
)
def test_parse_invalid_locus(raw):
    with pytest.raises(LocusInputError):
        parse_locus_input(raw)
