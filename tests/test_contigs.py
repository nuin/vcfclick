"""storage.contigs: the single definition of every contig alias set."""

from __future__ import annotations

import pytest

from storage import contigs

# (name, x, y, par, mito)  -- one row per alias the program recognises.
ALIAS_TABLE = [
    ("chrX", True, False, False, False),
    ("X", True, False, False, False),
    ("chr23", True, False, False, False),
    ("23", True, False, False, False),
    ("chrY", False, True, False, False),
    ("Y", False, True, False, False),
    ("chr24", False, True, False, False),
    ("24", False, True, False, False),
    ("chrXY", False, False, True, False),
    ("XY", False, False, True, False),
    ("chr25", False, False, True, False),
    ("25", False, False, True, False),
    ("PAR", False, False, True, False),
    ("chrPAR", False, False, True, False),
    ("chrM", False, False, False, True),
    ("M", False, False, False, True),
    ("MT", False, False, False, True),
    ("chrMT", False, False, False, True),
    ("chr26", False, False, False, True),
    ("26", False, False, False, True),
]


@pytest.mark.parametrize("name,x,y,par,mito", ALIAS_TABLE)
def test_alias_sets(name, x, y, par, mito):
    assert contigs.is_x(name) is x
    assert contigs.is_y(name) is y
    assert contigs.is_par_label(name) is par
    assert contigs.is_mito(name) is mito
    # every alias is non-autosomal, and exactly one class claims it
    assert contigs.is_sex_or_mt(name)
    assert not contigs.is_autosome(name)
    assert sum([x, y, par, mito]) == 1
    assert name in contigs.NON_AUTOSOMAL


def test_table_covers_every_declared_alias():
    """A new alias added to a set must be added to ALIAS_TABLE too."""
    declared = {
        *contigs.X_NAMES,
        *contigs.Y_NAMES,
        *contigs.PAR_NAMES,
        *contigs.MITO_NAMES,
    }
    assert declared == {row[0] for row in ALIAS_TABLE}
    assert contigs.NON_AUTOSOMAL == declared
    # the four sets are disjoint
    assert len(declared) == len(contigs.X_NAMES) + len(contigs.Y_NAMES) + len(
        contigs.PAR_NAMES
    ) + len(contigs.MITO_NAMES)


@pytest.mark.parametrize(
    "name", [str(n) for n in range(1, 23)] + [f"chr{n}" for n in range(1, 23)]
)
def test_autosomes(name):
    assert contigs.is_autosome(name)
    assert not contigs.is_sex_or_mt(name)
    assert name in contigs.autosome_names()


@pytest.mark.parametrize(
    "name",
    [
        "0", "chr0", "01", "chr01", "27", "chr27", "100",
        "CHR1", "Chr1", "chrx", "x", "mt",
        "chrUn_KI270302v1", "chr1_KI270706v1_random", "chrX_KI270880v1_alt",
        "chrEBV", "hs37d5", "GL000192.1", "HLA-A*01:01:01:01", "NC_000001.11",
        "", "chr",
    ],
)  # fmt: skip
def test_not_autosomes_and_not_aliases(name):
    """Unplaced/alt/decoy contigs and mis-cased names match nothing: they
    are not autosomes, and they are not X, Y, PAR or mitochondrial either."""
    assert not contigs.is_autosome(name)
    assert not contigs.is_sex_or_mt(name)


def test_autosome_names_for_sql():
    names = contigs.autosome_names()
    assert len(names) == 44 and len(set(names)) == 44
    assert names[:2] == ("1", "2") and names[22] == "chr1" and names[-1] == "chr22"
    assert all(contigs.is_autosome(n) for n in names)
