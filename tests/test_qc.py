"""Tests for `vcfclick db qc` — per-sample QC metrics + chrX sex check.

The chrX fixture (tests/fixtures/qc_sex.vcf.gz) has three samples: M1 with
a male hemizygous pattern (1 het / 24 hom-alt), F1 and SWAP with a female
pattern (~half het). The pedigree declares SWAP male, so its inferred
female sex is a flagged mismatch. ALTs alternate transition/transversion
so Ti/Tv is a real ratio.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCFCLICK_BIN = shutil.which("vcfclick") or str(REPO / ".venv" / "bin" / "vcfclick")
FIX = Path(__file__).parent / "fixtures"
QC_VCF = FIX / "qc_sex.vcf.gz"
QC_PED = FIX / "qc_sex.ped"


def _vc(home: Path, *args: str):
    env = os.environ.copy()
    env["VCFCLICK_HOME"] = str(home)
    env.pop("VCFCLICK_DB_NAME", None)
    r = subprocess.run(
        [VCFCLICK_BIN, *args], cwd=REPO, env=env, capture_output=True, text=True
    )
    assert r.returncode == 0, f"{' '.join(args)} failed:\n{r.stderr}"
    return r


def _setup(home: Path, *, with_ped: bool = True) -> None:
    _vc(home, "db", "create", "q")
    _vc(
        home,
        "db",
        "ingest",
        "q",
        str(QC_VCF),
        "--cohort",
        "c",
        "--ingest-id",
        "i1",
        "--serial",
    )
    if with_ped:
        _vc(home, "db", "ped", "q", str(QC_PED))


def _qc(home: Path) -> dict:
    out = _vc(home, "db", "qc", "q", "--format", "json").stdout
    return {s["sample_id"]: s for s in json.loads(out)}


def test_qc_het_hom_and_titv(vcfclick_home):
    """het / hom-alt counts and the Ti/Tv ratio come straight from the
    stored non-reference calls (13 transitions, 12 transversions)."""
    _setup(vcfclick_home)
    qc = _qc(vcfclick_home)
    assert qc["M1"]["het"] == 1 and qc["M1"]["hom_alt"] == 24
    assert qc["F1"]["het"] == 12 and qc["F1"]["hom_alt"] == 13
    assert qc["F1"]["ti_tv"] == "1.08"


def test_qc_sex_inference(vcfclick_home):
    """chrX het fraction infers sex: M1 hemizygous → male, F1 ~half het →
    female."""
    _setup(vcfclick_home)
    qc = _qc(vcfclick_home)
    assert qc["M1"]["inferred_sex"] == "male"
    assert qc["F1"]["inferred_sex"] == "female"


def test_qc_flags_pedigree_sex_mismatch(vcfclick_home):
    """SWAP has a female chrX pattern but the pedigree declares it male →
    sex_mismatch, the signal for a sample swap. M1/F1 agree, no flag."""
    _setup(vcfclick_home)
    qc = _qc(vcfclick_home)
    assert qc["SWAP"]["inferred_sex"] == "female"
    assert qc["SWAP"]["pedigree_sex"] == "male"
    assert qc["SWAP"]["sex_mismatch"] is True
    assert qc["M1"]["sex_mismatch"] is False
    assert qc["F1"]["sex_mismatch"] is False


def test_qc_without_pedigree_has_no_declared_sex(vcfclick_home):
    """No pedigree loaded → inference still runs, but nothing to flag."""
    _setup(vcfclick_home, with_ped=False)
    qc = _qc(vcfclick_home)
    assert qc["M1"]["inferred_sex"] == "male"
    assert qc["M1"]["pedigree_sex"] is None
    assert all(not s["sex_mismatch"] for s in qc.values())
