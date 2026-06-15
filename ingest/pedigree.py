"""Load family relationships from a PED/FAM file.

PED is the standard plink/linkage pedigree format: whitespace-delimited,
six columns, one row per individual:

    family_id  individual_id  paternal_id  maternal_id  sex  phenotype

  * paternal_id / maternal_id are '0' for a founder (no parent in file).
  * sex: 1 = male, 2 = female, other = unknown.
  * phenotype (affected): 1 = unaffected, 2 = affected, 0/-9 = unknown.

Lines beginning with '#' are skipped. We normalise sex and affected to
readable strings ('male'/'female', 'affected'/'unaffected') so SQL
filters read naturally, matching vcfclick's filter='PASS' style.

v1 assumes a JOINT-called trio: all members share one ingest_id, so
father_id/mother_id resolve to sample_ids within the same ingest_id.
"""

from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)

_SEX = {"1": "male", "2": "female"}
_AFFECTED = {"1": "unaffected", "2": "affected"}


def parse_ped(path: str | Path) -> list[dict]:
    """Parse a PED file into rows ready for the pedigree table.

    Each row is a dict with keys: sample_id, family_id, father_id,
    mother_id, sex, affected. sex/affected are normalised strings or
    None. ingest_id is NOT set here — the caller fills it (the PED's
    sample_ids are resolved against an ingested cohort).
    """
    rows: list[dict] = []
    seen: set[str] = set()
    for lineno, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(
                f"{path}:{lineno}: PED needs 6 whitespace-delimited columns "
                f"(family individual paternal maternal sex phenotype), got "
                f"{len(fields)}: {line!r}"
            )
        family_id, sample_id, father_id, mother_id, sex_code, pheno_code = fields[:6]
        if sample_id in seen:
            raise ValueError(f"{path}:{lineno}: duplicate individual id {sample_id!r}")
        seen.add(sample_id)
        rows.append(
            {
                "sample_id": sample_id,
                "family_id": family_id,
                "father_id": father_id,
                "mother_id": mother_id,
                "sex": _SEX.get(sex_code),
                "affected": _AFFECTED.get(pheno_code),
            }
        )
    if not rows:
        raise ValueError(f"{path}: no pedigree rows found")
    return rows


def load_pedigree(ingest_id: str, ped_path: str | Path) -> int:
    """Load a PED file into the active DB's pedigree table under
    `ingest_id`. Replaces any prior pedigree rows for that ingest_id
    (idempotent re-load). Returns the number of individuals loaded.

    Validates that every individual and every named parent is a real
    sample in the cohort under this ingest_id — a PED referencing
    sample_ids that aren't in the data is almost always a mistake.
    """
    from ingest._arrow import PEDIGREE_ARROW_SCHEMA
    from storage import (
        delete_where_sql,
        get_session,
        insert_via_parquet,
        sql_quote_str,
        validate_ingest_id,
    )

    validate_ingest_id(ingest_id)
    rows = parse_ped(ped_path)

    sess = get_session()
    # Which sample_ids actually exist in this cohort?
    raw = (
        sess.query(
            f"SELECT DISTINCT sample_id FROM samples WHERE ingest_id = "
            f"{sql_quote_str(ingest_id)} FORMAT TabSeparated"
        )
        .bytes()
        .decode()
    )
    known = {s for s in raw.splitlines() if s.strip()}
    if not known:
        raise ValueError(
            f"no samples found under ingest_id={ingest_id!r}; ingest the "
            f"VCF before loading its pedigree."
        )

    for r in rows:
        if r["sample_id"] not in known:
            raise ValueError(
                f"PED individual {r['sample_id']!r} is not a sample under "
                f"ingest_id={ingest_id!r}. Known samples: {sorted(known)}"
            )
        for parent_key in ("father_id", "mother_id"):
            pid = r[parent_key]
            if pid not in ("0", "") and pid not in known:
                raise ValueError(
                    f"PED lists {pid!r} as {parent_key} of {r['sample_id']!r} "
                    f"but it is not a sample under ingest_id={ingest_id!r}."
                )

    # Idempotent replace: clear prior pedigree for this ingest_id first.
    sess.query(delete_where_sql("pedigree", f"ingest_id = '{ingest_id}'"))

    insert_via_parquet(
        "pedigree",
        PEDIGREE_ARROW_SCHEMA,
        [{"ingest_id": ingest_id, **r} for r in rows],
    )
    log.info("[ped] loaded %d individuals under ingest_id=%s", len(rows), ingest_id)
    return len(rows)
