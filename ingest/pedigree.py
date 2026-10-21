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

All members share one ingest_id, so father_id/mother_id resolve to
sample_ids within the same ingest_id. A pedigree may hold several families
and several generations. Everything is checked before the table is touched
(an unknown sample, a self-parent, a cycle, a parent in another family), so
a refused load leaves the stored pedigree as it was.
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


def resulting_rows(
    rows: list[dict], existing: list[dict] | None, *, merge: bool
) -> list[dict]:
    """The pedigree a load will leave behind. Without `merge` it is `rows`
    alone (the whole pedigree of the ingest_id is replaced). With `merge` it
    is `rows` plus every stored row of a family that `rows` does not name.

    Raises ValueError when `merge` would move a sample out of the family it
    is stored in: --merge replaces whole families, so a sample that is
    named under a different family id is an error, not a silent move.
    """
    if not merge:
        return list(rows)
    if existing is None:
        raise ValueError("a merge needs the stored pedigree rows (existing=...)")
    incoming = {r["family_id"] for r in rows}
    named = {r["sample_id"]: r["family_id"] for r in rows}
    for e in existing:
        fam = named.get(e["sample_id"])
        if fam is not None and e["family_id"] not in incoming:
            raise ValueError(
                f"{e['sample_id']!r} is already in family {e['family_id']!r}; "
                f"it cannot also join {fam!r}. A sample belongs to one "
                "family; --merge replaces whole families."
            )
    return [e for e in existing if e["family_id"] not in incoming] + list(rows)


def load_pedigree(
    ingest_id: str,
    ped_path: str | Path,
    *,
    merge: bool = False,
    existing: list[dict] | None = None,
) -> int:
    """Load a PED file into the active DB's pedigree table under
    `ingest_id` (see `load_rows`). Returns the number of individuals
    loaded."""
    return load_rows(ingest_id, parse_ped(ped_path), merge=merge, existing=existing)


def load_rows(
    ingest_id: str,
    rows: list[dict],
    *,
    merge: bool = False,
    existing: list[dict] | None = None,
) -> int:
    """Write pedigree `rows` (the shape `parse_ped` returns) into the active
    DB under `ingest_id`. Returns the number of individuals loaded.

    Without `merge`, the prior pedigree of the ingest_id is replaced
    (idempotent re-load). With `merge`, only the families named in `rows`
    are replaced and every other family stays; `existing` must then hold
    the stored rows (cli/pedigree_load.pedigree_rows).

    Everything is validated BEFORE anything is deleted, so a refused load
    changes nothing: every individual and every named parent must be a
    real sample under this ingest_id, and the resulting pedigree must have
    no self-parent, no cycle and no parent in another family. A father
    recorded as female (or a mother as male) is only warned about.
    """
    from ingest._arrow import PEDIGREE_ARROW_SCHEMA
    from storage import (
        delete_where_sql,
        get_session,
        insert_via_parquet,
        sql_quote_str,
        validate_ingest_id,
    )
    from storage.pedigree_graph import Pedigree
    from storage.sql import str_list

    validate_ingest_id(ingest_id)
    if not rows:
        raise ValueError("no pedigree rows to load")

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

    result = Pedigree.from_rows(resulting_rows(rows, existing, merge=merge))
    errors = result.errors()
    if errors:
        raise ValueError(
            "the pedigree has errors, so nothing was changed: "
            + "; ".join(q.message for q in errors)
        )
    families = {r["family_id"] for r in rows}
    for q in result.validate():
        if q.code == "parent_sex" and q.family_id in families:
            log.warning("warning: %s", q.message)

    # Replace: clear the prior pedigree (or only the families being
    # written) first.
    where = f"ingest_id = {sql_quote_str(ingest_id)}"
    if merge:
        where += f" AND family_id IN ({str_list(sorted(families))})"
    sess.query(delete_where_sql("pedigree", where))

    insert_via_parquet(
        "pedigree",
        PEDIGREE_ARROW_SCHEMA,
        [{"ingest_id": ingest_id, **r} for r in rows],
    )
    log.info("[ped] loaded %d individuals under ingest_id=%s", len(rows), ingest_id)
    return len(rows)
