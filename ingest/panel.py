"""Load a sample -> population panel into the `populations` table.

Accepted out of the box: the 1000 Genomes panel format, a tab-separated
file with a header

    sample  pop  super_pop  gender

(the published `integrated_call_samples_*.panel` files carry trailing
empty tab fields, which are ignored). Any other delimited file works too:
tab, comma (`.csv` or a comma-only header) or whitespace separated, with
a header row. Columns are found by name — `sample`/`sample_id`/`IID`,
`pop`/`population`, `super_pop`/`super_population`, `gender`/`sex` —
or named explicitly by the caller.

`gender`/`sex` values are normalised to 'male' / 'female' / NULL
(male, m, 1 and female, f, 2, case-insensitive; anything else is NULL),
matching the pedigree table's convention.

A panel is keyed by sample_id only, but vcfclick identifies a sample by
(ingest_id, sample_id). Without an explicit ingest_id the panel applies
to every ingestion that contains the sample. Re-loading replaces the
rows for the affected (ingest_id, sample_id) pairs, so it is idempotent.
"""

from __future__ import annotations

import csv
import io
import logging
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

_SAMPLE_NAMES = ("sample", "sample_id", "sampleid", "iid", "id", "sample_name")
_POP_NAMES = ("pop", "population")
_SUPER_NAMES = ("super_pop", "super_population", "superpopulation", "superpop")
_SEX_NAMES = ("gender", "sex")

_MALE = {"male", "m", "1"}
_FEMALE = {"female", "f", "2"}

# Deletes are issued per chunk of sample ids so the IN list stays small.
_DELETE_CHUNK = 500


class PanelError(ValueError):
    """A malformed panel file or an unusable column choice."""


def normalise_sex(value: str | None) -> str | None:
    """'male' / 'female' / None from the usual panel and PED spellings."""
    if value is None:
        return None
    v = value.strip().lower()
    if v in _MALE:
        return "male"
    if v in _FEMALE:
        return "female"
    return None


def _sniff_delimiter(path: Path, header: str) -> str | None:
    """Tab if the header has one, comma for CSV, else None (whitespace)."""
    if "\t" in header:
        return "\t"
    if path.suffix.lower() == ".csv" or "," in header:
        return ","
    return None


def _split(line: str, delimiter: str | None) -> list[str]:
    if delimiter is None:
        return line.split()
    return next(csv.reader(io.StringIO(line), delimiter=delimiter))


def _find_column(
    header: list[str], explicit: str | None, candidates: tuple[str, ...], what: str
) -> int | None:
    lowered = [h.strip().lower() for h in header]
    if explicit is not None:
        try:
            return lowered.index(explicit.strip().lower())
        except ValueError:
            raise PanelError(
                f"{what} column {explicit!r} is not in the panel header {header}"
            ) from None
    for name in candidates:
        if name in lowered:
            return lowered.index(name)
    return None


@dataclass(frozen=True)
class PanelColumns:
    """Explicit column names; None means auto-detect by header name."""

    sample: str | None = None
    population: str | None = None
    super_population: str | None = None
    sex: str | None = None


def parse_panel(path: str | Path, columns: PanelColumns = PanelColumns()) -> list[dict]:
    """Parse a panel file into rows with keys sample_id, population,
    super_population, sex. Raises PanelError on a malformed file."""
    path = Path(path)
    lines = [
        ln
        for ln in path.read_text().splitlines()
        if ln.strip() and not ln.lstrip().startswith("##")
    ]
    if not lines:
        raise PanelError(f"{path}: empty panel file")
    header_line = lines[0].lstrip("#")
    delimiter = _sniff_delimiter(path, header_line)
    header = [h.strip() for h in _split(header_line, delimiter)]
    while header and not header[-1]:
        header.pop()

    i_sample = _find_column(header, columns.sample, _SAMPLE_NAMES, "sample")
    i_pop = _find_column(header, columns.population, _POP_NAMES, "population")
    if i_sample is None or i_pop is None:
        raise PanelError(
            f"{path}: cannot find the sample and population columns in header "
            f"{header}; name them with --sample-col / --pop-col"
        )
    i_super = _find_column(
        header, columns.super_population, _SUPER_NAMES, "super-population"
    )
    i_sex = _find_column(header, columns.sex, _SEX_NAMES, "sex")

    rows: dict[str, dict] = {}
    for lineno, line in enumerate(lines[1:], start=2):
        fields = [f.strip() for f in _split(line, delimiter)]

        def get(i: int | None, fields: list[str] = fields) -> str | None:
            if i is None or i >= len(fields):
                return None
            return fields[i] or None

        sample, pop = get(i_sample), get(i_pop)
        if not sample:
            raise PanelError(f"{path}:{lineno}: empty sample id: {line!r}")
        if not pop:
            raise PanelError(
                f"{path}:{lineno}: sample {sample!r} has no population: {line!r}"
            )
        row = {
            "sample_id": sample,
            "population": pop,
            "super_population": get(i_super),
            "sex": normalise_sex(get(i_sex)),
        }
        prior = rows.get(sample)
        if prior is not None and prior != row:
            raise PanelError(
                f"{path}:{lineno}: sample {sample!r} listed twice with "
                f"different labels ({prior} vs {row})"
            )
        rows[sample] = row
    if not rows:
        raise PanelError(f"{path}: no panel rows found")
    return list(rows.values())


@dataclass
class PanelReport:
    """What a panel load did, for the CLI to print."""

    loaded: int = 0
    by_ingest: dict[str, int] = field(default_factory=dict)
    panel_not_in_db: list[str] = field(default_factory=list)
    db_not_in_panel: list[tuple[str, str]] = field(default_factory=list)


def _db_samples(sess, ingest_id: str | None) -> list[tuple[str, str]]:
    import json

    from storage import sql_quote_str

    where = f" WHERE ingest_id = {sql_quote_str(ingest_id)}" if ingest_id else ""
    raw = (
        sess.query(
            "SELECT DISTINCT ingest_id, sample_id FROM samples"
            f"{where} ORDER BY ingest_id, sample_id",
            "JSONCompact",
        )
        .bytes()
        .decode()
    )
    return [(r[0], r[1]) for r in json.loads(raw)["data"]]


def load_panel(rows: list[dict], ingest_id: str | None = None) -> PanelReport:
    """Write `rows` (from parse_panel) to the active DB's populations table.

    Applies to `ingest_id` only when given, else to every ingestion that
    contains each sample. Replaces prior rows for the affected
    (ingest_id, sample_id) pairs.
    """
    from ingest._arrow import POPULATIONS_ARROW_SCHEMA
    from storage import (
        delete_where_sql,
        get_session,
        insert_via_parquet,
        sql_quote_str,
        upgrade_schema,
        validate_ingest_id,
    )

    if ingest_id is not None:
        validate_ingest_id(ingest_id)
    upgrade_schema()  # an older database may not have `populations` yet
    sess = get_session()

    in_db = _db_samples(sess, ingest_id)
    if not in_db:
        scope = f"ingest_id={ingest_id!r}" if ingest_id else "this database"
        raise PanelError(f"no samples found under {scope}; ingest a VCF first")

    panel = {r["sample_id"]: r for r in rows}
    matched = [(ing, s) for ing, s in in_db if s in panel]
    report = PanelReport(
        panel_not_in_db=sorted(set(panel) - {s for _, s in in_db}),
        db_not_in_panel=[(ing, s) for ing, s in in_db if s not in panel],
    )
    if not matched:
        return report

    by_ingest: dict[str, list[str]] = {}
    for ing, s in matched:
        by_ingest.setdefault(ing, []).append(s)
    for ing, samples in by_ingest.items():
        for start in range(0, len(samples), _DELETE_CHUNK):
            chunk = samples[start : start + _DELETE_CHUNK]
            in_list = ", ".join(sql_quote_str(s) for s in chunk)
            sess.query(
                delete_where_sql(
                    "populations",
                    f"ingest_id = {sql_quote_str(ing)} AND sample_id IN ({in_list})",
                )
            )
    insert_via_parquet(
        "populations",
        POPULATIONS_ARROW_SCHEMA,
        [{"ingest_id": ing, **panel[s]} for ing, s in matched],
    )
    report.loaded = len(matched)
    report.by_ingest = {ing: len(s) for ing, s in by_ingest.items()}
    log.info("[panel] loaded %d sample labels", report.loaded)
    return report
