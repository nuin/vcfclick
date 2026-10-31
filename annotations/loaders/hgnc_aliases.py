"""Load the HGNC complete set into the optional `gene_aliases` table.

    vcfclick annotations load-hgnc hgnc_complete_set.txt

Optional. Without it a gene list resolves by current symbol (any case) and
Ensembl gene id; with it, previous symbols, aliases, HGNC ids (`HGNC:1234`)
and Entrez ids resolve to the current symbol too (cli/scope.py). Panel files
from labs routinely carry retired symbols.

The file is the HGNC "complete set" TSV (download it yourself from
https://www.genenames.org/download/archive/; nothing is fetched here). Columns
are found by their header name, so a custom download with a subset of
columns works as long as `symbol` is there. Only rows whose `status` is
`Approved` are read when the column exists. `prev_symbol` and `alias_symbol`
hold several values separated by `|`, sometimes in double quotes.

The table lives in the annotation store (DuckDB), never the genotype
database.
"""

from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path

from annotations.db import get_connection

log = logging.getLogger(__name__)

# header -> the kind of token it holds
_LIST_COLUMNS = {"prev_symbol": "previous", "alias_symbol": "alias"}
_ID_COLUMNS = {"hgnc_id": "hgnc_id", "entrez_id": "entrez"}


def _values(cell: str | None) -> list[str]:
    if not cell:
        return []
    return [
        v.strip().strip('"') for v in cell.strip().strip('"').split("|") if v.strip()
    ]


def iter_aliases(path: Path):
    """Yield (alias, kind, symbol) for every approved gene in the TSV."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "symbol" not in reader.fieldnames:
            raise ValueError(f"{path.name}: not an HGNC TSV (no `symbol` column)")
        for row in reader:
            symbol = (row.get("symbol") or "").strip()
            if not symbol:
                continue
            status = row.get("status")
            if status is not None and status.strip() != "Approved":
                continue
            for column, kind in _ID_COLUMNS.items():
                value = (row.get(column) or "").strip()
                if value:
                    yield value, kind, symbol
            for column, kind in _LIST_COLUMNS.items():
                for value in _values(row.get(column)):
                    if value.upper() != symbol.upper():
                        yield value, kind, symbol


def load(path: Path, replace: bool = True) -> dict[str, int]:
    """Fill `gene_aliases` from an HGNC TSV. Returns {kind: rows written}."""
    rows = list(dict.fromkeys(iter_aliases(path)))
    if not rows:
        raise ValueError(f"{path.name}: no approved genes found")
    conn = get_connection()
    try:
        if replace:
            conn.execute("DELETE FROM gene_aliases")
        conn.executemany(
            "INSERT INTO gene_aliases (alias, kind, symbol) VALUES (?, ?, ?)", rows
        )
    finally:
        conn.close()
    counts: dict[str, int] = {}
    for _, kind, _ in rows:
        counts[kind] = counts.get(kind, 0) + 1
    log.info("[hgnc] loaded %s aliases into gene_aliases", f"{len(rows):,}")
    return counts
