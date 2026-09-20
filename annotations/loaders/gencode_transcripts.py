"""Load the transcript / exon / CDS hierarchy from a GENCODE GFF3.

The same file `gencode_genes` already downloads and caches: that loader keeps
only `gene` rows, this one keeps `transcript`, `exon`, and `CDS`. No new data
source, no second download.

MANE Select is read from the transcript's `tag=` attribute (comma-separated),
which is what makes `canonical_transcript()` meaningful — it avoids spurious
hits on rare isoforms.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from annotations.db import get_connection

log = logging.getLogger(__name__)

_WANTED = {"transcript", "exon", "CDS"}


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_attributes(attrs: str) -> dict[str, str]:
    """GFF3 attribute column: semicolon-delimited key=value pairs."""
    out: dict[str, str] = {}
    for pair in attrs.split(";"):
        if "=" in pair:
            k, v = pair.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def _exon_number(attrs: dict[str, str], fallback: int) -> int:
    raw = attrs.get("exon_number")
    if raw and raw.isdigit():
        return int(raw)
    return fallback


def iter_features(gff_path: Path):
    """Yield (kind, row) for transcript/exon/CDS features.

    `row` is the tuple in the order its table's INSERT expects.
    """
    seen_exons: dict[str, int] = {}
    seen_cds: dict[str, int] = {}
    with _open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _src, kind, start, end, _sc, strand, phase, attrs_str = fields[:9]
            if kind not in _WANTED:
                continue
            attrs = parse_attributes(attrs_str)
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            start_i, end_i = int(start), int(end)

            if kind == "transcript":
                tags = set(attrs.get("tag", "").split(","))
                yield (
                    "transcript",
                    (
                        tid,
                        attrs.get("gene_name", ""),
                        chrom,
                        start_i,
                        end_i,
                        strand,
                        attrs.get("transcript_type"),
                        "MANE_Select" in tags,
                        "MANE_Plus_Clinical" in tags,
                    ),
                )
            elif kind == "exon":
                n = seen_exons.get(tid, 0) + 1
                seen_exons[tid] = n
                yield ("exon", (tid, _exon_number(attrs, n), chrom, start_i, end_i))
            else:  # CDS
                n = seen_cds.get(tid, 0) + 1
                seen_cds[tid] = n
                yield (
                    "cds",
                    (
                        tid,
                        _exon_number(attrs, n),
                        chrom,
                        start_i,
                        end_i,
                        int(phase) if phase.isdigit() else None,
                    ),
                )


def load(gff_path: Path | None = None, replace: bool = True) -> dict[str, int]:
    """Populate transcripts/exons/cds from a GENCODE GFF3.

    Returns `{"transcripts": n, "exons": n, "cds": n}`. Pass `replace=False`
    to keep prior rows.
    """
    if gff_path is None:
        from annotations.loaders.gencode_genes import download_gencode

        gff_path = download_gencode()

    conn = get_connection()
    if replace:
        for t in ("cds", "exons", "transcripts"):
            conn.execute(f"DELETE FROM {t}")

    batches: dict[str, list[tuple]] = {"transcript": [], "exon": [], "cds": []}
    for kind, row in iter_features(Path(gff_path)):
        batches[kind].append(row)

    if not batches["transcript"]:
        raise RuntimeError(f"No transcript features parsed from {gff_path}")

    conn.executemany(
        "INSERT INTO transcripts (transcript_id, gene_symbol, chrom, start_pos, "
        "end_pos, strand, biotype, is_canonical, is_mane_plus) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        batches["transcript"],
    )
    conn.executemany(
        "INSERT INTO exons (transcript_id, exon_number, chrom, start_pos, end_pos) "
        "VALUES (?, ?, ?, ?, ?)",
        batches["exon"],
    )
    conn.executemany(
        "INSERT INTO cds (transcript_id, exon_number, chrom, start_pos, end_pos, "
        "phase) VALUES (?, ?, ?, ?, ?, ?)",
        batches["cds"],
    )
    counts = {
        "transcripts": len(batches["transcript"]),
        "exons": len(batches["exon"]),
        "cds": len(batches["cds"]),
    }
    log.info(
        "[gencode] loaded %(transcripts)s transcripts, %(exons)s exons, "
        "%(cds)s CDS rows",
        {k: f"{v:,}" for k, v in counts.items()},
    )
    return counts


# Library module — invoke via `vcfclick annotations load-transcripts`.
