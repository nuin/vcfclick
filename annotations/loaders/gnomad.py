"""gnomAD sites VCF → DuckDB gnomad_af table.

gnomAD population allele frequencies for rarity filtering. Unlike
ClinVar, the gnomAD release is far too large to bundle or auto-download
(tens of GB per chromosome), so this loader takes a VCF the caller
supplies — a region slice, an `af-only` file, or the full per-chromosome
sites VCF. It reads `AF` (overall) and `AF_grpmax` (the highest
genetic-ancestry-group AF, i.e. popmax) per allele.

The gnomAD sites VCFs are already decomposed (one ALT per record). A
small region can be pulled directly with tabix-over-HTTPS, e.g.:

    tabix https://storage.googleapis.com/gcp-public-data--gnomad/release/\\
4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr7.vcf.bgz chr7:117480000-117510000

Loaded via `vcfclick annotations load-gnomad <vcf>`; the MCP
`gnomad_lookup` tool and `db trio --gnomad-max-af` resolve against it.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pyarrow as pa
from cyvcf2 import VCF

from annotations.db import get_connection

log = logging.getLogger(__name__)


def _af(info, key: str) -> float | None:
    """Read a per-allele AF INFO value as a float, or None. cyvcf2 returns
    a scalar for a decomposed record's Number=A field; tolerate a tuple."""
    v = info.get(key)
    if v is None:
        return None
    if isinstance(v, (list, tuple)):
        v = v[0] if v else None
    try:
        return None if v is None else float(v)
    except (TypeError, ValueError):
        return None


def iter_gnomad(vcf_path: Path):
    """Stream a gnomAD sites VCF, yielding (chrom, pos, ref, alt, af,
    af_grpmax). Multi-allelic records are skipped — gnomAD sites VCFs are
    decomposed, so a multi-allelic line means a non-standard input."""
    vcf = VCF(str(vcf_path))
    for rec in vcf:
        if not rec.ALT or len(rec.ALT) != 1:
            continue
        info = dict(rec.INFO)
        yield (
            rec.CHROM,
            int(rec.POS),
            rec.REF,
            rec.ALT[0],
            _af(info, "AF"),
            _af(info, "AF_grpmax"),
        )


def _to_arrow(rows: list) -> pa.Table:
    cols = list(zip(*rows)) if rows else [[]] * 6
    return pa.table(
        {
            "chrom": pa.array(cols[0], type=pa.string()),
            "pos": pa.array(cols[1], type=pa.uint32()),
            "ref": pa.array(cols[2], type=pa.string()),
            "alt": pa.array(cols[3], type=pa.string()),
            "af": pa.array(cols[4], type=pa.float64()),
            "af_grpmax": pa.array(cols[5], type=pa.float64()),
        }
    )


def load(vcf_path: Path, replace: bool = False) -> int:
    """Populate the gnomad_af table from a gnomAD sites VCF.

    `replace=True` clears the table first; the default appends, so several
    per-chromosome slices can be loaded incrementally. Returns the number
    of rows written. Duplicate (chrom, pos, ref, alt) keys keep the first
    occurrence.
    """
    seen: set[tuple[str, int, str, str]] = set()
    rows: list[tuple] = []
    for row in iter_gnomad(Path(vcf_path)):
        key = (row[0], row[1], row[2], row[3])
        if key in seen:
            continue
        seen.add(key)
        rows.append(row)

    if not rows:
        raise RuntimeError(f"No gnomAD variants parsed from {vcf_path}")

    conn = get_connection()
    if replace:
        conn.execute("DELETE FROM gnomad_af")

    arrow_table = _to_arrow(rows)
    conn.register("gnomad_arrow", arrow_table)
    try:
        # Append, replacing any colliding key already loaded.
        conn.execute(
            "DELETE FROM gnomad_af WHERE (chrom, pos, ref, alt) IN "
            "(SELECT chrom, pos, ref, alt FROM gnomad_arrow)"
        )
        conn.execute(
            "INSERT INTO gnomad_af (chrom, pos, ref, alt, af, af_grpmax) "
            "SELECT chrom, pos, ref, alt, af, af_grpmax FROM gnomad_arrow"
        )
    finally:
        conn.unregister("gnomad_arrow")

    log.info("[gnomad] loaded %s allele frequencies into gnomad_af", f"{len(rows):,}")
    return len(rows)
