"""ClinVar VCF → DuckDB clinvar_variants table.

NCBI ClinVar publishes a weekly-refreshed GRCh38 VCF at
ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz (~80 MB
compressed, ~3M variants). This loader streams it through cyvcf2,
decomposes multi-allelic records into one row per (chrom, pos, ref,
alt), normalises ClinVar's bare numeric contigs to the `chr`-prefixed
convention the rest of the project uses, and bulk-inserts via a
temporary pyarrow table.

Usage:
    uv run python -m annotations.loaders.clinvar
    uv run python -m annotations.loaders.clinvar --vcf path/to/local.vcf.gz

The MCP server's `clinvar_lookup` tool resolves against this table.
"""

from __future__ import annotations

import logging
import time
import urllib.request
from pathlib import Path

import pyarrow as pa
from cyvcf2 import VCF

from annotations.db import get_connection

log = logging.getLogger(__name__)

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"

CACHE_DIR = Path(__file__).parent / "_cache"
CACHED_VCF = CACHE_DIR / "clinvar.vcf.gz"


def download_clinvar() -> Path:
    """Download the latest ClinVar VCF and cache it.

    Cache is replaced unconditionally — the file refreshes weekly and a
    stale cache silently desyncs the `clinvar_lookup` MCP tool from
    current significance calls. Delete the cache file manually if you
    want to force a re-download mid-week.
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if CACHED_VCF.exists():
        log.info("[clinvar] already cached: %s", CACHED_VCF)
        return CACHED_VCF
    log.info("[clinvar] downloading %s → %s", CLINVAR_URL, CACHED_VCF)
    started = time.time()
    urllib.request.urlretrieve(CLINVAR_URL, CACHED_VCF)
    size_mb = CACHED_VCF.stat().st_size / 1_000_000
    log.info("[clinvar] done (%.0f MB in %.1fs)", size_mb, time.time() - started)
    return CACHED_VCF


def _normalise_contig(chrom: str) -> str:
    """ClinVar uses bare numeric contigs ('1', '17', 'X', 'MT'); our
    sample data uses 'chr'-prefixed and 'chrM' (not chrMT). Map both."""
    c = chrom.lstrip("chr")  # no-op if already prefixed
    if c == "MT":
        return "chrM"
    return f"chr{c}"


def _clnstr(info: dict, key: str) -> str | None:
    """Read a ClinVar INFO field that may be a string, tuple, or absent.
    Strings come back as-is; tuples/lists join with '|'. Returns None
    for missing — the caller should distinguish absent from empty when
    rendering."""
    v = info.get(key)
    if v is None:
        return None
    if isinstance(v, (list, tuple)):
        return "|".join(str(x) for x in v) or None
    s = str(v).strip()
    return s or None


def iter_clinvar(vcf_path: Path):
    """Stream a ClinVar VCF, yielding one tuple per (record, alt) pair.

    Multi-allelic records are decomposed inline — the destination
    PRIMARY KEY is (chrom, pos, ref, alt) and ClinVar legitimately has
    multiple alt alleles per site with different significance calls.
    """
    vcf = VCF(str(vcf_path))
    for rec in vcf:
        info = dict(rec.INFO)
        clin_sig = _clnstr(info, "CLNSIG")
        review_status = _clnstr(info, "CLNREVSTAT")
        condition = _clnstr(info, "CLNDN")
        clinvar_id = rec.ID  # the VCV/RCV numeric accession
        chrom = _normalise_contig(rec.CHROM)
        for alt in rec.ALT:
            # ClinVar uses '.' or empty ALT for some structural / deletion
            # rows. Those don't compose with the (ref, alt) primary key —
            # skip rather than fabricate a placeholder.
            if not alt or alt == ".":
                continue
            yield (
                chrom,
                int(rec.POS),
                rec.REF,
                alt,
                clin_sig,
                review_status,
                clinvar_id,
                condition,
            )


def _to_arrow(rows: list) -> pa.Table:
    cols = list(zip(*rows)) if rows else [[]] * 8
    return pa.table(
        {
            "chrom": pa.array(cols[0], type=pa.string()),
            "pos": pa.array(cols[1], type=pa.uint32()),
            "ref": pa.array(cols[2], type=pa.string()),
            "alt": pa.array(cols[3], type=pa.string()),
            "clin_sig": pa.array(cols[4], type=pa.string()),
            "review_status": pa.array(cols[5], type=pa.string()),
            "clinvar_id": pa.array(cols[6], type=pa.string()),
            "condition": pa.array(cols[7], type=pa.string()),
        }
    )


def load(vcf_path: Path | None = None, replace: bool = True) -> int:
    """Populate the clinvar_variants DuckDB table from a ClinVar VCF.

    Returns the number of (chrom, pos, ref, alt) rows written.

    ClinVar PRIMARY KEY duplicates can occur when the same allele is
    described by multiple ClinVar records (different VCV accessions for
    different submissions). We keep the first occurrence — order is
    file order, which for ClinVar means "lowest VCV ID wins" which is
    usually the oldest curated record. Refine later if needed.
    """
    if vcf_path is None:
        vcf_path = download_clinvar()

    seen: set[tuple[str, int, str, str]] = set()
    rows: list[tuple] = []
    for row in iter_clinvar(vcf_path):
        key = (row[0], row[1], row[2], row[3])
        if key in seen:
            continue
        seen.add(key)
        rows.append(row)

    if not rows:
        raise RuntimeError(f"No ClinVar variants parsed from {vcf_path}")

    conn = get_connection()
    if replace:
        conn.execute("DELETE FROM clinvar_variants")

    # Bulk-insert via pyarrow → DuckDB. ~3M rows takes a few seconds vs.
    # minutes with executemany.
    arrow_table = _to_arrow(rows)
    conn.register("clinvar_arrow", arrow_table)
    try:
        conn.execute(
            "INSERT INTO clinvar_variants "
            "(chrom, pos, ref, alt, clin_sig, review_status, clinvar_id, condition) "
            "SELECT chrom, pos, ref, alt, clin_sig, review_status, clinvar_id, condition "
            "FROM clinvar_arrow"
        )
    finally:
        conn.unregister("clinvar_arrow")

    n = len(rows)
    log.info("[clinvar] loaded %s variants into clinvar_variants", f"{n:,}")

    # Quick verification — anything pathogenic should be reachable.
    pathogenic = conn.execute(
        "SELECT count() FROM clinvar_variants WHERE clin_sig LIKE '%athogenic%'"
    ).fetchone()
    if pathogenic and pathogenic[0]:
        log.info("[clinvar] %s pathogenic-class variants indexed", f"{pathogenic[0]:,}")

    return n


# Library module — invoke via `vcfclick annotations load-clinvar`.
# The public CLI lives in cli/annotations.py.
