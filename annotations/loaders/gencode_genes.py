"""GENCODE GFF3 → DuckDB refseq_genes table.

GENCODE is used (not RefSeq's NCBI GFF) because:
  - It's GRCh38-native with 'chr'-prefixed contig names — matches our
    sample-data convention without any contig-name remapping.
  - Single download, clean attribute format, ubiquitous in research
    bioinformatics.
  - Same gene symbols as HGNC (which is what bioinformaticians type).

Only `gene` feature rows are extracted — we need name + coordinates.
Transcripts/exons/CDS land in `annotations/transcripts.py` in Phase 2
when that depth of annotation matters.

Usage:
    uv run python -m annotations.loaders.gencode_genes
    uv run python -m annotations.loaders.gencode_genes --gff path/to/local.gff3.gz
"""

from __future__ import annotations

import gzip
import logging
import time
import urllib.request
from pathlib import Path

from annotations.db import get_connection

log = logging.getLogger(__name__)

GENCODE_VERSION = "45"
GENCODE_URL = (
    f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    f"release_{GENCODE_VERSION}/gencode.v{GENCODE_VERSION}.annotation.gff3.gz"
)

CACHE_DIR = Path(__file__).parent / "_cache"
CACHED_GFF = CACHE_DIR / f"gencode.v{GENCODE_VERSION}.annotation.gff3.gz"


def download_gencode() -> Path:
    """Download the GENCODE GFF3 once and cache it. Returns the local path."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if CACHED_GFF.exists():
        log.info("[gencode] already cached: %s", CACHED_GFF)
        return CACHED_GFF
    log.info("[gencode] downloading %s → %s", GENCODE_URL, CACHED_GFF)
    started = time.time()
    urllib.request.urlretrieve(GENCODE_URL, CACHED_GFF)
    size_mb = CACHED_GFF.stat().st_size / 1_000_000
    log.info("[gencode] done (%.0f MB in %.1fs)", size_mb, time.time() - started)
    return CACHED_GFF


def parse_attributes(attrs: str) -> dict[str, str]:
    """GFF3 attribute column: semicolon-delimited key=value pairs."""
    out = {}
    for pair in attrs.split(";"):
        if "=" in pair:
            k, v = pair.split("=", 1)
            out[k.strip()] = v.strip()
    return out


# Primary GRCh38 contigs only: chr1-22, X, Y, M. Drops alt-locus
# (`_alt`), patch (`_fix`), random (`_random`), unplaced (`chrUn_`)
# scaffolds so a gene symbol resolves to exactly one canonical
# coordinate range. Research-bioinformatics queries almost never want
# alt-locus coords; a Phase 2 schema change can lift this when needed.
PRIMARY_CONTIGS = frozenset(
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
)


def iter_genes(gff_path: Path):
    """Stream a GFF3 file, yielding (gene_symbol, chrom, start, end, strand,
    refseq_id, description) tuples for each gene on a primary contig.

    Duplicate gene_symbols on the same primary contig are also possible
    (rare; usually PAR genes like XG that appear on both chrX and chrY).
    We keep the first occurrence and drop the rest; PAR resolution is a
    Phase 2 concern.
    """
    opener = gzip.open if gff_path.suffix == ".gz" else open
    seen: set[str] = set()
    with opener(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand, _, attrs_str = fields
            if chrom not in PRIMARY_CONTIGS:
                continue
            attrs = parse_attributes(attrs_str)
            # gene_name is the HGNC symbol (e.g. 'BRCA1'). gene_id is the
            # Ensembl ID (e.g. 'ENSG00000012048.24'). We prefer the name.
            symbol = attrs.get("gene_name")
            if not symbol or symbol in seen:
                continue
            seen.add(symbol)
            yield (
                symbol,
                chrom,
                int(start),
                int(end),
                strand,
                attrs.get("gene_id"),
                attrs.get("gene_type"),
            )


def load(gff_path: Path | None = None, replace: bool = True) -> int:
    """Populate the refseq_genes DuckDB table from a GENCODE GFF3 file.

    Returns the number of gene rows written. Pass `replace=False` to
    keep prior rows (DuckDB will error on duplicate primary keys).
    """
    if gff_path is None:
        gff_path = download_gencode()

    conn = get_connection()
    if replace:
        conn.execute("DELETE FROM refseq_genes")

    rows = list(iter_genes(gff_path))
    if not rows:
        raise RuntimeError(f"No gene features parsed from {gff_path}")

    conn.executemany(
        "INSERT INTO refseq_genes "
        "(gene_symbol, chrom, start_pos, end_pos, strand, refseq_id, description) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    log.info("[gencode] loaded %s genes into refseq_genes", f"{len(rows):,}")

    # Quick verification — the demo gene should be reachable.
    brca1 = conn.execute(
        "SELECT chrom, start_pos, end_pos FROM refseq_genes WHERE gene_symbol = 'BRCA1'"
    ).fetchone()
    if brca1:
        log.info("[gencode] BRCA1 → %s:%s-%s", brca1[0], brca1[1], brca1[2])
    else:
        log.warning("[gencode] WARNING: BRCA1 not found in loaded data")

    return len(rows)


# Library module — invoke via `vcfclick annotations load`.
# The public CLI lives in cli/main.py.
