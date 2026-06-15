"""Merge per-sample VCFs into one joint multi-sample VCF.

The motivating case is trios: a proband + father + mother are usually
delivered as three separate single-sample VCFs, but vcfclick's trio
analysis needs them as ONE joint cohort (one ingest_id) so a single
query can compare the three genotypes at a site. `merge_vcfs` produces
that joint VCF, ready for `vcfclick db ingest`.

This wraps `bcftools merge` rather than reimplementing it. A correct,
memory-efficient VCF merger IS bcftools merge — years of edge-case
handling (overlapping indels, FORMAT reconciliation, missing-genotype
fill). bcftools is already vcfclick's prerequisite for the
`bcftools norm -m -` decomposition step, so this adds no new dependency
category; it just orchestrates it with vcfclick-friendly defaults and
validation.

Scope / honesty: a sample absent from a site that another sample calls
becomes `./.` (missing) in the merged output, NOT `0/0`. A variant-only
single-sample VCF does not assert hom-reference at sites it omits, so
the merge cannot invent confident hom-ref calls. This is the same
limitation that makes trio de-novo detection "candidate" rather than
rigorous (see cli.db trio).
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)

BGZF_MAGIC = b"\x1f\x8b\x08\x04"


class MergeError(RuntimeError):
    """Raised for any precondition failure or bcftools error during merge."""


def _require_bcftools() -> str:
    path = shutil.which("bcftools")
    if not path:
        raise MergeError(
            "bcftools not found on PATH. vcfclick merge wraps `bcftools merge`; "
            "install bcftools (e.g. `conda install -c bioconda bcftools` or "
            "`brew install bcftools`) and retry."
        )
    return path


def _is_bgzipped(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(4) == BGZF_MAGIC


def _has_index(path: Path) -> bool:
    return (
        path.with_suffix(path.suffix + ".tbi").exists()
        or path.with_suffix(path.suffix + ".csi").exists()
    )


def _ensure_indexed(path: Path) -> None:
    """bcftools merge needs each input bgzipped + indexed. Index in place
    if the index is missing; error clearly if the file isn't bgzipped."""
    if not _is_bgzipped(path):
        raise MergeError(
            f"{path} is not bgzip-compressed. Compress it first: "
            f"`bgzip {path}` (and re-pass the .gz), or pass an already "
            f"bgzipped + indexed VCF."
        )
    if not _has_index(path):
        log.info("[merge] indexing %s", path)
        bcftools = _require_bcftools()
        proc = subprocess.run(
            [bcftools, "index", "--tbi", str(path)],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise MergeError(f"failed to index {path}: {proc.stderr.strip()}")


def merge_vcfs(
    inputs: list[str | Path],
    output: str | Path,
    *,
    multiallelic: str = "none",
    index_output: bool = True,
) -> Path:
    """Merge `inputs` (>=2 VCFs) into a joint bgzipped VCF at `output`.

    `multiallelic` is passed to `bcftools merge -m`; default "none" keeps
    records decomposed (one ALT per row) so the merged VCF stays
    ingest-ready — vcfclick rejects multi-allelic sites on ingest.

    Returns the output Path.
    """
    bcftools = _require_bcftools()

    in_paths = [Path(p) for p in inputs]
    if len(in_paths) < 2:
        raise MergeError("merge needs at least two input VCFs.")
    for p in in_paths:
        if not p.exists():
            raise MergeError(f"input not found: {p}")

    # Detect overlapping sample names up front — bcftools merge errors on
    # duplicate samples, but a vcfclick-level message is clearer.
    from cyvcf2 import VCF

    seen: dict[str, Path] = {}
    for p in in_paths:
        for s in VCF(str(p)).samples:
            if s in seen:
                raise MergeError(
                    f"sample {s!r} appears in both {seen[s]} and {p}. "
                    f"Inputs must have disjoint sample names."
                )
            seen[s] = p

    for p in in_paths:
        _ensure_indexed(p)

    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        bcftools,
        "merge",
        "-m",
        multiallelic,
        "-Oz",
        "-o",
        str(out_path),
        *[str(p) for p in in_paths],
    ]
    log.info("[merge] %d inputs → %s", len(in_paths), out_path)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise MergeError(f"bcftools merge failed: {proc.stderr.strip()}")

    if index_output:
        idx = subprocess.run(
            [bcftools, "index", "--tbi", str(out_path)],
            capture_output=True,
            text=True,
        )
        if idx.returncode != 0:
            # Non-fatal: the merged VCF is valid even without an index.
            log.warning("[merge] could not index output: %s", idx.stderr.strip())

    log.info("[merge] wrote %s (%d samples)", out_path, len(seen))
    return out_path
