"""Combine multiple VCF call sets into one — the GATK3 CombineVariants
functionality that GATK4 dropped and never fully replaced.

Unlike `vcfclick merge` (which wraps `bcftools merge` to join *disjoint*
samples into a multi-sample VCF), `combine` merges *call sets* that may
share samples — e.g. the same cohort called by two different callers,
or pre/post a filter — and:

  * unions the sites across all inputs;
  * annotates each output record with `set=` showing which inputs
    contain it ("Intersection" when all do);
  * resolves a sample present in multiple inputs by PRIORITY — the
    genotype is taken from the highest-priority input (input order =
    priority) that has a non-missing call for that sample;
  * can keep only sites present in at least N inputs (consensus) — the
    "variants present in all / a fraction of the call sets" feature
    GATK4 specifically lost.

There is no bcftools equivalent for this, so it is implemented
natively: read each input with cyvcf2, union by (chrom, pos, ref, alt),
and write a fresh VCF. A plain `.vcf` output is fully native; a `.gz`
output is written then bgzip + tabix-indexed (via htslib) so it is BGZF,
not plain gzip — the format the rest of vcfclick (and region-parallel
ingest) assumes. Output carries GT + the `set=` provenance, plus the
GQ/DP/AD FORMAT fields (the ones trio quality gates read) carried from
the same priority-source record that supplied each genotype. Inputs must
be on the same reference and, by default, decomposed (one ALT per record)
like every vcfclick input — or pass `reference=` to split and left-align
multi-allelic input internally.

Opt-in GATK CombineVariants parity (defaults preserve the behaviour above):
`pass_only` counts only PASS calls toward the consensus filter (naming
filtered inputs `filterIn<name>`); `count_by="site"` counts by position
instead of by allele; `reference=` normalizes internally; `carry_info`
carries QUAL/FILTER/INFO from the priority input.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from ingest.combine_reader import AlleleNormalizer, CombinedInputs
from ingest.combine_records import CombineError as CombineError
from ingest.combine_writer import _write_combined

log = logging.getLogger(__name__)


def _require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise CombineError(
            f"{name} not found on PATH. Writing a .gz output requires "
            f"htslib's {name} (so the result is BGZF + tabix-indexable like "
            f"every other vcfclick VCF). Install htslib "
            f"(`brew install htslib` / `conda install -c bioconda htslib`), "
            f"or write a plain .vcf output instead."
        )
    return path


def _default_name(path: Path, used: set[str]) -> str:
    base = path.name
    for suffix in (".vcf.gz", ".vcf.bgz", ".vcf"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    name = base or "set"
    n = name
    i = 1
    while n in used:
        i += 1
        n = f"{name}.{i}"
    used.add(n)
    return n


def _validate_inputs(inputs, names, min_callsets, count_by, reference, atomize):
    in_paths = [Path(p) for p in inputs]
    if len(in_paths) < 2:
        raise CombineError("combine needs at least two input VCFs.")
    for p in in_paths:
        if not p.exists():
            raise CombineError(f"input not found: {p}")
    if atomize and reference is None:
        raise CombineError(
            "--atomize requires --reference (it follows left-alignment)."
        )
    if count_by not in ("allele", "site"):
        raise CombineError(f"--count-by must be 'allele' or 'site', got {count_by!r}.")
    if min_callsets < 1 or min_callsets > len(in_paths):
        raise CombineError(
            f"--min-callsets must be between 1 and {len(in_paths)} (the number "
            f"of inputs), got {min_callsets}."
        )

    used: set[str] = set()
    set_names = names or [_default_name(p, used) for p in in_paths]
    if len(set_names) != len(in_paths):
        raise CombineError("number of --name values must match number of inputs.")

    return in_paths, set_names


def combine_vcfs(
    inputs: list[str | Path],
    output: str | Path,
    *,
    names: list[str] | None = None,
    min_callsets: int = 1,
    pass_only: bool = False,
    count_by: str = "allele",
    carry_info: bool = False,
    reference: str | Path | None = None,
    atomize: bool = False,
) -> Path:
    """Combine `inputs` (>=2 VCFs, priority = input order) into `output`.

    `names` overrides the per-input set names used in the `set=` field
    (default: derived from filenames). `min_callsets` keeps only sites
    present in at least that many inputs.

    Returns the output Path.
    """
    in_paths, set_names = _validate_inputs(
        inputs, names, min_callsets, count_by, reference, atomize
    )
    combined = CombinedInputs(carry_info)
    normalizer = AlleleNormalizer(reference, atomize, combined.union.contigs)
    combined.read(in_paths, normalizer)
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # A compressed target needs BGZF and tabix, rather than plain gzip.
    gz = str(out_path).endswith(".gz")
    plain_path = out_path.with_suffix("") if gz else out_path
    n_kept = _write_combined(
        plain_path,
        combined.samples,
        set_names,
        combined.union,
        min_callsets,
        pass_only,
        count_by,
        carry_info,
        combined.info_headers,
        combined.filter_headers,
    )
    if gz:
        _bgzip_and_index(plain_path, out_path)
    log.info(
        "[combine] %d inputs → %s (%d sites, %d samples)",
        len(in_paths),
        out_path,
        n_kept,
        len(combined.samples),
    )
    return out_path


def _bgzip_and_index(plain_path: Path, gz_path: Path) -> None:
    """Compress `plain_path` to BGZF at `gz_path` and build a tabix index.

    `bgzip <file>` replaces `file` with `file.gz`; since plain_path is
    gz_path without the .gz suffix, the result lands exactly at gz_path.
    """
    bgzip = _require_tool("bgzip")
    tabix = _require_tool("tabix")
    proc = subprocess.run(
        [bgzip, "-f", str(plain_path)], capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise CombineError(f"bgzip failed: {proc.stderr.strip()}")
    proc = subprocess.run(
        [tabix, "-f", "-p", "vcf", str(gz_path)], capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise CombineError(f"tabix index failed: {proc.stderr.strip()}")
