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
    genotype_merge: str = "prioritize",
    merge_info_with_max_ac: bool = False,
    gatk3: bool = False,
    filtered_are_uncalled: bool = False,
    priority: list[str] | None = None,
) -> Path:
    """Combine `inputs` (>=2 VCFs, priority = input order) into `output`.

    `names` overrides the per-input set names used in the `set=` field
    (default: derived from filenames). `min_callsets` keeps only sites
    present in at least that many inputs.

    `gatk3=True` runs the GATK 3.8-1 emulation in `ingest.combine_gatk3`
    (no normalisation, GATK's own counting, QUAL/FILTER/INFO/header rules).
    `genotype_merge="unsorted"` keeps, in the normalising engine, the first
    input's record for a sample even when its call is `./.` (GATK's UNSORTED
    and PRIORITIZE-in-input-order). `merge_info_with_max_ac` takes INFO,
    QUAL and FILTER from the input with the highest INFO AC for the allele.

    Returns the output Path.
    """
    if genotype_merge not in ("prioritize", "unsorted"):
        raise CombineError(
            f"--genotype-merge must be 'prioritize' or 'unsorted', got {genotype_merge!r}."
        )
    if gatk3:
        return _combine_gatk3(
            inputs,
            output,
            names,
            min_callsets,
            reference,
            genotype_merge,
            merge_info_with_max_ac,
            filtered_are_uncalled,
            priority,
            conflicting=dict(
                pass_only=pass_only,
                carry_info=carry_info,
                atomize=atomize,
                count_by=count_by != "allele",
            ),
        )
    if filtered_are_uncalled or priority:
        raise CombineError("--filtered-are-uncalled and --priority need --gatk3.")
    in_paths, set_names = _validate_inputs(
        inputs, names, min_callsets, count_by, reference, atomize
    )
    carry_info = carry_info or merge_info_with_max_ac
    combined = CombinedInputs(
        carry_info, genotype_merge == "unsorted", merge_info_with_max_ac
    )
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


def _combine_gatk3(
    inputs,
    output,
    names,
    min_callsets,
    reference,
    genotype_merge,
    merge_info_with_max_ac,
    filtered_are_uncalled,
    priority,
    conflicting,
) -> Path:
    from ingest.combine_gatk3 import combine_gatk3

    bad = [k.replace("_", "-") for k, v in conflicting.items() if v]
    if bad:
        raise CombineError(
            "--gatk3 does not take "
            + ", ".join("--count-by site" if b == "count-by" else "--" + b for b in bad)
            + "."
        )
    in_paths = [Path(p) for p in inputs]
    if len(in_paths) < 2:
        raise CombineError("combine needs at least two input VCFs.")
    for p in in_paths:
        if not p.exists():
            raise CombineError(f"input not found: {p}")
    if reference is not None and not Path(reference).exists():
        raise CombineError(f"reference not found: {reference}")
    used: set[str] = set()
    set_names = list(names) if names else [_default_name(p, used) for p in in_paths]
    if len(set_names) != len(in_paths):
        raise CombineError("number of --name values must match number of inputs.")
    if len(set(set_names)) != len(set_names):
        raise CombineError("input names must be unique.")
    if min_callsets < 1:
        raise CombineError("--min-callsets must be at least 1.")
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    gz = str(out_path).endswith(".gz")
    plain_path = out_path.with_suffix("") if gz else out_path
    # GATK only sorts by priority under PRIORITIZE; UNSORTED keeps input order.
    prio = None
    if genotype_merge == "prioritize" and priority:
        prio = list(priority)
    elif priority:
        raise CombineError("--priority needs --genotype-merge prioritize.")
    n = combine_gatk3(
        in_paths,
        plain_path,
        set_names,
        min_callsets=min_callsets,
        merge_info_with_max_ac=merge_info_with_max_ac,
        filtered_are_uncalled=filtered_are_uncalled,
        priority=prio,
        reference=reference,
    )
    if gz:
        _bgzip_and_index(plain_path, out_path)
    log.info("[combine] gatk3 %d inputs → %s (%d records)", len(in_paths), out_path, n)
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
