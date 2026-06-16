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
and write a fresh VCF. v1 output carries GT + the `set=` provenance;
FORMAT passthrough (GQ/DP/AD from the priority source) is a future
refinement. Inputs must be on the same reference and decomposed
(one ALT per record), like every vcfclick input.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

log = logging.getLogger(__name__)


class CombineError(RuntimeError):
    """Raised for precondition failures during combine."""


def _gt_str(genotype) -> str | None:
    """Format a cyvcf2 genotype [a, b, phased] as 'a/b' / 'a|b', with
    '.' for missing alleles. Returns None if the call is fully missing
    (so PRIORITIZE skips it in favour of a lower-priority real call)."""
    if not genotype:
        return None
    *alleles, phased = genotype
    if all(a < 0 for a in alleles):
        return None  # ./.  — no information
    sep = "|" if phased else "/"
    return sep.join("." if a < 0 else str(a) for a in alleles)


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


def combine_vcfs(
    inputs: list[str | Path],
    output: str | Path,
    *,
    names: list[str] | None = None,
    min_callsets: int = 1,
) -> Path:
    """Combine `inputs` (>=2 VCFs, priority = input order) into `output`.

    `names` overrides the per-input set names used in the `set=` field
    (default: derived from filenames). `min_callsets` keeps only sites
    present in at least that many inputs.

    Returns the output Path.
    """
    from cyvcf2 import VCF

    in_paths = [Path(p) for p in inputs]
    if len(in_paths) < 2:
        raise CombineError("combine needs at least two input VCFs.")
    for p in in_paths:
        if not p.exists():
            raise CombineError(f"input not found: {p}")
    if min_callsets < 1 or min_callsets > len(in_paths):
        raise CombineError(
            f"--min-callsets must be between 1 and {len(in_paths)} (the number "
            f"of inputs), got {min_callsets}."
        )

    used: set[str] = set()
    set_names = names or [_default_name(p, used) for p in in_paths]
    if len(set_names) != len(in_paths):
        raise CombineError("number of --name values must match number of inputs.")

    # Union state, keyed by (chrom, pos, ref, alt).
    #   site_gts[key][sample]  = prioritised GT string (highest-priority
    #                            input with a non-missing call wins)
    #   site_inputs[key]       = set of input indices that have the site
    #   site_order[key]        = (contig_index, pos) for sorting
    site_gts: dict[tuple, dict[str, str]] = {}
    site_inputs: dict[tuple, set[int]] = {}
    contig_order: dict[str, int] = {}
    all_samples: list[str] = []
    seen_samples: set[str] = set()

    for idx, path in enumerate(in_paths):
        vcf = VCF(str(path))
        samples = list(vcf.samples)
        # Output sample order = first appearance across inputs (priority).
        for s in samples:
            if s not in seen_samples:
                seen_samples.add(s)
                all_samples.append(s)

        for variant in vcf:
            if variant.ALT and len(variant.ALT) != 1:
                raise CombineError(
                    f"{path} has a multi-allelic site at "
                    f"{variant.CHROM}:{variant.POS}. Decompose first: "
                    f"bcftools norm -m -."
                )
            alt = variant.ALT[0] if variant.ALT else "."
            key = (variant.CHROM, variant.POS, variant.REF, alt)
            contig_order.setdefault(variant.CHROM, len(contig_order))

            site_inputs.setdefault(key, set()).add(idx)
            gts = site_gts.setdefault(key, {})
            genotypes = variant.genotypes
            for s_i, sample in enumerate(samples):
                if sample in gts:
                    continue  # higher-priority input already filled it
                g = _gt_str(genotypes[s_i]) if s_i < len(genotypes) else None
                if g is not None:
                    gts[sample] = g

    # Write the combined VCF.
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_kept = _write_combined(
        out_path,
        all_samples,
        set_names,
        site_gts,
        site_inputs,
        contig_order,
        min_callsets,
    )
    log.info(
        "[combine] %d inputs → %s (%d sites, %d samples)",
        len(in_paths),
        out_path,
        n_kept,
        len(all_samples),
    )
    return out_path


def _set_field(present: set[int], n_inputs: int, set_names: list[str]) -> str:
    """The `set=` value: 'Intersection' if in all inputs, else the
    dash-joined names of the inputs that contain the site (priority
    order)."""
    if len(present) == n_inputs:
        return "Intersection"
    return "-".join(set_names[i] for i in sorted(present))


def _write_combined(
    out_path: Path,
    samples: list[str],
    set_names: list[str],
    site_gts: dict,
    site_inputs: dict,
    contig_order: dict,
    min_callsets: int,
) -> int:
    n_inputs = len(set_names)
    ordered_keys = sorted(
        (k for k in site_inputs if len(site_inputs[k]) >= min_callsets),
        key=lambda k: (contig_order.get(k[0], 0), k[1], k[2], k[3]),
    )

    header = [
        "##fileformat=VCFv4.3",
        '##source=vcfclick combine',
        (
            '##INFO=<ID=set,Number=1,Type=String,Description='
            '"Source call sets; Intersection when present in all inputs">'
        ),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    for contig in sorted(contig_order, key=lambda c: contig_order[c]):
        header.append(f"##contig=<ID={contig}>")
    header.append(
        "#" + "\t".join(
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + samples
        )
    )

    n = 0
    opener = gzip.open if str(out_path).endswith(".gz") else open
    with opener(out_path, "wt") as fh:
        fh.write("\n".join(header) + "\n")
        for key in ordered_keys:
            chrom, pos, ref, alt = key
            present = site_inputs[key]
            info = f"set={_set_field(present, n_inputs, set_names)}"
            gts = site_gts[key]
            cells = [gts.get(s, "./.") for s in samples]
            fh.write(
                "\t".join(
                    [chrom, str(pos), ".", ref, alt, ".", ".", info, "GT", *cells]
                )
                + "\n"
            )
            n += 1
    return n
