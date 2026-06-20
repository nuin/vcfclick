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
be on the same reference and decomposed (one ALT per record), like every
vcfclick input.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import NamedTuple

log = logging.getLogger(__name__)


class CombineError(RuntimeError):
    """Raised for precondition failures during combine."""


# FORMAT fields carried through from the priority source — exactly the
# ones trio quality gates read (gq, dp, ad_ref/ad_alt). Output FORMAT is
# GT plus whichever of these actually appear in some input.
_PASSTHROUGH = ("GQ", "DP", "AD")

_FORMAT_HEADERS = {
    "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "GQ": '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
    "AD": (
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths '
        'for the ref and alt alleles in the order listed">'
    ),
}


class _Union(NamedTuple):
    """The accumulated cross-input union, keyed by (chrom, pos, ref, alt).

    gts[key][sample] — the prioritized per-sample cell: a dict with "GT"
    (the highest-priority input with a non-missing call wins) plus the
    GQ/DP/AD tokens from that same source record (None where absent).
    inputs[key] — set of input indices that have the site (drives set= and
    --min-callsets). contigs[contig] — reference (header) order, for
    coordinate-sorting the output. fields — which _PASSTHROUGH FORMAT
    fields appeared in any input, so the output FORMAT lists only those.
    """

    gts: dict
    inputs: dict
    contigs: dict
    fields: set


def _format_arr(variant, field):
    """A FORMAT field as a per-sample cyvcf2 array, or None if absent."""
    try:
        return variant.format(field)
    except KeyError:
        return None


def _scalar_token(arr, i: int) -> str | None:
    """Sample i's scalar FORMAT value (GQ/DP) as a VCF token, or None."""
    if arr is None:
        return None
    try:
        v = arr[i]
        if hasattr(v, "__len__") and not isinstance(v, (str, bytes)):
            v = v[0]
        v = int(v)
    except (IndexError, TypeError, ValueError):
        return None
    return str(v) if v >= 0 else None


def _ad_token(arr, i: int) -> str | None:
    """Sample i's AD as 'ref,alt', or None if absent/missing. Output
    records are always biallelic, so AD is Number=R = exactly two values.
    A source AD with any other length (e.g. an improperly decomposed input
    still carrying original multi-allelic depths) is dropped rather than
    written as an uninterpretable cell for a single-ALT record."""
    if arr is None:
        return None
    try:
        vals = [int(x) for x in arr[i]]
    except (IndexError, TypeError, ValueError):
        return None
    if len(vals) != 2 or any(x < 0 for x in vals):
        return None
    return ",".join(str(x) for x in vals)


def _sample_cell(gt: str, variant, fmt_arrs: dict, i: int) -> dict:
    """Build sample i's output cell from the source record: its GT plus
    the GQ/DP/AD tokens, so passed-through quality travels with the
    genotype it describes."""
    cell = {"GT": gt}
    for f in _PASSTHROUGH:
        cell[f] = (
            _ad_token(fmt_arrs[f], i) if f == "AD" else _scalar_token(fmt_arrs[f], i)
        )
    return cell


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

    union = _Union(gts={}, inputs={}, contigs={}, fields=set())
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
        # Seed contig order from the header sequence dictionary (reference
        # order), not from whichever variant happens to appear first — an
        # input that starts on chr2 must not push chr2 ahead of chr1.
        for contig in vcf.seqnames:
            union.contigs.setdefault(contig, len(union.contigs))

        for variant in vcf:
            if variant.ALT and len(variant.ALT) != 1:
                raise CombineError(
                    f"{path} has a multi-allelic site at "
                    f"{variant.CHROM}:{variant.POS}. Decompose first: "
                    f"bcftools norm -m -."
                )
            alt = variant.ALT[0] if variant.ALT else "."
            key = (variant.CHROM, variant.POS, variant.REF, alt)
            union.contigs.setdefault(variant.CHROM, len(union.contigs))

            union.inputs.setdefault(key, set()).add(idx)
            gts = union.gts.setdefault(key, {})
            genotypes = variant.genotypes
            fmt_arrs = {f: _format_arr(variant, f) for f in _PASSTHROUGH}
            for f, arr in fmt_arrs.items():
                if arr is not None:
                    union.fields.add(f)
            for s_i, sample in enumerate(samples):
                if sample in gts:
                    continue  # higher-priority input already filled it
                g = _gt_str(genotypes[s_i]) if s_i < len(genotypes) else None
                if g is not None:
                    gts[sample] = _sample_cell(g, variant, fmt_arrs, s_i)

    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Write plain VCF text first. For a .gz target, write the uncompressed
    # body next to it, then bgzip + tabix so the result is BGZF (not plain
    # gzip) and tabix-indexable — the format every other vcfclick VCF and
    # the region-parallel ingest path assume.
    gz = str(out_path).endswith(".gz")
    plain_path = out_path.with_suffix("") if gz else out_path
    n_kept = _write_combined(plain_path, all_samples, set_names, union, min_callsets)
    if gz:
        _bgzip_and_index(plain_path, out_path)
    log.info(
        "[combine] %d inputs → %s (%d sites, %d samples)",
        len(in_paths),
        out_path,
        n_kept,
        len(all_samples),
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
    union: _Union,
    min_callsets: int,
) -> int:
    n_inputs = len(set_names)
    ordered_keys = sorted(
        (k for k in union.inputs if len(union.inputs[k]) >= min_callsets),
        key=lambda k: (union.contigs.get(k[0], 0), k[1], k[2], k[3]),
    )
    # Output FORMAT = GT plus whichever passthrough fields any input had.
    out_fields = ["GT"] + [f for f in _PASSTHROUGH if f in union.fields]
    format_col = ":".join(out_fields)
    missing_cell = ":".join(["./." if f == "GT" else "." for f in out_fields])

    header = [
        "##fileformat=VCFv4.3",
        "##source=vcfclick combine",
        (
            "##INFO=<ID=set,Number=1,Type=String,Description="
            '"Source call sets; Intersection when present in all inputs">'
        ),
    ]
    header += [_FORMAT_HEADERS[f] for f in out_fields]
    for contig in sorted(union.contigs, key=lambda c: union.contigs[c]):
        header.append(f"##contig=<ID={contig}>")
    header.append(
        "#"
        + "\t".join(
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + samples
        )
    )

    n = 0
    with open(out_path, "w") as fh:
        fh.write("\n".join(header) + "\n")
        for key in ordered_keys:
            chrom, pos, ref, alt = key
            present = union.inputs[key]
            info = f"set={_set_field(present, n_inputs, set_names)}"
            gts = union.gts[key]
            cells = [
                _render_cell(gts.get(s), out_fields, missing_cell) for s in samples
            ]
            fh.write(
                "\t".join(
                    [chrom, str(pos), ".", ref, alt, ".", ".", info, format_col, *cells]
                )
                + "\n"
            )
            n += 1
    return n


def _render_cell(cell: dict | None, out_fields: list[str], missing: str) -> str:
    """A sample's FORMAT cell text: its GT plus each output field token
    ('.' where that field was absent in the source). A sample with no
    record at this site is the all-missing cell."""
    if cell is None:
        return missing
    return ":".join(cell.get(f) or "." for f in out_fields)
