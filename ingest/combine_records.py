"""Shared VCF cell and header handling for call-set combination."""

from __future__ import annotations

import re
from typing import NamedTuple


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
    passes: dict  # key -> set of input idx whose call at that allele is PASS
    site_inputs: dict  # (chrom,pos) -> set of input idx with any record there
    site_passes: dict  # (chrom,pos) -> set of input idx with a PASS record there
    meta: dict  # key -> {"qual","filter","info"} from the priority input (carry_info)
    contigs: dict
    fields: set


def _hdr_lines_by_id(raw_header: str, tag: str) -> dict[str, str]:
    """The `##<tag>=<ID=...>` header lines from a raw header, keyed by ID —
    so carried INFO/FILTER fields keep their definitions in the output."""
    out: dict[str, str] = {}
    for line in raw_header.splitlines():
        if line.startswith(f"##{tag}=<"):
            m = re.search(r"ID=([^,>]+)", line)
            if m:
                out.setdefault(m.group(1), line)
    return out


def _hdr_fields(line: str) -> tuple[str, str | None]:
    """(Number, Type) declared by an INFO header line, for compatibility checks."""
    num = re.search(r"Number=([^,>]+)", line)
    typ = re.search(r"Type=([^,>]+)", line)
    return (num.group(1) if num else ".", typ.group(1) if typ else None)


def _project_info(info: str, alt_index: int | None, numbers: dict[str, str]) -> str:
    """Project a record's INFO onto a single (split) allele: drop the input's own
    `set` (combine recomputes it), and — when the record was split (`alt_index`
    set) — subset `Number=A` fields to that alt and `Number=R` fields to
    [ref, alt]. `Number=G` fields are dropped (they need genotype context to
    project). Returns `.` when nothing remains."""
    if info in ("", "."):
        return "."
    out: list[str] = []
    for field in info.split(";"):
        if not field:
            continue
        key, _, val = field.partition("=")
        if key == "set":
            continue
        if not val or alt_index is None:
            out.append(field)
            continue
        num = numbers.get(key, ".")
        if num == "A":
            parts = val.split(",")
            if alt_index - 1 < len(parts):
                val = parts[alt_index - 1]
        elif num == "R":
            parts = val.split(",")
            if alt_index < len(parts):
                val = f"{parts[0]},{parts[alt_index]}"
        elif num == "G":
            continue  # genotype-cardinality: cannot project per-allele
        out.append(f"{key}={val}")
    return ";".join(out) if out else "."


def _is_pass(variant) -> bool:
    """A record counts as PASS when FILTER is PASS, `.`, or empty — cyvcf2
    reports all of those as None or the literal string (§6: `.` is PASS)."""
    f = variant.FILTER
    return f is None or f in ("PASS", ".", "")


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


def _ad_token(arr, i: int, alt_index: int | None = None) -> str | None:
    """Sample i's AD as 'ref,alt' for the output biallelic record, or None if
    absent/missing. When `--reference` split a multi-allelic record, `alt_index`
    (1-based) selects [ref, alt_j] from the source AD [ref, alt1, ...]; otherwise
    the source AD must already be exactly two values."""
    if arr is None:
        return None
    try:
        vals = [int(x) for x in arr[i]]
    except (IndexError, TypeError, ValueError):
        return None
    if alt_index is not None:
        if alt_index < len(vals) and vals[0] >= 0 and vals[alt_index] >= 0:
            return f"{vals[0]},{vals[alt_index]}"
        return None
    if len(vals) != 2 or any(x < 0 for x in vals):
        return None
    return ",".join(str(x) for x in vals)


def _sample_cell(gt: str, variant, fmt_arrs: dict, i: int, alt_index=None) -> dict:
    """Build sample i's output cell from the source record: its GT plus
    the GQ/DP/AD tokens, so passed-through quality travels with the
    genotype it describes. `alt_index` subsets AD when a record was split."""
    cell = {"GT": gt}
    for f in _PASSTHROUGH:
        cell[f] = (
            _ad_token(fmt_arrs[f], i, alt_index)
            if f == "AD"
            else _scalar_token(fmt_arrs[f], i)
        )
    return cell


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


def _remap_gt(genotype, alt_index: int) -> str | None:
    """Remap a cyvcf2 genotype to the biallelic locus for `alt_index` (1-based):
    the target alt becomes 1, every other allele 0, missing stays `.`. Returns
    None for a fully-missing call. Used when `--reference` splits multi-allelics."""
    if not genotype:
        return None
    *alleles, phased = genotype
    if all(a < 0 for a in alleles):
        return None
    sep = "|" if phased else "/"
    return sep.join("." if a < 0 else ("1" if a == alt_index else "0") for a in alleles)
