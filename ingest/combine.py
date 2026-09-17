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
import re
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
    pass_only: bool = False,
    count_by: str = "allele",
    carry_info: bool = False,
    reference: str | Path | None = None,
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

    union = _Union(
        gts={},
        inputs={},
        passes={},
        site_inputs={},
        site_passes={},
        meta={},
        contigs={},
        fields=set(),
    )
    all_samples: list[str] = []
    seen_samples: set[str] = set()
    info_hdrs: dict[str, str] = {}
    filter_hdrs: dict[str, str] = {}

    # --reference: split multi-allelics + left-align/trim internally (the
    # `bcftools norm -m - -f` equivalent) instead of refusing. Reuses the
    # benchmark normalizer, which needs pyfaidx.
    ref_fetch = None
    left_align = None
    if reference is not None:
        if not Path(reference).exists():
            raise CombineError(f"reference not found: {reference}")
        try:
            from benchmark.normalize import left_align as _left_align
            from benchmark.reference import Reference
        except ImportError as e:  # pragma: no cover - environment guard
            raise CombineError(
                "--reference requires pyfaidx; install 'vcfclick[benchmark]'."
            ) from e
        try:
            ref_obj = Reference(reference)
        except ImportError as e:
            raise CombineError(
                "--reference requires pyfaidx; install 'vcfclick[benchmark]'."
            ) from e
        ref_fetch, left_align = ref_obj.fetch, _left_align
        # Seed contig order from the reference .fai (§6: inputs may lack
        # ##contig headers; the reference is the authoritative order).
        fai = Path(f"{reference}.fai")
        if fai.exists():
            for line in fai.read_text().splitlines():
                if line.strip():
                    union.contigs.setdefault(line.split("\t")[0], len(union.contigs))

    for idx, path in enumerate(in_paths):
        vcf = VCF(str(path))
        samples = list(vcf.samples)
        if carry_info:
            for _id, ln in _hdr_lines_by_id(vcf.raw_header, "INFO").items():
                info_hdrs.setdefault(_id, ln)
            for _id, ln in _hdr_lines_by_id(vcf.raw_header, "FILTER").items():
                filter_hdrs.setdefault(_id, ln)
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
            alts = list(variant.ALT)
            # One (key, alt_index) per output allele. Without --reference the
            # record must already be biallelic (alt_index None → GT used as-is);
            # with --reference each ALT is split and left-aligned.
            if ref_fetch is None:
                if len(alts) > 1:
                    raise CombineError(
                        f"{path} has a multi-allelic site at "
                        f"{variant.CHROM}:{variant.POS}. Decompose first "
                        f"(bcftools norm -m -), or pass --reference to split "
                        f"internally."
                    )
                alt = alts[0] if alts else "."
                allele_rows = [((variant.CHROM, variant.POS, variant.REF, alt), None)]
            else:
                allele_rows = []
                for j, alt in enumerate(alts, start=1):
                    npos, nref, nalt = left_align(
                        ref_fetch, variant.CHROM, variant.POS, variant.REF, alt
                    )
                    allele_rows.append(((variant.CHROM, npos, nref, nalt), j))

            genotypes = variant.genotypes
            fmt_arrs = {f: _format_arr(variant, f) for f in _PASSTHROUGH}
            for f, arr in fmt_arrs.items():
                if arr is not None:
                    union.fields.add(f)
            record_pass = _is_pass(variant)

            for key, alt_index in allele_rows:
                union.contigs.setdefault(key[0], len(union.contigs))
                spos = (key[0], key[1])
                union.inputs.setdefault(key, set()).add(idx)
                union.site_inputs.setdefault(spos, set()).add(idx)
                if record_pass:
                    union.passes.setdefault(key, set()).add(idx)
                    union.site_passes.setdefault(spos, set()).add(idx)
                if carry_info and key not in union.meta:
                    # First (highest-priority) input to call this allele supplies
                    # its QUAL/FILTER/INFO, taken verbatim from the source record.
                    cols = str(variant).rstrip("\n").split("\t")
                    union.meta[key] = {
                        "qual": cols[5],
                        "filter": cols[6],
                        "info": cols[7],
                    }
                gts = union.gts.setdefault(key, {})
                for s_i, sample in enumerate(samples):
                    if sample in gts:
                        continue  # higher-priority input already filled it
                    if s_i >= len(genotypes):
                        continue
                    g = (
                        _gt_str(genotypes[s_i])
                        if alt_index is None
                        else _remap_gt(genotypes[s_i], alt_index)
                    )
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
    n_kept = _write_combined(
        plain_path,
        all_samples,
        set_names,
        union,
        min_callsets,
        pass_only,
        count_by,
        carry_info,
        info_hdrs,
        filter_hdrs,
    )
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


def _set_field(
    present: set[int],
    passes: set[int],
    n_inputs: int,
    set_names: list[str],
    pass_only: bool,
) -> str:
    """The `set=` value. Default: 'Intersection' if in all inputs, else the
    dash-joined names of the inputs that contain the site (priority order).
    With `pass_only`, a present-but-filtered input is named `filterIn<name>`
    (GATK's convention) and 'Intersection' requires all inputs present AND
    all PASS."""
    if not pass_only:
        if len(present) == n_inputs:
            return "Intersection"
        return "-".join(set_names[i] for i in sorted(present))
    if len(present) == n_inputs and len(passes) == n_inputs:
        return "Intersection"
    return "-".join(
        set_names[i] if i in passes else f"filterIn{set_names[i]}"
        for i in sorted(present)
    )


def _write_combined(
    out_path: Path,
    samples: list[str],
    set_names: list[str],
    union: _Union,
    min_callsets: int,
    pass_only: bool = False,
    count_by: str = "allele",
    carry_info: bool = False,
    info_hdrs: dict | None = None,
    filter_hdrs: dict | None = None,
) -> int:
    n_inputs = len(set_names)
    by_site = count_by == "site"

    def _present(key: tuple) -> set[int]:
        # inputs backing a kept allele: the position's inputs in site mode,
        # the exact-allele inputs otherwise.
        return union.site_inputs[(key[0], key[1])] if by_site else union.inputs[key]

    def _passset(key: tuple) -> set[int]:
        if by_site:
            return union.site_passes.get((key[0], key[1]), set())
        return union.passes.get(key, set())

    def _counting(key: tuple) -> set[int]:
        # inputs that count toward --min-callsets: PASS-only when requested.
        return _passset(key) if pass_only else _present(key)

    ordered_keys = sorted(
        (k for k in union.inputs if len(_counting(k)) >= min_callsets),
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
    if carry_info:
        # Keep the definitions for the carried INFO/FILTER fields (not `set`,
        # already declared above).
        header += [ln for _id, ln in (info_hdrs or {}).items() if _id != "set"]
        header += list((filter_hdrs or {}).values())
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
            present = _present(key)
            passes = _passset(key)
            setval = _set_field(present, passes, n_inputs, set_names, pass_only)
            qual, flt = ".", "."
            info = f"set={setval}"
            if carry_info:
                m = union.meta.get(key, {"qual": ".", "filter": ".", "info": "."})
                qual = m["qual"]
                flt = m["filter"] if m["filter"] not in ("", ".") else "PASS"
                parts = [f"set={setval}"]
                if m["info"] not in ("", "."):
                    parts.append(m["info"])
                info = ";".join(parts)
            gts = union.gts[key]
            cells = [
                _render_cell(gts.get(s), out_fields, missing_cell) for s in samples
            ]
            fh.write(
                "\t".join(
                    [
                        chrom,
                        str(pos),
                        ".",
                        ref,
                        alt,
                        qual,
                        flt,
                        info,
                        format_col,
                        *cells,
                    ]
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
