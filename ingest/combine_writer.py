"""Render combined call sets as VCF records and provenance headers."""

from __future__ import annotations

from pathlib import Path

from ingest.combine_records import _FORMAT_HEADERS, _PASSTHROUGH, _Union


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


def _combined_header(samples, union, out_fields, carry_info, info_hdrs, filter_hdrs):
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

    return header


def _record_metadata(union, key, setval, carry_info):
    qual, flt = ".", "."
    info = f"set={setval}"
    if carry_info:
        m = union.meta.get(key, {"qual": ".", "filter": ".", "info": "."})
        qual = m["qual"] or "."
        flt = m["filter"] or "."  # verbatim: '.' (not applied) stays '.'
        parts = [f"set={setval}"]
        if m["info"] not in ("", "."):
            parts.append(m["info"])
        info = ";".join(parts)
    return qual, flt, info


def _record_line(key, union, setval, samples, out_fields, carry_info):
    chrom, pos, ref, alt = key
    qual, flt, info = _record_metadata(union, key, setval, carry_info)
    missing = ":".join("./." if f == "GT" else "." for f in out_fields)
    gts = union.gts[key]
    cells = [_render_cell(gts.get(s), out_fields, missing) for s in samples]
    return (
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
                ":".join(out_fields),
                *cells,
            ]
        )
        + "\n"
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

    header = _combined_header(
        samples, union, out_fields, carry_info, info_hdrs, filter_hdrs
    )
    n = 0
    with open(out_path, "w") as fh:
        fh.write("\n".join(header) + "\n")
        for key in ordered_keys:
            setval = _set_field(
                _present(key), _passset(key), n_inputs, set_names, pass_only
            )
            fh.write(_record_line(key, union, setval, samples, out_fields, carry_info))
            n += 1
    return n


def _render_cell(cell: dict | None, out_fields: list[str], missing: str) -> str:
    """A sample's FORMAT cell text: its GT plus each output field token
    ('.' where that field was absent in the source). A sample with no
    record at this site is the all-missing cell."""
    if cell is None:
        return missing
    return ":".join(cell.get(f) or "." for f in out_fields)
