"""Build table rows from cyvcf2 records."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ingest._arrow import GENOTYPES_COLUMNS, VARIANTS_COLUMNS
from ingest.routing import FORMAT_PAIR, FORMAT_SCALAR, FORMAT_TRIPLE
from ingest.routing import INFO_FLAG, INFO_PAIR, INFO_SCALAR


def _stringify(value) -> str:
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    return str(value)


def _format_arr(variant, field):
    """Read a FORMAT field as a per-sample array via cyvcf2, or None."""
    try:
        arr = variant.format(field)
    except KeyError:
        return None
    return arr


def _cell(arr, i):
    """Extract sample `i`'s value from a FORMAT array, unwrapping (n,1) shapes."""
    v = arr[i]
    if hasattr(v, "__len__") and not isinstance(v, (str, bytes)) and len(v) == 1:
        v = v[0]
    return v


def _scalar_int(arr, i):
    if arr is None:
        return None
    try:
        v = int(_cell(arr, i))
    except (IndexError, TypeError, ValueError):
        return None
    return v if v >= 0 else None


def _pair_int(arr, i):
    if arr is None:
        return (None, None)
    try:
        v = arr[i]
        if len(v) < 2:
            return (None, None)
        a, b = int(v[0]), int(v[1])
    except (IndexError, TypeError, ValueError):
        return (None, None)
    return (a if a >= 0 else None, b if b >= 0 else None)


def _triple_int(arr, i):
    if arr is None:
        return (None, None, None)
    try:
        v = arr[i]
        if len(v) < 3:
            return (None, None, None)
        a, b, c = int(v[0]), int(v[1]), int(v[2])
    except (IndexError, TypeError, ValueError):
        return (None, None, None)
    return (
        a if a >= 0 else None,
        b if b >= 0 else None,
        c if c >= 0 else None,
    )


def _triple_float(arr, i):
    if arr is None:
        return (None, None, None)
    try:
        v = arr[i]
        if len(v) < 3:
            return (None, None, None)
        a, b, c = float(v[0]), float(v[1]), float(v[2])
    except (IndexError, TypeError, ValueError):
        return (None, None, None)
    return (
        None if a != a else a,
        None if b != b else b,
        None if c != c else c,
    )


@dataclass(frozen=True)
class CalledCounts:
    """Per-site called-genotype accounting for one decomposed record.

    `n_called`  samples whose GT is fully called (no missing allele),
    `an_called` called alleles (a haploid call counts 1, a diploid 2, a
                partially missing ./1 counts its one called allele),
    `ac_called` ALT alleles among the called ones (haploid `1` counts 1),
    `missing`   indices of samples whose GT is fully missing (./. or .).

    The three counts are None when the record carries no GT to count
    (sites-only VCF, or GT absent from FORMAT): "not recorded", not 0.
    """

    n_called: int | None
    an_called: int | None
    ac_called: int | None
    missing: tuple[int, ...] = ()


NOT_RECORDED = CalledCounts(None, None, None, ())

# A staging batch is flushed after this many genotype + missing-call rows
# as well as after the per-batch site count. On a large cohort one site can
# carry tens of thousands of rows (missing calls, or 0/0 under
# --keep-reference), so the site count alone does not bound memory.
MAX_BATCH_ROWS = 2_000_000


def called_counts(variant, n_samples: int) -> CalledCounts:
    """Count called samples/alleles from cyvcf2's per-sample allele array.

    Reads `variant.genotype.array()` (one row per sample: allele indices,
    then the phase flag; -1 = missing allele, -2 = padding for a lower
    ploidy) rather than `gt_types`, which collapses ploidy and partial
    missingness. Single-ALT records only, so any allele index > 0 is ALT.
    """
    if n_samples == 0:
        return NOT_RECORDED
    try:
        genotype = variant.genotype
    except Exception:  # cyvcf2 raises when GT is absent from the record
        return NOT_RECORDED
    if genotype is None:
        return NOT_RECORDED
    alleles = genotype.array()[:, :-1]
    called = alleles >= 0
    any_called = called.any(axis=1)
    any_missing = (alleles == -1).any(axis=1)
    return CalledCounts(
        n_called=int(np.count_nonzero(any_called & ~any_missing)),
        an_called=int(np.count_nonzero(called)),
        ac_called=int(np.count_nonzero(alleles > 0)),
        missing=tuple(int(i) for i in np.flatnonzero(~any_called)),
    )


def build_missing_rows(
    variant, samples: list[str], counts: CalledCounts, ingest_id: str
) -> list[list]:
    """`missing_genotypes` rows: one per sample whose GT is fully missing."""
    return [
        [ingest_id, variant.CHROM, variant.POS, variant.REF, variant.ALT[0], samples[i]]
        for i in counts.missing
    ]


def build_variant_row(
    variant, ingest_id: str, counts: CalledCounts = NOT_RECORDED
) -> list:
    """One row for the variants table. Caller has verified len(ALT) == 1."""
    info = dict(variant.INFO)

    row: dict = {col: None for col in VARIANTS_COLUMNS}
    row["ingest_id"] = ingest_id
    row["chrom"] = variant.CHROM
    row["pos"] = variant.POS
    row["ref"] = variant.REF
    row["alt"] = variant.ALT[0]
    row["vcf_id"] = variant.ID
    row["qual"] = variant.QUAL
    row["filter"] = variant.FILTER

    extra: dict[str, str] = {}
    for field, value in info.items():
        if field in INFO_SCALAR:
            row[INFO_SCALAR[field]] = value
        elif field in INFO_PAIR:
            ref_col, alt_col = INFO_PAIR[field]
            if isinstance(value, (list, tuple)) and len(value) >= 2:
                row[ref_col] = value[0]
                row[alt_col] = value[1]
        elif field in INFO_FLAG:
            row[INFO_FLAG[field]] = 1 if value else 0
        else:
            if value is not None:
                extra[field] = _stringify(value)

    for flag_col in INFO_FLAG.values():
        if row[flag_col] is None:
            row[flag_col] = 0

    row["info_extra"] = extra
    row["n_called"] = counts.n_called
    row["an_called"] = counts.an_called
    row["ac_called"] = counts.ac_called
    return [row[c] for c in VARIANTS_COLUMNS]


# cyvcf2 gt_types → stored `gt`. Default (sparse) keeps only
# non-reference calls: 1=HET→1, 3=HOM_ALT→2. HOM_REF (0) and the
# missing/no-call UNKNOWN (2) are dropped (absence is the signal).
GT_ENCODE = {1: 1, 3: 2}

# keep_reference also stores confident HOM_REF as gt=0, so trio
# analysis can prove a parent is genuinely 0/0 at a site (vs simply
# absent = no-call). The missing/no-call UNKNOWN (gt_type 2) is STILL
# dropped — a ./. asserts nothing, so it must stay distinguishable
# from a stored 0/0 by its absence.
GT_ENCODE_KEEP_REF = {0: 0, 1: 1, 3: 2}


def build_genotype_rows(
    variant,
    samples: list[str],
    extra_format_fields: list[str],
    ingest_id: str,
    keep_reference: bool = False,
) -> list[list]:
    encode = GT_ENCODE_KEEP_REF if keep_reference else GT_ENCODE
    gt_types = variant.gt_types
    gt_phases = variant.gt_phases
    scalar_arrs = {f: _format_arr(variant, f) for f in FORMAT_SCALAR}
    pair_arrs = {f: _format_arr(variant, f) for f in FORMAT_PAIR}
    triple_arrs = {f: _format_arr(variant, f) for f in FORMAT_TRIPLE}

    extra_arrays: dict[str, object] = {}
    for f in extra_format_fields:
        arr = _format_arr(variant, f)
        if arr is not None:
            extra_arrays[f] = arr

    rows = []
    for i, sample_id in enumerate(samples):
        encoded = encode.get(int(gt_types[i]))
        if encoded is None:
            continue

        row: dict = {col: None for col in GENOTYPES_COLUMNS}
        row["ingest_id"] = ingest_id
        row["chrom"] = variant.CHROM
        row["pos"] = variant.POS
        row["ref"] = variant.REF
        row["alt"] = variant.ALT[0]
        row["sample_id"] = sample_id
        row["gt"] = encoded
        row["phased"] = 1 if gt_phases[i] else 0

        for src, col in FORMAT_SCALAR.items():
            if col != "ft":
                row[col] = _scalar_int(scalar_arrs[src], i)

        ft_arr = scalar_arrs.get("FT")
        if ft_arr is not None:
            try:
                v = _cell(ft_arr, i)
                if isinstance(v, bytes):
                    v = v.decode()
                s = str(v).strip() if v is not None else ""
                if s and s != ".":
                    row["ft"] = s
            except (IndexError, TypeError):
                pass

        for src, (ref_col, alt_col) in FORMAT_PAIR.items():
            row[ref_col], row[alt_col] = _pair_int(pair_arrs[src], i)

        for src, (c1, c2, c3) in FORMAT_TRIPLE.items():
            if src == "GL":
                row[c1], row[c2], row[c3] = _triple_float(triple_arrs[src], i)
            else:
                row[c1], row[c2], row[c3] = _triple_int(triple_arrs[src], i)

        extra: dict[str, str] = {}
        for f, arr in extra_arrays.items():
            try:
                extra[f] = _stringify(arr[i])
            except (IndexError, TypeError):
                continue
        row["format_extra"] = extra

        rows.append([row[c] for c in GENOTYPES_COLUMNS])

    return rows


def build_record_rows(
    variant,
    samples: list[str],
    extra_format_fields: list[str],
    ingest_id: str,
    keep_reference: bool = False,
    record_missing: bool = True,
) -> tuple[list, list[list], list[list]]:
    """All rows one decomposed record contributes: (variant row,
    genotype rows, missing-genotype rows). Shared by the serial and
    parallel loaders so both stay byte-identical."""
    counts = called_counts(variant, len(samples))
    variant_row = build_variant_row(variant, ingest_id, counts)
    genotype_rows = build_genotype_rows(
        variant, samples, extra_format_fields, ingest_id, keep_reference
    )
    missing_rows = (
        build_missing_rows(variant, samples, counts, ingest_id)
        if record_missing
        else []
    )
    return variant_row, genotype_rows, missing_rows
