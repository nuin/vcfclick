"""Build table rows from cyvcf2 records."""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pyarrow as pa

from ingest._arrow import GENOTYPES_ARROW_SCHEMA, GENOTYPES_COLUMNS, VARIANTS_COLUMNS
from ingest.routing import FORMAT_PAIR, FORMAT_SCALAR, FORMAT_TRIPLE
from ingest.routing import INFO_FLAG, INFO_PAIR, INFO_SCALAR

log = logging.getLogger(__name__)


# htslib sentinels in cyvcf2's int32 FORMAT arrays: a missing value and
# the padding after a vector shorter than the header's Number.
_INT_MISSING = -2147483648
_INT_END_OF_VECTOR = -2147483647
# The float sentinels are NaNs with these payloads (bcf_float_missing /
# bcf_float_vector_end), told apart from a real NaN by their bits.
_FLOAT_MISSING_BITS = 0x7F800001
_FLOAT_END_OF_VECTOR_BITS = 0x7F800002


def _stringify_numpy(arr: np.ndarray) -> str:
    """One FORMAT value from a cyvcf2 array as the VCF would spell it:
    elements joined by ',', '.' for a missing element, end-of-vector
    padding dropped. A value that is entirely missing is '.'."""
    if arr.ndim == 0:
        arr = arr.reshape(1)
    out: list[str] = []
    if arr.dtype.kind == "f":
        bits = arr.astype(np.float32).view(np.uint32)
        for x, b in zip(arr, bits):
            if b == _FLOAT_END_OF_VECTOR_BITS:
                break
            out.append("." if (b == _FLOAT_MISSING_BITS or x != x) else str(x))
    elif arr.dtype.kind in "iu":
        for x in arr:
            if x == _INT_END_OF_VECTOR:
                break
            out.append("." if x == _INT_MISSING else str(x))
    else:
        out = [x.decode() if isinstance(x, bytes) else str(x) for x in arr]
    if not out or all(s == "." for s in out):
        return "."
    return ",".join(out)


def _stringify(value) -> str:
    if isinstance(value, np.ndarray):
        return _stringify_numpy(value)
    if isinstance(value, np.generic):
        return _stringify_numpy(np.asarray(value))
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


# Upper bound of each unsigned genotypes column, from the Arrow schema:
# 65,535 for the UInt16 columns (gq, dp, ad_*, pl_*, mq, pq) and
# 4,294,967,295 for ps. A value above the bound (mitochondrial or
# amplicon depth above 65,535) is clamped to it rather than aborting the
# ingest with an ArrowInvalid, and the clamp is logged once per process.
_COLUMN_MAX: dict[str, int] = {
    f.name: 2**f.type.bit_width - 1
    for f in GENOTYPES_ARROW_SCHEMA
    if pa.types.is_unsigned_integer(f.type)
}
_UINT16_MAX = 65535
_clamp_warned = False


def _bounded(v: int, max_value: int) -> int | None:
    global _clamp_warned
    if v < 0:
        return None
    if v > max_value:
        if not _clamp_warned:
            _clamp_warned = True
            log.warning(
                "[ingest] FORMAT value %d exceeds the column range; "
                "clamping such values to %d",
                v,
                max_value,
            )
        return max_value
    return v


def _scalar_int(arr, i, max_value: int = _UINT16_MAX):
    if arr is None:
        return None
    try:
        v = int(_cell(arr, i))
    except (IndexError, TypeError, ValueError):
        return None
    return _bounded(v, max_value)


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
    return (_bounded(a, _UINT16_MAX), _bounded(b, _UINT16_MAX))


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
        _bounded(a, _UINT16_MAX),
        _bounded(b, _UINT16_MAX),
        _bounded(c, _UINT16_MAX),
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
# 2,000,000 rows measured 1.5 to 1.75 GB per process (each of the default
# 4 workers holds one batch); 250,000 keeps a process near 200 MB.
MAX_BATCH_ROWS = 250_000


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


def require_single_alt(variant, vcf_path: str) -> None:
    """Raise ValueError unless the record has exactly one ALT allele.

    The loaders expect `bcftools norm -m -` output. A record with two or
    more ALTs is told to re-normalise; a record with none (ALT '.', as
    bcftools mpileup without -v or some joint-call exports write) cannot
    be fixed by norm, so it is told which records to drop instead.
    """
    n_alt = len(variant.ALT)
    if n_alt == 1:
        return
    site = f"{variant.CHROM}:{variant.POS}"
    if n_alt == 0:
        raise ValueError(
            f"Record at {site} has no ALT allele (ALT is '.'); vcfclick "
            f"stores variant sites only. Drop such records with: "
            f"bcftools view -e 'ALT=\".\"' {vcf_path} -Oz -o out.vcf.gz"
        )
    raise ValueError(
        f"Multi-allelic site at {site} ({n_alt} ALTs). Re-normalise with: "
        f"bcftools norm -m - {vcf_path} | bgzip > out.vcf.gz"
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
            # A header may declare a routed field with Number=4 (GATK's
            # SB) or Number=. and carry several values; a scalar column
            # cannot hold them, so the tuple goes to info_extra instead.
            if isinstance(value, (list, tuple)):
                if len(value) == 1:
                    row[INFO_SCALAR[field]] = value[0]
                elif value:
                    extra[field] = _stringify(value)
            else:
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
# from a stored 0/0 by its absence. A partial call ('./0', '0/.') is
# also HOM_REF to cyvcf2 but proves nothing about the missing allele, so
# build_genotype_rows stores no gt=0 row for it; './1' stays a het (gt=1).
GT_ENCODE_KEEP_REF = {0: 0, 1: 1, 3: 2}


def _partial_calls(variant):
    """Per-sample flag: True where the GT has a missing allele (-1)."""
    try:
        return (variant.genotype.array()[:, :-1] == -1).any(axis=1)
    except Exception:
        return None


def build_genotype_rows(
    variant,
    samples: list[str],
    extra_format_fields: list[str],
    ingest_id: str,
    keep_reference: bool = False,
) -> list[list]:
    encode = GT_ENCODE_KEEP_REF if keep_reference else GT_ENCODE
    try:
        gt_types = variant.gt_types
        gt_phases = variant.gt_phases
    except Exception:
        # cyvcf2 raises when GT is absent from the record (FORMAT without
        # GT, or a sites-only record among called ones). called_counts
        # reports the same record as NOT_RECORDED; there is nothing to
        # store per sample.
        return []
    partial = _partial_calls(variant) if keep_reference else None
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
        if encoded == 0 and partial is not None and partial[i]:
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
                row[col] = _scalar_int(
                    scalar_arrs[src], i, _COLUMN_MAX.get(col, _UINT16_MAX)
                )

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
