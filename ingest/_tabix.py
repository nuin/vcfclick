"""Minimal tabix (.tbi) index reader for density-aware region splitting.

The tabix linear index stores, for every 16Kb position bucket on each
contig, the virtual file offset of the first record overlapping that
bucket. Walking the differences between successive offsets gives a
direct proxy for compressed-byte density per 16Kb bucket, which is in
turn a direct proxy for variant density — exactly what the parallel
splitter needs to balance worker load.

This lets us skip the ~28s cyvcf2 pre-pass and read density from a
~40KB index file in milliseconds.

Format reference:
  https://samtools.github.io/hts-specs/tabix.pdf  (section 4)
"""

from __future__ import annotations

import gzip
import struct
from pathlib import Path


TBI_MAGIC = b"TBI\x01"


def _read(buf: bytes, pos: int, fmt: str) -> tuple:
    """struct.unpack a single value from buf at pos. Returns (value, new_pos)."""
    size = struct.calcsize(fmt)
    return struct.unpack_from(fmt, buf, pos)[0], pos + size


def parse_tbi(tbi_path: Path) -> dict[str, list[int]]:
    """Read a tabix index. Returns {contig_name: linear_index_offsets}.

    `linear_index_offsets` is a list of uint64 virtual file offsets, one
    per 16Kb position bucket on that contig. A virtual offset is
    `(bgzf_block_offset << 16) | uncompressed_offset_within_block`.

    Empty buckets share the offset of the next non-empty bucket
    (htslib convention) — so `offsets[i+1] - offsets[i]` reflects the
    actual compressed-byte cost of bucket i.
    """
    with gzip.open(tbi_path, "rb") as f:
        buf = f.read()

    if buf[:4] != TBI_MAGIC:
        raise ValueError(f"{tbi_path}: bad magic {buf[:4]!r}, expected {TBI_MAGIC!r}")

    pos = 4
    (n_ref, pos) = _read(buf, pos, "<i")
    # 7 header ints we don't need for density computation: format,
    # col_seq, col_beg, col_end, meta, skip, l_nm.
    for _ in range(6):
        (_, pos) = _read(buf, pos, "<i")
    (l_nm, pos) = _read(buf, pos, "<i")

    names_blob = buf[pos : pos + l_nm]
    pos += l_nm
    # Sequence names: concatenated null-terminated strings. Filter empty
    # trailing slot from the final NUL.
    names = [n.decode() for n in names_blob.split(b"\x00") if n]

    contigs: dict[str, list[int]] = {}
    for ref_idx in range(n_ref):
        name = names[ref_idx] if ref_idx < len(names) else f"_unknown_{ref_idx}"

        # Skip the binning index — we only need the linear index for
        # density. (The binning index supports point queries; the linear
        # index is the per-16Kb-bucket structure.)
        (n_bin, pos) = _read(buf, pos, "<i")
        for _ in range(n_bin):
            (_, pos) = _read(buf, pos, "<I")  # bin id
            (n_chunk, pos) = _read(buf, pos, "<i")
            pos += n_chunk * 16  # each chunk is two uint64s

        (n_intv, pos) = _read(buf, pos, "<i")
        offsets = list(struct.unpack_from(f"<{n_intv}Q", buf, pos))
        pos += n_intv * 8

        contigs[name] = offsets

    return contigs


def variant_density(
    linear_offsets: list[int],
    position_bucket_size: int = 100_000,
) -> dict[int, int]:
    """Roll up 16Kb-resolution offsets to `position_bucket_size` buckets.

    Returns {bucket_index: byte_cost} for non-empty buckets only —
    same shape as the cyvcf2-derived `bucket_counts` the splitter
    consumes, just with bytes as the count proxy instead of records.
    The splitter only needs RELATIVE counts to balance, so the units
    don't matter.
    """
    LINEAR_BUCKET_BP = 16384  # tabix convention

    # Per-16Kb byte deltas. The upper 48 bits of the virtual offset
    # are the BGZF block offset; intra-block uoffset noise washes out
    # at the position-bucket scale.
    deltas = []
    for i in range(len(linear_offsets) - 1):
        a = linear_offsets[i] >> 16
        b = linear_offsets[i + 1] >> 16
        deltas.append(max(b - a, 0))
    if linear_offsets:
        # Final 16Kb bucket has no successor to compute byte_cost
        # against, but it is in the index because at least one variant
        # falls inside it (tabix never writes trailing empty buckets).
        # Use a placeholder cost of 1 so it appears in the rolled-up
        # density map — otherwise the splitter's last region cuts off
        # before this bucket and any variants in it are silently
        # dropped from the parallel ingest.
        deltas.append(1)

    counts: dict[int, int] = {}
    for linear_idx, byte_cost in enumerate(deltas):
        if byte_cost == 0:
            continue
        pos_bucket = (linear_idx * LINEAR_BUCKET_BP) // position_bucket_size
        counts[pos_bucket] = counts.get(pos_bucket, 0) + byte_cost
    return counts


def split_via_tbi(
    vcf_path: str | Path,
    n_workers: int,
    bucket_size: int = 100_000,
) -> list[str] | None:
    """Read the .tbi sibling index and produce balanced regions.

    Returns None if the index doesn't exist; the caller should fall
    back to the cyvcf2 pre-pass in that case.
    """
    tbi_path = Path(str(vcf_path) + ".tbi")
    if not tbi_path.exists():
        return None

    from ingest.parallel_split import _split_contig_balanced

    contig_offsets = parse_tbi(tbi_path)
    regions: list[str] = []
    for contig, offsets in contig_offsets.items():
        if not offsets:
            continue
        density = variant_density(offsets, bucket_size)
        if not density:
            continue
        regions.extend(_split_contig_balanced(contig, density, n_workers, bucket_size))
    return regions
