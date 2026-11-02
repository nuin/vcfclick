"""Minimal tabix (.tbi) index reader for density-aware region splitting.

The tabix linear index stores, for every 16Kb position bucket on each
contig, the virtual file offset of the first record overlapping that
bucket. Walking the differences between successive offsets gives a
direct proxy for compressed-byte density per 16Kb bucket, which is in
turn a direct proxy for variant density — exactly what the parallel
splitter needs to balance worker load.

This lets us skip the ~28s cyvcf2 pre-pass and read density from a
~40KB index file in milliseconds.

The index also carries, per contig, a pseudo-bin (37450) whose second
chunk is (n_mapped, n_unmapped): the number of records htslib indexed
on that contig. `tbi_record_count` sums it, which gives the parallel
ingest a record count that does not depend on the regions it planned.

Format reference:
  https://samtools.github.io/hts-specs/tabix.pdf  (section 4)
"""

from __future__ import annotations

import gzip
import struct
from pathlib import Path


TBI_MAGIC = b"TBI\x01"
LINEAR_BUCKET_BP = 16384  # tabix convention
# htslib's META_BIN for the 5-level, 14-bit-minimum-shift tabix scheme:
# n_bins + 1 with n_bins = ((1 << 18) - 1) / 7 = 37449.
META_BIN = 37450


def _read(buf: bytes, pos: int, fmt: str) -> tuple:
    """struct.unpack a single value from buf at pos. Returns (value, new_pos)."""
    size = struct.calcsize(fmt)
    return struct.unpack_from(fmt, buf, pos)[0], pos + size


def _parse_tbi(tbi_path: Path) -> tuple[dict[str, list[int]], dict[str, int | None]]:
    """Read a tabix index.

    Returns ({contig: linear_index_offsets}, {contig: n_records}). The
    record count is None for a contig whose index has no pseudo-bin
    (an index written by something other than htslib).
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
    records: dict[str, int | None] = {}
    for ref_idx in range(n_ref):
        name = names[ref_idx] if ref_idx < len(names) else f"_unknown_{ref_idx}"

        # Walk the binning index only for the pseudo-bin's record count;
        # the linear index is the per-16Kb-bucket structure density
        # needs.
        n_records: int | None = None
        (n_bin, pos) = _read(buf, pos, "<i")
        for _ in range(n_bin):
            (bin_id, pos) = _read(buf, pos, "<I")
            (n_chunk, pos) = _read(buf, pos, "<i")
            if bin_id == META_BIN and n_chunk == 2:
                # chunk 0: (off_beg, off_end); chunk 1: (n_mapped, n_unmapped)
                n_mapped, n_unmapped = struct.unpack_from("<2Q", buf, pos + 16)
                n_records = n_mapped + n_unmapped
            pos += n_chunk * 16  # each chunk is two uint64s

        (n_intv, pos) = _read(buf, pos, "<i")
        offsets = list(struct.unpack_from(f"<{n_intv}Q", buf, pos))
        pos += n_intv * 8

        contigs[name] = offsets
        records[name] = n_records

    return contigs, records


def parse_tbi(tbi_path: Path) -> dict[str, list[int]]:
    """Read a tabix index. Returns {contig_name: linear_index_offsets}.

    `linear_index_offsets` is a list of uint64 virtual file offsets, one
    per 16Kb position bucket on that contig. A virtual offset is
    `(bgzf_block_offset << 16) | uncompressed_offset_within_block`.

    Empty buckets share the offset of the next non-empty bucket
    (htslib convention) — so `offsets[i+1] - offsets[i]` reflects the
    actual compressed-byte cost of bucket i.
    """
    return _parse_tbi(tbi_path)[0]


def tbi_record_count(vcf_path: str | Path) -> int | None:
    """Number of records the .tbi sibling index says the VCF holds.

    None when there is no index or when any contig lacks the htslib
    pseudo-bin that carries the count. A header-only VCF has an index
    with no contigs and a count of 0.
    """
    tbi_path = Path(str(vcf_path) + ".tbi")
    if not tbi_path.exists():
        return None
    _, records = _parse_tbi(tbi_path)
    if any(n is None for n in records.values()):
        return None
    return sum(records.values())


def variant_density(
    linear_offsets: list[int],
    position_bucket_size: int = 100_000,
) -> dict[int, int]:
    """Roll up 16Kb-resolution offsets to `position_bucket_size` buckets.

    Returns {bucket_index: byte_cost} for every position bucket from 0
    through the one holding the last base of the last linear window —
    same shape as the cyvcf2-derived `bucket_counts` the splitter
    consumes, just with bytes as the count proxy instead of records.
    The splitter only needs RELATIVE counts to balance, so the units
    don't matter.

    Every linear window costs at least 1. Records in one BGZF block
    share a block offset, so a small or sparse file has a byte delta of
    0 for every window; dropping those windows used to leave only the
    trailing one, and the planned region started at the last bucket.
    A window that straddles a position-bucket boundary is credited to
    the bucket of its first base, and the bucket of its last base is
    kept in the map (at cost 0 if nothing else lands there) so the
    planned regions reach the last record.
    """
    if not linear_offsets:
        return {}

    # Per-16Kb byte deltas. The upper 48 bits of the virtual offset
    # are the BGZF block offset; intra-block uoffset noise washes out
    # at the position-bucket scale. The final window has no successor
    # (tabix never writes trailing empty windows, so it holds a record).
    costs = []
    for i in range(len(linear_offsets) - 1):
        a = linear_offsets[i] >> 16
        b = linear_offsets[i + 1] >> 16
        costs.append(max(b - a, 1))
    costs.append(1)

    counts: dict[int, int] = {}
    for linear_idx, byte_cost in enumerate(costs):
        first_bucket = (linear_idx * LINEAR_BUCKET_BP) // position_bucket_size
        last_bucket = ((linear_idx + 1) * LINEAR_BUCKET_BP - 1) // position_bucket_size
        for bucket in range(first_bucket, last_bucket + 1):
            counts.setdefault(bucket, 0)
        counts[first_bucket] += byte_cost
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
