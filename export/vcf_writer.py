"""BGZF compression and final FORMAT projection for exported VCF files."""

from __future__ import annotations

import struct
import sys
import zlib
from pathlib import Path
from typing import IO

_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


class BgzfWriter:
    """Minimal BGZF writer (the blocked gzip htslib and tabix read)."""

    _MAX_INPUT = 65280  # leaves room for the block header in the 64 KiB limit

    def __init__(self, fh: IO[bytes]):
        self._fh = fh
        self._buf = bytearray()

    def write(self, text: str) -> None:
        self._buf += text.encode()
        while len(self._buf) >= self._MAX_INPUT:
            self._block(bytes(self._buf[: self._MAX_INPUT]))
            del self._buf[: self._MAX_INPUT]

    def _block(self, data: bytes) -> None:
        comp = zlib.compressobj(6, zlib.DEFLATED, -15)
        cdata = comp.compress(data) + comp.flush()
        if len(cdata) > 65536 - 26:  # incompressible: split and retry
            half = len(data) // 2
            self._block(data[:half])
            self._block(data[half:])
            return
        bsize = 18 + len(cdata) + 8 - 1
        header = (
            b"\x1f\x8b\x08\x04"
            + b"\x00\x00\x00\x00"
            + b"\x00\xff"
            + struct.pack("<H", 6)
            + b"BC"
            + struct.pack("<HH", 2, bsize)
        )
        self._fh.write(
            header
            + cdata
            + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
        )

    def close(self) -> None:
        if self._buf:
            self._block(bytes(self._buf))
            self._buf.clear()
        self._fh.write(_BGZF_EOF)
        self._fh.close()


def write_vcf(
    out: str, header: str, tmp, placeholder: str, used_fmt: list[str]
) -> None:
    keep = [i for i, f in enumerate(("GT", "GQ", "DP", "AD", "FT")) if f in used_fmt]
    fmt = ":".join(used_fmt)

    def finish(line: str) -> str:
        if placeholder not in line:
            return line
        fixed, _, rest = line.partition(placeholder)
        out_cols = []
        for cell in rest.rstrip("\n").lstrip("\t").split("\t"):
            if "\x01" in cell:
                parts = cell.split("\x01")
                cell = ":".join(parts[i] for i in keep)
            out_cols.append(cell)
        return fixed + fmt + "\t" + "\t".join(out_cols) + "\n"

    tmp.seek(0)
    if out == "-":
        sys.stdout.write(header)
        for line in tmp:
            sys.stdout.write(finish(line))
        sys.stdout.flush()
        return
    path = Path(out)
    if path.suffix == ".gz":
        w = BgzfWriter(path.open("wb"))
        w.write(header)
        for line in tmp:
            w.write(finish(line))
        w.close()
    else:
        with path.open("w", encoding="utf-8") as f:
            f.write(header)
            for line in tmp:
                f.write(finish(line))
