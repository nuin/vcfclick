"""Download a reference file into the loader cache, safely.

The file is written to `<name>.part` and renamed into place only once every
byte the server announced has arrived, so an interrupted download never
leaves a partial file where the next run would reuse it. A cache older than
`max_age_days` is fetched again (ClinVar refreshes weekly); `None` keeps a
versioned file (one GENCODE release) for good.
"""

from __future__ import annotations

import logging
import os
import shutil
import time
import urllib.request
from pathlib import Path

log = logging.getLogger(__name__)

TIMEOUT_S = 60  # per socket operation, not for the whole download
_CHUNK = 1 << 20


def is_fresh(path: Path, max_age_days: float | None) -> bool:
    if not path.exists():
        return False
    if max_age_days is None:
        return True
    return time.time() - path.stat().st_mtime < max_age_days * 86400


def fetch(url: str, dest: Path, *, max_age_days: float | None, tag: str) -> Path:
    """`dest`, downloaded from `url` unless a fresh copy is cached."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if is_fresh(dest, max_age_days):
        log.info("[%s] cached: %s", tag, dest)
        return dest
    part = dest.with_name(dest.name + ".part")
    log.info("[%s] downloading %s → %s", tag, url, dest)
    started = time.time()
    try:
        with urllib.request.urlopen(url, timeout=TIMEOUT_S) as resp:
            expected = resp.headers.get("Content-Length")
            with open(part, "wb") as fh:
                shutil.copyfileobj(resp, fh, _CHUNK)
        got = part.stat().st_size
        if expected is not None and got != int(expected):
            raise OSError(
                f"[{tag}] incomplete download of {url}: {got:,} of "
                f"{int(expected):,} bytes"
            )
        os.replace(part, dest)
    finally:
        part.unlink(missing_ok=True)
    log.info(
        "[%s] done (%.0f MB in %.1fs)",
        tag,
        dest.stat().st_size / 1_000_000,
        time.time() - started,
    )
    return dest
