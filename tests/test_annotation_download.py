"""The loader download cache: never a partial file, ClinVar refreshed weekly.

`urllib.request.urlopen` is replaced by a fake; nothing touches the network.
"""

from __future__ import annotations

import io
import os
import time

import pytest


class _Resp(io.BytesIO):
    def __init__(self, body: bytes, length: int | None):
        super().__init__(body)
        self.headers = {} if length is None else {"Content-Length": str(length)}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


@pytest.fixture
def served(monkeypatch):
    """(calls, set_response): every urlopen call is recorded and answered
    with the body and Content-Length set last."""
    import annotations.loaders._download as dl

    calls: list[str] = []
    state = {"body": b"new", "length": 3}

    def urlopen(url, timeout=None):
        assert timeout, "a download must have a timeout"
        calls.append(url)
        return _Resp(state["body"], state["length"])

    monkeypatch.setattr(dl.urllib.request, "urlopen", urlopen)
    return calls, state


def _clinvar_cache(tmp_path, monkeypatch):
    import annotations.loaders.clinvar as cv

    cached = tmp_path / "clinvar.vcf.gz"
    monkeypatch.setattr(cv, "CACHED_VCF", cached)
    return cv, cached


def test_a_short_download_leaves_no_cache(tmp_path, monkeypatch, served):
    calls, state = served
    cv, cached = _clinvar_cache(tmp_path, monkeypatch)
    state.update(body=b"abc", length=10)  # the connection dropped at 3 bytes
    with pytest.raises(OSError, match="incomplete"):
        cv.download_clinvar()
    assert not cached.exists()
    assert list(tmp_path.iterdir()) == []  # no .part left either
    state.update(body=b"abcdefghij", length=10)
    assert cv.download_clinvar().read_bytes() == b"abcdefghij"
    assert len(calls) == 2


def test_a_failed_download_keeps_the_old_cache(tmp_path, monkeypatch, served):
    _, state = served
    cv, cached = _clinvar_cache(tmp_path, monkeypatch)
    cached.write_bytes(b"old")
    old = time.time() - 30 * 86400
    os.utime(cached, (old, old))
    state.update(body=b"x", length=99)
    with pytest.raises(OSError):
        cv.download_clinvar()
    assert cached.read_bytes() == b"old"


def test_clinvar_cache_is_refreshed_when_stale(tmp_path, monkeypatch, served):
    calls, _ = served
    cv, cached = _clinvar_cache(tmp_path, monkeypatch)
    cached.write_bytes(b"old")
    assert cv.download_clinvar().read_bytes() == b"old"  # fresh: reused
    assert calls == []
    stale = time.time() - (cv.CACHE_MAX_AGE_DAYS + 1) * 86400
    os.utime(cached, (stale, stale))
    assert cv.download_clinvar().read_bytes() == b"new"
    assert len(calls) == 1


def test_gencode_cache_is_kept(tmp_path, monkeypatch, served):
    import annotations.loaders.gencode_genes as gg

    calls, _ = served
    cached = tmp_path / "gencode.gff3.gz"
    monkeypatch.setattr(gg, "CACHED_GFF", cached)
    assert gg.download_gencode().read_bytes() == b"new"
    old = time.time() - 400 * 86400
    os.utime(cached, (old, old))
    assert gg.download_gencode().read_bytes() == b"new"
    assert len(calls) == 1
