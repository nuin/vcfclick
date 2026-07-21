"""Skip the entire benchmark test suite when its optional deps are absent.

The `benchmark` extra (pyfaidx, numpy) is optional; without it these tests
can't import their fixtures. Skip the directory cleanly instead of erroring at
collection, so a bare `pytest` (no `--extra benchmark`) still succeeds.
"""

try:
    import numpy  # noqa: F401
    import pyfaidx  # noqa: F401
except ImportError:
    collect_ignore_glob = ["*"]
