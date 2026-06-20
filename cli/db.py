"""Register `vcfclick db ...` subcommands.

`cli.main` imports this module for side effects. The command
implementations live in focused sibling modules so each file stays
small enough for the repository quality gate.
"""

from __future__ import annotations

import importlib

for _module in (
    "cli.db_basic",
    "cli.db_bundle",
    "cli.db_diff",
    "cli.db_batch",
    "cli.db_stats",
    "cli.db_trio",
    "cli.db_qc",
):
    importlib.import_module(_module)


# Compatibility re-exports for tests and any local callers that imported
# these helpers from cli.db before the command modules were split.
from cli.db_batch import _derive_ingest_id  # noqa: E402,F401
from cli.db_diff import _quote_str  # noqa: E402,F401

__all__ = ["_derive_ingest_id", "_quote_str"]
