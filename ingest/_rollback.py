"""The failure path shared by the serial, parallel and Parquet loaders."""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)


def handle_ingest_failure(commit_started: bool, ingest_id: str, tag: str) -> None:
    """Undo a failed ingest's engine writes, if any were made.

    Rolls back only once the commit phase began: a failure while parsing
    wrote nothing to the engine, and a rollback then would delete the
    prior data under the same ingest_id, turning a failed re-ingest into
    a wipe. A rollback that itself fails is logged and swallowed so the
    caller re-raises the original error. `tag` names the loader in the
    log ("ingest", "parallel-ingest", "parquet-ingest").
    """
    import storage

    if commit_started:
        log.warning("[%s] failed mid-commit — rolling back %s", tag, ingest_id)
        try:
            storage.rollback_ingest(ingest_id)
        except Exception as rb_err:
            log.error("[%s] rollback FAILED: %s", tag, rb_err)
    else:
        log.warning(
            "[%s] failed during parse — no engine writes occurred, "
            "prior data under ingest_id=%s preserved",
            tag,
            ingest_id,
        )
