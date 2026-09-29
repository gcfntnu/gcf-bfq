"""Delivery orchestration shared by the daemon and explicit operator retries."""

import logging

from bcl2fastq_pipeline import notifications
from bcl2fastq_pipeline.state import ExecutionLeaseError

log = logging.getLogger(__name__)


def deliver_pending(  # noqa: PLR0913
    cfg, store, run_id, *, kind=None, retry=False, retry_uncertain=False
):
    """Attempt eligible intents while the caller holds the run execution lease.

    Contain *notification* errors, including persistence failures after sending.
    Processing failures are handled by the pipeline, outside this function.
    """
    successful = True
    try:
        entries = store.read(run_id).get("delivery_notifications", [])
        selected = [entry for entry in entries if kind is None or entry["kind"] == kind]
        for entry in selected:
            if not retry and entry["status"] in {"failed", "uncertain"}:
                continue
            if entry["status"] in {"sent", "superseded"}:
                continue
            try:
                sent = store.deliver_notification(
                    run_id,
                    entry["id"],
                    lambda current: notifications.send_notification(cfg, current),
                    retry=retry,
                    retry_uncertain=retry_uncertain,
                )
                if sent:
                    log.info("Notification %s sent for %s", entry["id"], run_id)
                else:
                    successful = False
                    current = next(
                        item
                        for item in store.read(run_id)["delivery_notifications"]
                        if item["id"] == entry["id"]
                    )
                    log.warning(
                        "Notification %s for %s: %s (%s). Processing state is unchanged; "
                        "use flowcell-manager retry-notifications %s",
                        entry["id"],
                        run_id,
                        current["status"],
                        current.get("last_error"),
                        run_id,
                    )
            except Exception:
                successful = False
                log.exception(
                    "Notification %s for %s could not be recorded/delivered; "
                    "processing state is unchanged. Inspect notification status before retrying.",
                    entry["id"],
                    run_id,
                )
    except Exception:
        successful = False
        log.exception("Unable to load notifications for %s; processing state is unchanged", run_id)
    return successful


def recover_pending(cfg, store):
    """Recover never-attempted intents, not failed or uncertain deliveries.

    One automatic attempt per intent. Explicit retry is required after failure;
    a previous interrupted claim is classified by the store without resending.
    """
    try:
        states = store.list_states()
    except Exception:
        log.exception("Unable to scan pending notifications")
        return
    for state in states:
        if not any(
            entry["status"] in {"pending", "sending"}
            for entry in state.get("delivery_notifications", [])
        ):
            continue
        try:
            with store.execution_lease(state["run_id"]):
                deliver_pending(cfg, store, state["run_id"])
        except ExecutionLeaseError:
            continue
        except Exception:
            log.exception("Unable to recover notifications for %s", state["run_id"])
