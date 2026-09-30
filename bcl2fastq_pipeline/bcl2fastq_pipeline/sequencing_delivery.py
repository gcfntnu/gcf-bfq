"""Early QC recovery, independent of workflow execution and SMTP delivery."""

import logging
import time

from pathlib import Path

from bcl2fastq_pipeline import notification_delivery, notifications, sequencing_qc
from bcl2fastq_pipeline.state import StateConflictError

log = logging.getLogger(__name__)


def record_conversion(cfg, store, tool, version, *, run_time=""):
    context = notifications.make_payload(cfg, "sequencing", projects=[], run_time=run_time)
    # Keep the effective sheet from this conversion, even if an operator later
    # refreshes analysis inputs before recovering a failed QC report.
    sheet = cfg.run.sample_sheet or cfg.output_path / "SampleSheet.csv"
    try:
        context["sequencing_sample_sheet"] = sheet.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeError) as error:
        context["sequencing_sample_sheet_error"] = str(error)
    return store.record_demultiplexing_execution(cfg.run.run_id, tool, version, context)


def ensure_report(cfg, store, run_id, *, retry=False, run_time=""):
    """Generate once per successful BCL conversion while holding its execution lease.

    Report-tool failures are saved and do not fail downstream processing. State
    write errors propagate: we cannot safely send without a durable intent.
    """
    state = store.read(run_id)
    qc = state.get("sequencing_qc")
    if not qc or state["status"] == "archived":
        raise StateConflictError(
            "No retained early-QC execution. Legacy runs are not automatically notified; "
            "a new demultiplexing execution creates an early-QC event."
        )
    if qc["status"] == "completed":
        report = Path(qc["result"]["report_path"])
        if report.is_file() and report.stat().st_size:
            return True
    if not retry and qc["status"] != "pending":
        return False
    store.start_sequencing_qc(run_id)
    started = time.monotonic()
    try:
        saved = notifications._saved_config(cfg, qc["context"])
        if "sequencing_sample_sheet" in qc["context"]:
            sheet = saved.output_path / "Stats" / "sequencing_qc_samplesheet.csv"
            sheet.parent.mkdir(parents=True, exist_ok=True)
            sheet.write_text(qc["context"]["sequencing_sample_sheet"], encoding="utf-8")
            saved.run.sample_sheet = sheet
        else:
            # Use statistics alone when capture failed, never a possibly edited
            # analysis sheet. The report will flag planned samples unavailable.
            saved.run.sample_sheet = None
            log.warning(
                "Conversion SampleSheet unavailable: %s",
                qc["context"].get("sequencing_sample_sheet_error", "not captured"),
            )
        result = sequencing_qc.generate(saved, tool=qc["tool"])
    except Exception as error:
        store.finish_sequencing_qc(run_id, error=f"{type(error).__name__}: {error}")
        log.exception(
            "Sequencing QC unavailable for %s; processing continues. Recover with "
            "flowcell-manager retry-sequencing-qc %s",
            run_id,
            run_id,
        )
        return False
    duration = time.monotonic() - started
    payload = dict(qc["context"])
    payload.pop("sequencing_sample_sheet", None)
    payload.update(
        projects=result["projects"],
        sequencing_qc=result,
        run_time=str(run_time or qc["context"].get("run_time") or "Unavailable (report recovery)"),
    )
    store.finish_sequencing_qc(run_id, result=result, duration=duration, payload=payload)
    return True


def report_and_notify(cfg, store, run_id, *, retry=False, run_time=""):
    if ensure_report(cfg, store, run_id, retry=retry, run_time=run_time):
        notification_delivery.deliver_pending(cfg, store, run_id, kind="sequencing")
    return store.read(run_id)["sequencing_qc"]
