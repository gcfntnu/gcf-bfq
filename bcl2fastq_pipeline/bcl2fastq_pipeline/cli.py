#!/usr/bin/env python3
import datetime
import importlib
import logging
import os
import signal
import sys
import time

from pathlib import Path
from threading import Event

import urllib3

import bcl2fastq_pipeline.afterFastq
import bcl2fastq_pipeline.findFlowCells
import bcl2fastq_pipeline.makeFastq
import bcl2fastq_pipeline.misc

from bcl2fastq_pipeline import (
    analysis_snapshots,
    notification_delivery,
    notifications,
    preflight,
    workflow_config,
)
from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.state import (
    ExecutionLeaseError,
    FlowcellStateStore,
    StateConflictError,
    output_entries,
    resolve_output_path,
)

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
gotHUP = Event()


def breakSleep(signo, _frame):
    gotHUP.set()


def sleep(cfg):
    gotHUP.wait(timeout=float(cfg.static.system["sleeptime"]) * 60 * 60)
    gotHUP.clear()


def setup_logging(verbosity: int = 1) -> None:
    level = logging.DEBUG if verbosity > 1 else logging.INFO
    fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    logging.basicConfig(level=level, format=fmt, datefmt="%Y-%m-%d %H:%M:%S")


def report_run_error(cfg, log, message, store=None, stage=None):
    """Log, persist, and state-track a run failure before clearing its context."""
    error_info = sys.exc_info()
    log.exception(message)
    report_path = None
    try:
        report_path = bcl2fastq_pipeline.misc.write_error_report(error_info, message)
    except Exception:
        log.exception("Unable to write the flowcell error report")

    signature = bcl2fastq_pipeline.misc.error_failure_signature(stage, error_info, message)
    recorded = False
    if store is not None and stage is not None and cfg.run.run_id:
        try:
            store.fail_stage(
                cfg.run.run_id,
                stage,
                summary=message,
                report_path=report_path,
                failure_signature=signature,
            )
            recorded = True
        except Exception:
            log.exception("Unable to record the flowcell failure in state")
    try:
        if report_path is not None and recorded:
            bcl2fastq_pipeline.misc.send_error_report(
                cfg, report_path, stage, error_info, store, signature
            )
    except Exception:
        log.exception("Unable to deliver the flowcell error email; original failure retained")
    finally:
        cfg.run.reset()
    return report_path


def _run_reporting(cfg, start_time):
    # Email-only metrics and composition happen after this stage is committed.
    return bcl2fastq_pipeline.afterFastq.reporting_steps()


def _prepare_workflow(cfg, store, state):
    """Use completed analysis metadata for downstream-only retries."""
    stage = state["current_stage"]
    workflow = state["stages"]["analysis"]["metadata"].get("workflow")
    if stage in ("reporting", "finalization") and workflow is None:
        projects = state["projects"] or bcl2fastq_pipeline.afterFastq.get_project_names(
            bcl2fastq_pipeline.afterFastq.get_project_dirs(cfg)
        )
        workflow = workflow_config.workflow_from_projects(cfg, projects)

        def record_workflow(current):
            current["stages"]["analysis"]["metadata"]["workflow"] = workflow
            return current

        store.mutate(cfg.run.run_id, record_workflow)
    workflow_config.prepare_execution(cfg, stage, workflow=workflow)


def _run_state_backed_flowcell(cfg, store, log, *, prepare=False):
    run_id = cfg.run.run_id
    with store.execution_lease(run_id):
        analysis_snapshots.recover(store, run_id)
        if prepare:
            current = store.read(run_id)
            if current["status"] != "queued":
                raise StateConflictError(f"Run {run_id} is no longer queued; inspect its state")
            output = resolve_output_path(current["output_path"], run_id, cfg)
            output_entries(output, allow_missing=current["current_stage"] == "demultiplexing")
            if current["output_path"] != str(output):
                store.write({**current, "output_path": str(output)})
            bcl2fastq_pipeline.findFlowCells.newFlowCell()
            if not cfg.run.run_id:
                return
        state = store.begin_attempt(run_id, cfg=cfg)
        first_stage = state["current_stage"]
        start_time = datetime.datetime.now()
        processing_started = time.monotonic()
        notification_seconds = 0.0

        if not bcl2fastq_pipeline.misc.enoughFreeSpace():
            raise RuntimeError("Insufficient free space!")

        stages = ("demultiplexing", "analysis", "reporting", "finalization")
        start_index = stages.index(first_stage)

        for index, stage in enumerate(stages[start_index:], start=start_index):
            current = store.read(run_id)
            if current["stages"][stage]["status"] == "queued":
                store.start_stage(run_id, stage)

            try:
                if stage in {"demultiplexing", "analysis"}:
                    preflight.run_preflight(cfg, store, stage)
                    if prepare and not cfg.run.custom:
                        raise RuntimeError("BFQ SampleSheet is missing usable [CustomOptions]")
                if prepare and index == start_index:
                    _prepare_workflow(cfg, store, current)
                if index == start_index and stage != "demultiplexing":
                    log.info("Checking FASTQ manifests before %s: %s", stage, run_id)
                    bcl2fastq_pipeline.afterFastq.md5sum_worker(cfg)
                if stage == "demultiplexing":
                    log.info("Starting demultiplexing: %s", run_id)
                    tool, version = bcl2fastq_pipeline.makeFastq.bcl2fq()
                    bcl2fastq_pipeline.makeFastq.rename_fastqs()
                    bcl2fastq_pipeline.afterFastq.md5sum_worker(cfg, force=True)
                    store.complete_stage(
                        run_id,
                        stage,
                        {"tool": tool, "version": version},
                    )
                elif stage == "analysis":
                    log.info("Starting analysis: %s", run_id)
                    workdirs = bcl2fastq_pipeline.afterFastq.analysis_steps()
                    store.complete_stage(
                        run_id, stage, {"workflow": cfg.run.pipeline, "workdirs": workdirs or {}}
                    )
                elif stage == "reporting":
                    log.info("Starting reporting: %s", run_id)
                    projects = _run_reporting(cfg, start_time)
                    if projects is None:
                        projects = current["projects"]
                    run_time = datetime.timedelta(
                        seconds=time.monotonic() - processing_started - notification_seconds
                    )
                    payload = notifications.make_payload(
                        cfg, "processed", projects=projects, run_time=str(run_time)
                    )
                    store.complete_stage(
                        run_id, stage, notification={"kind": "processed", "payload": payload}
                    )
                elif stage == "finalization":
                    log.info("Starting finalization: %s", run_id)
                    before_finalize = time.monotonic()
                    bcl2fastq_pipeline.afterFastq.finalize()
                    finalize_time = datetime.timedelta(seconds=time.monotonic() - before_finalize)
                    run_time = datetime.timedelta(
                        seconds=time.monotonic() - processing_started - notification_seconds
                    )
                    projects = bcl2fastq_pipeline.findFlowCells.markFinished()
                    payload = notifications.make_payload(
                        cfg,
                        "finalized",
                        projects=projects,
                        finalize_time=str(finalize_time),
                        run_time=str(run_time),
                    )
                    snapshots = analysis_snapshots.prepare(store.read(run_id), projects)
                    store.complete_run(
                        run_id,
                        projects,
                        notification={"kind": "finalized", "payload": payload},
                        snapshots=snapshots,
                    )
            except Exception as error:
                try:
                    analysis_snapshots.recover(store, run_id)
                except Exception:
                    log.exception("Snapshot recovery deferred for %s", run_id)
                report_run_error(
                    cfg,
                    log,
                    f"Got an error during {stage}: {error}",
                    store=store,
                    stage=stage,
                )
                return

            if stage == "finalization":
                try:
                    analysis_snapshots.recover(store, run_id)
                except Exception:
                    log.exception(
                        "Analysis snapshot publication pending for %s; completed delivery "
                        "and staged snapshot preserved. BFQ will retry publication on its next scan.",
                        run_id,
                    )

            # Delivery errors (including state I/O after SMTP) must never enter
            # the processing failure path above. The intent is already durable.
            if stage in {"reporting", "finalization"}:
                before_notification = time.monotonic()
                notification_delivery.deliver_pending(cfg, store, run_id)
                notification_seconds += time.monotonic() - before_notification

            # complete_run finalizes the final stage itself.
            if index == len(stages) - 1:
                break

        log.info("bfq finished processing for %s", cfg.output_path)
        cfg.run.reset()


def candidate_flowcells(cfg, store):
    """Combine instrument discoveries with explicitly queued state records."""
    completion_files = {
        "SN7001334": "ImageAnalysis_Netcopy_complete.txt",
        "NB501038": "RunCompletionStatus.xml",
        "M026575": "ImageAnalysis_Netcopy_complete.txt",
        "M03942": "ImageAnalysis_Netcopy_complete.txt",
        "M05617": "ImageAnalysis_Netcopy_complete.txt",
        "M71102": "ImageAnalysis_Netcopy_complete.txt",
        "K00251": "SequencingComplete.txt",
        "A01990": "CopyComplete.txt",
        "MN00686": "CopyComplete.txt",
    }
    candidates = {}
    for base in (cfg.static.paths.nova_base_dir, cfg.static.paths.ekista_base_dir):
        for machine, finish_file in completion_files.items():
            for completion in base.glob(f"*_{machine}_*/{finish_file}"):
                candidates[completion.parent.name] = completion.parent

    # Explicitly queued state remains discoverable even if its completion marker
    # is no longer visible. JSON is authoritative, so its source path wins.
    for state in store.list_states():
        if state["status"] in {"queued", "running"}:
            candidates[state["run_id"]] = Path(state["source_path"])

    return [candidates[run_id] for run_id in sorted(candidates)]


def main():
    signal.signal(signal.SIGHUP, breakSleep)

    verbosity = 2 if os.environ.get("BFQ_DEBUG", None) else 1
    setup_logging(verbosity)
    log = logging.getLogger("bfq")
    log.info("Starting bcl2fastq pipeline")

    PipelineConfig.load("/config/bcl2fastq.ini")

    while True:
        importlib.reload(bcl2fastq_pipeline.findFlowCells)
        importlib.reload(bcl2fastq_pipeline.makeFastq)
        importlib.reload(bcl2fastq_pipeline.afterFastq)
        importlib.reload(bcl2fastq_pipeline.misc)

        cfg = PipelineConfig.get()
        if not cfg:
            log.error("Unable to read configfile")
            sys.exit(1)

        store = FlowcellStateStore(cfg.static.paths.manager_dir)
        try:
            store.ensure_writable()
        except Exception:
            log.exception("Flowcell state directory is unavailable; refusing to process")
            sleep(cfg)
            continue

        analysis_snapshots.recover_pending(store)
        notification_delivery.recover_pending(cfg, store)

        for flowcell_path in candidate_flowcells(cfg, store):
            cfg.run.begin(flowcell_path, cfg.static.paths)
            log.debug("Initiate %s", flowcell_path)
            try:
                if bcl2fastq_pipeline.findFlowCells.flowCellProcessed():
                    log.debug("Already processed or not queued: %s", flowcell_path)
                    cfg.run.reset()
                    continue
            except Exception:
                log.exception("Flowcell discovery failed for %s", flowcell_path)
                cfg.run.reset()
                continue

            try:
                _run_state_backed_flowcell(cfg, store, log, prepare=True)
            except ExecutionLeaseError:
                log.info("Skipping active flowcell %s", cfg.run.run_id)
                cfg.run.reset()
            except StateConflictError as error:
                log.error("Cannot prepare %s: %s", cfg.run.run_id, error)
                cfg.run.reset()
            except Exception as error:
                state = store.read(cfg.run.run_id)
                stage = state["current_stage"]
                report_run_error(
                    cfg,
                    log,
                    f"Got an unexpected error during {stage}: {error}",
                    store=store,
                    stage=stage,
                )

        sleep(cfg)


if __name__ == "__main__":
    main()
