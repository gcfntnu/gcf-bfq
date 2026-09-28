#!/usr/bin/env python3
import datetime
import importlib
import logging
import os
import signal
import sys

from pathlib import Path
from threading import Event

import urllib3

import bcl2fastq_pipeline.afterFastq
import bcl2fastq_pipeline.findFlowCells
import bcl2fastq_pipeline.makeFastq
import bcl2fastq_pipeline.misc

from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.state import ExecutionLeaseError, FlowcellStateStore

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
    message = bcl2fastq_pipeline.afterFastq.reporting_steps()
    message += bcl2fastq_pipeline.misc.getFCmetricsImproved()
    run_time = datetime.datetime.now() - start_time

    retry_email = False
    try:
        bcl2fastq_pipeline.misc.finishedEmail(message, run_time)
    except Exception:
        if cfg.run.libprep.startswith(("10X Genomics Chromium Single Cell", "Parse Biosciences")):
            retry_email = True
        else:
            raise

    if retry_email:
        logging.getLogger("bfq").info("Retry completion email without extra html")
        bcl2fastq_pipeline.misc.finishedEmail(message, run_time, False)


def _run_state_backed_flowcell(cfg, store, log):
    run_id = cfg.run.run_id
    with store.execution_lease(run_id):
        state = store.begin_attempt(run_id, cfg=cfg)
        first_stage = state["current_stage"]
        start_time = datetime.datetime.now()

        if not bcl2fastq_pipeline.misc.enoughFreeSpace():
            raise RuntimeError("Insufficient free space!")

        stages = ("demultiplexing", "analysis", "reporting", "finalization")
        start_index = stages.index(first_stage)

        for index, stage in enumerate(stages[start_index:], start=start_index):
            current = store.read(run_id)
            if current["stages"][stage]["status"] == "queued":
                store.start_stage(run_id, stage)

            try:
                if stage == "demultiplexing":
                    log.info("Starting demultiplexing: %s", run_id)
                    tool, version = bcl2fastq_pipeline.makeFastq.bcl2fq()
                    bcl2fastq_pipeline.makeFastq.rename_fastqs()
                    store.complete_stage(
                        run_id,
                        stage,
                        {"tool": tool, "version": version},
                    )
                elif stage == "analysis":
                    log.info("Starting analysis: %s", run_id)
                    bcl2fastq_pipeline.afterFastq.analysis_steps()
                    store.complete_stage(run_id, stage)
                elif stage == "reporting":
                    log.info("Starting reporting: %s", run_id)
                    _run_reporting(cfg, start_time)
                    store.complete_stage(run_id, stage)
                elif stage == "finalization":
                    log.info("Starting finalization: %s", run_id)
                    before_finalize = datetime.datetime.now()
                    bcl2fastq_pipeline.afterFastq.finalize()
                    finalize_time = datetime.datetime.now() - before_finalize
                    run_time = datetime.datetime.now() - start_time
                    bcl2fastq_pipeline.misc.finalizedEmail("", finalize_time, run_time)
                    projects = bcl2fastq_pipeline.findFlowCells.markFinished()
                    store.complete_run(run_id, projects)
            except Exception as error:
                report_run_error(
                    cfg,
                    log,
                    f"Got an error during {stage}: {error}",
                    store=store,
                    stage=stage,
                )
                return

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

            bcl2fastq_pipeline.findFlowCells.newFlowCell()
            if not cfg.run.run_id:
                continue
            cfg.run.set_pipeline_from_yaml(
                os.environ.get("BFQ_LIBPREP_CONFIG", "/opt/gcf-workflows/libprep.config")
            )

            try:
                _run_state_backed_flowcell(cfg, store, log)
            except ExecutionLeaseError:
                log.info("Skipping active flowcell %s", cfg.run.run_id)
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
