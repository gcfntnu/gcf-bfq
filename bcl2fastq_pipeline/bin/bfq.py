#!/usr/bin/env python3
import datetime
import importlib
import logging
import os
import signal
import sys

from threading import Event

import bcl2fastq_pipeline.afterFastq
import bcl2fastq_pipeline.findFlowCells
import bcl2fastq_pipeline.makeFastq
import bcl2fastq_pipeline.misc
import urllib3

from bcl2fastq_pipeline.config import PipelineConfig

# Disable excess warning messages if we disable SSL checks
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
gotHUP = Event()


def breakSleep(signo, _frame):
    gotHUP.set()


def sleep(cfg):
    gotHUP.wait(timeout=float(cfg.static.system["sleeptime"]) * 60 * 60)
    gotHUP.clear()


def setup_logging(verbosity: int = 1) -> None:
    """Initialize global logging configuration."""
    level = logging.DEBUG if verbosity > 1 else logging.INFO
    fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    logging.basicConfig(level=level, format=fmt, datefmt="%Y-%m-%d %H:%M:%S")


signal.signal(signal.SIGHUP, breakSleep)

verbosity = 2 if os.environ.get("BFQ_DEBUG", None) else 1
setup_logging(verbosity)
log = logging.getLogger("bfq")
log.info("Starting bcl2fastq pipeline")

PipelineConfig.load("/config/bcl2fastq.ini")

while True:
    # Reimport to allow reloading a new version
    importlib.reload(bcl2fastq_pipeline.findFlowCells)
    importlib.reload(bcl2fastq_pipeline.makeFastq)
    importlib.reload(bcl2fastq_pipeline.afterFastq)
    importlib.reload(bcl2fastq_pipeline.misc)

    # Read the config file
    cfg = PipelineConfig.get()
    if not cfg:
        # There's no recovering from this!
        log.error("Unable to read configfile")
        sys.exit(1)

    in_pths = [cfg.static.paths.nova_base_dir, cfg.static.paths.ekista_base_dir]
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
    dirs = list()
    # Get the next flow cell to process, or sleep
    for pth in in_pths:
        for machine, fin_file in completion_files.items():
            dirs += list(pth.glob(f"*_{machine}_*/{fin_file}"))

    for d in sorted(dirs):
        cfg.run.begin(d.parent, cfg.static.paths)
        log.debug(f"Initiate {d}")
        if bcl2fastq_pipeline.findFlowCells.flowCellProcessed():
            log.debug(f"Already processed {d}")
            cfg.run.reset()
            continue

        bcl2fastq_pipeline.findFlowCells.newFlowCell()
        if not cfg.run.run_id:
            continue
        # Ensure we have sufficient space
        if not bcl2fastq_pipeline.misc.enoughFreeSpace():
            log.error("Insufficient free space!")
            bcl2fastq_pipeline.misc.errorEmail(sys.exc_info(), "Insufficient free space!")
            cfg.run.reset()
            break

        startTime = datetime.datetime.now()

        # Make the fastq files, if not already done
        if not (cfg.output_path / "bcl.done").exists():
            try:
                log.info(f"Starting demultiplexing: {cfg.run.run_id}")
                bcl_done = bcl2fastq_pipeline.makeFastq.bcl2fq()
                (cfg.output_path / "bcl.done").write_text("\t".join(bcl_done))
            except Exception as e:
                cfg.run.reset()
                log.exception("Got an error in bcl2fq")
                bcl2fastq_pipeline.misc.errorEmail(sys.exc_info(), f"Got an error in bcl2fq: {e}")
                continue
        else:
            log.info(f"Demultiplexing already done for {cfg.output_path}")

        if not (cfg.output_path / "files.renamed").exists():
            try:
                log.info("Renaming files")
                bcl2fastq_pipeline.makeFastq.rename_fastqs()
                (cfg.output_path / "files.renamed").write_text("")
            except Exception as e:
                cfg.run.reset()
                log.exception("Got an error in rename_fastqs")
                bcl2fastq_pipeline.misc.errorEmail(
                    sys.exc_info(), f"Got an error in rename_fastqs: {e}"
                )
                continue

        # Run post-processing steps
        try:
            log.info("Starting post-processing")
            message = bcl2fastq_pipeline.afterFastq.postMakeSteps()
        except Exception as e:
            cfg.run.reset()
            log.exception("Got an error during postMakeSteps")
            bcl2fastq_pipeline.misc.errorEmail(
                sys.exc_info(), f"Got an error during postMakeSteps: {e}"
            )
            continue

        # Get more statistics and create PDFs
        try:
            message += bcl2fastq_pipeline.misc.getFCmetricsImproved()
        except Exception as e:
            cfg.run.reset()
            log.exception("Got an error during getFCmetrics")
            bcl2fastq_pipeline.misc.errorEmail(
                sys.exc_info(), f"Got an error during getFCmetrics: {e}"
            )
            continue
        endTime = datetime.datetime.now()
        runTime = endTime - startTime

        # Email finished message
        retry_email = None
        try:
            bcl2fastq_pipeline.misc.finishedEmail(message, runTime)
        except Exception as e:
            if cfg.run.libprep.startswith(
                ("10X Genomics Chromium Single Cell", "Parse Biosciences")
            ):
                retry_email = True
                log.info("Got an error during finishedEmail().")
            else:
                cfg.run.reset()
                log.exception("Got an error in finishedEmail")
                bcl2fastq_pipeline.misc.errorEmail(
                    sys.exc_info(), f"Got an error during finishedEmail(): {e}"
                )
                continue

        if retry_email:
            try:
                log.info("Retry without extra html")
                extra_html = False
                bcl2fastq_pipeline.misc.finishedEmail(message, runTime, extra_html)
            except Exception as e:
                cfg.run.reset()
                log.exception("Retry failed. Got an error in finishedEmail")
                bcl2fastq_pipeline.misc.errorEmail(
                    sys.exc_info(), f"Got an error during finishedEmail(): {e}"
                )
                continue

        # Finalize
        try:
            bcl2fastq_pipeline.afterFastq.finalize()
        except Exception as e:
            cfg.run.reset()
            log.exception("Got an error during finalize!")
            bcl2fastq_pipeline.misc.errorEmail(
                sys.exc_info(), f"Got an error during finalize(): {e}"
            )
            continue
        finalizeTime = datetime.datetime.now() - endTime
        runTime += finalizeTime
        try:
            bcl2fastq_pipeline.misc.finalizedEmail("", finalizeTime, runTime)
        except Exception as e:
            cfg.run.reset()
            log.exception("Got an error during finishedEmail")
            bcl2fastq_pipeline.misc.errorEmail(
                sys.exc_info(), f"Got an error during finishedEmail(): {e}"
            )
            continue
        # Mark the flow cell as having been processed
        bcl2fastq_pipeline.findFlowCells.markFinished()
        log.info(f"bfq finished processing for {cfg.output_path}")
        cfg.run.reset()

    # done processing, no more flowcells in queue
    sleep(cfg)
