#!/usr/bin/env python3
import datetime
import importlib
import signal
import sys
import syslog

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


signal.signal(signal.SIGHUP, breakSleep)

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
        sys.exit("Error: couldn't read the config file!")

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

    for d in dirs:
        cfg.run.begin(d.parent, cfg.static.paths)

        if bcl2fastq_pipeline.findFlowCells.flowCellProcessed():
            cfg.run.reset()
            continue

        bcl2fastq_pipeline.findFlowCells.newFlowCell()
        if not cfg.run.run_id:
            continue
        # Ensure we have sufficient space
        if not bcl2fastq_pipeline.misc.enoughFreeSpace():
            syslog.syslog("Error: insufficient free space!\n")
            bcl2fastq_pipeline.misc.errorEmail(sys.exc_info(), "Error: insufficient free space!")
            cfg.run.reset()
            break

        startTime = datetime.datetime.now()

        # Make the fastq files, if not already done
        if not (cfg.output_path / "bcl.done").exists():
            try:
                bcl_done = bcl2fastq_pipeline.makeFastq.bcl2fq()
                (cfg.output_path / "bcl.done").write_text("\t".join(bcl_done))
            except Exception as e:
                print(e)
                syslog.syslog("Got an error in bcl2fq\n")
                bcl2fastq_pipeline.misc.errorEmail(sys.exc_info(), "Got an error in bcl2fq")
                continue

        if not (cfg.output_path / "files.renamed").exists():
            try:
                bcl2fastq_pipeline.makeFastq.rename_fastqs()
                (cfg.output_path / "files.renamed").write_text("")
            except Exception as e:
                print(e)
                cfg.run.reset()
                syslog.syslog("Got an error in rename_fastqs\n")
                bcl2fastq_pipeline.misc.errorEmail(sys.exc_info(), "Got an error in rename_fastqs")
                continue

        # Run post-processing steps
        try:
            message = bcl2fastq_pipeline.afterFastq.postMakeSteps()
        except Exception as e:
            print(e)
            cfg.run.reset()
            syslog.syslog(f"Got an error during postMakeSteps:\n {e}")
            bcl2fastq_pipeline.misc.errorEmail(sys.exc_info(), "Got an error during postMakeSteps")
            continue

        # Get more statistics and create PDFs
        try:
            message += bcl2fastq_pipeline.misc.getFCmetricsImproved()
        except Exception as e:
            print(e)
            cfg.run.reset()
            syslog.syslog("Got an error during getFCmetrics\n")
            bcl2fastq_pipeline.misc.errorEmail(sys.exc_info(), "Got an error during getFCmetrics")
            continue
        endTime = datetime.datetime.now()
        runTime = endTime - startTime

        # Email finished message
        try:
            bcl2fastq_pipeline.misc.finishedEmail(message, runTime)
        except Exception as e:
            print(e)
            cfg.run.reset()
            syslog.syslog("Couldn't send the finished email! Quiting")
            bcl2fastq_pipeline.misc.errorEmail(
                sys.exc_info(), "Got an error during finishedEmail()"
            )
            continue

        # Finalize
        try:
            bcl2fastq_pipeline.afterFastq.finalize()
        except Exception as e:
            print(e)
            cfg.run.reset()
            syslog.syslog("Got an error during finalize!\n")
            bcl2fastq_pipeline.misc.errorEmail(sys.exc_info(), str(e))
            continue
        finalizeTime = datetime.datetime.now() - endTime
        runTime += finalizeTime
        try:
            bcl2fastq_pipeline.misc.finalizedEmail("", finalizeTime, runTime)
        except Exception as e:
            print(e)
            cfg.run.reset()
            syslog.syslog("Couldn't send the finalize email! Quiting")
            bcl2fastq_pipeline.misc.errorEmail(
                sys.exc_info(), "Got an error during finishedEmail()"
            )
            continue
        # Mark the flow cell as having been processed
        bcl2fastq_pipeline.findFlowCells.markFinished()
        cfg.run.reset()

    # done processing, no more flowcells in queue
    sleep(cfg)
