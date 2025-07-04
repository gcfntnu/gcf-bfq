#!/usr/bin/env python3
import sys
import os
import datetime
import time
import glob
import syslog
import bcl2fastq_pipeline.getConfig
import bcl2fastq_pipeline.findFlowCells
import bcl2fastq_pipeline.makeFastq
import bcl2fastq_pipeline.afterFastq
import bcl2fastq_pipeline.misc
import importlib
import signal
from threading import Event
import urllib3

# Disable excess warning messages if we disable SSL checks
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
gotHUP = Event()

def breakSleep(signo, _frame):
    gotHUP.set()

def sleep(config) :
    gotHUP.wait(timeout=float(config['Options']['sleepTime'])*60*60)
    gotHUP.clear()

signal.signal(signal.SIGHUP, breakSleep)

while True:
    #Reimport to allow reloading a new version
    importlib.reload(bcl2fastq_pipeline.getConfig)
    importlib.reload(bcl2fastq_pipeline.findFlowCells)
    importlib.reload(bcl2fastq_pipeline.makeFastq)
    importlib.reload(bcl2fastq_pipeline.afterFastq)
    importlib.reload(bcl2fastq_pipeline.misc)

    #Read the config file
    config = bcl2fastq_pipeline.getConfig.getConfig()
    if(config is None) :
        #There's no recovering from this!
        sys.exit("Error: couldn't read the config file!")


    in_pths = [config.get("Paths","nova_baseDir"), config.get("Paths","ekista_baseDir")]
    completion_files = {
        'SN7001334': 'ImageAnalysis_Netcopy_complete.txt',
        'NB501038': 'RunCompletionStatus.xml',
        'M026575': 'ImageAnalysis_Netcopy_complete.txt',
        'M026575': 'ImageAnalysis_Netcopy_complete.txt',
        'M03942': 'ImageAnalysis_Netcopy_complete.txt',
        'M05617': 'ImageAnalysis_Netcopy_complete.txt',
        'M71102': 'ImageAnalysis_Netcopy_complete.txt',
        'K00251': 'SequencingComplete.txt',
        'A01990': 'CopyComplete.txt',
        'MN00686': 'CopyComplete.txt',
    }
    dirs = list()
    #Get the next flow cell to process, or sleep
    for pth in in_pths:
        for machine, fin_file in completion_files.items():
            dirs.extend(glob.glob(f"{pth}/*_{machine}_*/{fin_file}"))

    for d in dirs :
        flow_path = os.path.dirname(d)
        config.set('Options','runID', os.path.basename(flow_path))
        config.set('Paths', 'baseDir', os.path.dirname(flow_path))

        if bcl2fastq_pipeline.findFlowCells.flowCellProcessed(config):
            continue

        config = bcl2fastq_pipeline.findFlowCells.newFlowCell(config)
        if(config.get('Options','runID') == ""):
            continue
        #Ensure we have sufficient space
        if(bcl2fastq_pipeline.misc.enoughFreeSpace(config) == False) :
            syslog.syslog("Error: insufficient free space!\n")
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Error: insufficient free space!")
            break

        startTime=datetime.datetime.now()

        #Make the fastq files, if not already done
        if not os.path.exists("{}/{}/bcl.done".format(config["Paths"]["outputDir"], config["Options"]["runID"])):
            try:
                bcl_done = bcl2fastq_pipeline.makeFastq.bcl2fq(config)
                with open("{}/{}/bcl.done".format(config["Paths"]["outputDir"], config["Options"]["runID"]), "w") as fh:
                    fh.write("\t".join(bcl_done))
            except :
                syslog.syslog("Got an error in bcl2fq\n")
                bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error in bcl2fq")
                continue

        if not os.path.exists("{}/{}/files.renamed".format(config["Paths"]["outputDir"], config["Options"]["runID"])):
            try:
                bcl2fastq_pipeline.makeFastq.fixNames(config)
                open("{}/{}/files.renamed".format(config["Paths"]["outputDir"], config["Options"]["runID"]), "w").close()
            except :
                syslog.syslog("Got an error in fixNames\n")
                bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error in fixNames")
                continue


        #Run post-processing steps
        try :
            message = bcl2fastq_pipeline.afterFastq.postMakeSteps(config)
        except Exception as e:
            syslog.syslog("Got an error during postMakeSteps:\n {}".format(e))
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during postMakeSteps")
            continue

        #Get more statistics and create PDFs
        try :
            #message += "\n\n"+bcl2fastq_pipeline.misc.parseConversionStats(config)
            message += bcl2fastq_pipeline.misc.getFCmetricsImproved(config)
        except :
            syslog.syslog("Got an error during parseConversionStats\n")
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during parseConversionStats")
            continue
        endTime = datetime.datetime.now()
        runTime = endTime-startTime

        #Email finished message
        try :
            bcl2fastq_pipeline.misc.finishedEmail(config, message, runTime)
        except :
            #Unrecoverable error
            syslog.syslog("Couldn't send the finished email! Quiting")
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during finishedEmail()")
            continue

        #Finalize
        try:
            bcl2fastq_pipeline.afterFastq.finalize(config)
        except Exception as e:
            syslog.syslog("Got an error during finalize!\n")
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), str(e))
            continue
        finalizeTime = datetime.datetime.now()-endTime
        runTime += finalizeTime
        try:
            bcl2fastq_pipeline.misc.finalizedEmail(config, "", finalizeTime, runTime)
        except:
            #Unrecoverable error
            syslog.syslog("Couldn't send the finalize email! Quiting")
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during finishedEmail()")
            continue
        #Mark the flow cell as having been processed
        bcl2fastq_pipeline.findFlowCells.markFinished(config)

    #done processing, no more flowcells in queue
    sleep(config)

