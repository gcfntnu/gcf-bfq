'''
This file contains functions required to actually convert the bcl files to fastq
'''
import multiprocessing as mp
import subprocess
import os
import sys
import shutil
import glob
import syslog
import csv
import codecs
import tempfile
import xml.etree.ElementTree as ET
import re
from distutils.dir_util import copy_tree
from reportlab.lib import colors, utils
from reportlab.platypus import BaseDocTemplate, Table, Preformatted, Paragraph, Spacer, Image, Frame, NextPageTemplate, PageTemplate, TableStyle
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
import bcl2fastq_pipeline.afterFastq as bfq_afq
import pandas as pd
import yaml
from  Bio.Seq import Seq


MKFASTQ_10X = {
    '10X Genomics Visium Spatial Gene Expression Slide & Reagents Kit': 'cellranger_spatial_mkfastq',
    '10X Genomics Chromium Next GEM Single Cell ATAC Library & Gel Bead Kit v1.1': 'cellranger_atac_mkfastq',
    '10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3': 'cellranger_mkfastq'
}


def fixNames(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    if "10X Genomics" in config.get("Options","Libprep"):
        return

    names = glob.glob("%s/%s%s/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes))
    names.extend(glob.glob("%s/%s%s/*/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes)))
    for fname in names:
        if "_001.fastq.gz" in fname:
            fnew = fname.replace("_001.fastq.gz", ".fastq.gz")
            fnew = re.sub("_S[0-9]+","",fnew) 
            syslog.syslog("Moving %s to %s\n" % (fname, fnew))
            shutil.move(fname, fnew)


def bcl2fq(config) :
    '''
    takes things from /dont_touch_this/solexa_runs/XXX/Data/Intensities/BaseCalls
    and writes most output into config.outputDir/XXX, where XXX is the run ID.
    '''
    lanes = config.get("Options", "lanes")
    if lanes != '':
        lanes = '_lanes{}'.format(lanes)

    #Make the output directories
    os.makedirs("%s/%s%s" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes), exist_ok=True)
    #Make log directory
    os.makedirs("%s" % (os.path.join(config.get("Paths","logDir"),os.path.dirname(config.get("Options","runID")))), exist_ok=True)
    os.makedirs(os.path.join(config.get("Paths","outputDir"),config.get("Options","runID"),'InterOp'),exist_ok=True)
    copy_tree(
        os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get("Options","runID"),'InterOp'),
        os.path.join(config.get("Paths","outputDir"),config.get("Options","runID"),'InterOp')
        )
    old_wd = os.getcwd()
    os.chdir(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')))

    if "10X Genomics" in config.get("Options","Libprep"):
        #TODO: --interop-dir not supported for cellranger
        cmd = "{cellranger_cmd} --output-dir={output_dir} --sample-sheet={sample_sheet} --run={run_dir} {cellranger_options}".format(
                cellranger_cmd = config.get("cellranger",MKFASTQ_10X[config.get("Options","Libprep")]),
                output_dir = "{}/{}".format(
                    config.get("Paths","outputDir"),
                    config.get("Options","runID")
                    ),
                sample_sheet = config.get("Options","sampleSheet"),
                run_dir = "{}/{}/data/{}".format(
                    config.get("Paths","baseDir"),
                    config.get("Options","sequencer"),
                    config.get("Options","runID")
                    ),
                cellranger_options = config.get("cellranger","cellranger_mkfastq_options")
                )
    else:
        cmd = "%s %s --sample-sheet %s -o %s/%s%s -R %s/%s/data/%s --interop-dir %s/%s/InterOp" % (
            config.get("bcl2fastq","bcl2fastq"),
            config.get("bcl2fastq","bcl2fastq_options"),
            config.get("Options","sampleSheet"),
            config.get("Paths","outputDir"),
            config.get("Options","runID"),
            lanes,
            config.get("Paths","baseDir"),
            config.get("Options","sequencer"),
            config.get("Options","runID"),
            config.get("Paths","outputDir"),
            config.get("Options","runID"),
        )
    try:
        syslog.syslog("[bcl2fq] Running: %s\n" % cmd)
        with open("%s/%s%s.log" % (config.get("Paths","logDir"), config.get("Options","runID"), lanes), "w") as logOut:
            subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True)
    except:
        if "10X Genomics" not in config.get("Options", "Libprep"):
            with open("%s/%s%s.log" % (config.get("Paths","logDir"), config.get("Options","runID"), lanes), "r") as logOut:
                log_content = logOut.read()
            if "<bcl2fastq::layout::BarcodeCollisionError>" in log_content:
                cmd += " --barcode-mismatches 0 "
                with open("%s/%s%s.log" % (config.get("Paths","logDir"), config.get("Options","runID"), lanes), "w") as logOut:
                    syslog.syslog("[bcl2fq] Retrying with --barcode-mismatches 0 : %s\n" % cmd)
                    subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True)
    #retrieve mkfastq logs
    if "10X Genomics" in config.get("Options","Libprep"):
        project_dirs = bfq_afq.get_project_dirs(config)
        pnames = bfq_afq.get_project_names(project_dirs)
        for p in pnames:
            os.makedirs(os.path.join(config.get("Paths", "outputDir"), config.get("Options","runID"), "QC_{}".format(p), "10X_mkfastq"), exist_ok=True)
            cmd = "cp {} {}".format(
                os.path.join(config.get("Paths", "outputDir"), config.get("Options","runID"), config.get("Options","runID").split("_")[-1][1:], "outs", "qc_summary.json"),
                os.path.join(config.get("Paths", "outputDir"), config.get("Options","runID"), "QC_{}".format(p), "10X_mkfastq", "qc_summary.json")

            )
            subprocess.check_call(cmd, shell=True)

    logOut.close()
    os.chdir(old_wd)


