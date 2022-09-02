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
import bcl2fastq_pipeline.afterFastq as bfq_afq
import pandas as pd
import yaml


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

def fix_stats_json(stats_fn):
	stats = list()
	barcodes = False

	with open(stats_fn) as fh:
		for l in fh.readlines():
			if not barcodes:
				stats.append(l)
				if "UnknownBarcodes" in l:
					barcodes = True
			else:
				if not any(s in l for s in ["{", "}", "Barcodes", "Lane"]):
					l = l.replace("\n", ",\n")
				elif l.startswith("  }"):
					l = l.replace("\n", ",\n")
					stats[-2] = stats[-2].replace(",\n", "\n")
				stats.append(l)

	stats[-2] = stats[-2].replace(",\n", "\n")
	stats[-3] = stats[-3].replace(",\n", "\n")
	with open(stats_fn, "w+") as fh:
		fh.writelines(stats)



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
    force_bcl2fastq = os.environ.get("FORCE_BCL2FASTQ", None)

    if "10X Genomics" in config.get("Options","Libprep"):
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
    elif force_bcl2fastq:
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
    else:
        cmd = "bcl-convert --force --bcl-input-directory {in_dir} --output-directory {out_dir} --sample-sheet {samplesheet} --bcl-sampleproject-subdirectories true --no-lane-splitting true --output-legacy-stats true".format(
            in_dir = os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),"data",config.get("Options","runID")),
            out_dir = os.path.join(config.get("Paths","outputDir"), config.get("Options","runID")),
            samplesheet = os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"), "SampleSheet.csv")
            )
    try:
        syslog.syslog("[convert bcl] Running: %s\n" % cmd)
        with open("%s/%s%s.log" % (config.get("Paths","logDir"), config.get("Options","runID"), lanes), "w") as logOut:
            subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True)
    except:
        if "10X Genomics" not in config.get("Options", "Libprep") and force_bcl2fastq:
            with open("%s/%s%s.log" % (config.get("Paths","logDir"), config.get("Options","runID"), lanes), "r") as logOut:
                log_content = logOut.read()
            if "<bcl2fastq::layout::BarcodeCollisionError>" in log_content:
                cmd += " --barcode-mismatches 0 "
                with open("%s/%s%s.log" % (config.get("Paths","logDir"), config.get("Options","runID"), lanes), "w") as logOut:
                    syslog.syslog("[bcl2fq] Retrying with --barcode-mismatches 0 : %s\n" % cmd)
                    subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True)

    if os.path.exists(os.path.join("Reports", "legacy", "Stats")):
        cmd = "ln -sr {} {}".format(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"), "Reports", "legacy", "Stats"), os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"), "Stats"))
        subprocess.check_call(cmd, shell=True)
		fix_stats_json(os.path.join("Stats", "Stats.json"))

    logOut.close()
    os.chdir(old_wd)


