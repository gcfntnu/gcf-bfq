'''
This file contains functions required to actually convert the bcl files to fastq
'''
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

def determineMask(config):
    '''
    If there's already a mask set in the config file then return it.

    Otherwise:
     1. Check for RunInfo.xml
     2. Parse each <Read> child, adding it to a list.
     3. Join the list by commas
     4. If there's no mask then return nothing
    '''
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes2 = []
        for l in lanes.split("_"):
            lanes2.append("s_{}".format(l))
        lanes = "--tiles {}".format(",".join(lanes2))

    mask = config.get("Options", "index_mask")
    bcLens = [int(x) for x in config.get("Options","bcLen").split(",")]
    bcNum = 0
    if mask != "":
        return "--use-bases-mask {} {}".format(mask, lanes)
    elif os.path.isfile("{}/{}/RunInfo.xml".format(config.get("Paths","baseDir"),config.get("Options","runID"))):
        xml = ET.parse("{}/{}/RunInfo.xml".format(config.get("Paths","baseDir"),config.get("Options","runID")))
        root = xml.getroot()[0][3]
        l = []
        for read in root.findall("Read"):
            if read.get("IsIndexedRead") == "N":
                l.append("Y*")
            else:
                nc = int(read.get("NumCycles"))
                if nc > bcLens[bcNum]:
                    if bcLens[bcNum] > 0:
                        l.append("I{}{}".format(bcLens[bcNum], "n" * (nc - bcLens[bcNum])))
                    else:
                        l.append("{}".format("n" * nc))
                else:
                    l.append("I{}".format(bcLens[bcNum]))
                bcNum += 1
        if len(l) > 0:
            return "--use-bases-mask {} {}".format(",".join(l), lanes)
    return lanes


def fixNames(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    if config.get("Options","singleCell") == "1":
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

    if config.get("Options","singleCell") == "1":
        #TODO: --interop-dir not supported for cellranger
        cmd = "{cellranger_cmd} --output-dir={output_dir} --sample-sheet={sample_sheet} --run={run_dir} {cellranger_options}".format(
                cellranger_cmd = config.get("cellranger","cellranger_mkfastq"),
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
    syslog.syslog("[bcl2fq] Running: %s\n" % cmd)
    logOut = open("%s/%s%s.log" % (config.get("Paths","logDir"), config.get("Options","runID"), lanes), "w")
    subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True)
    logOut.close()
    os.chdir(old_wd)

def getOffSpecies(fname) :
    total = 0
    species=[]
    ohol=[]
    mhol=[]
    i = 0
    maxi = 0
    for line in csv.reader(open(fname, "r"), dialect="excel-tab") :
        if(len(line) == 0) :
            break
        if(line[0].startswith("#")) :
            continue
        if(line[0].startswith("Library")) :
            continue
        if(line[0].startswith("PhiX") or line[0].startswith("Adapters") or line[0].startswith("Vectors") or line[0].startswith("rRNA")):
            continue
        species.append(line[0])
        ohol.append(float(line[5]))

        if(ohol[maxi] < ohol[i]) :
            maxi = i
        i += 1

    off = 0
    for i in range(len(ohol)) :
        if(i != maxi) :
            off += ohol[i]
    return off


