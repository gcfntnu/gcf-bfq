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


def demultiplex_16s_its(config):
    old_dir = os.getcwd()
    tmp_dir = os.path.join(os.environ["TMPDIR"],config.get("Options","runID"))
    os.makedirs(tmp_dir, exist_ok=True)
    os.chdir(tmp_dir)
    project_dirs = bfq_afq.get_project_dirs(config)
    pnames = bfq_afq.get_project_names(project_dirs)
    for p in pnames:
        os.makedirs(p, exist_ok=True)
        os.chdir(p)
        os.makedirs("log", exist_ok=True)
        os.makedirs("rev_comp_log", exist_ok=True)
        os.makedirs(os.path.join(config.get("Paths", "outputDir"), config.get("Options","runID"), "QC_{}".format(p), "cutadapt"), exist_ok=True)
        primers = get_16s_its_primers(config)
        r1 = glob.glob(os.path.join(config.get("Paths","outputDir"),config.get("Options","runID"),p,"*R1.fastq.gz"))
        for f in r1:
            cutadapt_worker(config, primers, f)

        """
        p = mp.Pool(4)
        p.map(cutadapt_worker, r1)
        p.close()
        p.join()
        """
        #cutadapt logs in try except
        cmd = "/opt/conda/bin/python /opt/qiaseq/qiaseq_region_summary.py log/*_qiaseq_demultiplex.log > {odir}/qiaseq_regions_mqc.yaml".format(
            odir = os.path.join(config.get("Paths", "outputDir"), config.get("Options","runID"), "QC_{}".format(p), "cutadapt")
        )
        try:
            subprocess.check_call(cmd, shell=True)
        except Exception as e:
            err = "Got an error in qiaseq_region_summary.py: {}\n".format(e)
            syslog.syslog(err)
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), err)

        #keep original files in utputfolder/runid/raw_fastq_pnr
        cmd = "mv {src} {dst}".format(
            src = os.path.join(config.get("Paths", "outputDir"), config.get("Options","runID"), p),
            dst = os.path.join(config.get("Paths", "outputDir"), config.get("Options","runID"), "raw_fastq_{}".format(p)),
        )
        print(cmd)
        subprocess.check_call(cmd, shell=True)

        #move demultiplexed files to utputfolder/runid/pnr
        cmd = "cp -rv {src} {dst} && rm -rf {src}/* ".format(
            dst = os.path.join(config.get("Paths", "outputDir"), config.get("Options","runID"), p),
            src = os.path.join(tmp_dir,p)
        )
        print(cmd)
        subprocess.check_call(cmd, shell=True)

        os.chdir(tmp_dir)

def get_16s_its_primers(config):
    with open("/opt/gcf-workflows/libprep.config","r") as libprepconf:
        libconf = yaml.load(libprepconf, Loader=yaml.Loader)
    print(os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID")))
    conf_libprep = config.get("Options", "Libprep") + (" PE" if bfq_afq.is_paired_end(os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"))) else " SE")
    l_conf = libconf.get(conf_libprep, {})
    """
    print(config.get("Options", "runID"))
    l_conf = libconf[config.get("Options","Libprep")]
    print(config.get("Options", "Libprep"))
    """
    forward_primers = {}
    reverse_primers = {}
    for region, primers in l_conf['db']['primers'].items():
        forward_primers[region] = primers.split('-')[0]
        reverse_primers[region] = primers.split('-')[1]

    return {'forward': forward_primers, 'reverse': reverse_primers}

def cutadapt_worker(config, primers, fname):
    sample = os.path.basename(fname).replace("_R1.fastq.gz","")
    forward = primers.get('forward', {})
    reverse = primers.get('reverse', {})

    if os.path.exists("log/{}_qiaseq_demultiplex.log".format(sample)):
        cmd = "rm log/{}_qiaseq_demultiplex.log".format(sample)
        subprocess.check_call(cmd, shell=True)

    rev_comp_logs = glob.glob("rev_comp_log/{}_*_rev_comp.log".format(sample))
    if rev_comp_logs:
        cmd = "rm {}".format(" ".join(rev_comp_logs))

    regions = ['unknown']
    for region, primer in forward.items():
        if os.path.exists("{}_unknown_R1.fastq".format(sample)):
            fname = "{}_unknown_R1.fastq".format(sample)
            cp_cmd = "mv -f {} input_{}".format(fname, fname)
            subprocess.check_call(cp_cmd, shell=True)
            cp_cmd = "mv -f {} input_{}".format(fname.replace("R1.fastq","R2.fastq"), fname.replace("R1.fastq","R2.fastq"))
            subprocess.check_call(cp_cmd, shell=True)
            sed_cmd = "sed -i -e s/:region=no_adapter//g input_{}".format(fname)
            subprocess.check_call(sed_cmd, shell=True)
            sed_cmd = "sed -i -e s/:region=no_adapter//g input_{}".format(fname.replace("R1.fastq","R2.fastq"))
            subprocess.check_call(sed_cmd, shell=True)

        overlap_fwd = len(primer.replace("N",""))
        overlap_rev = len(reverse[region].replace("N",""))

        cmd = "cutadapt -g \"{region}={fwd_primer};min_overlap={overlap_fwd};max_error_rate=0.25\" -G \"{region}={rev_primer};min_overlap={overlap_rev};max_error_rate=0.3\" --pair-adapters --untrimmed-output {unknown_r1} --untrimmed-paired-output {unknown_r2} --suffix ':region={{name}}' -o {sample}_{{name}}_R1.fastq -p {sample}_{{name}}_R2.fastq --minimum-length 20 {r1} {r2} >> log/{sample}_qiaseq_demultiplex.log".format(
            sample = sample,
            unknown_r1 = "{}_unknown_R1.fastq".format(sample),
            unknown_r2 = "{}_unknown_R2.fastq".format(sample),
            region = region,
            fwd_primer = primer,
            rev_primer = reverse[region],
            overlap_fwd = overlap_fwd,
            overlap_rev = overlap_rev,
            r1 = ("input_" + fname) if "unknown_R1.fastq" in fname else fname,
            r2 = ("input_" + fname.replace("R1.fastq", "R2.fastq")) if "unknown_R1.fastq" in fname else fname.replace("R1.fastq", "R2.fastq"),
            )
        subprocess.check_call(cmd, shell=True)

        #SECOND PASS OF CUtADAPT TO REMOVE REV COMP PRIMERS IN SHORT AMPLICONS
        if os.path.exists(f"{sample}_{region}_R1.fastq"):
            rev_comp_r = str(Seq(reverse[region]).reverse_complement())
            rev_comp_f = str(Seq(primer).reverse_complement())
            cmd = "cutadapt -a {rev_comp_r} -A {rev_comp_f} --pair-adapters -o tmp_{sample}_{region}_R1.fastq -p tmp_{sample}_{region}_R2.fastq --minimum-length 20 {sample}_{region}_R1.fastq {sample}_{region}_R2.fastq > rev_comp_log/{sample}_{region}_rev_comp.log".format(
                sample = sample,
                region = region,
                rev_comp_r = rev_comp_r,
                rev_comp_f = rev_comp_f,
                )
            subprocess.check_call(cmd, shell=True)

            cmd = f"mv tmp_{sample}_{region}_R1.fastq {sample}_{region}_R1.fastq"
            subprocess.check_call(cmd, shell=True)
            cmd = f"mv tmp_{sample}_{region}_R2.fastq {sample}_{region}_R2.fastq"
            subprocess.check_call(cmd, shell=True)


        regions.append(region)

    if glob.glob("input_{}_*fastq".format(sample)):
        rm_cmd = "rm input_{}_*fastq".format(sample)
        subprocess.check_call(rm_cmd, shell=True)

    R1_comb = ["{}_{}_R1.fastq".format(sample, r) for r in regions]
    R1 = [r1 for r1 in R1_comb if os.path.exists(r1)]
    r1 = " ".join(R1)
    r2 = r1.replace("R1.fastq", "R2.fastq")

    #cat and compress
    cmd = "cat {r1} | pigz -6 -p 8 > {sample}_R1.fastq.gz".format(r1 = r1, sample = sample)
    subprocess.check_call(cmd, shell=True)
    cmd = "cat {r2} | pigz -6 -p 8 > {sample}_R2.fastq.gz".format(r2 = r2, sample = sample)
    subprocess.check_call(cmd, shell=True)

    #unlink r1 and r2 from above
    cmd = "rm -f {}".format(r1)
    subprocess.check_call(cmd, shell=True)
    cmd = "rm -f {}".format(r2)
    subprocess.check_call(cmd, shell=True)
    print("finished sample: {}".format(sample))


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


