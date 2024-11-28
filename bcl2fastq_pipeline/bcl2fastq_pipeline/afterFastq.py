'''
This file includes code that actually runs FastQC and any other tools after the fastq files have actually been made. This uses a pool of workers to process each request.
'''
import multiprocessing as mp
import glob
import sys
import subprocess
import os
import os.path
import shlex
import shutil
import xml.etree.ElementTree as ET
import syslog
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as img
import yaml
import json
import re
import datetime as dt
import flowcell_manager.flowcell_manager as fm
import configmaker.configmaker as cm
import pandas as pd
from configmaker.configmaker import SNAKEFILE_TEMPLATE
from configmaker.configmaker import SEQUENCERS

localConfig = None

"""
SEQUENCERS = {
        'NB501038' : 'NextSeq 500',
        'SN7001334' : 'HiSeq 2500',
        'K00251' : 'HiSeq 4000',
        'M02675' : 'MiSeq NTNU',
        'M03942' : 'MiSeq StOlav',
        'M05617' : 'MiSeq SINTEF',
        'A01990' : 'NovaSeq NTNU'
        }
"""
SEQUENCER_OUTPUTFOLDER = {
        'NB501038' : 'nextseq',
        'SN7001334' : 'hiseq2500',
        'K00251' : 'hiseq',
        'M02675' : 'miseq',
        'M03942' : 'miseq',
        'M05617' : 'miseq',
        'M71102' : 'miseq',
        'MN00686' : 'miniseq',
        'A01990' : 'novaseq'
}


def get_gcf_name(fname):
    for name in fname.split('/'):
        if re.match(r'GCF-[0-9]{4}-[0-9]{3,}',name):
            return re.search(r'GCF-[0-9]{4}-[0-9]{3,}',name)[0]
    raise Exception('Unable to determine GCF project number for filename {}\n'.format(fname))

def toDirs(files) :
    s = set()
    for f in files :
        d = os.path.dirname(f)
        if d.split('/')[-1].startswith("GCF-"):
            s.add(d)
        else:
            s.add(d[:d.rfind('/')])
    return s

def get_sequencer(run_id):
    return SEQUENCERS.get(run_id.split('_')[1],'Sequencer could not be automatically determined.')

def get_sequencer_outputfolder(run_id):
    return SEQUENCER_OUTPUTFOLDER.get(run_id.split('_')[1],'Sequencer could not be automatically determined.')

def get_read_geometry(run_dir):
    stats_file = open('{}/Stats/Stats.json'.format(run_dir),'r')
    stats_json = json.load(stats_file)
    lane_info = stats_json['ReadInfosForLanes'][0].get('ReadInfos', None)
    if not lane_info:
        return 'Read geometry could not be automatically determined.'
    R1 = None
    R2 = None
    for read in lane_info:
        if read['IsIndexedRead'] == True:
            continue
        elif read['Number'] == 1:
            R1 = int(read['NumCycles'])
        elif read['Number'] == 2:
            R2 = int(read['NumCycles'])
    if R1 and R2:
        return 'Paired end - forward read length (R1): {}, reverse read length (R2): {}'.format(R1,R2)
    elif R1 and not R2:
        return 'Single end - read length (R1): {}'.format(R1)
    elif not R1 and not R2:
        return 'Read geometry could not be automatically determined.'

def md5sum_worker(config):
    global localConfig
    config = localConfig
    old_wd = os.getcwd()
    os.chdir(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')))
    project_dirs = get_project_dirs(config)
    pnames = get_project_names(project_dirs)
    for p in pnames:
        if os.path.exists('md5sum_{}_fastq.txt'.format(p)):
            continue
        cmd = "find {} -type f -name '*.fastq.gz' | parallel -j 5 md5sum > {}".format(
            p,
            'md5sum_{}_fastq.txt'.format(p)
        )
        syslog.syslog("[md5sum_worker] Processing %s\n" % os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),p))
        subprocess.check_call(cmd, shell=True)
    os.chdir(old_wd)

def md5sum_archive(archive):
    a = archive.replace(".7za","")
    if (not os.path.exists('md5sum_{}_archive.txt'.format(a))) or (os.path.exists('md5sum_{}_archive.txt'.format(a)) and (os.path.getmtime('{}.7za'.format(a)) > os.path.getmtime('md5sum_{}_archive.txt'.format(a)))) :
        cmd = "md5sum {a}.7za > md5sum_{a}_archive.txt".format(a=a)
        subprocess.check_call(cmd, shell=True)

def md5sum_archive_worker(config):
    global localConfig
    config = localConfig
    old_wd = os.getcwd()
    os.chdir(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')))
    archives = glob.glob("*.7za")
    p = mp.Pool()
    p.map(md5sum_archive, archives)
    p.close()
    os.chdir(old_wd)



def multiqc_stats(project_dirs) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats'))

    shutil.copyfile(
        os.path.join(config.get("Paths","baseDir"),config.get('Options','runID'),'RunInfo.xml'),
        os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'RunInfo.xml'),
        )
    #Illumina sequencer update - RunParameters.xml -> runParameters.xml
    if os.path.isfile(os.path.join(config.get("Paths","baseDir"),config.get('Options','runID'),'RunParameters.xml')):
        shutil.copyfile(
            os.path.join(config.get("Paths","baseDir"),config.get('Options','runID'),'RunParameters.xml'),
            os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'RunParameters.xml'),
            )
    else:
        shutil.copyfile(
            os.path.join(config.get("Paths","baseDir"),config.get('Options','runID'),'runParameters.xml'),
            os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'RunParameters.xml'),
            )

    #Illumina interop
    cmd = "interop_summary {} --csv=1 > {}".format(
            os.path.join(config.get("Paths","baseDir"),config.get("Options","runID")),
            os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats','interop_summary.csv'),
        )
    syslog.syslog("[multiqc_worker] Interop summary on %s\n" % os.path.join(config.get("Paths","baseDir"),config.get("Options","runID")))
    subprocess.check_call(cmd,shell=True)
    
    cmd = "interop_index-summary {} --csv=1 > {}".format(
            os.path.join(config.get("Paths","outputDir"),config.get("Options","runID")),
            os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats','interop_index-summary.csv'),
        )
    syslog.syslog("[multiqc_worker] Interop index summary on %s\n" % os.path.join(config.get("Paths","baseDir"),config.get("Options","runID")))
    subprocess.check_call(cmd,shell=True)

    run_dir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'))

    in_confs = glob.glob(os.path.join(run_dir, ".multiqc_config*.yaml"))
    samples_custom_data = dict()
    for c in in_confs:
        with open(c, 'r') as c_fh:
            mqc_conf = yaml.load(c_fh,Loader=yaml.FullLoader)
        samples_custom_data.update(mqc_conf["custom_data"]["general_statistics"]["data"])

    #use one of the existing multiqc_config.yaml as template
    with open(in_confs[0], 'r') as in_conf_fh:
        mqc_conf = yaml.load(in_conf_fh,Loader=yaml.FullLoader)
    pnames = get_project_names(project_dirs)
    pnames = ', '.join(pnames)
    mqc_conf['title'] = pnames
    mqc_conf['intro_text'] = "This report is generated for projects run at Genomics Core Facility, NTNU, Trondheim. The results are reported per sample."
    mqc_conf["custom_data"]["general_statistics"]["data"] = samples_custom_data

    conf_name = os.path.join(run_dir, "Stats", ".multiqc_config.yaml")
    with open(conf_name, "w+") as out_conf_fh:
        yaml.dump(mqc_conf,out_conf_fh)

    modules = "-m interop "
    FORCE_BCL2FASTQ = os.environ.get("FORCE_BCL2FASTQ", None)
    modules += "-m bclconvert " if not FORCE_BCL2FASTQ else "-m bcl2fastq"
    #stats_dir = "Reports" if not FORCE_BCL2FASTQ else "Stats"
    stats_dir = "Stats"

    cmd = "{multiqc_cmd} {multiqc_opts} --config {conf} {flow_dir}/{stats_dir} --filename {flow_dir}/Stats/sequencer_stats_{pname}.html {modules}".format(
            multiqc_cmd = config.get("MultiQC", "multiqc_command"),
            multiqc_opts = config.get("MultiQC", "multiqc_options"),
            conf = conf_name,
            flow_dir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
            stats_dir = stats_dir,
            pname = pnames.replace(", ","_"),
            modules = modules
            )
    syslog.syslog("[multiqc_worker] Processing %s\n" % os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats'))
 
    if os.environ.get("BFQ_TEST", None) and not FORCE_BCL2FASTQ:
        if not os.path.exists(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'), stats_dir, "Demultiplex_Stats.csv")):
            print("BFQ-TEST: Testflowcell was generated with bcl2fastq but environment is configured for bcl-convert. Using bcl2fastq paths and mqc modules.")
            cmd = cmd.replace("Reports", "Stats")
            cmd = cmd.replace("-m bclconvert", "-m bcl2fastq")
            print(cmd)

    subprocess.check_call(cmd, shell=True)

    os.chdir(oldWd)

def archive_worker(config):
    project_dirs = get_project_dirs(config)
    pnames = get_project_names(project_dirs)
    run_date = config.get('Options','runID').split("_")[0]

    for p in pnames:
        """
        Archive fastq
        """
        if os.path.exists(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'{}.7za'.format(p))):
            os.remove(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'{}.7za'.format(p)))
        pw = None
        if config.get("Options","SensitiveData") == "1":
            pw = subprocess.check_output("xkcdpass -n 5 -d '-' -v '[a-z]'",shell=True).decode().strip('\n')
            with open(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),"encryption.{}".format(p)),'w') as pwfile:
                pwfile.write('{}\n'.format(pw))
        opts = "-p{}".format(pw) if pw else ""
        flowdir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'))
        report_dir = os.path.join(flowdir, "Reports") if os.path.exists(os.path.join(flowdir, "Reports")) else ""

        cmd = "7za a {opts} {flowdir}/{pnr}_{date}.7za {flowdir}/{pnr}/ {flowdir}/Stats {report_dir} {flowdir}/Undetermined*.fastq.gz {flowdir}/{pnr}_samplesheet.tsv {flowdir}/SampleSheet.csv {flowdir}/Sample-Submission-Form.xlsx {flowdir}/md5sum_{pnr}_fastq.txt ".format(
                opts = opts,
                date = run_date,
                flowdir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
                report_dir = report_dir,
                pnr = p,
            )
        if "10X Genomics" in config.get("Options","Libprep"):
            cmd += " {}".format(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),config.get("Options","runID").split("_")[-1][1:]))
        syslog.syslog("[archive_worker] Zipping %s\n" % os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'{}_{}.7za'.format(p, run_date)))
        subprocess.check_call(cmd, shell=True)
        """
        Archive pipeline output
        """
        if os.path.exists(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'QC_{}.7za'.format(p))):
            os.remove(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'QC_{}.7za'.format(p)))
        pw = None
        if config.get("Options","SensitiveData") == "1":
            pw = subprocess.check_output("xkcdpass -n 5 -d '-' -v '[a-z]'",shell=True).decode().strip('\n')
            with open(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),"encryption.QC_{}".format(p)),'w') as pwfile:
                pwfile.write('{}\n'.format(pw))

        opts = "-p{}".format(pw) if pw else ""
        flowdir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'))
        pipeline = config.get("Options","pipeline")
        qc_dir = os.path.join(os.environ["TMPDIR"], "{}_{}".format(p, run_date), "data", "tmp", pipeline, "bfq" )

        cmd = f"7z a -l {opts} {flowdir}/QC_{p}_{run_date}.7za {qc_dir} "

        """
        cmd = "7z a -l {opts} {flowdir}/QC_{pnr}_{date}.7za {qc_dir} ".format(
                opts = opts,
                flowdir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
                pnr = p,
                date = run_date,
                qc_dir = qc_dir
            )
        """
        subprocess.check_call(cmd, shell=True)


def get_project_names(dirs):
    gcf = set()
    for d in dirs:
        for catalog in d.split('/'):
            if catalog.startswith('GCF-'):
                gcf.add(catalog)
    return gcf

def get_project_dirs(config):
    projectDirs = glob.glob("%s/%s/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    projectDirs.extend(glob.glob("%s/%s/*/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID"))))
    return toDirs(projectDirs)

def post_workflow(project_id, base_dir, pipeline):
    analysis_dir = os.path.join(os.environ["TMPDIR"], "{}_{}".format(project_id, os.path.basename(base_dir).split("_")[0]))
    os.makedirs(os.path.join(base_dir, "QC_{}".format(project_id), "bfq"), exist_ok=True)
    cmd = "rsync -rvLp {}/ {}".format(
        os.path.join(analysis_dir, "data", "tmp", pipeline, "bfq"),
        os.path.join(base_dir, "QC_{}".format(project_id), "bfq"),
    )
    subprocess.check_call(cmd, shell=True)

    #Copy sample info
    cmd = "cp {} {}".format(
        os.path.join(analysis_dir, "data", "tmp", "sample_info.tsv"),
        os.path.join(base_dir, "{}_samplesheet.tsv".format(project_id))
    )
    subprocess.check_call(cmd, shell=True)

    return True

def full_align(config):
    old_wd = os.getcwd()
    libprep = config.get("Options","Libprep")
    pipeline = config.get("Options","pipeline")

    base_dir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'))
    os.chdir(os.environ["TMPDIR"])
    project_names = get_project_names(get_project_dirs(config))
    run_date = os.path.basename(base_dir).split("_")[0]
    for p in project_names:
        analysis_dir = os.path.join(os.environ["TMPDIR"],"{}_{}".format(p, run_date))
        os.makedirs(analysis_dir, exist_ok=True)
        os.makedirs(os.path.join(analysis_dir,"src"), exist_ok=True)
        os.makedirs(os.path.join(analysis_dir,"data"), exist_ok=True)

        os.chdir(analysis_dir)

        #copy snakemake pipeline
        cmd = "rm -rf {dst} && cp -r {src} {dst}".format(
            src = "/opt/gcf-workflows",
            dst = os.path.join(analysis_dir, "src", "gcf-workflows")
        )
        subprocess.check_call(cmd,shell=True)

        #create config.yaml
        cmd = "/opt/conda/bin/python /opt/conda/bin/configmaker.py {runfolder} -p {project} --libkit '{lib}' --machine '{machine}' {create_fastq}".format(
            runfolder = base_dir,
            project = p,
            lib = libprep,
            samplesheet = os.path.join(base_dir,"SampleSheet.csv"),
            sample_sub = os.path.join(base_dir,"Sample-Submission-Form.xlsx"),
            machine = get_sequencer(config.get("Options","runID")),
            create_fastq = " --skip-create-fastq-dir" if os.path.exists("data/raw/fastq") else ""
        )
        subprocess.check_call(cmd,shell=True)

        #write template Snakefile for pipeline
        #with open(os.path.join(analysis_dir,"Snakefile"), "w") as sn:
        #    sn.write(SNAKEFILE_TEMPLATE.format(workflow=pipeline))

        #run snakemake pipeline
        cmd = "snakemake --use-singularity --singularity-prefix $SINGULARITY_CACHEDIR --cores 32 --verbose -p multiqc_report"
        subprocess.check_call(cmd,shell=True)

        #copy report
        cmd = "cp {src} {dst}".format(
            src = os.path.join("data", "tmp", pipeline, "bfq", "multiqc_{}.html".format(p)),
            dst = os.path.join(base_dir, "multiqc_{}_{}.html".format(p, run_date))
            )
        subprocess.check_call(cmd,shell=True)

        #if additional html reports exists (single cell), copy
        extra_html = glob.glob(os.path.join(analysis_dir, "data", "tmp", pipeline, "bfq", "summaries", "all_samples*.html"))
        if extra_html:
            src = extra_html[0]
            dst = os.path.join(base_dir, f"all_samples_web_summary_{p}_{run_date}.html")
            cmd = f"cp {src} {dst}"
            subprocess.check_call(cmd,shell=True)

        #copy mqc_config
        cmd = "cp {src} {dst}".format(
            src = os.path.join("data", "tmp", pipeline, "bfq", ".multiqc_config.yaml"),
            dst = os.path.join(base_dir, ".multiqc_config_{}.yaml".format(p))
            )
        subprocess.check_call(cmd,shell=True)


        """
        #push workflow bfq output to project directory
        post_workflow(p, base_dir, pipeline)
        """
        """
        #push analysis folder
        analysis_export_dir = os.path.join(config.get("Paths","analysisDir"),"{}_{}".format(p,config.get("Options","runID").split("_")[0]))

        cmd = "rm -rf {dst} && cp -rv {src}/ {dst} ".format(
            src = analysis_dir,
            dst = analysis_export_dir,
        )
        subprocess.check_call(cmd,shell=True)

        #touch bfq_all to avoid rerunning pipelines from scratch
        if not config.get("Options", "libprep").startswith("QIAseq 16S ITS Region Panels"):
            os.chdir(analysis_export_dir)
            #cmd = "find . -type d -exec chmod a+rwx {} \; && find . -type f -exec chmod a+rw {} \; && snakemake --touch -j1 multiqc_report "
            cmd = "snakemake --touch --cores 1 multiqc_report "
            subprocess.check_call(cmd,shell=True)
        """

    os.chdir(old_wd)
    open(os.path.join(config["Paths"]["outputDir"], config["Options"]["runID"],"analysis.made"), "w").close()
    return True


#All steps that should be run after `make` go here
def postMakeSteps(config) :
    '''
    Current steps are:
      1) Run FastQC on each fastq.gz file
      2) Run md5sum on the files in each project directory
    Other steps could easily be added to follow those. Note that this function
    will try to use a pool of threads. The size of the pool is set by config.postMakeThreads
    '''

    projectDirs = get_project_dirs(config)

    global localConfig
    localConfig = config

    with open("/opt/gcf-workflows/libprep.config", "r") as lc_f:
        libprep_config = yaml.load(lc_f,Loader=yaml.FullLoader)

    if config.get("Options", "libprep") + " SE" in libprep_config.keys():
        config["Options"]["pipeline"] = libprep_config[config.get("Options", "libprep") + ' SE']['workflow']
    elif config.get("Options", "libprep") + " PE" in libprep_config.keys():
        config["Options"]["pipeline"] = libprep_config[config.get("Options", "libprep") + ' PE']['workflow']
    else:
        print('failed to identify pipeline from libprep name, using default workflow')
        config["Options"]["pipeline"] = "default"

    #md5sum fastqs
    md5sum_worker(config)

    if not os.path.exists(os.path.join(config["Paths"]["outputDir"], config["Options"]["runID"],"analysis.made")):
        full_align(config)
    
    # multiqc_stats
    multiqc_stats(projectDirs)


    #disk usage
    (tot,used,free) = shutil.disk_usage(config.get("Paths","outputDir"))
    tot /= 1024*1024*1024 #Convert to gigs
    used /= 1024*1024*1024
    free /= 1024*1024*1024

    message = "Current free space for output: %i of %i GiB (%5.2f%%)\n<br>" % (
        free,tot,100*free/tot)

    (tot,used,free) = shutil.disk_usage(config.get("Paths","baseDir"))
    tot /= 1024*1024*1024 #Convert to gigs
    used /= 1024*1024*1024
    free /= 1024*1024*1024


    message += "Current free space for instruments: %i of %i GiB (%5.2f%%)\n<br>\n<br>" % (
        free,tot,100*free/tot)

    #save configfile to flowcell
    with open(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),'bcl2fastq.ini'), 'w+') as configfile:
        config.write(configfile)


    return(message)

def finalize(config):
    #zip arhive
    archive_worker(config)
    #md5sum archive
    md5sum_archive_worker(config)

    return None
