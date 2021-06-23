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

localConfig = None

SEQUENCERS = {
        'NB501038' : 'NextSeq 500',
        'SN7001334' : 'HiSeq 2500',
        'K00251' : 'HiSeq 4000',
        'M02675' : 'MiSeq NTNU',
        'M03942' : 'MiSeq StOlav',
        'M05617' : 'MiSeq SINTEF'
        }

SEQUENCER_OUTPUTFOLDER = {
        'NB501038' : 'nextseq',
        'SN7001334' : 'hiseq2500',
        'K00251' : 'hiseq',
        'M02675' : 'miseq',
        'M03942' : 'miseq',
        'M05617' : 'miseq'
}

QC_PLACEMENT = {
    'External_ID': 0,
    'Sample_Biosource': 10,
    'Sample_Group': 20,
    'Customer_Comment': 30,
    'Fragment_Length': 35,
    '260/230': 40,
    '260/280': 50,
    'Concentration': 50,
    'RIN': 60
}

PIPELINE_MULTIQC_MODULES = {
    'rnaseq': ["fastq_screen","star","picard","fastp","fastqc_rnaseq","custom_content"],
    'microbiome': ["fastq_screen","star","picard","fastp","fastqc_rnaseq","custom_content", "qiime2"],
    'singlecell': ["fastq_screen", "cellranger", "starsolo", "fastp","fastqc_rnaseq","custom_content", "cellranger_count"],
    'smallrna': ["fastq_screen","star","picard","fastp","fastqc_rnaseq", "unitas", "custom_content"],
    'default': ["fastq_screen","fastp","fastqc_rnaseq","custom_content"],
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

def is_paired_end(run_dir):
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
    return R1 and R2

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

def get_sequencer(run_id):
    return SEQUENCERS.get(run_id.split('_')[1],'Sequencer could not be automatically determined.')

def get_sequencer_outputfolder(run_id):
    return SEQUENCER_OUTPUTFOLDER.get(run_id.split('_')[1],'Sequencer could not be automatically determined.')

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

def md5sum_archive_worker(config):
    global localConfig
    config = localConfig
    old_wd = os.getcwd()
    os.chdir(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')))
    project_dirs = get_project_dirs(config)
    pnames = get_project_names(project_dirs)
    for p in pnames:
        if (not os.path.exists('md5sum_{}_archive.txt'.format(p))) or (os.path.exists('md5sum_{}_archive.txt'.format(p)) and (os.path.getmtime('{}.7za'.format(p)) > os.path.getmtime('md5sum_{}_archive.txt'.format(p)))) :
            cmd = "md5sum {p}.7za > md5sum_{p}_archive.txt".format(p=p)
            syslog.syslog("[md5sum_worker] Processing %s\n" % os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')))
            subprocess.check_call(cmd, shell=True)
    os.chdir(old_wd)

def md5sum_instrument_worker(config):
    global localConfig
    config = localConfig
    old_wd = os.getcwd()
    os.chdir(os.path.join(config.get('Paths','archiveInstr'), config.get('Options','runID')))

    if os.path.exists('md5sum_{}.txt'.format(config.get('Options','runID'))):
        os.chdir(old_wd)
        return
    cmd = "md5sum {r}.7za > md5sum_{r}.txt".format(r=config.get('Options','runID'))
    syslog.syslog("[md5sum_worker] Processing %s\n" % os.path.join(config.get('Paths','archiveInstr'), config.get('Options','runID')))
    subprocess.check_call(cmd, shell=True)

    os.chdir(old_wd)

def set_mqc_conf_header(config, mqc_conf, seq_stats=False):
    odir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'))
    read_geometry = get_read_geometry(odir)
    contact = config.get('MultiQC','report_contact')
    sequencer = get_sequencer(config.get('Options','runID'))
    prepkit = config.get('Options','Libprep')
    organism = config.get('Options','Organism')

    if not seq_stats:
        pipeline = config.get("Options","pipeline")
        with open("/opt/mqc_headers/{}_mqc_header.txt".format(pipeline),"r") as header_file:
            header_text = header_file.read()
        mqc_conf['intro_text'] = header_text.format(pname=mqc_conf['title'])
        software = get_software_versions(config)
        format_software = '<br/>'.join(["<strong>Software versions</strong>"] + software.split("\n"))
        mqc_conf['intro_text'] = '<br/><br/>'.join([mqc_conf['intro_text'],format_software])

    report_header = [
    {'Contact E-mail': contact},
    {'Sequencing Platform': sequencer},
    {'Read Geometry': read_geometry},
    ]
    if not organism in ['N/A', 'n/a', '', ' ']:
        report_header.append({'Organism': organism})
    if not prepkit in ['N/A', 'n/a', '', ' ']:
        report_header.append({'Lib prep kit': prepkit})

    mqc_conf['report_header_info'] = report_header

    if read_geometry.startswith('Single end'):
        mqc_conf['extra_fn_clean_exts'].append('_R1')

    if os.path.exists(os.path.join(odir,'{}_samplesheet.tsv'.format(mqc_conf['title']))) and "10X Genomics" not in config.get("Options","Libprep"):
        s_df = pd.read_csv(os.path.join(odir,'{}_samplesheet.tsv'.format(mqc_conf['title'])),sep='\t')
        s_df.index = s_df['Sample_ID']

        max_260_230 = float(s_df['260/230'].max()) if '260/230' in s_df.columns else 3
        max_260_280 = float(s_df['260/280'].max()) if '260/280' in s_df.columns else 3

        MAX = {
            'RIN': 10,
            '260/230': max_260_230,
            '260/280': max_260_280
            }

        COL_SCALE = {
            'RIN': 'OrRd',
            '260/230': 'PuBu',
            '260/280': 'BuGn'
        }

        s_df.drop(['Sample_ID'], axis=1,inplace=True)
        s_df.dropna(how='all', axis=1, inplace=True)
        s_df = s_df.round(2)
        s_dict = s_df.to_dict(orient='index')

        pconfig = {}
        for col in list(s_df.columns.values):
            pconfig[col] = {'format': '{}', 'min': 0, 'placement': QC_PLACEMENT[col]}
            pconfig[col]['max'] = MAX.get(col,None)
            pconfig[col]['scale'] = COL_SCALE.get(col,False)

        data = s_dict

        general_statistics = {
            'plot_type': 'generalstats',
            'pconfig': [pconfig],
            'data': data
        }
        custom_data = {'general_statistics': general_statistics}

        mqc_conf['custom_data'] = custom_data



    return mqc_conf

def get_pipeline_multiqc_modules(config):
    modules = PIPELINE_MULTIQC_MODULES.get(config.get("Options","pipeline"),None)
    return ('-m ' + ' -m '.join(modules)) if modules else '-e fastqc -e qiime2'

def multiqc_worker(d) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)
    #os.chdir("{}/{}".format(config.get('Paths','outputDir'), config.get('Options','runID')))
    dname = d.split("/")
    pname = dname[-1]

    conf_name = "{}/{}/QC_{}/.multiqc_config.yaml".format(config.get('Paths','outputDir'), config.get('Options','runID'),pname)
    pipeline = config.get("Options","pipeline")
    in_conf = open("/config/multiqc_config" + ("-" + pipeline if pipeline else "") + ".yaml","r")
    out_conf = open(conf_name,"w+")
    mqc_conf = yaml.load(in_conf,Loader=yaml.FullLoader)

    mqc_conf['title'] = pname
    mqc_conf = set_mqc_conf_header(config,mqc_conf)

    yaml.dump(mqc_conf,out_conf)
    in_conf.close()
    out_conf.close()

    cmd = "{multiqc_cmd} {multiqc_opts} --config {conf} {flow_dir}/QC_{pname} {flow_dir}/{pname} {modules} --filename {flow_dir}/QC_{pname}/multiqc_{pname}.html".format(
            multiqc_cmd = config.get("MultiQC", "multiqc_command"),
            multiqc_opts = config.get("MultiQC", "multiqc_options"),
            conf = conf_name,
            flow_dir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
            pname=pname,
            modules = get_pipeline_multiqc_modules(config)
            )
    syslog.syslog("[multiqc_worker] Processing %s\n" % d)
    subprocess.check_call(cmd, shell=True)

    os.chdir(oldWd)

def multiqc_stats(project_dirs) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats'))

    shutil.copyfile(
        os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get('Options','runID'),'RunInfo.xml'),
        os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'RunInfo.xml'),
        )
    #Illumina sequencer update - RunParameters.xml -> runParameters.xml
    if os.path.isfile(os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get('Options','runID'),'RunParameters.xml')):
        shutil.copyfile(
            os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get('Options','runID'),'RunParameters.xml'),
            os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'RunParameters.xml'),
            )
    else:
        shutil.copyfile(
            os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get('Options','runID'),'runParameters.xml'),
            os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'RunParameters.xml'),
            )

    #Illumina interop
    cmd = "interop_summary {} --csv=1 > {}".format(
            os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get("Options","runID")),
            os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats','interop_summary.csv'),
        )
    syslog.syslog("[multiqc_worker] Interop summary on %s\n" % os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get("Options","runID")))
    subprocess.check_call(cmd,shell=True)
    
    cmd = "interop_index-summary {} --csv=1 > {}".format(
            os.path.join(config.get("Paths","outputDir"),config.get("Options","runID")),
            os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats','interop_index-summary.csv'),
        )
    syslog.syslog("[multiqc_worker] Interop index summary on %s\n" % os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get("Options","runID")))
    subprocess.check_call(cmd,shell=True)

    conf_name = "{}/{}/Stats/.multiqc_config.yaml".format(config.get('Paths','outputDir'), config.get('Options','runID'))
    in_conf = open("/config/multiqc_config.yaml","r")
    out_conf = open(conf_name,"w+")
    mqc_conf = yaml.load(in_conf,Loader=yaml.FullLoader)

    pnames = get_project_names(project_dirs)
    pnames = ', '.join(pnames)
    mqc_conf['title'] = pnames

    mqc_conf = set_mqc_conf_header(config,mqc_conf,seq_stats=True)

    yaml.dump(mqc_conf,out_conf)
    in_conf.close()
    out_conf.close()

    cmd = "{multiqc_cmd} {multiqc_opts} --config {conf} {flow_dir}/Stats --filename {flow_dir}/Stats/sequencer_stats_{pname}.html -e qiime2".format(
            multiqc_cmd = config.get("MultiQC", "multiqc_command"),
            multiqc_opts = config.get("MultiQC", "multiqc_options"),
            conf = conf_name,
            flow_dir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
            pname = pnames.replace(", ","_")
            )
    syslog.syslog("[multiqc_worker] Processing %s\n" % os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats'))
    subprocess.check_call(cmd, shell=True)
    os.chdir(oldWd)

def archive_worker(config):
    project_dirs = get_project_dirs(config)
    pnames = get_project_names(project_dirs)

    for p in pnames:
        if os.path.exists(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'{}.7za'.format(p))):
            os.remove(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'{}.7za'.format(p)))
        pw = None
        if config.get("Options","SensitiveData") == "1":
            pw = subprocess.check_output("xkcdpass -n 5 -d '-' -v '[a-z]'",shell=True).decode().strip('\n')
            with open(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),"encryption.{}".format(p)),'w') as pwfile:
                pwfile.write('{}\n'.format(pw))
        opts = "-p{}".format(pw) if pw else ""

        cmd = "7za a {opts} {flowdir}/{pnr}.7za {flowdir}/{pnr}/ {flowdir}/QC_{pnr} {flowdir}/Stats {flowdir}/Undetermined*.fastq.gz {flowdir}/{pnr}_samplesheet.tsv {flowdir}/SampleSheet.csv {flowdir}/Sample-Submission-Form.xlsx {flowdir}/md5sum_{pnr}_fastq.txt {flowdir}/software.versions ".format(
                opts = opts,
                flowdir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
                pnr = p,
            )
        if "10X Genomics" in config.get("Options","Libprep"):
            cmd += " {}".format(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),config.get("Options","runID").split("_")[-1][1:]))
        syslog.syslog("[archive_worker] Zipping %s\n" % os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'{}.7za'.format(p)))
        subprocess.check_call(cmd, shell=True)

def instrument_archive_worker(config):
    if not os.path.exists(os.path.join(config.get('Paths','archiveInstr'), config.get('Options','runID'))):
        os.makedirs(os.path.join(config.get('Paths','archiveInstr'), config.get('Options','runID')),exist_ok=True)
    if os.path.exists(os.path.join(config.get('Paths','archiveInstr'), config.get('Options','runID'), '{}.7za'.format(config.get('Options','runID')))):
        return
    pw = None
    if config.get("Options","SensitiveData") == "1":
        pw = subprocess.check_output("xkcdpass -n 5 -d '-' -v '[a-z]'",shell=True).decode().strip('\n')
        with open(os.path.join(config.get('Paths','archiveInstr'), config.get('Options','runID'),"encryption.{}".format(config.get('Options','runID'))),'w') as pwfile:
            pwfile.write('{}\n'.format(pw))
    opts = "-p{}".format(pw) if pw else ""
    seq_out = get_sequencer_outputfolder(config.get('Options','runID'))
    cmd = "7za a {opts} {arch_dir}/{fnm}.7za {instr}".format(
            opts = opts,
            arch_dir = os.path.join(config.get('Paths','archiveInstr'), config.get('Options','runID')),
            fnm = config.get('Options','runID'),
            instr = os.path.join(config.get('Paths','baseDir'), seq_out,"data",config.get('Options','runID'))
        )
    syslog.syslog("[instrument_archive_worker] Zipping instruments." )
    subprocess.check_call(cmd, shell=True)
    

def samplesheet_worker(config,project_dirs):
    """
    TODO:
    Import configmaker
    Get samplesheet data as dict from 'get_data_from_samplesheet(fh)'
    If sample submission form - use 'merge_samples_with_submission_form(ssub,sample_dict)'
    Keep selected columns and write customer samplesheet
    Use data structure to add data to general statistics in multiqc config

    """
    #TODO: Account for multiple projects
    project_names = get_project_names(project_dirs)
    for pid in project_names:
        with open(config.get("Options","sampleSheet"),'r') as ss:
            #sample_df, _ = cm.get_data_from_samplesheet(ss)
            sample_df, _ = cm.get_project_samples_from_samplesheet(ss,[os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'))] , [pid])
        sample_ids = sample_df['Sample_ID']
        #project_dirs = cm.inspect_dirs([os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'))])
        sample_dict = cm.find_samples(sample_df,[os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'), pid)])

        keep_cols = ['Sample_ID']
        try:
            if not config.get("Options","sampleSubForm") == "":
                ssub_d = {config.get("Options","sampleSubForm"): open(config.get("Options","sampleSubForm"),'rb')}
                #TODO: get message from merge (check intersection between sample sheet and sample-sub-form and attach message to email
                sample_dict = cm.merge_samples_with_submission_form(ssub_d, sample_dict)
                keep_cols.extend(['External_ID', 'Sample_Group','Sample_Biosource','Customer_Comment', 'Fragment_Length','RIN', '260/280', '260/230', 'Concentration'])
                try:
                    sample_df = pd.DataFrame.from_dict(sample_dict,orient='index')[keep_cols]
                except Exception as e:
                    #TODO: attach message to error email
                    sample_df = pd.DataFrame.from_dict(sample_dict,orient='index')[['Sample_ID']]
            else:
                sample_df = pd.DataFrame.from_dict(sample_dict,orient='index')[keep_cols]
        except Exception as e:
            #TODO: attach message to error email
            sample_df = sample_ids
        sample_df.to_csv(
            os.path.join(config.get("Paths","outputDir"),config.get("Options","runID"),"{}_samplesheet.tsv".format(pid)),
            index=False,
            sep='\t',
            header=True
            )



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
    projectDirs = [x for x in projectDirs if "raw_fastq_GCF" not in x]
    return toDirs(projectDirs)

def get_software_versions(config):
    versions = {}
    versions['GCF-NTNU bcl2fastq pipeline'] = subprocess.check_output("cat /bfq_version/bfq.version",shell=True).rstrip()
    if "10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit" in config.get("Options","Libprep"):
        versions["cellranger"] = subprocess.check_output("cellranger mkfastq --version",stderr=subprocess.STDOUT,shell=True).split(b'\n')[1].split(b' ')[-1].rstrip().strip(b'(').strip(b')').split(b'-')[-1]
    elif "10X Genomics Chromium Next GEM Single Cell ATAC Library & Gel Bead Kit" in config.get("Options","Libprep"):
        versions["cellranger-atac"] = subprocess.check_output("cellranger-atac mkfastq --version",stderr=subprocess.STDOUT,shell=True).split(b'\n')[0].split(b' ')[-1].rstrip().strip(b'(').strip(b')')
    elif "10X Genomics Visium Spatial Gene Expression Slide & Reagents Kit" in config.get("Options","Libprep"):
        versions["spaceranger"] = subprocess.check_output("spaceranger mkfastq --version",stderr=subprocess.STDOUT,shell=True).split(b'\n')[1].split(b' ')[-1].rstrip().strip(b'(').strip(b')').split(b'-')[-1]
    else:
        versions['bcl2fastq'] = subprocess.check_output("bcl2fastq --version",stderr=subprocess.STDOUT,shell=True).split(b'\n')[1].split(b' ')[1].rstrip()
    pipeline = config.get("Options","pipeline")
    branch = subprocess.check_output("cd /opt/gcf-workflows && git branch",stderr=subprocess.STDOUT,shell=True).split(b' ')[-1].rstrip()
    commit = subprocess.check_output("cd /opt/gcf-workflows && git log",stderr=subprocess.STDOUT,shell=True).split(b'\n')[0].split(b' ')[1]

    versions["Analysis pipeline"] = "github.com/gcfntnu/gcf-workflows/tree/{} commit {}".format(branch.decode(), commit.decode()).encode()

    software = '\n'.join('{}: {}'.format(key,val.decode()) for (key,val) in versions.items())
    with open(os.path.join(config.get("Paths","outputDir"),config.get("Options","runID"),"software.versions"),'w+') as sf:
        sf.write(software)
    return software

def post_workflow(project_id, base_dir, pipeline):
    analysis_dir = os.path.join(os.environ["TMPDIR"], "analysis_{}_{}".format(project_id, os.path.basename(base_dir).split("_")[0]))
    os.makedirs(os.path.join(base_dir, "QC_{}".format(project_id), "bfq"), exist_ok=True)
    cmd = "rsync -rvLp {}/ {}".format(
        os.path.join(analysis_dir, "data", "tmp", pipeline, "bfq"),
        os.path.join(base_dir, "QC_{}".format(project_id), "bfq"),
    )
    subprocess.check_call(cmd, shell=True)

    return None

def full_align(config):
    old_wd = os.getcwd()
    libprep = config.get("Options","Libprep")
    pipeline = config.get("Options","pipeline")

    base_dir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'))
    os.chdir(os.environ["TMPDIR"])
    project_names = get_project_names(get_project_dirs(config))
    for p in project_names:
        analysis_dir = os.path.join(os.environ["TMPDIR"],"analysis_{}_{}".format(p,os.path.basename(base_dir).split("_")[0]))
        os.makedirs(analysis_dir, exist_ok=True)
        os.makedirs(os.path.join(analysis_dir,"src"), exist_ok=True)
        os.makedirs(os.path.join(analysis_dir,"data"), exist_ok=True)

        os.chdir(analysis_dir)

        #create config.yaml
        cmd = "/opt/conda/bin/python /opt/conda/bin/configmaker.py {runfolder} -p {project} --libkit '{lib}' --machine '{machine}' {create_fastq}".format(
            runfolder = base_dir,
            project = p,
            lib = libprep,
            samplesheet = os.path.join(base_dir,"SampleSheet.csv"),
            sample_sub = os.path.join(base_dir,"Sample-Submission-Form.xlsx"),
            machine = get_sequencer(config.get("Options","runID")),
            create_fastq = " --create-fastq-dir" if not os.path.exists("data/raw/fastq") else ""
        )
        subprocess.check_call(cmd,shell=True)

        #copy snakemake pipeline
        cmd = "rm -rf {dst} && cp -r {src} {dst}".format(
            src = "/opt/gcf-workflows",
            dst = os.path.join(analysis_dir, "src", "gcf-workflows")
        )
        subprocess.check_call(cmd,shell=True)

        #write template Snakefile for pipeline
        with open(os.path.join(analysis_dir,"Snakefile"), "w") as sn:
            sn.write(SNAKEFILE_TEMPLATE.format(workflow=pipeline))

        #run snakemake pipeline
        cmd = "snakemake --reason --use-singularity --singularity-prefix $SINGULARITY_CACHEDIR -j32 -p bfq_all"
        subprocess.check_call(cmd,shell=True)

        #push workflow bfq output to project directory
        post_workflow(p, base_dir, pipeline)

        #push analysis folder
        analysis_export_dir = os.path.join(config.get("Paths","analysisDir"),"{}_{}".format(p,config.get("Options","runID").split("_")[0]))

        cmd = "rm -rf {dst} && cp -rv {src}/ {dst} ".format(
            src = analysis_dir,
            dst = analysis_export_dir,
        )
        subprocess.check_call(cmd,shell=True)

        #touch bfq_all to avoid rerunning pipelines from scratch
        os.chdir(analysis_export_dir)
        cmd = "snakemake --touch -j1 bfq_all"
        subprocess.check_call(cmd,shell=True)

        #chmod outputdir
        #cmd = "chmod -R 775 {}".format(analysis_export_dir)
        #subprocess.check_call(cmd,shell=True)

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
        libprep_config = yaml.load(lc_f)

    if config.get("Options", "libprep") + " SE" in libprep_config.keys():
        config["Options"]["pipeline"] = libprep_config[config.get("Options", "libprep") + ' SE']['workflow']
    elif config.get("Options", "libprep") + " PE" in libprep_config.keys():
        config["Options"]["pipeline"] = libprep_config[config.get("Options", "libprep") + ' PE']['workflow']
    else:
        print('failed to identify pipeline from libprep name, using default workflow')
        config["Options"]["pipeline"] = "default"



    #customer_samplesheet
    samplesheet_worker(config,projectDirs)

    for d in get_project_dirs(config):
        project_name = os.path.basename(d)
        os.makedirs(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"), "QC_{}".format(project_name)), exist_ok=True)

    if not os.path.exists(os.path.join(config["Paths"]["outputDir"], config["Options"]["runID"],"analysis.made")):
        full_align(config)

    # multiqc
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(multiqc_worker, projectDirs)
    p.close()
    p.join()

    # multiqc_stats
    multiqc_stats(projectDirs)

    #disk usage
    (tot,used,free) = shutil.disk_usage(config.get("Paths","outputDir"))
    tot /= 1024*1024*1024 #Convert to gigs
    used /= 1024*1024*1024
    free /= 1024*1024*1024

    #Undetermined indices
    #undeter = parserDemultiplexStats(config)

    message = "Current free space for output: %i of %i GiB (%5.2f%%)\n<br>" % (
        free,tot,100*free/tot)

    (tot,used,free) = shutil.disk_usage(config.get("Paths","baseDir"))
    tot /= 1024*1024*1024 #Convert to gigs
    used /= 1024*1024*1024
    free /= 1024*1024*1024


    message += "Current free space for instruments: %i of %i GiB (%5.2f%%)\n<br>\n<br>" % (
        free,tot,100*free/tot)

    #message += undeter

    #save configfile to flowcell
    with open(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),'bcl2fastq.ini'), 'w+') as configfile:
        config.write(configfile)


    return(message)

def finalize(config):
    #md5sum fastqs
    md5sum_worker(config)
    #zip arhive
    archive_worker(config)
    #md5sum archive
    md5sum_archive_worker(config)
    #archive instruments
    instrument_archive_worker(config)
    #md5sum instrument
    md5sum_instrument_worker(config)

    return None
