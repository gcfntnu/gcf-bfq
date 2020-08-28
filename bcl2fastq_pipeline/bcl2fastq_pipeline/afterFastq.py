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

PIPELINE_MAP = {
    'Lexogen SENSE Total RNA-Seq Library Prep Kit (w/RiboCop rRNA Depletion Kit V1.2)': 'rna-seq',
    'Lexogen SENSE mRNA-Seq Library Prep Kit V2': 'rna-seq',
    'Illumina TruSeq Stranded Total RNA Library Prep (Human/Mouse/Rat)': 'rna-seq',
    'Illumina TruSeq Stranded Total RNA Library Prep (Globin)': 'rna-seq',
    'Illumina TruSeq Stranded mRNA Library Prep': 'rna-seq',
    'QIAseq 16S ITS Region Panels': 'microbiome',
    '16S Metagenomic Sequencing Library Prep': 'microbiome',
    'ITS Low Input GCF Custom': 'microbiome',
    'ITS2 Low Input GCF Custom': 'microbiome',
    '10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3': 'single-cell'
}

PIPELINE_MULTIQC_MODULES = {
    'rna-seq': ["fastq_screen","star","picard","fastp","fastqc_rnaseq","custom_content"],
    'microbiome': ["fastq_screen","star","picard","fastp","fastqc_rnaseq","custom_content", "qiime2"],
    'single-cell': ["fastq_screen","star", "picard", "cellranger", "starsolo", "fastp","fastqc_rnaseq","custom_content"],
}

PIPELINE_ORGANISMS = {
    'rna-seq': ['homo_sapiens', 'mus_musculus', 'rattus_norvegicus', 'salmo_salar'],
    'microbiome': ['N/A'],
    'single-cell': ['homo_sapiens', 'mus_musculus', 'rattus_norvegicus'],
}


def bgzip_worker(fname) :
    global localConfig
    config = localConfig
    cmd = "%s -r %s" % (
        config.get("bgzip","bgzip_command"),
        fname)
    syslog.syslog("[bgzip_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

def clumpify_worker(fname):
    global localConfig
    config = localConfig

    if "10X Genomics" in config.get("Options","Libprep"):
        return

    #We use R1 files to construct R2 filenames for paired end
    if "R2.fastq.gz" in fname:
        return

    if os.path.exists(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"clumpify.done")):
        return

    r1 = fname
    r2 = fname.replace("R1.fastq.gz","R2.fastq.gz") if os.path.exists(fname.replace("R1.fastq.gz","R2.fastq.gz")) else None

    if r2:
        out_r1 = r1.replace(".fastq.gz","_clumped.fastq.gz") 
        out_r2 = r2.replace(".fastq.gz","_clumped.fastq.gz")
        cmd = "{clump_cmd} {clump_opts} in1={in1} in2={in2} out1={out1} out2={out2} rcomp=f rename=f overwrite=true".format(
                clump_cmd = config.get("clumpify","clumpify_cmd"),
                clump_opts = config.get("clumpify", "clumpify_opts"),
                in1 = r1,
                in2 = r2,
                out1 = out_r1, #yes, we want to overwrite the files
                out2 = out_r2
                )
    else:
        out_r1 = r1.replace(".fastq.gz","_clumped.fastq.gz") 
        cmd = "{clump_cmd} {clump_opts} in={in1} out={out1} rename=f overwrite=true".format(
                clump_cmd = config.get("clumpify","clumpify_cmd"),
                clump_opts = config.get("clumpify", "clumpify_opts"),
                in1 = r1,
                out1 = out_r1 
                )
    syslog.syslog("[clumpify_worker] Processing %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)
    os.rename(out_r1,r1)
    if r2:
        os.rename(out_r2,r2)

def RemoveHumanReads_worker(fname):
    global localConfig
    config = localConfig

    #We use R1 files to construct R2 filenames for paired end
    if "R2.fastq.gz" in fname:
        return

    if os.path.exists(fname.replace(os.path.basename(fname),"contaminated/{}".format(os.path.basename(fname)))):
        return

    r2 = fname.replace("R1.fastq.gz","R2.fastq.gz") if os.path.exists(fname.replace("R1.fastq.gz","R2.fastq.gz")) else None

    masked_path = config.get("MaskedGenomes","HGDir")
    if not masked_path:
        raise Exception("HGDir not set in config bcl2fastq.ini!\n")

    if not os.path.exists(masked_path):
        raise Exception("HGDir {} does not exist!\n".format(masked_path))
    
    os.makedirs("{}/contaminated/".format(os.path.dirname(fname)),exist_ok=True)

    if r2:
        cont_out = fname.replace(os.path.basename(fname),"contaminated/{}".format(os.path.basename(fname).replace("R1.fastq.gz","interleaved.fastq.gz")))
        cmd = "{bbmap_cmd} {bbmap_opts} path={masked_path} in={infile} in2={infile2} outu={clean_out} outm={contaminated_out}".format(
                bbmap_cmd = config.get("MaskedGenomes","bbmap_cmd"),
                bbmap_opts = config.get("MaskedGenomes","bbmap_opts"),
                masked_path = config.get("MaskedGenomes","HGDir"),
                infile = fname,
                infile2 = r2,
                clean_out = fname.replace("R1.fastq.gz","interleaved.fastq.gz"),
                contaminated_out = cont_out 
                )
    else:   
        cmd = "{bbmap_cmd} {bbmap_opts} path={masked_path} in={infile} outu={clean_out} outm={contaminated_out}".format(
                bbmap_cmd = config.get("MaskedGenomes","bbmap_cmd"),
                bbmap_opts = config.get("MaskedGenomes","bbmap_opts"),
                masked_path = config.get("MaskedGenomes","HGDir"),
                infile = fname,
                clean_out = fname,
                contaminated_out = fname.replace(os.path.basename(fname),"contaminated/{}".format(os.path.basename(fname)))
                )

    syslog.syslog("[RemoveHumanReads_worker] Processing %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)
    # clean up
    os.remove(fname)
    if r2:
        os.remove(r2)
        #If paired end, split interleaved files to R1 and R2
        cmd = "{rename_cmd} {rename_opts} in={interleaved} out1={out_r1} out2={out_r2}".format(
                rename_cmd = "rename.sh",
                rename_opts = "renamebymapping=t",
                interleaved = fname.replace("R1.fastq.gz","interleaved.fastq.gz"),
                out_r1 = fname,
                out_r2 = r2
                )
        syslog.syslog("[RemoveHumanReads_worker] De-interleaving %s\n" % cmd)
        subprocess.check_call(cmd, shell=True)
        os.remove(fname.replace("R1.fastq.gz","interleaved.fastq.gz"))
        #De-interleave contaminated file
        cmd = "{rename_cmd} {rename_opts} in={interleaved} out1={out_r1} out2={out_r2}".format(
                rename_cmd = "rename.sh",
                rename_opts = "renamebymapping=t",
                interleaved = cont_out,
                out_r1 = cont_out.replace("interleaved.fastq.gz","R1.fastq.gz"),
                out_r2 = cont_out.replace("interleaved.fastq.gz","R2.fastq.gz")
                )
        syslog.syslog("[RemoveHumanReads_worker] De-interleaving %s\n" % cmd)
        subprocess.check_call(cmd, shell=True)
        os.remove(cont_out)


def fastq_screen_worker(fname) :
    global localConfig
    config = localConfig

    os.chdir(os.path.dirname(fname))

    project_nr = get_gcf_name(fname)

    #Skip read 1 when single cell
    if "10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3" == config.get("Options","Libprep") and ("R1.fastq" in fname or "R1_001.fastq" in fname):
        return

    ofile="{}/{}/QC_{}/fastq_screen/{}".format(
            os.path.join(config.get("Paths","outputDir")),
            config.get("Options","runID"),
            project_nr,
            os.path.basename(fname)
            )

    if os.path.exists(ofile.replace(".fastq.gz","_screen.html")) :
        return

    os.makedirs(os.path.dirname(ofile), exist_ok=True)

    #fastq_screen
    cmd = "%s %s --outdir '%s' '%s'" % (
        config.get("fastq_screen", "fastq_screen_command"),
        config.get("fastq_screen", "fastq_screen_options"),
        os.path.dirname(ofile),
        fname)
    syslog.syslog("[fastq_screen_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

    #Unlink/rename
    #os.unlink(ofile)
    #os.rename(ofile.replace(".fastq","_screen.png"), fname.replace("_R1_001.fastq.gz", "_R1_001_screen.png"))


def suprDUPr_worker(fname) :
    global localConfig
    config = localConfig

    if "10X Genomics" in config.get("Options","Libprep"):
        return

    ofname = "{}/filtered/{}".format(os.path.dirname(fname),os.path.basename(fname).replace(".fastq.gz","_filtered.fastq.gz"))

    # WHEN ./filterfq is working cmd will be 
    #  "{cmd} {opts} {infile} | filterfq {infile} | tee {ofile}"
    cmd = "{cmd} {opts} {infile} > {ofile}".format(
          cmd = config.get("suprDUPr","suprdupr_command"),
          opts = config.get("suprDUPr","suprdupr_options"),
          infile = fname,
          ofile = ofname
          )

    # Skip if the output exists
    if os.path.exists(ofname):
        return

    os.makedirs(os.path.dirname(ofname), exist_ok=True)

    syslog.syslog("[suprDUPr_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

def get_gcf_name(fname):
    for name in fname.split('/'):
        if re.match(r'GCF-[0-9]{4}-[0-9]{3,}',name):
            return re.search(r'GCF-[0-9]{4}-[0-9]{3,}',name)[0]
    raise Exception('Unable to determine GCF project number for filename {}\n'.format(fname))



def FastQC_worker(fname) :
    global localConfig
    config = localConfig
    lanes = config.get("Options", "lanes")

    #projectName = fname.split("/")[-3]
    projectName = get_gcf_name(fname)

    # Skip if the output exists
    fastqc_fname = glob.glob("{}/{}/QC_{}/FASTQC/*/{}".format(
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          projectName,
          os.path.basename(fname).replace(".fastq.gz","_fastqc.zip")
          ))

    fastqc_fname.extend(glob.glob("{}/{}/QC_{}/FASTQC/{}".format(
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          projectName,
          os.path.basename(fname).replace(".fastq.gz","_fastqc.zip")
          )))

    if fastqc_fname:
        return

    os.makedirs("%s/%s/QC_%s/FASTQC" % (
            config.get("Paths","outputDir"),
            config.get("Options","runID"),
            projectName), exist_ok=True
            )

    cmd = "%s %s -o %s/%s/QC_%s/FASTQC %s" % (
          config.get("FastQC","fastqc_command"),
          config.get("FastQC","fastqc_options"),
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          projectName,
          fname
          )

    syslog.syslog("[FastQC_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

def fastp_worker(fname):
    global localConfig
    config = localConfig

    if config.get("Options","Libprep") in ["10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3", "10X Genomics Visium Spatial Gene Expression Slide & Reagents Kit"]:
        in2 = ""
    elif config.get("Options","Libprep") == "10X Genomics Chromium Next GEM Single Cell ATAC Library & Gel Bead Kit v1.1":
        if "R3_001.fastq.gz" in fname:
            return
        in2 = ("--in2=" + fname.replace("R1_001.fastq.gz","R3_001.fastq.gz")) if os.path.exists(fname.replace("R1_001.fastq.gz","R3_001.fastq.gz")) else ""
    else:
        #We use R1 files to construct R2 filenames for paired end
        if "R2.fastq.gz" in fname:
            return
        in2 = ("--in2=" + fname.replace("R1.fastq.gz","R2.fastq.gz")) if os.path.exists(fname.replace("R1.fastq.gz","R2.fastq.gz")) else ""

    project_name = get_gcf_name(fname)
    if os.path.exists(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"QC_{}".format(project_name),"fastp","{}.fastp.json".format(os.path.basename(fname).replace("_R1.fastq.gz","")))):
        return
    elif os.path.exists(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"QC_{}".format(project_name),"fastp","{}.fastp.json".format(os.path.basename(fname).replace("_R1_001.fastq.gz","")))):
        return

    if "10X Genomics" in config.get("Options","Libprep"):
        json_out = os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"QC_{}".format(project_name),"fastp","{}.fastp.json".format(os.path.basename(fname).replace("_R1_001.fastq.gz","")))
        html_out = os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"QC_{}".format(project_name),"fastp","{}.fastp.html".format(os.path.basename(fname).replace("_R1_001.fastq.gz","")))
    else:
        json_out = os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"QC_{}".format(project_name),"fastp","{}.fastp.json".format(os.path.basename(fname).replace("_R1.fastq.gz","")))
        html_out = os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"QC_{}".format(project_name),"fastp","{}.fastp.html".format(os.path.basename(fname).replace("_R1.fastq.gz","")))


    tmpdir = os.path.join(os.environ["TMPDIR"], "{}_{}".format(project_name, config.get("Options","runID")))
    os.makedirs(tmpdir, exist_ok=True)

    with open("/opt/gcfdb/libprep.config","r") as libprepconf:
        libconf = yaml.load(libprepconf)

    l_conf = libconf.get(config.get("Options","Libprep"), {})
    if l_conf:
        if in2:
            #Paired End
            if 'paired_end' in l_conf.keys():
                adapter1 = l_conf.get('paired_end').get('adapter','auto')
                adapter2 = l_conf.get('paired_end').get('adapter2','auto')
            else:
                adapter1 = 'auto'
                adapter2 = 'auto'
        else:
            #Single End
            adapter2 = ''
            if 'single_end' in l_conf.keys():
                adapter1 = l_conf.get('single_end').get('adapter','auto')
            elif 'adapter' in l_conf.keys():
                adapter1 = l_conf.get('adapter','auto')
            else:
                adapter1 = 'auto'
    else:
        adapter1 = 'auto'
        adapter2 = 'auto' if in2 else ''

    adapter = "--adapter_sequence=" + adapter1
    if adapter2:
        adapter += " --adapter_sequence_r2=" + adapter2 + (" --detect_adapter_for_pe " if adapter2 == 'auto' else '')

    out1 = os.path.join(tmpdir ,os.path.basename(fname))
    out2 = ("--out2=" + os.path.join(tmpdir, os.path.basename(in2))) if in2 else ""

    cmd = "{fastp} --in1={in1} {in2} --out1={out1} {out2} --json={json} --html={html} {adapter} --thread={thread}".format(
            fastp = config.get("fastp","command"),
            in1 = fname,
            in2 = in2,
            out1 = out1,
            out2 = out2,
            json = json_out,
            html = html_out,
            adapter = adapter,
            thread = config.get("fastp","threads")
            )
    syslog.syslog("[fastp_worker] Processing %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)


def toDirs(files) :
    s = set()
    for f in files :
        d = os.path.dirname(f)
        if d.split('/')[-1].startswith("GCF-"):
            s.add(d)
        else:
            s.add(d[:d.rfind('/')])
    return s


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
        pipeline = PIPELINE_MAP.get(config.get("Options","Libprep"),"default")
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

        """
        data = {}
        if read_geometry.startswith('Paired end') and not seq_stats and not PIPELINE_MAP.get(config.get("Options","Libprep"),'') == 'rna-seq':
            for k,v in s_dict.items():
                data['{}_R1'.format(k)] = v
                data['{}_R2'.format(k)] = v
        else:
            data = s_dict
        """
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
    modules = PIPELINE_MULTIQC_MODULES.get(PIPELINE_MAP.get(config.get("Options","Libprep"),None),None)
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
    pipeline = PIPELINE_MAP.get(config.get("Options","Libprep"),None)
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
        cmd = "7za a {opts} {flowdir}/{pnr}.7za {flowdir}/{pnr}/ {flowdir}/QC_{pnr} {flowdir}/Stats {flowdir}/Undetermined*.fastq.gz {flowdir}/{pnr}_samplesheet.tsv {flowdir}/SampleSheet.csv {flowdir}/Sample-Submission-Form.xlsx {flowdir}/md5sum_{pnr}_fastq.txt {flowdir}/software.versions".format(
                opts = opts,
                flowdir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
                pnr = p
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


def clumpify_mark_done(config):
    global localConfig
    config = localConfig
    open(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"clumpify.done"),"w+").close()

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

def get_software_versions(config):
    versions = {}
    versions['GCF-NTNU bcl2fastq pipeline'] = subprocess.check_output("cat /bfq_version/bfq.version",shell=True).rstrip()
    if "10X Genomics" in config.get("Options","Libprep"):
        subprocess.check_output("cellranger mkfastq --version",stderr=subprocess.STDOUT,shell=True).split(b'\n')[1].split(b' ')[-1].rstrip().strip(b'(').strip(b')')
    else:
        versions['bcl2fastq'] = subprocess.check_output("bcl2fastq --version",stderr=subprocess.STDOUT,shell=True).split(b'\n')[1].split(b' ')[1].rstrip()
        versions['fastq_screen'] = subprocess.check_output("fastq_screen --version",shell=True).split(b' ')[-1].rstrip()
    versions['FastQC'] = subprocess.check_output("fastqc --version",shell=True).split(b' ')[-1].rstrip()
    versions['clumpify/bbmap'] = subprocess.check_output("clumpify.sh --version",stderr=subprocess.STDOUT,shell=True).split(b'\n')[1].split(b' ')[-1].rstrip()
    versions['multiqc'] = subprocess.check_output("multiqc --version",stderr=subprocess.STDOUT,shell=True).split(b' ')[-1].rstrip()
    software = '\n'.join('{}: {}'.format(key,val.decode()) for (key,val) in versions.items())
    with open(os.path.join(config.get("Paths","outputDir"),config.get("Options","runID"),"software.versions"),'w+') as sf:
        sf.write(software)
    return software

def post_rna_seq(var_d):
    p = var_d['p']
    base_dir = var_d['base_dir']

    #move expressions and bam
    analysis_dir = os.path.join(os.environ["TMPDIR"], "analysis_{}_{}".format(p,os.path.basename(base_dir).split("_")[0]))
    os.makedirs(os.path.join(base_dir, "QC_{}".format(p), "bfq"), exist_ok=True)
    cmd = "rsync -rvLp {}/ {}".format(
        os.path.join(analysis_dir, "data", "tmp", "rnaseq", "bfq"),
        os.path.join(base_dir, "QC_{}".format(p), "bfq"),
    )
    subprocess.check_call(cmd, shell=True)

    #move logs
    cmd = "rsync -rvLp {}/ {}".format(
        os.path.join(analysis_dir,"logs"),
        os.path.join(base_dir, "QC_{}".format(p),"logs"),
    )
    subprocess.check_call(cmd, shell=True)
    return None

def post_microbiome(var_d):
    p = var_d['p']
    base_dir = var_d['base_dir']

    #move expressions and bam
    analysis_dir = os.path.join(os.environ["TMPDIR"], "analysis_{}_{}".format(p,os.path.basename(base_dir).split("_")[0]))
    os.makedirs(os.path.join(base_dir, "QC_{}".format(p), "bfq"), exist_ok=True)
    cmd = "rsync -rvLp {}/ {}".format(
        os.path.join(analysis_dir, "data", "tmp", "microbiome", "bfq"),
        os.path.join(base_dir, "QC_{}".format(p), "bfq"),
    )
    subprocess.check_call(cmd, shell=True)

    return None

def post_single_cell(var_d):
    p = var_d['p']
    base_dir = var_d['base_dir']
    #move expressions and bam
    analysis_dir = os.path.join(os.environ["TMPDIR"], "analysis_{}_{}".format(p,os.path.basename(base_dir).split("_")[0]))
    os.makedirs(os.path.join(base_dir, "QC_{}".format(p), "bfq"), exist_ok=True)
    cmd = "rsync -rvLp {}/ {}".format(
        os.path.join(analysis_dir, "data", "tmp", "singlecell", "bfq"),
        os.path.join(base_dir, "QC_{}".format(p), "bfq"),
    )
    subprocess.check_call(cmd, shell=True)

    return None

POST_PIPELINE_MAP = {
    'rna-seq': post_rna_seq,
    'microbiome': post_microbiome,
    'single-cell': post_single_cell,
}

def full_align(config):
    old_wd = os.getcwd()
    libprep = config.get("Options","Libprep")
    pipeline = PIPELINE_MAP.get(libprep,None)
    if not pipeline:
        return False
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
            src = os.path.join("/opt",pipeline),
            dst = os.path.join(analysis_dir,"src",pipeline)
        )
        subprocess.check_call(cmd,shell=True)

        #copy template Snakefile for pipeline
        shutil.copyfile(
            os.path.join("/opt","snakefiles","Snakefile-{}".format(pipeline)),
            os.path.join(analysis_dir,"Snakefile"),
            )

        if pipeline in ['microbiome', 'single-cell']:
            os.makedirs(os.path.join(analysis_dir, "data", "tmp"), exist_ok=True)
            cmd = "cp {} {}".format(
                os.path.join(base_dir, "{}_samplesheet.tsv".format(p)),
                os.path.join(analysis_dir, "data", "tmp", "sample_info.tsv")
            )
            subprocess.check_call(cmd,shell=True)

        #run snakemake pipeline
        cmd = "snakemake --reason --use-singularity --singularity-prefix $SINGULARITY_CACHEDIR -j32 -p bfq_all"
        subprocess.check_call(cmd,shell=True)

        POST_PIPELINE_MAP[pipeline]({'p': p, 'base_dir': base_dir})

        #push analysis folder
        analysis_export_dir = os.path.join(config.get("Paths","analysisDir"),"{}_{}".format(p,config.get("Options","runID").split("_")[0]))

        cmd = "cp -rv {src}/ {dst} ".format(
            src = analysis_dir,
            dst = analysis_export_dir,
        )
        subprocess.check_call(cmd,shell=True)
    os.chdir(old_wd)
    open(os.path.join(config["Paths"]["outputDir"], config["Options"]["runID"],"analysis.made"), "w").close()
    return True

def get_tmp(config, fname):
    project_name = get_gcf_name(fname)
    tmpdir = os.path.join(os.environ["TMPDIR"], "{}_{}".format(project_name, config.get("Options","runID")))
    if os.path.exists(os.path.join(tmpdir, os.path.basename(fname))):
        return os.path.join(tmpdir, os.path.basename(fname))

def get_tmp_sample_files(config, sample_files):
    return [get_tmp(config, x) for x in sample_files]

def clean_up_tmp_sample_files(tmp_sample_files):
    for x in tmp_sample_files:
        os.unlink(x)

def concat_sc_rna_seq_input_files(config):
    return None

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

    if config.get("Options", "Libprep") == "10X Genomics Chromium Next GEM Single Cell ATAC Library & Gel Bead Kit v1.1":
        sampleFiles = glob.glob("{}/{}/*/*R[13].fastq.gz".format(config.get("Paths","outputDir"),config.get("Options","runID")))
        sampleFiles.extend(glob.glob("{}/{}/*/*/*R[13].fastq.gz".format(config.get("Paths","outputDir"),config.get("Options","runID"))))
        sampleFiles.extend(glob.glob("{}/{}/*/*R[13]_001.fastq.gz".format(config.get("Paths","outputDir"),config.get("Options","runID"))))
        sampleFiles.extend(glob.glob("{}/{}/*/*/*R[13]_001.fastq.gz".format(config.get("Paths","outputDir"),config.get("Options","runID"))))
    elif config.get("Options", "Libprep") == "10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3":
        sampleFiles = []
        for d in get_project_dirs(config):
            project_name = os.path.basename(d)
            tmp_raw = os.path.join(os.environ["TMPDIR"], "{}_{}".format(project_name, config.get("Options", "runID")), "raw")
            os.makedirs(tmp_raw, exist_ok=True)
            for sample in os.listdir(d):
                s_files = sorted(glob.glob(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"), project_name, sample, "{}_*_R2_001.fastq.gz".format(sample))))
                tmp_sample = os.path.join(tmp_raw, "{}.fastq.gz".format(sample))
                cmd = "cat {infiles} > {tmp_sample}".format(
                    infiles = " ".join(s_files),
                    tmp_sample = tmp_sample,
                )
                subprocess.check_call(cmd, shell=True)
                sampleFiles.append(tmp_sample)
    elif config.get("Options", "Libprep") ==  "10X Genomics Visium Spatial Gene Expression Slide & Reagents Kit":
        sampleFiles = []
        for d in get_project_dirs(config):
            project_name = os.path.basename(d)
            tmp_raw = os.path.join(os.environ["TMPDIR"], "{}_{}".format(project_name, config.get("Options", "runID")), "raw")
            os.makedirs(tmp_raw, exist_ok=True)
            sample_names = set()
            dna_reads = glob.glob(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"), project_name, "*_R2_001.fastq.gz"))
            for sample in dna_reads:
                fn_format = re.search("_S[0-9]+_L[0-9]{3}_[RI][123]_001.fastq.gz", sample)
                sample = sample.replace(fn_format.group(0),"")
                sample_names.add(os.path.basename(sample))
            for sample in sample_names:
                s_files = sorted(glob.glob(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"), project_name, "{}_*_R2_001.fastq.gz".format(sample))))
                tmp_sample = os.path.join(tmp_raw, "{}.fastq.gz".format(sample))
                cmd = "cat {infiles} > {tmp_sample}".format(
                    infiles = " ".join(s_files),
                    tmp_sample = tmp_sample,
                )
                subprocess.check_call(cmd, shell=True)
                sampleFiles.append(tmp_sample)

    else:
        sampleFiles = glob.glob("{}/{}/*/*R[12].fastq.gz".format(config.get("Paths","outputDir"),config.get("Options","runID")))
        sampleFiles.extend(glob.glob("{}/{}/*/*/*R[12].fastq.gz".format(config.get("Paths","outputDir"),config.get("Options","runID"))))
        sampleFiles.extend(glob.glob("{}/{}/*/*R[12]_001.fastq.gz".format(config.get("Paths","outputDir"),config.get("Options","runID"))))
        sampleFiles.extend(glob.glob("{}/{}/*/*/*R[12]_001.fastq.gz".format(config.get("Paths","outputDir"),config.get("Options","runID"))))

    sampleFiles = [sf for sf in sampleFiles if not 'contaminated' in sf]
    sampleFiles = [sf for sf in sampleFiles if not 'filtered' in sf]

    global localConfig
    localConfig = config

    #Decontaminate with masked genome
    if config.get("Options","RemoveHumanReads") == "1":
        p = mp.Pool(int(1))
        p.map(RemoveHumanReads_worker, sampleFiles)
        p.close()
        p.join()

    if not os.path.exists(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"qc.done")):
        #fastp
        for d in get_project_dirs(config):
            project_name = os.path.basename(d)
            os.makedirs(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"QC_{}".format(project_name),"fastp"), exist_ok=True)

        p = mp.Pool(int(int(config.get("system","threads"))/int(config.get("fastp","threads"))))
        p.map(fastp_worker, sampleFiles)
        p.close()
        p.join()
        tmp_sample_files = get_tmp_sample_files(config, sampleFiles)

        if config.get("Options", "Libprep") in ["10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3", "10X Genomics Visium Spatial Gene Expression Slide & Reagents Kit"]:
            for d in get_project_dirs(config):
                project_name = os.path.basename(d)
                tmp_raw = os.path.join(os.environ["TMPDIR"], "{}_{}".format(project_name, config.get("Options", "runID")), "raw")
                cmd = "rm -rf {}".format(tmp_raw)
                subprocess.check_call(cmd, shell=True)
        #FastQC
        p = mp.Pool(int(config.get("Options","fastqcThreads")))
        p.map(FastQC_worker, tmp_sample_files)
        p.close()
        p.join()
        #fastq_screen
        p = mp.Pool(int(config.get("Options", "fastqScreenThreads")))
        p.map(fastq_screen_worker, tmp_sample_files)
        p.close()
        p.join()
        if tmp_sample_files:
            clean_up_tmp_sample_files(tmp_sample_files)
        open(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"qc.done"),"w").close()
    #customer_samplesheet
    samplesheet_worker(config,projectDirs)
    if config.get("Options","Organism") in PIPELINE_ORGANISMS.get(PIPELINE_MAP.get(config.get("Options","Libprep"),None),[]) and not os.path.exists(os.path.join(config["Paths"]["outputDir"], config["Options"]["runID"],"analysis.made")):
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
