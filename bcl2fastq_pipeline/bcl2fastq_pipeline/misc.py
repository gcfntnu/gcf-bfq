"""
Misc. functions
"""
import pandas as pd
import tempfile as tmp
import configparser
import shutil
import smtplib
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.utils import COMMASPACE, formatdate
import xml.etree.ElementTree as ET
from reportlab.lib import colors, utils
from reportlab.platypus import BaseDocTemplate, Table, Preformatted, Paragraph, Spacer, Image, Frame, NextPageTemplate, PageTemplate, TableStyle, PageBreak, ListFlowable, ListItem
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4, landscape
from time import strftime
from reportlab.pdfgen import canvas
import csv
import sys
import glob
import pathlib
import os
import os.path
import syslog
import stat
import codecs
import requests
import json
import configmaker.configmaker as cm
from bcl2fastq_pipeline.afterFastq import get_project_dirs, get_project_names, get_sequencer, get_read_geometry

style = """
<style>
table {
  border-collapse: collapse;
}
th {
  padding: 4px;
}
td {
  text-align:center;
}
table, th, td {
  border: 1px solid black;
}
</style>
"""


def getSampleID(sampleTuple, project, lane, sampleName) :
    if(sampleTuple is None) :
        return " "
    for item in sampleTuple :
        if(sampleName == item[1] and
            lane == item[2] and
            project == item[3]) :
            return item[0]
    return " "


def getFCmetrics(root) :
    barcode = root[0][0] #Sample "all", barcode "all"
    message = "Lane\t# Clusters (% pass)\t% Bases >=Q30\tAve. base qual.\n"
    for lane in barcode.findall("Lane") :
        message += "Lane %s" % lane.get("number")
        clusterCount = 0
        clusterCountPass = 0
        baseYield = [0,0]
        baseYieldQ30 = [0,0]
        QualSum = [0,0]
        rlens=[0,0]
        for tile in lane :
            clusterCount += int(tile[0][0].text)
            clusterCountPass += int(tile[1][0].text)
            #Yield
            baseYield[0] += int(tile[1][1][0].text)
            if(len(tile[1]) == 3) :
                baseYield[1] += int(tile[1][2][0].text)
            #YieldQ30
            baseYieldQ30[0] += int(tile[1][1][1].text)
            if(len(tile[1]) == 3) :
                baseYieldQ30[1] += int(tile[1][2][1].text)
            #QualSum
            QualSum[0] += int(tile[1][1][2].text)
            if(len(tile[1]) == 3) :
                QualSum[1] += int(tile[1][2][2].text)
        #Number of clusters (%passing filter)
        try:
            message += "\t%s (%5.2f%%)" % ("{:,}".format(clusterCount).replace(","," "),100*clusterCountPass/clusterCount)
        except:
            message += "\t%s (NA)" % ("{:,}".format(clusterCount).replace(","," "))
        #%bases above Q30
        if(baseYield[1] > 0) :
            try:
                message += "\t%5.2f%%/%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]),
                    100*(baseYieldQ30[1]/baseYield[1]))
            except:
                message += "\tNA/NA"
        else :
            try:
                message += "\t%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]))
            except:
                message += "\tNA"
        #Average base quality
        if(baseYield[1] > 0) :
            try:
                message += "\t%4.1f/%4.1f\n" % (QualSum[0]/float(baseYield[0]),
                    QualSum[1]/float(baseYield[1]))
            except:
                message += "\tNA/NA\n"
        else :
            try:
                message += "\t%4.1f\n" % (QualSum[0]/float(baseYield[0]))
            except:
                message += "\tNA\n"

    return message

def getFCmetricsImproved(config):
    message = ""
    try:
        with open(os.path.join(config.get("Paths","outputDir"),config.get("Options","runID"),"Stats","interop_summary.csv"),"r") as fh:
            header = False
            while not header:
                line = fh.readline()
                if line.startswith("\n"):
                    line = fh.readline()
                    header = True
            lines = fh.readlines()
    except Exception as e:
        return "Not able to generate table for flowcell metrics."

    read_start = []
    for i,l in enumerate(lines):
        if l.startswith("Read"):
            read_start.append(i)
        elif l.startswith("Extracted"):
            read_start.append(i)
            break
    dfs = []
    for i in range(len(read_start)-1):
        if lines[read_start[i]].endswith("(I)\n"):
            continue
        tmpfh = tmp.NamedTemporaryFile(mode="w+")
        tmpfh.writelines(lines[read_start[i]+1:read_start[i+1]])
        tmpfh.seek(0)
        df = pd.read_csv(tmpfh)
        tmpfh.close()
        df = df[["Lane","Surface","Density","Reads","Cluster PF","Aligned","%>=Q30"]]
        df = df[df["Surface"]=='-']
        df = df.drop(columns=["Surface"])

        df["Density"] = [float(v.split(" ")[0]) for v in df["Density"]]
        df["Cluster PF"] = [float(v.split(" ")[0]) for v in df["Cluster PF"]]
        df["Aligned"] = [float(v.split(" ")[0]) for v in df["Aligned"]]

        df = df.round(2)

        mapper = {"Cluster PF": "% Cluster PF", "Reads": " Total Reads (M)", "Aligned": "% PhiX"}
        df = df.rename(columns=mapper)
        dfs.append(df)


    undeter = parserDemultiplexStats(config)
    if len(dfs) > 1:
        dfs[0]["R2 %>=Q30"] = dfs[1]["%>=Q30"]
        dfs[0] = dfs[0].rename(columns={"%>=Q30": "R1 %>=Q30"})
    dfs[0] = dfs[0].join(undeter.set_index("Lane"),on="Lane")
    #message += "\n<br>\n<br><strong>{} metrics </strong>\n<br>".format(lines[read_start[i]].rstrip())
    message += "\n<br><strong>Flowcell metrics </strong>\n<br>".format(lines[read_start[i]].rstrip())
    message += dfs[0].to_html(index=False,classes="border-collapse: collapse",border=1,justify="center",col_space=12)
    return message

def parseSampleSheetMetrics(config):
    project_dirs = get_project_dirs(config)
    project_names = get_project_names(project_dirs)
    sample_sheet_metrics = dict()
    msg = "<strong>Sample sheet info</strong>\n"
    for pid in project_names:
        with open(config.get("Options","sampleSheet"),'r') as ss:
            sample_df, _ = cm.get_project_samples_from_samplesheet(ss,[os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'))] , [pid])
            sample_dict = cm.find_samples(sample_df,[os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'), pid)])
        msg += "<strong>{}</strong>: Found {} samples in samplesheet.\n".format(pid, len(sample_df))

    sample_sub_form_d = cm.sample_submission_form_parser(config.get("Options","sampleSubForm"))
    msg += "\nFound {} samples in sample submission form.\n".format(len(sample_sub_form_d))
    ssub_df = pd.DataFrame.from_dict(sample_sub_form_d, orient="index")
    if "Sample_Group" in ssub_df:
        if ssub_df["Sample_Group"].notnull().all():
            msg += "Sample_Group has {} unique values: {}.\n".format(len(ssub_df["Sample_Group"].unique()), ", ".join(ssub_df["Sample_Group"].unique().astype(str)))
        elif ssub_df["Sample_Group"].isnull().all():
            msg += "Sample_Group has not been provided.\n"
        else:
            n_missing = ssub_df["Sample_Group"].isnull().values.sum()
            n_groups = len(ssub_df["Sample_Group"].dropna().unique())
            msg += "Sample_Group has {} unique values: {}.\n".format(n_groups, ", ".join(ssub_df["Sample_Group"].dropna().unique().astype(str)))
            msg += "Missing Sample_Group for {} sample{}.\n".format(str(n_missing), "s" if n_missing > 1 else "")
    else:
        msg += "Sample_Group has not been provided.\n"

    return msg

def parserDemultiplexStats(config) :
    '''
    Parse DemultiplexingStats.xml under outputDir/Stats/ to get the
    number/percent of undetermined indices.

    In particular, we extract the BarcodeCount values from Project "default"
    Sample "all" and Project "all" Sample "all", as the former gives the total
    undetermined and the later simply the total clusters
    '''
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    totals = [0,0,0,0,0,0,0,0]
    undetermined = [0,0,0,0,0,0,0,0]
    tree = ET.parse("%s/%s%s/Stats/DemultiplexingStats.xml" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
    root = tree.getroot()
    for child in root[0].findall("Project") :
        if(child.get("name") == "default") :
            break
    for sample in child.findall("Sample") :
        if(sample.get("name") == "all") :
            break
    child = sample[0] #Get inside Barcode
    for lane in child.findall("Lane") :
        lnum = int(lane.get("number"))
        undetermined[lnum-1] += int(lane[0].text)

    for child in root[0].findall("Project") :
        if(child.get("name") == "all") :
            break
    for sample in child.findall("Sample") :
        if(sample.get("name") == "all") :
            break
    child = sample[0] #Get Inside Barcode
    for lane in child.findall("Lane") :
        lnum = int(lane.get("number"))
        totals[lnum-1] += int(lane[0].text)

    out = ""
    lanes = []
    undeter = []
    for i in range(8) :
        if(totals[i] == 0) :
            continue
        #Bad hack with "{:,}".format(val).replace(","," ") for separator, but avoids using locale. The "right" locale would also yield unwanted results (comma as separatpr)
        """
        out += "Lane %i: %s of %s reads/pairs had undetermined indices (%5.2f%%)\n<br>" % (
            i+1,"{:,}".format(undetermined[i]).replace(","," "),"{:,}".format(totals[i]).replace(","," "),100*undetermined[i]/totals[i])
        """
        lanes.append(i+1)
        undeter.append(100*undetermined[i]/totals[i])
        #out_d.append({"Lane": i+1, "Undetermined": 100*undetermined[i]/totals[i]})
    return pd.DataFrame.from_dict({'Lane': lanes, "% Undetermined": undeter}).round(2)


def parseConversionStats(config) :
    """
    Parse ConversionStats.xml, producing:
     1) A PDF file for each project
     2) A message that will be included in the email message
    """
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    try :
        tree = ET.parse("%s/%s%s/Stats/ConversionStats.xml" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
        root = tree.getroot()[0] #We only ever have a single flow cell
    except :
        return None
    metrics = None
    #Per-project PDF files
    for project in root.findall("Project") :
        if(project.get("name") == "default") :
            continue
        if(project.get("name") == "all") :
            metrics = getFCmetrics(project)
    return metrics

def enoughFreeSpace(config) :
    """
    Ensure that outputDir has at least minSpace gigs
    """
    (tot,used,free) = shutil.disk_usage(config.get("Paths","outputDir"))
    free /= 1024*1024*1024
    if(free >= float(config.get("Options","minSpace"))) :
        return True
    return False

def errorEmail(config, errTuple, msg) :
    msg = msg + "\nError type: %s\nError value: %s\n%s\n" % (errTuple[0], errTuple[1], errTuple[2])
    """
    msg['Subject'] = "[bcl2fastq_pipeline] Error"
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","errorTo")
    """
    with open(os.path.join(config.get("Paths", "reportDir"),'{}.error'.format(config.get("Options","runID"))),'w') as report:
        report.write(msg)
    """
    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()
    """

def finishedEmail(config, msg, runTime) :
    lanes = config.get("Options", "lanes")

    projects = get_project_names(get_project_dirs(config))

    message = "<strong>Short summary for {}. </strong>\n\n".format(", ".join(projects))
    message += "<strong>User: {} </strong>\n".format(config.get("Options","User")) if config.get("Options","User") != "N/A" else ""
    message += "Flow cell: %s \n" % (config.get("Options","runID"))
    message += "Sequencer: {} \n".format(get_sequencer(os.path.join(config.get("Paths","baseDir"),config.get("Options","runID"))))
    message += "Read geometry: {} \n\n".format(get_read_geometry(os.path.join(config.get("Paths","outputDir"),config.get("Options","runID"))))
    message += "bcl2fastq_pipeline run time: %s \n" % runTime
    #message += "Data transfer: %s\n" % transferTime
    message = message.replace("\n","\n<br>")
    message += msg

    sample_sheet_metrics = parseSampleSheetMetrics(config)
    sample_sheet_metrics = sample_sheet_metrics.replace("\n", "\n<br>")

    message = "<html>\n<body>\n<head>\n" + style + "\n</head>\n" + message + "<br>"  + sample_sheet_metrics  + "\n</body>\n</html>"

    odir = os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"))

    #with open(os.path.join(config.get("Paths", "reportDir"),'{}.report'.format(config.get("Options","runID"))),'w') as report:
    #    report.write(msg)
    msg = MIMEMultipart()
    msg['Subject'] = "[bcl2fastq_pipeline] {} processed".format(", ".join(projects))
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","finishedTo")
    msg['Date'] = formatdate(localtime=True)

    msg.attach(MIMEText(message,'html'))

    for p in projects:
        with open(os.path.join(odir,"QC_{pnr}/multiqc_{pnr}.html".format(pnr=p)),"rb") as report:
            part = MIMEApplication(
                    report.read(),
                    report.name
                    )
        part['Content-Disposition'] = 'attachment; filename="multiqc_{}.html"'.format(p)
        msg.attach(part)


    with open(os.path.join(odir,"Stats/sequencer_stats_{}.html".format("_".join(projects))),"rb") as report:
        part = MIMEApplication(
                report.read(),
                report.name
                )
    part['Content-Disposition'] = 'attachment; filename="sequencer_stats_{}.html"'.format("_".join(projects))
    msg.attach(part)

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()

def finalizedEmail(config, msg, finalizeTime, runTime) :
    lanes = config.get("Options", "lanes")

    projects = get_project_names(get_project_dirs(config))

    message = "{} has been finalized and prepared for delivery.\n\n".format(", ".join(projects))
    message += "md5sum and 7zip runtime: %s\n" % finalizeTime
    message += "Total runtime for bcl2fastq_pipeline: %s\n" % runTime
    #message += "Data transfer: %s\n" % transferTime
    message += msg

    #message = "<html>\n<body>\n<head></head>\n" + message + "\n</body>\n</html>"


    odir = os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"))

    #with open(os.path.join(config.get("Paths", "reportDir"),'{}.report'.format(config.get("Options","runID"))),'w') as report:
    #    report.write(msg)
    msg = MIMEMultipart()
    msg['Subject'] = "[bcl2fastq_pipeline] {} finalized".format(", ".join(projects))
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","errorTo")
    msg['Date'] = formatdate(localtime=True)

    msg.attach(MIMEText(message))

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()
