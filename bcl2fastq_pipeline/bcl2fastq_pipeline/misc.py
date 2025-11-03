"""
Misc. functions
"""

import logging
import shutil
import smtplib
import tempfile as tmp
import xml.etree.ElementTree as ET

from argparse import Namespace
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import formatdate

import configmaker.configmaker as cm
import pandas as pd

from bcl2fastq_pipeline.afterFastq import (
    get_project_dirs,
    get_project_names,
    get_read_geometry,
    get_sequencer,
)
from bcl2fastq_pipeline.config import PipelineConfig

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

log = logging.getLogger(__name__)


def getSampleID(sampleTuple, project, lane, sampleName):
    if sampleTuple is None:
        return " "
    for item in sampleTuple:
        if sampleName == item[1] and lane == item[2] and project == item[3]:
            return item[0]
    return " "


def getFCmetricsImproved():
    cfg = PipelineConfig.get()
    message = ""
    try:
        with (cfg.output_path / "Stats" / "interop_summary.csv").open() as fh:
            header = False
            while not header:
                line = fh.readline()
                if line.startswith("\n"):
                    line = fh.readline()
                    header = True
            lines = fh.readlines()
    except Exception:
        return "Not able to generate table for flowcell metrics."

    read_start = []
    for i, line in enumerate(lines):
        if line.startswith("Read"):
            read_start.append(i)
        elif line.startswith("Extracted"):
            read_start.append(i)
            break
    dfs = []
    for i in range(len(read_start) - 1):
        if lines[read_start[i]].endswith("(I)\n"):
            continue
        tmpfh = tmp.NamedTemporaryFile(mode="w+")
        tmpfh.writelines(lines[read_start[i] + 1 : read_start[i + 1]])
        tmpfh.seek(0)
        df = pd.read_csv(tmpfh)
        tmpfh.close()
        df = df[["Lane", "Surface", "Density", "Reads", "Cluster PF", "Aligned", "%>=Q30"]]
        df = df[df["Surface"] == "-"]
        df = df.drop(columns=["Surface"])

        df["Density"] = [float(v.split(" ")[0]) for v in df["Density"]]
        df["Cluster PF"] = [float(v.split(" ")[0]) for v in df["Cluster PF"]]
        df["Aligned"] = [float(v.split(" ")[0]) for v in df["Aligned"]]

        df = df.round(2)

        mapper = {"Cluster PF": "% Cluster PF", "Reads": " Total Reads (M)", "Aligned": "% PhiX"}
        df = df.rename(columns=mapper)
        dfs.append(df)

    undeter = parserDemultiplexStats(cfg)
    if len(dfs) > 1:
        dfs[0]["R2 %>=Q30"] = dfs[1]["%>=Q30"]
        dfs[0] = dfs[0].rename(columns={"%>=Q30": "R1 %>=Q30"})
    dfs[0] = dfs[0].join(undeter.set_index("Lane"), on="Lane")
    message += "\n<br><strong>Flowcell metrics </strong>\n<br>"
    message += dfs[0].to_html(
        index=False, classes="border-collapse: collapse", border=1, justify="center", col_space=12
    )
    return message


def parseSampleSheetMetrics(cfg):
    project_dirs = get_project_dirs(cfg)
    project_names = get_project_names(project_dirs)
    msg = "<strong>Sample sheet info</strong>\n"
    for pid in project_names:
        args = Namespace(
            samplesheet=[cfg.run.sample_sheet],
            project_id=[pid],
        )
        sample_df, _, _ = cm.get_project_samples_from_samplesheet(args)
        msg += f"<strong>{pid}</strong>: Found {len(sample_df)} samples in samplesheet.\n"

    ssub_df, _ = cm.sample_submission_form_parser(cfg.run.sample_submission_form)
    msg += f"\nFound {len(ssub_df.index)} samples in sample submission form.\n"
    if "Sample_Group" in ssub_df:
        if ssub_df["Sample_Group"].notnull().all():
            unique = ssub_df["Sample_Group"].unique().astype(str)
            msg += f"Sample_Group has {len(unique)} unique values: {', '.join(unique)}.\n"
        elif ssub_df["Sample_Group"].isnull().all():
            msg += "Sample_Group has not been provided.\n"
        else:
            n_missing = ssub_df["Sample_Group"].isnull().values.sum()
            n_groups = len(ssub_df["Sample_Group"].dropna().unique())
            group_names = ", ".join(ssub_df["Sample_Group"].dropna().unique().astype(str))
            msg += f"Sample_Group has {n_groups} unique values: {group_names}.\n"
            grammar = "s" if n_missing > 1 else ""
            msg += f"Missing Sample_Group for {n_missing} sample{grammar}.\n"
    else:
        msg += "Sample_Group has not been provided.\n"

    return msg


def parserDemultiplexStats(cfg):
    """
    Parse DemultiplexingStats.xml under outputDir/Stats/ to get the
    number/percent of undetermined indices.

    In particular, we extract the BarcodeCount values from Project "default"
    Sample "all" and Project "all" Sample "all", as the former gives the total
    undetermined and the later simply the total clusters
    """

    totals = [0, 0, 0, 0, 0, 0, 0, 0]
    undetermined = [0, 0, 0, 0, 0, 0, 0, 0]
    tree = ET.parse(cfg.output_path / "Stats" / "DemultiplexingStats.xml")
    root = tree.getroot()
    for child in root[0].findall("Project"):
        if child.get("name") == "default":
            break
    for sample in child.findall("Sample"):
        if sample.get("name") == "all":
            break
    child = sample[0]  # Get inside Barcode
    for lane in child.findall("Lane"):
        lnum = int(lane.get("number"))
        undetermined[lnum - 1] += int(lane[0].text)

    for child in root[0].findall("Project"):
        if child.get("name") == "all":
            break
    for sample in child.findall("Sample"):
        if sample.get("name") == "all":
            break
    child = sample[0]  # Get Inside Barcode
    for lane in child.findall("Lane"):
        lnum = int(lane.get("number"))
        totals[lnum - 1] += int(lane[0].text)

    lanes = []
    undeter = []
    for i in range(8):
        if totals[i] == 0:
            continue
        lanes.append(i + 1)
        undeter.append(100 * undetermined[i] / totals[i])
        # out_d.append({"Lane": i+1, "Undetermined": 100*undetermined[i]/totals[i]})
    return pd.DataFrame.from_dict({"Lane": lanes, "% Undetermined": undeter}).round(2)


def enoughFreeSpace():
    """
    Ensure that outputDir has at least minSpace gigs
    """
    cfg = PipelineConfig.get()
    (tot, used, free) = shutil.disk_usage(cfg.static.paths.output_dir)
    free_gb = free / (1024**3)
    need = float(cfg.static.system["minspace"])
    log.debug(f"Free GiB in output_dir: {free_gb:.1f} (need ≥ {need:.1f})")
    return free_gb >= need


def errorEmail(errTuple, msg):
    cfg = PipelineConfig.get()
    msg = msg + f"\nError type: {errTuple[0]}\nError value: {errTuple[1]}\n{errTuple[2]}\n"
    (cfg.static.paths.report_dir / f"{cfg.run.run_id}.error").write_text(msg)


def finishedEmail(msg, runTime, extra_html=True):
    cfg = PipelineConfig.get()
    projects = get_project_names(get_project_dirs(cfg))

    message = f"<strong>Short summary for {', '.join(projects)}. </strong>\n\n"
    message += f"<strong>User: {cfg.run.user} </strong>\n" if cfg.run.user != "N/A" else ""
    message += f"Flow cell: {cfg.run.run_id} \n"
    seq = get_sequencer(cfg.run.run_id)
    message += f"Sequencer: {seq} \n"
    read_geo = get_read_geometry(cfg.output_path)
    message += f"Read geometry: {read_geo} \n\n"
    message += f"bcl2fastq_pipeline run time: {runTime} \n"
    # message += "Data transfer: %s\n" % transferTime
    message = message.replace("\n", "\n<br>")
    message += msg

    sample_sheet_metrics = parseSampleSheetMetrics(cfg)
    sample_sheet_metrics = sample_sheet_metrics.replace("\n", "\n<br>")

    message = (
        "<html>\n<body>\n<head>\n"
        + style
        + "\n</head>\n"
        + message
        + "<br>"
        + sample_sheet_metrics
        + "\n</body>\n</html>"
    )

    msg = MIMEMultipart()
    msg["Subject"] = f"[bcl2fastq_pipeline] {', '.join(projects)} processed"
    msg["From"] = cfg.static.email["from_address"]
    msg["To"] = cfg.static.email["finished_to"]
    msg["Date"] = formatdate(localtime=True)

    msg.attach(MIMEText(message, "html"))

    date = cfg.run.run_id.split("_")[0]

    for p in projects:
        with (cfg.output_path / f"multiqc_{p}_{date}.html").open("rb") as report:
            part = MIMEApplication(report.read(), report.name)
        part["Content-Disposition"] = f'attachment; filename="multiqc_{p}_{date}.html"'
        msg.attach(part)

        if (
            cfg.run.libprep.startswith(("10X Genomics Chromium Single Cell", "Parse Biosciences"))
            and extra_html
        ):
            f = cfg.output_path / f"all_samples_web_summary_{p}_{date}.html"
            if f.exists():
                fname = f.name
                with f.open("rb") as report:
                    part = MIMEApplication(report.read(), report.name)
                part["Content-Disposition"] = f'attachment; filename="{fname}"'
                msg.attach(part)
    project_str = "_".join(projects)
    with (cfg.output_path / "Stats" / f"sequencer_stats_{project_str}.html").open("rb") as report:
        part = MIMEApplication(report.read(), report.name)
    part["Content-Disposition"] = f'attachment; filename="sequencer_stats_{project_str}.html"'
    msg.attach(part)

    s = smtplib.SMTP(cfg.static.email["host"])
    s.send_message(msg)
    s.quit()


def finalizedEmail(msg, finalizeTime, runTime):
    cfg = PipelineConfig.get()

    projects = get_project_names(get_project_dirs(cfg))

    message = f"{', '.join(projects)} has been finalized and prepared for delivery.\n\n"
    message += f"md5sum and 7zip runtime: {finalizeTime}\n"
    message += f"Total runtime for bcl2fastq_pipeline: {runTime}\n"
    message += msg

    msg = MIMEMultipart()
    msg["Subject"] = f"[bcl2fastq_pipeline] {', '.join(projects)} finalized"
    msg["From"] = cfg.static.email["from_address"]
    msg["To"] = cfg.static.email["error_to"]
    msg["Date"] = formatdate(localtime=True)

    msg.attach(MIMEText(message))

    s = smtplib.SMTP(cfg.static.email["host"])
    s.send_message(msg)
    s.quit()
