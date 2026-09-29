"""
Misc. functions
"""

import hashlib
import json
import logging
import os
import shutil
import smtplib
import socket
import tempfile as tmp
import traceback
import xml.etree.ElementTree as ET

from argparse import Namespace
from datetime import UTC, datetime
from email.message import EmailMessage
from email.utils import formatdate, getaddresses

import configmaker.configmaker as cm
import pandas as pd

from bcl2fastq_pipeline.afterFastq import (
    get_project_dirs,
    get_project_names,
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


def getFCmetricsImproved(cfg=None):
    cfg = cfg if cfg is not None else PipelineConfig.get()
    message = ""
    try:
        with (cfg.output_path / "Stats" / "interop_summary.csv").open() as fh:
            header = False
            while not header:
                line = fh.readline()
                if not line:
                    raise ValueError("Missing table header in Stats/interop_summary.csv")
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


def parseSampleSheetMetrics(cfg, projects=None):
    project_names = projects if projects is not None else get_project_names(get_project_dirs(cfg))
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


def write_error_report(errTuple, msg):
    cfg = PipelineConfig.get()
    report_dir = cfg.static.paths.report_dir
    report_dir.mkdir(parents=True, exist_ok=True)

    if errTuple and errTuple[0] is not None:
        msg = f"{msg}\n\n{''.join(traceback.format_exception(*errTuple))}"
        command_output = getattr(errTuple[1], "output", None)
        if command_output:
            if isinstance(command_output, bytes):
                command_output = command_output.decode("utf-8", errors="replace")
            msg += f"\nCaptured command output (last 400 lines):\n{command_output}"

    report_path = report_dir / f"{cfg.run.run_id}.error"
    report_path.write_text(msg, encoding="utf-8")
    return report_path


def errorEmail(errTuple, msg):
    """Compatibility alias: write a report only; never send email."""
    return write_error_report(errTuple, msg)


def error_failure_signature(stage, error_info, message):
    """Identify failures without volatile report paths, timestamps or traceback lines."""
    error = error_info[1]
    error_type = type(error)
    details = [stage, f"{error_type.__module__}.{error_type.__qualname__}", str(error)]
    if error is None:
        details.append(message)
    output = getattr(error, "output", None)
    if isinstance(output, bytes):
        output = output.decode("utf-8", errors="replace")
    if output:
        details.append(output)
    return hashlib.sha256(json.dumps(details, ensure_ascii=False).encode("utf-8")).hexdigest()


def send_error_report(cfg, report_path, stage, error_info, store, signature):  # noqa: PLR0913
    """Deliver a saved report in production, using persistent duplicate protection."""
    log = logging.getLogger(__name__)
    if os.environ.get("BFQ_ENV") != "production":
        log.info("Error email suppressed for %s: BFQ_ENV is not production", cfg.run.run_id)
        return False
    if store is None:
        log.error("Error email suppressed: durable flowcell state is unavailable")
        return False
    recipients = list(
        dict.fromkeys(
            address
            for _name, address in getaddresses([cfg.static.email.get("error_to", "")])
            if address
        )
    )
    if not recipients:
        log.error("Error email suppressed for %s: error_to has no recipients", cfg.run.run_id)
        return False

    # Read and construct the message before claiming delivery. No report, no mail.
    report = report_path.read_text(encoding="utf-8")
    message = EmailMessage()
    message["Subject"] = f"[BFQ ERROR] {cfg.run.run_id} — {stage}"
    message["From"] = cfg.static.email["from_address"]
    message["To"] = ", ".join(recipients)
    message["Date"] = formatdate(localtime=True)
    error = error_info[1]
    message.set_content(
        f"Flowcell: {cfg.run.run_id}\nFailure stage: {stage}\n"
        f"Exception: {type(error).__name__}: {error}\n"
        f"Timestamp: {datetime.now(UTC).isoformat()}\n"
        f"Host: {socket.gethostname()}\nReport: {report_path.resolve()}\n"
    )
    message.add_attachment(report, subtype="plain", filename=report_path.name)

    def deliver():
        with smtplib.SMTP(cfg.static.email["host"], timeout=30) as smtp:
            refused = smtp.send_message(
                message, from_addr=cfg.static.email["from_address"], to_addrs=recipients
            )
            if refused:
                raise smtplib.SMTPRecipientsRefused(refused)

    sent = store.deliver_failure_notification(cfg.run.run_id, signature, deliver)
    log.info("Error email %s for %s", "sent" if sent else "duplicate suppressed", cfg.run.run_id)
    return sent
