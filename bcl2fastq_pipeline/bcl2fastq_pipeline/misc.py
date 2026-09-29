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

from datetime import UTC, datetime
from email.message import EmailMessage
from email.utils import formatdate, getaddresses
from html import escape

import pandas as pd

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
    """Render current planned metadata without relying on configmaker log files.

    Revalidate the actual bytes on every composition: an edited workbook must
    never inherit a previous success. Invalid inputs produce readable findings
    instead of breaking the notification that is meant to describe them.
    """
    from configmaker.validation import validate_inputs  # noqa: PLC0415

    result = validate_inputs(
        [cfg.run.sample_sheet] if cfg.run.sample_sheet else [],
        [cfg.run.sample_submission_form] if cfg.run.sample_submission_form else [],
    )
    summary = result.to_dict()["summary"]
    lines = ["<strong>Planned input samples (current metadata)</strong>"]
    selected = set(projects) if projects is not None else None
    for project in summary.get("projects", []):
        pid = project["project_id"]
        if selected is None or pid in selected:
            lines.append(
                f"<strong>{escape(str(pid))}</strong>: "
                f"{project['sample_count']} unique samples in SampleSheet."
            )
    lines.append(
        f"{summary.get('planned_sample_count', 0)} unique planned samples "
        f"across {summary.get('samplesheet_rows', 0)} SampleSheet rows."
    )
    lines.append(
        f"{summary.get('submission_sample_count', 0)} samples in effective submission metadata; "
        f"{summary.get('extra_submission_sample_count', 0)} additional submission samples allowed."
    )
    groups = summary.get("sample_groups", {})
    lines.append(
        "<strong>Sample_Group in effective submission metadata (including extras)</strong>"
    )
    values = groups.get("values", [])
    if values:
        lines.append(
            f"Sample_Group has {len(values)} unique values: "
            f"{escape(', '.join(str(value) for value in values))}."
        )
        missing = groups.get("missing_count", 0)
        if missing:
            lines.append(f"Missing Sample_Group for {missing} samples.")
    else:
        lines.append("Sample_Group has not been provided.")
    if not result.ok:
        lines.append("<strong>Current input validation failed; summary may be incomplete.</strong>")
    # Include the complete contextual findings, including warnings on valid pairs.
    lines.append(f"<pre>{escape(result.render_text())}</pre>")
    return "\n".join(lines) + "\n"


def analysisSampleMetrics(cfg, projects):
    """Report observed FASTQ samples independently of the current input plan."""
    lines = ["<strong>Samples discovered in FASTQs at analysis initialization</strong>"]
    for project in projects:
        report = cfg.output_path / f"configmaker-analysis-{project}.json"
        label = escape(str(project))
        try:
            summary = json.loads(report.read_text(encoding="utf-8"))
            if (
                not isinstance(summary, dict)
                or summary.get("kind") != "fastq_discovery"
                or summary.get("schema_version") != 1
            ):
                raise ValueError("Unsupported FASTQ discovery summary")
            count = summary["sample_count"]
            missing = summary.get("missing_sample_ids", [])
            if type(count) is not int or count < 0:
                raise ValueError("Invalid FASTQ discovery sample count")
            if not isinstance(missing, list) or not all(isinstance(sid, str) for sid in missing):
                raise ValueError("Invalid missing sample list in FASTQ discovery summary")
        except (OSError, ValueError, KeyError, TypeError) as error:
            log.info("FASTQ discovery summary unavailable for %s: %s", project, error)
            lines.append(f"{label}: discovery summary unavailable for this analysis.")
            continue
        lines.append(f"{label}: {count} samples discovered in FASTQs.")
        if missing:
            lines.append(f"Planned samples without FASTQs: {escape(', '.join(map(str, missing)))}.")
    return "\n".join(lines) + "\n"


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
        validation = getattr(errTuple[1], "validation_result", None)
        if validation is not None:
            msg += "\n\nComplete structured input preflight report:\n"
            msg += json.dumps(validation.to_dict(), ensure_ascii=False, indent=2)

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
    validation = getattr(error, "validation_result", None)
    if validation is not None:
        report = validation.to_dict()
        # Only stable content contributes; retry time and report path do not.
        details.append(
            {key: report.get(key) for key in ("validator", "inputs", "errors", "warnings", "info")}
        )
    return hashlib.sha256(
        json.dumps(details, ensure_ascii=False, sort_keys=True).encode("utf-8")
    ).hexdigest()


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
    validation = getattr(error, "validation_result", None)
    if validation is not None:
        # Attach the in-memory result, never a stale sidecar from an earlier run.
        message.add_attachment(
            json.dumps(validation.to_dict(), ensure_ascii=False, indent=2).encode("utf-8"),
            maintype="application",
            subtype="json",
            filename=f"{cfg.run.run_id}.input-preflight.json",
        )

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
