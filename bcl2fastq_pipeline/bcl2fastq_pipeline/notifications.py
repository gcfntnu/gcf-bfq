"""Compose completion mail from durable context and deliver with current settings.

State transitions and retry policy live in ``state``. This module performs one
delivery attempt only and never changes processing state or the global run context.
"""

from __future__ import annotations

import hashlib
import logging
import smtplib

from datetime import datetime
from email.headerregistry import HeaderRegistry
from email.message import EmailMessage
from email.policy import SMTP
from email.utils import format_datetime
from html import escape
from html.parser import HTMLParser
from pathlib import Path
from types import SimpleNamespace

from bcl2fastq_pipeline import email_attachments, processing_times
from bcl2fastq_pipeline.state import DeliveryUncertainError

log = logging.getLogger(__name__)
RECIPIENT_KEYS = {
    "sequencing": "finished_to",
    "processed": "finished_to",
    "finalized": "error_to",
}


class NotificationConfigError(ValueError):
    """An actionable error in the current email settings."""


def make_payload(  # noqa: PLR0913
    cfg,
    kind,
    *,
    message="",
    run_time="",
    finalize_time="",
    projects=None,
    sequencing_qc=None,
    analysis_qc=None,
    processing_timing=None,
):
    """Capture JSON-compatible run context without composing mail or reading inputs."""
    if kind not in RECIPIENT_KEYS:
        raise ValueError(f"Unsupported completion notification kind: {kind}")
    if kind == "sequencing" and sequencing_qc is None and projects is None:
        raise ValueError("Sequencing notification requires planned projects or saved QC")
    if projects is None and sequencing_qc is not None:
        projects = sequencing_qc["projects"]
    if projects is None:
        from bcl2fastq_pipeline import afterFastq  # noqa: PLC0415

        projects = afterFastq.get_project_names(afterFastq.get_project_dirs(cfg))
    run = cfg.run
    payload = {
        "email_version": 2,
        "run_id": run.run_id,
        "output_path": str(cfg.output_path),
        "flowcell_path": str(run.flowcell_path) if run.flowcell_path else None,
        "sample_sheet": str(run.sample_sheet) if run.sample_sheet else None,
        "sample_submission_form": (
            str(run.sample_submission_form) if run.sample_submission_form else None
        ),
        "user": run.user,
        "libprep": run.libprep,
        "custom": dict(run.custom),
        "projects": sorted(set(projects)),
        "message": str(message or ""),
        "run_time": str(run_time),
        "finalize_time": str(finalize_time),
    }
    if processing_timing is not None:
        payload["processing_timing"] = processing_timing
        payload["run_time"] = processing_times.format_duration(processing_timing["total_seconds"])
    if sequencing_qc is not None:
        payload["sequencing_qc"] = sequencing_qc
    if analysis_qc is not None:
        payload["analysis_qc"] = analysis_qc
    return payload


def _required(settings, key):
    value = settings.get(key)
    if not isinstance(value, str) or not value.strip():
        aliases = [
            name for name in settings if name.replace("_", "").lower() == key.replace("_", "")
        ]
        hint = f"; rename {aliases[0]!r} to {key!r}" if aliases else ""
        raise NotificationConfigError(f"[Email] {key} is required{hint}")
    if any(character in value for character in "\r\n\x00"):
        raise NotificationConfigError(f"[Email] {key} must not contain control characters")
    return value.strip()


def _addresses(settings, key):
    value = _required(settings, key)
    try:
        parsed = HeaderRegistry()("To", value)
        addresses = list(dict.fromkeys(address.addr_spec for address in parsed.addresses))
        valid = (
            not parsed.defects
            and addresses
            and all(
                address.username
                and address.domain
                and not any(character.isspace() for character in address.addr_spec)
                for address in parsed.addresses
            )
        )
    except (ValueError, IndexError) as error:
        raise NotificationConfigError(f"[Email] {key} contains an invalid email address") from error
    if not valid:
        raise NotificationConfigError(f"[Email] {key} must contain comma-separated email addresses")
    return value, addresses


def _settings(cfg, kind, *, warn_unknown=True):
    if kind not in RECIPIENT_KEYS:
        raise ValueError(f"Unsupported completion notification kind: {kind}")
    settings = cfg.static.email
    if warn_unknown:
        supported = {"host", "from_address", "max_message_bytes", *RECIPIENT_KEYS.values()}
        for key in sorted(set(settings) - supported):
            correction = next(
                (
                    name
                    for name in supported
                    if name.replace("_", "") == key.replace("_", "").lower()
                ),
                None,
            )
            hint = f"; use {correction!r}" if correction else ""
            log.warning("Unknown [Email] option %r is ignored%s", key, hint)
    host = _required(settings, "host")
    if any(character.isspace() for character in host):
        raise NotificationConfigError("[Email] host must be a hostname without whitespace")
    sender_header, senders = _addresses(settings, "from_address")
    if len(senders) != 1:
        raise NotificationConfigError("[Email] from_address must contain exactly one address")
    recipient_header, recipients = _addresses(settings, RECIPIENT_KEYS[kind])
    return host, sender_header, senders[0], recipient_header, recipients


def _saved_config(cfg, payload):
    """Provide parser context without changing the active daemon's run singleton."""
    run = SimpleNamespace(
        **{key: payload.get(key) for key in ("run_id", "user", "libprep", "custom")},
        **{
            key: Path(payload[key]) if payload.get(key) else None
            for key in ("flowcell_path", "sample_sheet", "sample_submission_form")
        },
    )
    return SimpleNamespace(static=cfg.static, run=run, output_path=Path(payload["output_path"]))


def _add_report(message, report_path):
    message.add_attachment(
        report_path.read_bytes(), maintype="text", subtype="html", filename=report_path.name
    )


def _legacy_processed_message(cfg, message, payload):
    from bcl2fastq_pipeline import afterFastq, misc  # noqa: PLC0415

    saved_cfg = _saved_config(cfg, payload)
    projects = payload["projects"]
    run_id = payload["run_id"]
    summary = [f"Short summary for {', '.join(projects)}."]
    if payload.get("user") not in (None, "", "N/A"):
        summary.append(f"User: {payload['user']}")
    summary.extend(
        (
            f"Flow cell: {run_id}",
            f"Sequencer: {afterFastq.get_sequencer(run_id)}",
            f"Read geometry: {afterFastq.get_read_geometry(saved_cfg.output_path)}",
            f"bcl2fastq_pipeline run time: {payload['run_time']}",
        )
    )
    text_summary = "\n".join(summary)
    message.set_content(text_summary)
    metadata = misc.parseSampleSheetMetrics(saved_cfg, projects=projects).replace("\n", "\n<br>")
    discovery = misc.analysisSampleMetrics(saved_cfg, projects).replace("\n", "\n<br>")
    metrics = misc.getFCmetricsImproved(saved_cfg)
    disk_usage = afterFastq._disk_usage_message(saved_cfg)
    message.add_alternative(
        "<html><head>"
        + misc.style
        + "</head><body>"
        + escape(text_summary).replace("\n", "\n<br>")
        + "<br>"
        + payload["message"]
        + metrics
        + "<br>"
        + metadata
        + "<br>"
        + discovery
        + "<br>"
        + disk_usage
        + "</body></html>",
        subtype="html",
    )
    for report in _analysis_reports(payload):
        _add_report(message, report.path)


class _PlainText(HTMLParser):
    """Keep the full operational summary accessible in plain-text mail clients."""

    def __init__(self):
        super().__init__()
        self.parts = []

    def handle_starttag(self, tag, attrs):
        if tag in {"br", "p", "tr", "div", "h2", "h3", "pre"}:
            self.parts.append("\n")
        elif tag in {"td", "th"}:
            self.parts.append("\t")

    def handle_endtag(self, tag):
        if tag in {"p", "tr", "div", "h2", "h3", "pre"}:
            self.parts.append("\n")

    def handle_data(self, data):
        self.parts.append(data)


def _plain_text(html):
    parser = _PlainText()
    parser.feed(html)
    return "\n".join(line.strip() for line in "".join(parser.parts).splitlines() if line.strip())


def _run_summary(payload):
    from bcl2fastq_pipeline import afterFastq  # noqa: PLC0415

    lines = [f"Projects: {', '.join(payload['projects']) or '(none)'}"]
    if payload.get("user") not in (None, "", "N/A"):
        lines.append(f"User: {payload['user']}")
    lines.extend(
        (
            f"Flow cell: {payload['run_id']}",
            f"Sequencer: {afterFastq.get_sequencer(payload['run_id'])}",
        )
    )
    if payload.get("processing_timing") is None:
        lines.append(f"bcl2fastq_pipeline elapsed time: {payload['run_time']}")
    return "\n".join(lines)


def _body(message, heading, summary, html_sections, text_sections=None):
    from bcl2fastq_pipeline import misc  # noqa: PLC0415

    if text_sections is None:
        text_sections = [_plain_text(section) for section in html_sections]
    message.set_content("\n\n".join([heading, summary, *text_sections]))
    message.add_alternative(
        "<html><head>"
        + misc.style
        + "</head><body>"
        + f"<h2>{escape(heading)}</h2>"
        + escape(summary).replace("\n", "<br>\n")
        + "<br><br>"
        + "<br><br>".join(html_sections)
        + "</body></html>",
        subtype="html",
    )


def _sequencing_message(cfg, message, payload):
    """Use only the saved early QC result and operational disk availability."""
    from bcl2fastq_pipeline import afterFastq  # noqa: PLC0415

    saved_cfg = _saved_config(cfg, payload)
    qc = payload["sequencing_qc"]
    disk = afterFastq._disk_usage_message(saved_cfg)
    explanation = (
        "Demultiplexing is complete. This report summarizes sequencing and index assignment. "
        "Analysis is a separate stage; no analysis results or automatic QC pass/fail decision "
        "are included. Review the attached sequencing report if intervention is needed."
    )
    html_sections = [escape(explanation), qc["summary_html"], disk]
    text_sections = [explanation, qc["summary_text"], _plain_text(disk)]
    if payload.get("message"):
        html_sections.append(payload["message"])
        text_sections.append(_plain_text(payload["message"]))
    _body(
        message,
        "Demultiplexing complete — sequencing QC",
        _run_summary(payload),
        html_sections,
        text_sections,
    )
    _add_report(message, Path(qc["report_path"]))


def _processed_message(cfg, message, payload):
    """Describe completed analysis; sequencing QC has its own earlier notification."""
    from bcl2fastq_pipeline import misc  # noqa: PLC0415

    saved_cfg = _saved_config(cfg, payload)
    projects = payload["projects"]
    html_sections = [
        "Analysis and project QC reports are complete. Review the MultiQC reports "
        "for detailed sample QC; this email does not assign a QC pass/fail status. "
        "Archiving and delivery preparation are reported separately after finalization."
    ]
    if payload.get("analysis_qc"):
        html_sections.append(payload["analysis_qc"]["summary_html"])
    else:
        html_sections.append(
            misc.analysisSampleMetrics(saved_cfg, projects).replace("\n", "<br>\n")
        )
    html_sections.append(
        misc.parseSampleSheetMetrics(saved_cfg, projects=projects).replace("\n", "<br>\n")
    )
    if payload.get("message"):
        html_sections.append(payload["message"])
    legacy_sequencing = payload.get("sequencing_qc")
    if legacy_sequencing:
        html_sections.extend(
            (
                "Sequencing QC is included here for this legacy run; "
                "no separate early notification was created.",
                legacy_sequencing["summary_html"],
            )
        )
    else:
        html_sections.append(
            "Sequencing yield, base quality, PhiX and index assignment are covered in the "
            "early sequencing report."
        )
    if payload.get("processing_timing") is not None:
        html_sections.append(
            escape(processing_times.format_summary(payload["processing_timing"])).replace(
                "\n", "<br>\n"
            )
        )
    _body(message, "Analysis complete — QC summary", _run_summary(payload), html_sections)
    for report in _analysis_reports(payload):
        _add_report(message, report.path)


def _analysis_reports(payload):
    output = Path(payload["output_path"])
    date = payload["run_id"].split("_")[0]
    projects = sorted(set(payload["projects"]))
    reports = []
    for project in projects:
        reports.append(
            email_attachments.Report(output / f"multiqc_{project}_{date}.html", "MultiQC report", 2)
        )
        if (payload.get("libprep") or "").startswith(
            ("10X Genomics Chromium Single Cell", "Parse Biosciences")
        ):
            report = output / f"all_samples_web_summary_{project}_{date}.html"
            if report.exists():
                reports.append(
                    email_attachments.Report(report, "additional single-cell HTML report", 0)
                )
    if payload.get("email_version", 1) < 2:
        report = output / "Stats" / f"sequencer_stats_{'_'.join(payload['projects'])}.html"
        reports.append(email_attachments.Report(report, "sequencing HTML report", 1))
    elif payload.get("sequencing_qc"):
        reports.append(
            email_attachments.Report(
                Path(payload["sequencing_qc"]["report_path"]), "sequencing HTML report", 1
            )
        )
    return reports


def _message_limit(cfg):
    value = cfg.static.email.get(
        "max_message_bytes", str(email_attachments.DEFAULT_MAX_MESSAGE_BYTES)
    )
    if not isinstance(value, str) or not value.strip().isascii() or not value.strip().isdecimal():
        raise NotificationConfigError(
            "[Email] max_message_bytes must be a positive integer in bytes"
        )
    limit = int(value)
    if limit <= 0:
        raise NotificationConfigError(
            "[Email] max_message_bytes must be a positive integer in bytes"
        )
    return limit


def _relay_limit(smtp, configured):
    # Negotiation is before DATA: failures here are safe to retry.
    smtp.ehlo_or_helo_if_needed()
    advertised = smtp.esmtp_features.get("size", "")
    limit = configured
    if isinstance(advertised, str) and advertised.isascii() and advertised.isdecimal():
        if int(advertised) > 0:
            limit = min(configured, int(advertised))
    elif isinstance(advertised, str) and advertised:
        log.warning("SMTP advertised an invalid SIZE limit; using configured budget")
    log.info(
        "Analysis email budget: configured=%d relay_SIZE=%s effective=%d bytes",
        configured,
        advertised if isinstance(advertised, str) else "unavailable",
        limit,
    )
    return limit


def build_message(cfg, entry, *, apply_size_limit=True):
    """Read reports and compose mail using saved run context and current recipients."""
    kind = entry["kind"]
    _host, sender_header, sender, recipient_header, recipients = _settings(cfg, kind)
    limit = _message_limit(cfg) if kind == "processed" else None
    payload = entry["payload"]
    projects = ", ".join(payload["projects"])
    international = any(not address.isascii() for address in (sender, *recipients))
    message = EmailMessage(policy=SMTP.clone(utf8=international))
    label = {
        "sequencing": "Demultiplexing complete — sequencing QC",
        "processed": "Analysis complete — QC summary",
    }.get(kind, kind)
    if kind == "processed" and payload.get("email_version", 1) < 2:
        label = kind
    message["Subject"] = f"[bcl2fastq_pipeline] {projects} {label}"
    message["From"] = sender_header
    message["To"] = recipient_header
    message["Date"] = format_datetime(datetime.fromisoformat(entry["created_at"]))
    identity = "\0".join((payload["run_id"], entry["id"], entry["created_at"]))
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()
    message["Message-ID"] = f"<{digest}@bfq.invalid>"
    if kind == "sequencing":
        _sequencing_message(cfg, message, payload)
    elif kind == "processed":
        if payload.get("email_version", 1) < 2:
            _legacy_processed_message(cfg, message, payload)
        else:
            _processed_message(cfg, message, payload)
    else:
        timing = (
            processing_times.format_summary(payload["processing_timing"]) + "\n"
            if payload.get("processing_timing") is not None
            else (
                f"md5sum and 7zip runtime: {payload['finalize_time']}\n"
                f"Total runtime for bcl2fastq_pipeline: {payload['run_time']}\n"
            )
        )
        message.set_content(
            f"{projects} has been finalized and prepared for delivery.\n\n"
            f"Flow cell: {payload['run_id']}\n"
            f"{timing}"
            f"{payload['message']}"
        )
    if kind == "processed" and apply_size_limit:
        message = email_attachments.limit_reports(message, _analysis_reports(payload), limit)
    return message


def send_notification(cfg, entry):
    """Attempt delivery once; signal ambiguous acceptance for an explicit retry."""
    message = build_message(cfg, entry, apply_size_limit=False)
    host, _sender_header, sender, _recipient_header, recipients = _settings(
        cfg, entry["kind"], warn_unknown=False
    )
    # Connection failures happen before sending and can safely be retried.
    smtp = smtplib.SMTP(host, timeout=30)
    try:
        if entry["kind"] == "processed":
            limit = _relay_limit(smtp, _message_limit(cfg))
            message = email_attachments.limit_reports(
                message, _analysis_reports(entry["payload"]), limit
            )
            log.info(
                "Analysis email serialized size=%d limit=%d bytes",
                email_attachments.message_size(message),
                limit,
            )
        try:
            refused = smtp.send_message(message, from_addr=sender, to_addrs=recipients)
        except (
            smtplib.SMTPRecipientsRefused,
            smtplib.SMTPSenderRefused,
            smtplib.SMTPDataError,
            smtplib.SMTPHeloError,
            smtplib.SMTPNotSupportedError,
        ) as error:
            if isinstance(error, smtplib.SMTPResponseException):
                log.error(
                    "SMTP rejected notification: code=%d response=%r",
                    error.smtp_code,
                    error.smtp_error,
                )
            raise
        except (OSError, smtplib.SMTPException) as error:
            raise DeliveryUncertainError(
                "SMTP delivery was interrupted; the relay may have accepted the message. "
                "An explicit retry may send a duplicate."
            ) from error
        if refused:
            if set(refused) >= set(recipients):
                raise smtplib.SMTPRecipientsRefused(refused)
            raise DeliveryUncertainError(
                "SMTP accepted the message for some recipients but rejected others "
                f"({', '.join(refused)}). An explicit retry may send duplicates."
            )
    finally:
        try:
            smtp.quit()
        except Exception:
            # DATA acceptance determines success. QUIT cannot undo delivery.
            log.warning("Could not close SMTP session cleanly", exc_info=True)
            try:
                smtp.close()
            except Exception:
                log.warning("Could not close SMTP connection", exc_info=True)
