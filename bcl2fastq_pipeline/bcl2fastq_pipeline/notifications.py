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
from email.utils import format_datetime
from html import escape
from pathlib import Path
from types import SimpleNamespace

from bcl2fastq_pipeline.state import DeliveryUncertainError

log = logging.getLogger(__name__)
RECIPIENT_KEYS = {"processed": "finished_to", "finalized": "error_to"}


class NotificationConfigError(ValueError):
    """An actionable error in the current email settings."""


def make_payload(  # noqa: PLR0913
    cfg, kind, *, message="", run_time="", finalize_time="", projects=None
):
    """Capture JSON-compatible run context without composing mail or reading inputs."""
    if kind not in RECIPIENT_KEYS:
        raise ValueError(f"Unsupported completion notification kind: {kind}")
    if projects is None:
        from bcl2fastq_pipeline import afterFastq  # noqa: PLC0415

        projects = afterFastq.get_project_names(afterFastq.get_project_dirs(cfg))
    run = cfg.run
    return {
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
        supported = {"host", "from_address", *RECIPIENT_KEYS.values()}
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


def _processed_message(cfg, message, payload):
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
    date = run_id.split("_")[0]
    for project in projects:
        _add_report(message, saved_cfg.output_path / f"multiqc_{project}_{date}.html")
        if (payload.get("libprep") or "").startswith(
            ("10X Genomics Chromium Single Cell", "Parse Biosciences")
        ):
            report = saved_cfg.output_path / f"all_samples_web_summary_{project}_{date}.html"
            if report.exists():
                _add_report(message, report)
    project_names = "_".join(projects)
    _add_report(message, saved_cfg.output_path / "Stats" / f"sequencer_stats_{project_names}.html")


def build_message(cfg, entry):
    """Read reports and compose mail using saved run context and current recipients."""
    kind = entry["kind"]
    _host, sender_header, _sender, recipient_header, _recipients = _settings(cfg, kind)
    payload = entry["payload"]
    projects = ", ".join(payload["projects"])
    message = EmailMessage()
    message["Subject"] = f"[bcl2fastq_pipeline] {projects} {kind}"
    message["From"] = sender_header
    message["To"] = recipient_header
    message["Date"] = format_datetime(datetime.fromisoformat(entry["created_at"]))
    identity = "\0".join((payload["run_id"], entry["id"], entry["created_at"]))
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()
    message["Message-ID"] = f"<{digest}@bfq.invalid>"
    if kind == "processed":
        _processed_message(cfg, message, payload)
    else:
        message.set_content(
            f"{projects} has been finalized and prepared for delivery.\n\n"
            f"Flow cell: {payload['run_id']}\n"
            f"md5sum and 7zip runtime: {payload['finalize_time']}\n"
            f"Total runtime for bcl2fastq_pipeline: {payload['run_time']}\n"
            f"{payload['message']}"
        )
    return message


def send_notification(cfg, entry):
    """Attempt delivery once; signal ambiguous acceptance for an explicit retry."""
    message = build_message(cfg, entry)
    host, _sender_header, sender, _recipient_header, recipients = _settings(
        cfg, entry["kind"], warn_unknown=False
    )
    # Connection failures happen before sending and can safely be retried.
    smtp = smtplib.SMTP(host, timeout=30)
    try:
        try:
            refused = smtp.send_message(message, from_addr=sender, to_addrs=recipients)
        except (
            smtplib.SMTPRecipientsRefused,
            smtplib.SMTPSenderRefused,
            smtplib.SMTPDataError,
            smtplib.SMTPHeloError,
            smtplib.SMTPNotSupportedError,
        ):
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
