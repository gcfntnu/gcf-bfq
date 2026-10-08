"""Budget analysis report attachments against the complete SMTP message size."""

from __future__ import annotations

import copy
import logging

from dataclasses import dataclass
from html import escape
from pathlib import Path

log = logging.getLogger(__name__)

# A BFQ fallback policy, not a verified production-relay limit. Decimal bytes.
DEFAULT_MAX_MESSAGE_BYTES = 20_000_000


class MessageSizeError(ValueError):
    """Even the summary and report-location notices exceed the message budget."""


@dataclass(frozen=True)
class Report:
    path: Path
    label: str
    # Lower values are omitted first: supplementary, legacy sequencing, MultiQC.
    priority: int


def message_size(message):
    """SMTP octets including CRLF, MIME encoding, headers and body (RFC 1870)."""
    return len(message.as_bytes(policy=message.policy.clone(linesep="\r\n")))


def limit_reports(message, reports, limit):
    """Remove lower-priority reports first; retain every omission in both bodies.

    Within a priority, omit larger encoded parts first, breaking ties by path.
    Every candidate is serialized again with its actual omission notices.
    Reports are never changed on disk. No delivery retry takes place here.
    """
    initial_size = message_size(message)
    if initial_size <= limit:
        return message
    attachments = list(message.iter_attachments())
    if len(attachments) != len(reports):
        raise ValueError("Analysis report attachment inventory does not match message")
    if not attachments:
        raise MessageSizeError(
            f"Analysis summary requires {initial_size} bytes, exceeding the email size limit "
            f"of {limit} bytes even without attachments"
        )
    base = copy.copy(message)
    base.set_payload(
        [copy.deepcopy(part) for part in message.iter_parts() if part not in attachments]
    )
    plain = base.get_body(preferencelist=("plain",)).get_content()
    html = base.get_body(preferencelist=("html",)).get_content()
    omitted = {}

    def render(selected, notices):
        candidate = copy.deepcopy(base)
        if notices:
            candidate.get_body(preferencelist=("plain",)).set_content(
                plain + "\nReport attachment notices:\n\n" + "\n\n".join(notices)
            )
            notice_html = "<h3>Report attachment notices</h3>" + "".join(
                f"<p>{escape(notice)}</p>" for notice in notices
            )
            candidate.get_body(preferencelist=("html",)).set_content(
                html.replace("</body>", notice_html + "</body>"), subtype="html"
            )
        for index in selected:
            candidate.attach(attachments[index])
        return candidate

    order = sorted(
        range(len(reports)),
        key=lambda index: (
            reports[index].priority,
            -message_size(attachments[index]),
            str(reports[index].path),
        ),
    )
    candidate = message
    current_size = initial_size
    for index in order:
        report = reports[index]
        individual_size = message_size(render([index], []))
        reason = (
            "it would exceed the email size limit even with no other report attachments"
            if individual_size > limit
            else "including it would exceed the combined email size limit"
        )
        omitted[index] = (
            f"The {report.label} was not attached because {reason} ({limit} bytes). "
            f"The report remains available at: {report.path}"
        )
        candidate = render(
            [index for index in range(len(reports)) if index not in omitted],
            [omitted[index] for index in sorted(omitted)],
        )
        next_size = message_size(candidate)
        log.warning(
            "Omitted %s: path=%s message_before=%d message_after=%d "
            "individual_message=%d limit=%d bytes; %s",
            report.label,
            report.path,
            current_size,
            next_size,
            individual_size,
            limit,
            reason,
        )
        current_size = next_size
        if current_size <= limit:
            return candidate
    raise MessageSizeError(
        f"Analysis summary and report-location notices require {current_size} bytes, "
        f"exceeding the email size limit of {limit} bytes even without attachments"
    )
