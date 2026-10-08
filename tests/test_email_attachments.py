"""Complete MIME budgets, deterministic priority and accessible omission notices."""

import copy

from email.message import EmailMessage
from email.policy import SMTP
from html import escape

import pytest

from bcl2fastq_pipeline.email_attachments import (
    MessageSizeError,
    Report,
    limit_reports,
    message_size,
)


def report_message(tmp_path, sizes):
    message = EmailMessage(policy=SMTP)
    message["Subject"] = "Analysis complete — QC summary"
    message["From"] = "bfq@example.org"
    message["To"] = "analyst@example.org"
    message.set_content("Analysis summary: read retention 95%.")
    message.add_alternative(
        "<html><body>Analysis summary: read retention 95%.</body></html>", subtype="html"
    )
    reports = []
    for name, size, priority in sizes:
        path = tmp_path / name
        path.write_bytes(b"x" * size)
        label = "MultiQC report" if priority == 2 else "additional single-cell HTML report"
        reports.append(Report(path, label, priority))
        message.add_attachment(path.read_bytes(), maintype="text", subtype="html", filename=name)
    # Freeze MIME boundaries so tests compare the exact wire message.
    message_size(message)
    return message, reports


def filenames(message):
    return [part.get_filename() for part in message.iter_attachments()]


def assert_notice(message, path, reason):
    for kind in ("plain", "html"):
        body = message.get_body(preferencelist=(kind,)).get_content()
        assert (str(path) if kind == "plain" else escape(str(path))) in body
        assert reason in body
        assert "report remains available at:" in body


def test_exact_boundary_counts_mime_and_crlf(tmp_path):
    message, reports = report_message(tmp_path, [("multiqc.html", 4_000, 2)])
    size = message_size(message)
    assert size > 4_000 * 4 / 3
    assert size > len(message.as_bytes(policy=SMTP.clone(linesep="\n")))
    assert filenames(limit_reports(message, reports, size)) == ["multiqc.html"]
    reduced = limit_reports(message, reports, size - 1)
    assert filenames(reduced) == []
    assert message_size(reduced) <= size - 1
    assert_notice(reduced, reports[0].path, "even with no other report attachments")


def test_supplementary_is_omitted_before_larger_multiqc(tmp_path):
    message, reports = report_message(
        tmp_path, [("multiqc.html", 10_000, 2), ("single-cell.html", 4_000, 0)]
    )
    limit = message_size(message) - 1_000
    reduced = limit_reports(message, reports, limit)
    assert filenames(reduced) == ["multiqc.html"]
    assert message_size(reduced) <= limit
    assert_notice(reduced, reports[1].path, "combined email size limit")
    for kind in ("plain", "html"):
        assert str(reports[0].path) not in reduced.get_body(preferencelist=(kind,)).get_content()


def test_all_reports_can_be_omitted_without_losing_summary_or_files(tmp_path, caplog):
    message, reports = report_message(
        tmp_path, [("multiqc.html", 12_000, 2), ("single-cell.html", 8_000, 0)]
    )
    before = {
        report.path: (report.path.read_bytes(), report.path.stat().st_mtime_ns)
        for report in reports
    }
    reduced = limit_reports(message, reports, 6_000)
    assert filenames(reduced) == []
    assert message_size(reduced) <= 6_000
    for report in reports:
        assert_notice(reduced, report.path, "even with no other report attachments")
        assert (report.path.read_bytes(), report.path.stat().st_mtime_ns) == before[report.path]
        assert str(report.path) in caplog.text
    for kind in ("plain", "html"):
        assert "read retention 95%" in reduced.get_body(preferencelist=(kind,)).get_content()
    assert "individual_message=" in caplog.text
    assert "limit=6000 bytes" in caplog.text


def test_multiple_projects_largest_first_and_path_tiebreak_are_repeatable(tmp_path):
    message, reports = report_message(
        tmp_path,
        [
            ("multiqc_b.html", 4_000, 2),
            ("multiqc_a.html", 4_000, 2),
            ("multiqc_large.html", 9_000, 2),
        ],
    )
    reduced = limit_reports(message, reports, 12_000)
    assert filenames(reduced) == ["multiqc_b.html"]
    assert_notice(reduced, reports[1].path, "combined email size limit")
    assert_notice(reduced, reports[2].path, "even with no other report attachments")
    assert message_size(reduced) <= 12_000
    again = limit_reports(copy.deepcopy(message), reports, 12_000)
    assert reduced.as_bytes() == again.as_bytes()


def test_final_notice_bytes_can_require_another_omission(tmp_path):
    message, reports = report_message(
        tmp_path, [("multiqc.html", 4_000, 2), ("single-cell.html", 4_000, 0)]
    )
    without_supplementary = copy.deepcopy(message)
    without_supplementary.set_payload(without_supplementary.get_payload()[:-1])
    # MultiQC alone fits, but the necessary notice for the supplementary report does not.
    limit = message_size(without_supplementary)
    reduced = limit_reports(message, reports, limit)
    assert filenames(reduced) == []
    assert message_size(reduced) <= limit
    for report in reports:
        assert_notice(reduced, report.path, "combined email size limit")


def test_notice_paths_are_escaped_in_html_and_preserved_in_plain_text(tmp_path):
    message, reports = report_message(tmp_path, [("multiqc_<sample>&ø.html", 8_000, 2)])
    reduced = limit_reports(message, reports, 5_000)
    assert_notice(reduced, reports[0].path, "even with no other report attachments")
    html = reduced.get_body(preferencelist=("html",)).get_content()
    assert "<sample>" not in html


def test_summary_too_large_is_an_actionable_failure(tmp_path):
    message, reports = report_message(tmp_path, [("multiqc.html", 4_000, 2)])
    with pytest.raises(MessageSizeError, match="even without attachments"):
        limit_reports(message, reports, 100)


def test_legacy_sequencing_report_is_omitted_before_multiqc(tmp_path):
    message, reports = report_message(
        tmp_path,
        [("multiqc.html", 8_000, 2), ("sequencing.html", 4_000, 1), ("single-cell.html", 4_000, 0)],
    )
    reduced = limit_reports(message, reports, 15_000)
    assert filenames(reduced) == ["multiqc.html"]
    assert message_size(reduced) <= 15_000
