"""Completion mail composition and SMTP failure boundaries, without real delivery."""

import json
import smtplib

from datetime import timedelta
from pathlib import Path
from smtplib import SMTP as SMTPClient
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from bcl2fastq_pipeline.config import RunContext
from bcl2fastq_pipeline.state import DeliveryUncertainError

from bcl2fastq_pipeline import afterFastq, email_attachments, misc, notifications


@pytest.fixture
def mail_cfg(tmp_path):
    output = tmp_path / "output" / "260918_MN00686_0026_A000HCMFHF"
    output.mkdir(parents=True)
    return SimpleNamespace(
        output_path=output,
        run=RunContext(
            run_id=output.name,
            flowcell_path=tmp_path / "instrument" / output.name,
            sample_sheet=output / "SampleSheet.csv",
            sample_submission_form=output / "Sample-Submission-Form.xlsx",
            user="Saved User",
            libprep="Illumina DNA Prep",
            custom={"User": "Saved User"},
        ),
        static=SimpleNamespace(
            email={
                "host": "relay.example.org",
                "from_address": "BFQ <bfq@example.org>",
                "finished_to": "Analyst <analyst@example.org>, other@example.org",
                "error_to": "delivery@example.org",
            },
            paths=SimpleNamespace(output_dir=output.parent),
        ),
    )


def entry_for(cfg, kind="finalized"):
    return {
        "id": f"{kind}:1",
        "kind": kind,
        "created_at": "2026-09-29T07:00:00+00:00",
        "payload": notifications.make_payload(
            cfg,
            kind,
            message="Saved operational summary",
            run_time=timedelta(hours=2),
            finalize_time=timedelta(minutes=5),
            projects=["GCF-2026-043"],
        ),
    }


@pytest.fixture
def smtp(monkeypatch):
    client = Mock()
    client.send_message.return_value = {}
    factory = Mock(return_value=client)
    monkeypatch.setattr(notifications.smtplib, "SMTP", factory)
    return factory, client


@pytest.fixture
def processed(mail_cfg, monkeypatch):
    entry = entry_for(mail_cfg, "processed")
    # Pending notifications from before the split retain their original contents.
    entry["payload"].pop("email_version")
    project = "GCF-2026-043"
    (mail_cfg.output_path / "Stats").mkdir(exist_ok=True)
    (mail_cfg.output_path / f"multiqc_{project}_260918.html").write_text("<html>Analysis</html>")
    (mail_cfg.output_path / "Stats" / f"sequencer_stats_{project}.html").write_text(
        "<html>Sequencer</html>"
    )
    monkeypatch.setattr(afterFastq, "get_read_geometry", lambda path: "2x150")
    monkeypatch.setattr(misc, "parseSampleSheetMetrics", lambda cfg, projects=None: "Sample groups")
    monkeypatch.setattr(misc, "getFCmetricsImproved", lambda cfg: "Flowcell metrics")
    monkeypatch.setattr(afterFastq, "_disk_usage_message", lambda cfg: "Disk summary")
    return entry


def test_payload_is_serializable_without_opening_inputs(mail_cfg, monkeypatch):
    monkeypatch.setattr(Path, "open", Mock(side_effect=AssertionError("Input read during capture")))
    entry = entry_for(mail_cfg)
    restored = json.loads(json.dumps(entry))
    assert restored["payload"]["run_time"] == "2:00:00"
    assert restored["payload"]["sample_sheet"] == str(mail_cfg.run.sample_sheet)


def test_finalized_mail_replays_saved_context_using_current_settings(mail_cfg):
    entry = entry_for(mail_cfg)
    mail_cfg.run = RunContext(run_id="different_active_flowcell")
    mail_cfg.static.email["error_to"] = "Corrected <corrected@example.org>"
    message = notifications.build_message(mail_cfg, entry)
    assert "GCF-2026-043 finalized" in message["Subject"]
    assert "corrected@example.org" in message["To"]
    assert "260918_MN00686_0026_A000HCMFHF" in message.get_content()
    assert "0:05:00" in message.get_content()
    assert "2:00:00" in message.get_content()
    assert "Saved operational summary" in message.get_content()
    assert mail_cfg.run.run_id == "different_active_flowcell"


def test_headers_are_stable_for_retries_but_distinct_for_new_events(mail_cfg):
    entry = entry_for(mail_cfg)
    first = notifications.build_message(mail_cfg, entry)
    second = notifications.build_message(mail_cfg, entry)
    assert first["Message-ID"] == second["Message-ID"]
    assert first["Date"] == second["Date"]
    entry["id"] = "finalized:2"
    assert notifications.build_message(mail_cfg, entry)["Message-ID"] != first["Message-ID"]


def test_processed_mail_preserves_metrics_and_attachments(mail_cfg, processed, monkeypatch):
    original_output = mail_cfg.output_path
    expected_sheet = mail_cfg.run.sample_sheet
    mail_cfg.run = RunContext(run_id="unrelated_flowcell")
    mail_cfg.output_path = Path("/unrelated/output")

    def sample_metrics(saved_cfg, projects=None):
        assert projects == ["GCF-2026-043"]
        assert saved_cfg.output_path == original_output
        assert saved_cfg.run.sample_sheet == expected_sheet
        assert saved_cfg.run.user == "Saved User"
        return "Sample groups"

    monkeypatch.setattr(misc, "parseSampleSheetMetrics", sample_metrics)
    message = notifications.build_message(mail_cfg, processed)
    assert message["Subject"] == "[bcl2fastq_pipeline] GCF-2026-043 processed"
    html = message.get_body(preferencelist=("html",)).get_content()
    for text in ("Sample groups", "Flowcell metrics", "Disk summary", "Saved User", "2x150"):
        assert text in html
    attachments = list(message.iter_attachments())
    assert [part.get_filename() for part in attachments] == [
        "multiqc_GCF-2026-043_260918.html",
        "sequencer_stats_GCF-2026-043.html",
    ]
    assert attachments[0].get_payload(decode=True) == b"<html>Analysis</html>"


def test_processed_mail_attaches_existing_single_cell_summary(mail_cfg, processed):
    processed["payload"]["libprep"] = "10X Genomics Chromium Single Cell 3prime"
    summary = mail_cfg.output_path / "all_samples_web_summary_GCF-2026-043_260918.html"
    summary.write_text("<html>Single cell</html>")
    message = notifications.build_message(mail_cfg, processed)
    assert summary.name in [part.get_filename() for part in message.iter_attachments()]


@pytest.mark.parametrize(
    "missing", ["multiqc_GCF-2026-043_260918.html", "Stats/sequencer_stats_GCF-2026-043.html"]
)
def test_attachment_failures_do_not_start_smtp(mail_cfg, processed, smtp, missing):
    (mail_cfg.output_path / missing).unlink()
    with pytest.raises(FileNotFoundError):
        notifications.send_notification(mail_cfg, processed)
    smtp[0].assert_not_called()


@pytest.mark.parametrize("source", ["parseSampleSheetMetrics", "getFCmetricsImproved"])
def test_metadata_composition_failure_does_not_start_smtp(
    mail_cfg, processed, smtp, monkeypatch, source
):
    monkeypatch.setattr(misc, source, Mock(side_effect=ValueError("malformed metadata")))
    with pytest.raises(ValueError, match="malformed metadata"):
        notifications.send_notification(mail_cfg, processed)
    smtp[0].assert_not_called()


def test_disk_summary_failure_does_not_start_smtp(mail_cfg, processed, smtp, monkeypatch):
    monkeypatch.setattr(
        afterFastq,
        "_disk_usage_message",
        Mock(side_effect=FileNotFoundError("instrument unmounted")),
    )
    with pytest.raises(FileNotFoundError, match="instrument unmounted"):
        notifications.send_notification(mail_cfg, processed)
    smtp[0].assert_not_called()


@pytest.mark.parametrize("contents", ["", "Incomplete InterOp output\nNo blank header separator\n"])
def test_malformed_interop_header_returns_unavailable_metrics(mail_cfg, contents):
    stats = mail_cfg.output_path / "Stats"
    stats.mkdir()
    (stats / "interop_summary.csv").write_text(contents)
    assert misc.getFCmetricsImproved(mail_cfg) == (
        "Not able to generate table for flowcell metrics."
    )


def test_missing_source_mount_does_not_prevent_disk_summary(mail_cfg):
    assert not mail_cfg.run.flowcell_path.parent.exists()
    summary = afterFastq._disk_usage_message(mail_cfg)
    assert "Current free space for output:" in summary
    assert "Current free space for instruments: unavailable" in summary


@pytest.mark.parametrize("alias", ["errorTo", "errorto"])
def test_missing_error_to_reports_misspelling_and_recovers(mail_cfg, smtp, alias):
    entry = entry_for(mail_cfg)
    mail_cfg.static.email[alias] = mail_cfg.static.email.pop("error_to")
    with pytest.raises(
        notifications.NotificationConfigError, match=f"rename '{alias}' to 'error_to'"
    ):
        notifications.send_notification(mail_cfg, entry)
    smtp[0].assert_not_called()
    mail_cfg.static.email["error_to"] = mail_cfg.static.email.pop(alias)
    notifications.send_notification(mail_cfg, entry)
    smtp[1].send_message.assert_called_once()


def test_unknown_settings_warn_without_exposing_values(mail_cfg, smtp, caplog):
    mail_cfg.static.email["errorTo"] = "private@example.org"
    mail_cfg.static.email["future_setting"] = "private-value"
    notifications.send_notification(mail_cfg, entry_for(mail_cfg))
    assert "Unknown [Email] option 'errorTo' is ignored; use 'error_to'" in caplog.text
    assert "Unknown [Email] option 'future_setting' is ignored" in caplog.text
    assert "private@example.org" not in caplog.text
    assert "private-value" not in caplog.text
    smtp[1].send_message.assert_called_once()


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("host", ""),
        ("host", "relay example.org"),
        ("host", "relay.example.org\nBcc: bad@example.org"),
        ("from_address", ""),
        ("from_address", "no-domain"),
        ("from_address", "first@example.org, second@example.org"),
        ("error_to", ""),
        ("error_to", "no-domain"),
        ("error_to", "person@"),
        ("error_to", "@example.org"),
        ("error_to", "valid@example.org, invalid"),
        ("error_to", "valid@example.org\nBcc: injected@example.org"),
    ],
)
def test_invalid_configuration_never_contacts_smtp(mail_cfg, smtp, key, value):
    mail_cfg.static.email[key] = value
    with pytest.raises(notifications.NotificationConfigError, match=key):
        notifications.send_notification(mail_cfg, entry_for(mail_cfg))
    smtp[0].assert_not_called()


def test_successful_send_has_bounded_timeout_and_explicit_envelope(mail_cfg, processed, smtp):
    notifications.send_notification(mail_cfg, processed)
    smtp[0].assert_called_once_with("relay.example.org", timeout=30)
    smtp[1].send_message.assert_called_once()
    arguments = smtp[1].send_message.call_args
    assert arguments.kwargs == {
        "from_addr": "bfq@example.org",
        "to_addrs": ["analyst@example.org", "other@example.org"],
    }
    smtp[1].quit.assert_called_once()


def test_duplicate_recipient_is_sent_once(mail_cfg, smtp):
    mail_cfg.static.email["error_to"] = "Person <person@example.org>, person@example.org"
    notifications.send_notification(mail_cfg, entry_for(mail_cfg))
    assert smtp[1].send_message.call_args.kwargs["to_addrs"] == ["person@example.org"]


@pytest.mark.parametrize(
    "error",
    [
        OSError("connection refused"),
        TimeoutError("connect timeout"),
        smtplib.SMTPConnectError(421, b"unavailable"),
    ],
)
def test_connection_failure_is_known_failure(mail_cfg, smtp, error):
    smtp[0].side_effect = error
    with pytest.raises(type(error)) as caught:
        notifications.send_notification(mail_cfg, entry_for(mail_cfg))
    assert caught.value is error
    smtp[1].send_message.assert_not_called()


@pytest.mark.parametrize(
    "error",
    [
        smtplib.SMTPRecipientsRefused({"delivery@example.org": (550, b"no mailbox")}),
        smtplib.SMTPSenderRefused(550, b"sender blocked", "bfq@example.org"),
        smtplib.SMTPDataError(554, b"data rejected"),
        smtplib.SMTPHeloError(550, b"hello rejected"),
        smtplib.SMTPNotSupportedError("SMTPUTF8 unavailable"),
    ],
)
def test_explicit_relay_rejection_is_known_failure(mail_cfg, smtp, error):
    smtp[1].send_message.side_effect = error
    with pytest.raises(type(error)) as caught:
        notifications.send_notification(mail_cfg, entry_for(mail_cfg))
    assert caught.value is error
    smtp[1].send_message.assert_called_once()
    smtp[1].quit.assert_called_once()


@pytest.mark.parametrize(
    "error",
    [
        TimeoutError("read timeout"),
        smtplib.SMTPServerDisconnected("lost connection"),
        ConnectionResetError("reset"),
    ],
)
def test_interrupted_send_is_uncertain_and_never_retried_inline(mail_cfg, smtp, error):
    smtp[1].send_message.side_effect = error
    with pytest.raises(DeliveryUncertainError, match="may have accepted"):
        notifications.send_notification(mail_cfg, entry_for(mail_cfg))
    smtp[0].assert_called_once()
    smtp[1].send_message.assert_called_once()
    smtp[1].quit.assert_called_once()


def test_partial_recipient_rejection_is_uncertain(mail_cfg, processed, smtp):
    smtp[1].send_message.return_value = {"other@example.org": (550, b"no mailbox")}
    with pytest.raises(DeliveryUncertainError, match="some recipients.*other@example.org"):
        notifications.send_notification(mail_cfg, processed)
    smtp[1].send_message.assert_called_once()


def test_complete_recipient_rejection_mapping_is_known_failure(mail_cfg, smtp):
    smtp[1].send_message.return_value = {"delivery@example.org": (550, b"no mailbox")}
    with pytest.raises(smtplib.SMTPRecipientsRefused):
        notifications.send_notification(mail_cfg, entry_for(mail_cfg))
    smtp[1].send_message.assert_called_once()


def test_quit_failure_does_not_change_successful_delivery(mail_cfg, smtp, caplog):
    smtp[1].quit.side_effect = smtplib.SMTPServerDisconnected("QUIT connection lost")
    notifications.send_notification(mail_cfg, entry_for(mail_cfg))
    smtp[1].send_message.assert_called_once()
    smtp[1].close.assert_called_once()
    assert "Could not close SMTP session cleanly" in caplog.text


def test_quit_failure_does_not_hide_original_rejection(mail_cfg, smtp):
    rejection = smtplib.SMTPDataError(554, b"data rejected")
    smtp[1].send_message.side_effect = rejection
    smtp[1].quit.side_effect = TimeoutError("QUIT timeout")
    with pytest.raises(smtplib.SMTPDataError) as caught:
        notifications.send_notification(mail_cfg, entry_for(mail_cfg))
    assert caught.value is rejection


def test_unsupported_notification_kind_fails_before_delivery(mail_cfg, smtp):
    entry = entry_for(mail_cfg)
    entry["kind"] = "early-qc"
    with pytest.raises(ValueError, match="Unsupported completion notification"):
        notifications.send_notification(mail_cfg, entry)
    smtp[0].assert_not_called()


@pytest.fixture
def sequencing(mail_cfg):
    report = mail_cfg.output_path / "sequencer_stats_GCF-2026-043_260918.html"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("<html>Early sequencing QC</html>")
    summary = {
        "projects": ["GCF-2026-043"],
        "report_path": str(report),
        "summary_text": "Planned samples: 3\nR1 Q30: 93.2%; Undetermined: 1.8%\nRead geometry: 2x150",
        "summary_html": "Planned samples: 3<br>R1 Q30: 93.2%; Undetermined: 1.8%<br>Read geometry: 2x150",
    }
    return {
        "id": "sequencing:1",
        "kind": "sequencing",
        "created_at": "2026-09-30T07:00:00+00:00",
        "payload": notifications.make_payload(
            mail_cfg, "sequencing", run_time=timedelta(minutes=12), sequencing_qc=summary
        ),
    }


def test_early_qc_requires_no_excel_analysis_or_stats_json(mail_cfg, sequencing, monkeypatch):
    forbidden = Mock(side_effect=AssertionError("Early email attempted analysis-dependent work"))
    monkeypatch.setattr(misc, "parseSampleSheetMetrics", forbidden)
    monkeypatch.setattr(misc, "analysisSampleMetrics", forbidden)
    monkeypatch.setattr(misc, "getFCmetricsImproved", forbidden)
    monkeypatch.setattr(afterFastq, "get_read_geometry", forbidden)
    monkeypatch.setattr(afterFastq, "get_project_dirs", forbidden)
    assert not mail_cfg.run.sample_sheet.exists()
    assert not mail_cfg.run.sample_submission_form.exists()
    assert not (mail_cfg.output_path / "Stats" / "Stats.json").exists()
    message = notifications.build_message(mail_cfg, sequencing)
    assert "Demultiplexing complete — sequencing QC" in message["Subject"]
    assert "analyst@example.org" in message["To"]
    for kind in ("plain", "html"):
        body = message.get_body(preferencelist=(kind,)).get_content()
        for value in (
            "Saved User",
            mail_cfg.run.run_id,
            "0:12:00",
            "Planned samples: 3",
            "93.2%",
            "1.8%",
            "2x150",
            "Current free space for output",
        ):
            assert value in body
    assert [part.get_filename() for part in message.iter_attachments()] == [
        "sequencer_stats_GCF-2026-043_260918.html"
    ]
    forbidden.assert_not_called()


def test_early_payload_does_not_discover_fastqs(mail_cfg, sequencing, monkeypatch):
    monkeypatch.setattr(
        afterFastq, "get_project_dirs", Mock(side_effect=AssertionError("FASTQ scan"))
    )
    monkeypatch.setattr(Path, "open", Mock(side_effect=AssertionError("Input read")))
    payload = notifications.make_payload(
        mail_cfg, "sequencing", sequencing_qc=sequencing["payload"]["sequencing_qc"]
    )
    assert json.loads(json.dumps(payload))["projects"] == ["GCF-2026-043"]


def test_early_missing_report_fails_before_smtp(mail_cfg, sequencing, smtp):
    Path(sequencing["payload"]["sequencing_qc"]["report_path"]).unlink()
    with pytest.raises(FileNotFoundError):
        notifications.send_notification(mail_cfg, sequencing)
    smtp[0].assert_not_called()


def test_early_delivery_uses_mock_smtp(mail_cfg, sequencing, smtp):
    notifications.send_notification(mail_cfg, sequencing)
    assert smtp[1].send_message.call_args.kwargs["to_addrs"] == [
        "analyst@example.org",
        "other@example.org",
    ]


def test_new_analysis_mail_keeps_metadata_and_only_analysis_attachments(
    mail_cfg, processed, monkeypatch
):
    processed["payload"]["email_version"] = 2
    processed["payload"]["analysis_qc"] = {
        "summary_text": "Read retention: 95.00%; 2 samples discovered; missing",
        "summary_html": "<strong>Read retention: 95.00%</strong><br>2 samples discovered<br>missing",
    }
    monkeypatch.setattr(
        misc,
        "analysisSampleMetrics",
        lambda cfg,
        projects: "2 samples discovered in FASTQs.\nPlanned samples without FASTQs: missing.",
    )
    forbidden = Mock(side_effect=AssertionError("Analysis email attempted old sequencing metrics"))
    monkeypatch.setattr(misc, "getFCmetricsImproved", forbidden)
    monkeypatch.setattr(afterFastq, "get_read_geometry", forbidden)
    (mail_cfg.output_path / "Stats" / "sequencer_stats_GCF-2026-043.html").unlink()
    message = notifications.build_message(mail_cfg, processed)
    assert "Analysis complete — QC summary" in message["Subject"]
    for kind in ("plain", "html"):
        body = message.get_body(preferencelist=(kind,)).get_content()
        for value in (
            "Sample groups",
            "95.00%",
            "2 samples discovered",
            "missing",
            "Saved User",
            "2:00:00",
        ):
            assert value in body
    assert [part.get_filename() for part in message.iter_attachments()] == [
        "multiqc_GCF-2026-043_260918.html"
    ]
    forbidden.assert_not_called()


def test_new_analysis_mail_escapes_user_and_retains_single_cell_report(mail_cfg, processed):
    processed["payload"].update(
        email_version=2, user="<script>User</script>", libprep="Parse Biosciences"
    )
    summary = mail_cfg.output_path / "all_samples_web_summary_GCF-2026-043_260918.html"
    summary.write_text("<html>Single cell</html>")
    message = notifications.build_message(mail_cfg, processed)
    html = message.get_body(preferencelist=("html",)).get_content()
    assert "&lt;script&gt;User&lt;/script&gt;" in html
    assert "<script>User</script>" not in html
    assert summary.name in [part.get_filename() for part in message.iter_attachments()]


def test_legacy_rerun_analysis_mail_includes_sequencing_qc(mail_cfg, processed, sequencing):
    processed["payload"].update(
        email_version=2,
        sequencing_qc=sequencing["payload"]["sequencing_qc"],
    )
    message = notifications.build_message(mail_cfg, processed)
    for kind in ("plain", "html"):
        body = message.get_body(preferencelist=(kind,)).get_content()
        assert "Sequencing QC is included here for this legacy run" in body
        assert "no separate early notification was created" in body
        assert "covered in the early sequencing report" not in body
        for value in ("93.2%", "1.8%", "2x150", "Sample groups"):
            assert value in body
    assert [part.get_filename() for part in message.iter_attachments()] == [
        "multiqc_GCF-2026-043_260918.html",
        "sequencer_stats_GCF-2026-043_260918.html",
    ]


def test_legacy_rerun_missing_sequencing_attachment_fails_before_smtp(
    mail_cfg, processed, sequencing, smtp
):
    processed["payload"].update(
        email_version=2,
        sequencing_qc=sequencing["payload"]["sequencing_qc"],
    )
    Path(processed["payload"]["sequencing_qc"]["report_path"]).unlink()
    with pytest.raises(FileNotFoundError):
        notifications.send_notification(mail_cfg, processed)
    smtp[0].assert_not_called()


@pytest.mark.parametrize("version", [1, 2])
def test_analysis_budget_applies_to_current_and_legacy_delivery(mail_cfg, processed, smtp, version):
    processed["payload"].update(email_version=version, libprep="Parse Biosciences")
    multiqc = mail_cfg.output_path / "multiqc_GCF-2026-043_260918.html"
    supplementary = mail_cfg.output_path / "all_samples_web_summary_GCF-2026-043_260918.html"
    multiqc.write_bytes(b"x" * 30_000)
    supplementary.write_bytes(b"y" * 20_000)
    mail_cfg.static.email["max_message_bytes"] = "15000"
    smtp[1].esmtp_features = {}
    before = {
        path: (path.read_bytes(), path.stat().st_mtime_ns) for path in (multiqc, supplementary)
    }
    notifications.send_notification(mail_cfg, processed)
    sent = smtp[1].send_message.call_args.args[0]
    assert email_attachments.message_size(sent) <= 15_000
    assert not any(part.get_filename() == multiqc.name for part in sent.iter_attachments())
    for kind in ("plain", "html"):
        body = sent.get_body(preferencelist=(kind,)).get_content()
        assert str(multiqc) in body
        assert str(supplementary) in body
        assert "even with no other report attachments" in body
        assert "Review the attached MultiQC" not in body
    for path, expected in before.items():
        assert (path.read_bytes(), path.stat().st_mtime_ns) == expected
    smtp[1].send_message.assert_called_once()


@pytest.mark.parametrize("advertised", ["12000", "0", "", "invalid", "-1", "999999999"])
def test_relay_size_uses_lower_positive_limit(mail_cfg, processed, smtp, advertised):
    processed["payload"]["email_version"] = 2
    report = mail_cfg.output_path / "multiqc_GCF-2026-043_260918.html"
    report.write_bytes(b"x" * 15_000)
    smtp[1].esmtp_features = {"size": advertised}
    mail_cfg.static.email["max_message_bytes"] = "25000"
    notifications.send_notification(mail_cfg, processed)
    sent = smtp[1].send_message.call_args.args[0]
    expected_limit = 12_000 if advertised == "12000" else 25_000
    assert email_attachments.message_size(sent) <= expected_limit
    assert bool(list(sent.iter_attachments())) == (advertised != "12000")
    smtp[1].ehlo_or_helo_if_needed.assert_called_once()
    smtp[1].quit.assert_called_once()


@pytest.mark.parametrize("value", ["0", "-1", "1.5", "20 MB", "", "NaN", "١٢٣"])
def test_invalid_message_budget_fails_before_smtp(mail_cfg, processed, smtp, value):
    mail_cfg.static.email["max_message_bytes"] = value
    with pytest.raises(notifications.NotificationConfigError, match="positive integer in bytes"):
        notifications.send_notification(mail_cfg, processed)
    smtp[0].assert_not_called()


def test_relay_budget_smaller_than_summary_fails_before_data(mail_cfg, processed, smtp):
    smtp[1].esmtp_features = {"size": "100"}
    with pytest.raises(email_attachments.MessageSizeError, match="even without attachments"):
        notifications.send_notification(mail_cfg, processed)
    smtp[1].send_message.assert_not_called()
    smtp[1].quit.assert_called_once()


def test_relay_negotiation_failure_is_known_and_closes_session(mail_cfg, processed, smtp):
    smtp[1].ehlo_or_helo_if_needed.side_effect = smtplib.SMTPServerDisconnected("EHLO lost")
    with pytest.raises(smtplib.SMTPServerDisconnected):
        notifications.send_notification(mail_cfg, processed)
    smtp[1].send_message.assert_not_called()
    smtp[1].quit.assert_called_once()


def test_international_envelope_budget_matches_smtp_utf8_serialization(mail_cfg, processed, smtp):
    mail_cfg.static.email["finished_to"] = "analyst@example.ø"
    smtp[1].esmtp_features = {}
    notifications.send_notification(mail_cfg, processed)
    sent = smtp[1].send_message.call_args.args[0]
    assert sent.policy.utf8
    assert email_attachments.message_size(sent) == len(
        sent.as_bytes(policy=sent.policy.clone(utf8=True, linesep="\r\n"))
    )


@pytest.mark.parametrize("international", [False, True])
@pytest.mark.parametrize("below_boundary", [False, True])
def test_budget_matches_real_send_message_wire_bytes(
    mail_cfg, processed, monkeypatch, international, below_boundary
):
    processed["payload"]["email_version"] = 2
    if international:
        mail_cfg.static.email["finished_to"] = "analyst@example.ø"
    report = mail_cfg.output_path / "multiqc_GCF-2026-043_260918.html"
    report.write_bytes(b"x" * 10_000)
    full = notifications.build_message(mail_cfg, processed)
    limit = email_attachments.message_size(full) - int(below_boundary)
    # Run the stdlib serializer with only SMTP transport/negotiation replaced.
    smtp = object.__new__(SMTPClient)
    smtp.ehlo_or_helo_if_needed = Mock()
    smtp.has_extn = Mock(return_value=True)
    smtp.esmtp_features = {"size": str(limit)}
    smtp.sendmail = Mock(return_value={})
    smtp.quit = Mock()
    monkeypatch.setattr(notifications.smtplib, "SMTP", Mock(return_value=smtp))
    notifications.send_notification(mail_cfg, processed)
    wire = smtp.sendmail.call_args.args[2]
    assert isinstance(wire, bytes)
    assert len(wire) <= limit
    assert len(wire) == limit if not below_boundary else len(wire) < limit
    assert (
        b"Content-Disposition: attachment" in wire
        if not below_boundary
        else b"Content-Disposition: attachment" not in wire
    )
    assert (b"not attached" in wire) == below_boundary
