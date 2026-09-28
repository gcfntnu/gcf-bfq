import logging
import subprocess
import sys

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore, new_state

from bcl2fastq_pipeline import cli, misc


def test_error_report_uses_flowcell_name_and_contains_traceback(tmp_path, monkeypatch):
    report_dir = tmp_path / "reports"
    cfg = SimpleNamespace(
        static=SimpleNamespace(paths=SimpleNamespace(report_dir=report_dir)),
        run=SimpleNamespace(run_id="260924_A01990_0221_TEST"),
    )
    monkeypatch.setattr(misc.PipelineConfig, "get", Mock(return_value=cfg))

    try:
        raise RuntimeError("workflow failed")
    except RuntimeError:
        report_path = misc.write_error_report(sys.exc_info(), "Got an error during postMakeSteps")

    assert report_path == report_dir / "260924_A01990_0221_TEST.error"
    report = report_path.read_text()
    assert "Got an error during postMakeSteps" in report
    assert "Traceback (most recent call last):" in report
    assert "RuntimeError: workflow failed" in report


def test_run_context_is_reset_after_error_report(monkeypatch):
    events = []
    cfg = SimpleNamespace(run=SimpleNamespace(reset=lambda: events.append("reset")))
    log = logging.getLogger("test-error-reporting")

    def write_report(error_info, message):
        assert error_info[0] is RuntimeError
        events.append("report")

    monkeypatch.setattr(misc, "write_error_report", write_report)

    try:
        raise RuntimeError("workflow failed")
    except RuntimeError:
        cli.report_run_error(cfg, log, "workflow failed")

    assert events == ["report", "reset"]


def test_error_report_includes_captured_command_output(tmp_path, monkeypatch):
    report_dir = tmp_path / "reports"
    cfg = SimpleNamespace(
        static=SimpleNamespace(paths=SimpleNamespace(report_dir=report_dir)),
        run=SimpleNamespace(run_id="260924_A01990_0221_TEST"),
    )
    monkeypatch.setattr(misc.PipelineConfig, "get", Mock(return_value=cfg))

    try:
        raise subprocess.CalledProcessError(
            1,
            ["snakemake", "multiqc_report"],
            output="/usr/bin/bash: wget: command not found\nError in rule ensembl_genome\n",
        )
    except subprocess.CalledProcessError:
        report_path = misc.write_error_report(sys.exc_info(), "Got an error during postMakeSteps")

    report = report_path.read_text()
    assert "Captured command output (last 400 lines):" in report
    assert "wget: command not found" in report
    assert "Error in rule ensembl_genome" in report


@pytest.fixture
def notification_run(tmp_path, monkeypatch):
    run_id = "260924_A01990_0221_TEST"
    cfg = SimpleNamespace(
        static=SimpleNamespace(
            paths=SimpleNamespace(report_dir=tmp_path / "reports"),
            email={
                "host": "relay.example.org",
                "from_address": "bfq@example.org",
                "error_to": "one@example.org, Two <two@example.org>, one@example.org",
            },
        ),
        run=SimpleNamespace(run_id=run_id, reset=Mock()),
    )
    store = FlowcellStateStore(tmp_path / "manager")
    store.create(
        new_state(
            run_id,
            tmp_path / "source",
            tmp_path / "output",
            origin="new",
            start_stage="demultiplexing",
        )
    )
    monkeypatch.setattr(misc.PipelineConfig, "get", Mock(return_value=cfg))
    smtp = Mock()
    smtp.return_value.__enter__ = Mock(return_value=smtp.return_value)
    smtp.return_value.__exit__ = Mock(return_value=False)
    smtp.return_value.send_message.return_value = {}
    monkeypatch.setattr(misc.smtplib, "SMTP", smtp)
    monkeypatch.setenv("BFQ_ENV", "production")
    return cfg, store, smtp


def fail_run(cfg, store, error=None, stage="analysis"):
    try:
        raise error or RuntimeError("workflow failed")
    except Exception:
        return cli.report_run_error(
            cfg, logging.getLogger("test"), "pipeline failure", store=store, stage=stage
        )


def test_production_sends_saved_report_to_all_recipients(notification_run):
    cfg, store, smtp = notification_run

    def check_saved(message, **kwargs):
        path = cfg.static.paths.report_dir / f"{cfg.run.run_id}.error"
        assert path.exists()
        assert store.read(cfg.run.run_id)["notification"]["attempted_at"]
        assert kwargs == {
            "from_addr": "bfq@example.org",
            "to_addrs": ["one@example.org", "two@example.org"],
        }
        assert cfg.run.run_id in message["Subject"]
        assert "analysis" in message["Subject"]
        body = message.get_body(preferencelist=("plain",)).get_content()
        for expected in [
            cfg.run.run_id,
            "analysis",
            "RuntimeError: workflow failed",
            "Timestamp:",
            "Host:",
            str(path),
        ]:
            assert expected in body
        (attachment,) = message.iter_attachments()
        assert attachment.get_content_type() == "text/plain"
        assert attachment.get_filename() == path.name
        assert attachment.get_content().rstrip() == path.read_text().rstrip()
        return {}

    smtp.return_value.send_message.side_effect = check_saved
    fail_run(cfg, store)
    smtp.assert_called_once_with("relay.example.org", timeout=30)
    assert store.read(cfg.run.run_id)["notification"]["notified"] is True
    cfg.run.reset.assert_called_once()


@pytest.mark.parametrize(
    "mode", [None, "test", "development", "invalid", "Production", "production ", ""]
)
def test_nonproduction_never_opens_smtp(notification_run, monkeypatch, mode):
    cfg, store, smtp = notification_run
    if mode is None:
        monkeypatch.delenv("BFQ_ENV")
    else:
        monkeypatch.setenv("BFQ_ENV", mode)
    report = fail_run(cfg, store)
    assert "RuntimeError: workflow failed" in report.read_text()
    smtp.assert_not_called()
    assert not store.read(cfg.run.run_id)["notification"].get("attempted_at")


def test_identical_failure_suppressed_across_rerun_and_store_reload(notification_run):
    cfg, store, smtp = notification_run
    fail_run(cfg, store)
    before = store.read(cfg.run.run_id)["notification"]
    store.set_preparing(cfg.run.run_id, "demultiplexing", reason="retry", refresh_inputs=False)
    store.queue(cfg.run.run_id, "demultiplexing")
    fresh_store = FlowcellStateStore(store.manager_dir)
    fail_run(cfg, fresh_store)
    assert smtp.call_count == 1
    assert fresh_store.read(cfg.run.run_id)["notification"] == before


@pytest.mark.parametrize(
    "error,stage",
    [
        (ValueError("workflow failed"), "analysis"),
        (RuntimeError("different failure"), "analysis"),
        (RuntimeError("workflow failed"), "reporting"),
    ],
)
def test_changed_failure_is_notified(notification_run, error, stage):
    cfg, store, smtp = notification_run
    fail_run(cfg, store)
    fail_run(cfg, store, error, stage)
    assert smtp.call_count == 2


def test_success_clears_notification_and_same_failure_can_notify_again(notification_run):
    cfg, store, smtp = notification_run
    fail_run(cfg, store)
    store.set_preparing(cfg.run.run_id, "demultiplexing", reason=None, refresh_inputs=False)
    store.queue(cfg.run.run_id, "demultiplexing")
    store.begin_attempt(cfg.run.run_id)
    for stage in ("demultiplexing", "analysis", "reporting", "finalization"):
        if stage != "demultiplexing":
            store.start_stage(cfg.run.run_id, stage)
        if stage != "finalization":
            store.complete_stage(cfg.run.run_id, stage)
    store.complete_run(cfg.run.run_id, [])
    assert store.read(cfg.run.run_id)["last_error"] is None
    fail_run(cfg, store)
    assert smtp.call_count == 2


@pytest.mark.parametrize("failure_point", ["connect", "send", "partial", "quit"])
def test_smtp_failure_retains_report_and_original_error(notification_run, caplog, failure_point):
    cfg, store, smtp = notification_run
    error = OSError("relay unavailable")
    if failure_point == "connect":
        smtp.side_effect = error
    elif failure_point == "send":
        smtp.return_value.send_message.side_effect = error
    elif failure_point == "partial":
        smtp.return_value.send_message.return_value = {"two@example.org": (550, b"refused")}
    else:
        smtp.return_value.__exit__.side_effect = error
    report = fail_run(cfg, store)
    assert "RuntimeError: workflow failed" in report.read_text()
    state = store.read(cfg.run.run_id)
    assert state["status"] == "failed"
    assert state["last_error"]["summary"] == "pipeline failure"
    assert not state["notification"]["notified"]
    assert state["notification"]["delivery_error"]
    assert "Unable to deliver the flowcell error email" in caplog.text
    cfg.run.reset.assert_called_once()
    fail_run(cfg, FlowcellStateStore(store.manager_dir))
    assert smtp.call_count == 1


def test_report_write_failure_still_records_failure_and_resets(notification_run, monkeypatch):
    cfg, store, smtp = notification_run
    monkeypatch.setattr(misc, "write_error_report", Mock(side_effect=OSError("disk full")))
    assert fail_run(cfg, store) is None
    assert store.read(cfg.run.run_id)["status"] == "failed"
    smtp.assert_not_called()
    cfg.run.reset.assert_called_once()


def test_state_write_failure_suppresses_email_and_resets(notification_run, monkeypatch):
    cfg, store, smtp = notification_run
    monkeypatch.setattr(store, "fail_stage", Mock(side_effect=OSError("state disk full")))
    assert fail_run(cfg, store).exists()
    smtp.assert_not_called()
    cfg.run.reset.assert_called_once()


def test_empty_recipients_do_not_claim_delivery(notification_run):
    cfg, store, smtp = notification_run
    cfg.static.email["error_to"] = " , "
    fail_run(cfg, store)
    smtp.assert_not_called()
    assert not store.read(cfg.run.run_id)["notification"].get("attempted_at")


def test_failure_signature_includes_command_diagnostics():
    first = subprocess.CalledProcessError(1, ["snakemake"], output="missing input")
    second = subprocess.CalledProcessError(1, ["snakemake"], output="disk full")
    assert misc.error_failure_signature(
        "analysis", (type(first), first, None), "x"
    ) != misc.error_failure_signature("analysis", (type(second), second, None), "x")


def test_suppressed_test_failure_can_notify_when_later_run_in_production(
    notification_run, monkeypatch
):
    cfg, store, smtp = notification_run
    monkeypatch.setenv("BFQ_ENV", "test")
    first = fail_run(cfg, store).read_text()
    monkeypatch.setenv("BFQ_ENV", "production")
    second = fail_run(cfg, store).read_text()
    assert first == second
    smtp.assert_called_once()


def test_concurrent_notification_attempts_send_once(notification_run):
    cfg, store, smtp = notification_run
    signature = "same failure"
    store.fail_stage(
        cfg.run.run_id, "analysis", summary="bad", report_path=None, failure_signature=signature
    )
    send = Mock()
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(
            pool.map(
                lambda _: FlowcellStateStore(store.manager_dir).deliver_failure_notification(
                    cfg.run.run_id, signature, send
                ),
                range(4),
            )
        )
    assert results.count(True) == 1
    send.assert_called_once()
    smtp.assert_not_called()


def test_crash_after_claim_does_not_resend(notification_run):
    cfg, store, smtp = notification_run
    signature = "same failure"
    store.fail_stage(
        cfg.run.run_id, "analysis", summary="bad", report_path=None, failure_signature=signature
    )
    with pytest.raises(KeyboardInterrupt):
        store.deliver_failure_notification(
            cfg.run.run_id, signature, Mock(side_effect=KeyboardInterrupt)
        )
    fresh = FlowcellStateStore(store.manager_dir)
    send = Mock()
    assert fresh.deliver_failure_notification(cfg.run.run_id, signature, send) is False
    send.assert_not_called()
    assert fresh.read(cfg.run.run_id)["notification"]["attempted_at"]
    assert fresh.read(cfg.run.run_id)["notification"]["notified"] is False
    smtp.assert_not_called()


def test_notification_claim_write_failure_never_opens_smtp(notification_run, monkeypatch):
    cfg, store, smtp = notification_run
    original = store._atomic_write_unlocked

    def fail_claim(state):
        if state["notification"].get("attempted_at"):
            raise OSError("cannot persist claim")
        return original(state)

    monkeypatch.setattr(store, "_atomic_write_unlocked", fail_claim)
    assert fail_run(cfg, store).exists()
    smtp.assert_not_called()
    cfg.run.reset.assert_called_once()


def test_stale_notification_cannot_send_after_failure_changes(notification_run):
    cfg, store, smtp = notification_run
    store.fail_stage(cfg.run.run_id, "analysis", summary="new failure", report_path=None)
    send = Mock()
    assert store.deliver_failure_notification(cfg.run.run_id, "old signature", send) is False
    send.assert_not_called()
    smtp.assert_not_called()


def test_images_set_explicit_error_notification_mode():
    root = Path(__file__).resolve().parents[1]
    assert "ENV BFQ_ENV=production" in (root / "dockerfile-prod").read_text()
    assert "ENV BFQ_ENV=test" in (root / "dockerfile-test").read_text()
