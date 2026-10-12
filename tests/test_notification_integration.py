"""Exercise notification recovery through the daemon and operator commands."""

import copy
import logging
import smtplib

from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from threading import Event
from unittest.mock import Mock

import flowcell_manager.flowcell_manager as manager
import pytest

from bcl2fastq_pipeline.config import PipelineConfig, RunContext
from bcl2fastq_pipeline.state import FlowcellStateStore, StateConflictError, new_state
from test_state_integration import RUN_ID, configured_bfq, write_fastq, write_inputs

from bcl2fastq_pipeline import (
    afterFastq,
    cli,
    email_attachments,
    findFlowCells,
    misc,
    notification_delivery,
    notifications,
)

PROJECT = "GCF-2026-001"


def prepare_pipeline(tmp_path, monkeypatch, *, start_stage="analysis"):
    cfg, source, output = configured_bfq(tmp_path)
    cfg.static.email.update(
        host="smtp.invalid",
        from_address="bfq@example.test",
        finished_to="finished@example.test",
        error_to="operator@example.test",
    )
    write_inputs(output)
    write_fastq(output / PROJECT / "sample_R1.fastq.gz")
    afterFastq.md5sum_worker(cfg)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep", "User": "test"},
        output / "SampleSheet.csv",
        output / "Sample-Submission-Form.xlsx",
    )
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(
        new_state(
            RUN_ID,
            source,
            output,
            origin="restored_legacy_fastq",
            start_stage=start_stage,
            cfg=cfg,
        )
    )
    calls = []

    def analyze():
        calls.append("analysis")
        (output / "analysis-config.yaml").write_text("workflow: tested\n")

    def report(*_args):
        calls.append("reporting")
        (output / f"multiqc_{PROJECT}_260918.html").write_text("<html>QC</html>")
        return [PROJECT]

    def finalize(**_kwargs):
        calls.append("finalization")
        (output / f"{PROJECT}_260918.7za").write_bytes(b"prepared delivery archive")
        (output / f"md5sum_{PROJECT}_archive.txt").write_text("archive checksum\n")

    monkeypatch.setattr(misc, "enoughFreeSpace", lambda: True)
    monkeypatch.setattr(afterFastq, "analysis_steps", analyze)
    monkeypatch.setattr(cli, "_run_reporting", report)
    monkeypatch.setattr(afterFastq, "finalize", finalize)
    monkeypatch.setattr(findFlowCells, "markFinished", lambda: [PROJECT])
    monkeypatch.setattr(
        misc, "write_error_report", lambda *_args: cfg.static.paths.report_dir / "run.error"
    )
    monkeypatch.setattr(misc, "send_error_report", Mock())
    return cfg, store, output, calls


def run_pipeline(cfg, store):
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("notification-integration"))


def products(output):
    """Assert both content and modification time survive a notification-only operation."""
    return {
        str(path.relative_to(output)): (path.read_bytes(), path.stat().st_mtime_ns)
        for path in output.rglob("*")
        if path.is_file()
    }


def processing_state(state):
    return {
        key: value
        for key, value in state.items()
        if key not in {"updated_at", "delivery_notifications"}
    }


def test_missing_error_to_does_not_fail_processing_and_corrected_config_retries(
    tmp_path, monkeypatch
):
    cfg, store, output, calls = prepare_pipeline(tmp_path, monkeypatch)
    cfg.static.email["errorTo"] = cfg.static.email.pop("error_to")
    smtp = Mock()
    smtp.send_message.return_value = {}
    monkeypatch.setattr(notifications.smtplib, "SMTP", Mock(return_value=smtp))
    monkeypatch.setattr(
        notifications,
        "_processed_message",
        lambda _cfg, message, _payload: message.set_content("Finished processing"),
    )

    run_pipeline(cfg, store)

    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert state["last_error"] is None
    assert all(stage["status"] == "completed" for stage in state["stages"].values())
    assert calls == ["analysis", "reporting", "finalization"]
    processed, finalized = state["delivery_notifications"]
    assert processed["status"] == "sent"
    assert finalized["status"] == "failed"
    assert "rename 'errorTo' to 'error_to'" in finalized["last_error"]
    assert smtp.send_message.call_count == 1
    assert not misc.send_error_report.called
    before = products(output)
    before_processing = processing_state(state)

    # A fresh process reads corrected settings; no live run context is required.
    updated_email = {**cfg.static.email, "error_to": "corrected@example.test"}
    del updated_email["errorTo"]
    fresh_cfg = PipelineConfig(static=replace(cfg.static, email=updated_email), run=RunContext())
    PipelineConfig._instance = fresh_cfg
    state = manager.retry_notifications(flowcell=RUN_ID)

    assert [entry["status"] for entry in state["delivery_notifications"]] == ["sent", "sent"]
    assert smtp.send_message.call_count == 2
    sent = smtp.send_message.call_args
    assert sent.kwargs["to_addrs"] == ["corrected@example.test"]
    assert RUN_ID in sent.args[0].get_content()
    assert processing_state(state) == before_processing
    assert products(output) == before
    assert calls == ["analysis", "reporting", "finalization"]
    assert not fresh_cfg.run.run_id
    manager.retry_notifications(flowcell=RUN_ID)
    notification_delivery.recover_pending(fresh_cfg, store)
    assert smtp.send_message.call_count == 2


@pytest.mark.parametrize("error", [ValueError("bad metadata"), FileNotFoundError("report missing")])
def test_composition_and_attachment_failures_retain_completion_and_saved_context(
    tmp_path, monkeypatch, error
):
    cfg, store, output, calls = prepare_pipeline(tmp_path, monkeypatch)
    monkeypatch.setattr(notifications, "build_message", Mock(side_effect=error))

    run_pipeline(cfg, store)

    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert {entry["status"] for entry in state["delivery_notifications"]} == {"failed"}
    assert all(str(error) in entry["last_error"] for entry in state["delivery_notifications"])
    assert calls == ["analysis", "reporting", "finalization"]
    before = products(output)
    before_processing = processing_state(state)
    expected_payloads = [
        copy.deepcopy(entry["payload"]) for entry in state["delivery_notifications"]
    ]
    assert all(payload["run_id"] == RUN_ID for payload in expected_payloads)
    assert expected_payloads[0]["sample_sheet"] == str(output / "SampleSheet.csv")
    assert expected_payloads[0]["projects"] == [PROJECT]

    # Retrying another flowcell must not hijack the daemon's active run context.
    cfg.run.begin(cfg.static.paths.nova_base_dir / "another-flowcell", cfg.static.paths)
    active_context = copy.deepcopy(cfg.run)
    delivered = []
    monkeypatch.setattr(
        notifications, "send_notification", lambda _cfg, entry: delivered.append(entry["payload"])
    )
    notification_delivery.recover_pending(cfg, store)
    assert delivered == []  # Failures are explicit retries, not repeated on every scan.
    retried = manager.retry_notifications(flowcell=RUN_ID)
    assert delivered == expected_payloads
    assert cfg.run == active_context
    assert products(output) == before
    assert processing_state(retried) == before_processing


@pytest.mark.parametrize("stage", ["reporting", "finalization"])
def test_real_processing_failure_is_not_converted_into_notification_failure(
    tmp_path, monkeypatch, stage
):
    cfg, store, output, _calls = prepare_pipeline(tmp_path, monkeypatch)
    sender = Mock()
    monkeypatch.setattr(notifications, "send_notification", sender)
    failure = Mock(side_effect=OSError(f"{stage} disk failure"))
    if stage == "reporting":
        monkeypatch.setattr(cli, "_run_reporting", failure)
    else:
        monkeypatch.setattr(afterFastq, "finalize", failure)

    run_pipeline(cfg, store)

    state = store.read(RUN_ID)
    assert state["status"] == "failed"
    assert state["current_stage"] == stage
    assert state["stages"][stage]["status"] == "failed"
    assert f"{stage} disk failure" in state["last_error"]["summary"]
    assert misc.send_error_report.call_count == 1
    assert not (output / f"{PROJECT}_260918.7za").exists()
    kinds = [entry["kind"] for entry in state["delivery_notifications"]]
    assert kinds == ([] if stage == "reporting" else ["processed"])
    assert sender.call_count == (0 if stage == "reporting" else 1)


def test_accepted_email_with_failed_success_write_needs_explicit_uncertain_retry(
    tmp_path, monkeypatch
):
    cfg, store, output, calls = prepare_pipeline(tmp_path, monkeypatch)
    real_write = store._atomic_write_unlocked
    injected = False

    def fail_after_acceptance(state):
        nonlocal injected
        if not injected and any(
            entry["kind"] == "finalized" and entry["status"] == "sent"
            for entry in state["delivery_notifications"]
        ):
            injected = True
            raise OSError("state disk unavailable after SMTP accepted")
        return real_write(state)

    sender = Mock()
    monkeypatch.setattr(notifications, "send_notification", sender)
    monkeypatch.setattr(store, "_atomic_write_unlocked", fail_after_acceptance)
    run_pipeline(cfg, store)

    assert injected
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert state["last_error"] is None
    assert state["delivery_notifications"][-1]["status"] == "sending"
    assert sender.call_count == 2
    before = products(output)
    before_processing = processing_state(state)
    restarted_store = FlowcellStateStore(cfg.static.paths.manager_dir)
    notification_delivery.recover_pending(cfg, restarted_store)
    assert restarted_store.read(RUN_ID)["delivery_notifications"][-1]["status"] == "uncertain"
    assert sender.call_count == 2
    with pytest.raises(StateConflictError, match="--retry-uncertain"):
        manager.retry_notifications(flowcell=RUN_ID, kind="finalized")
    assert sender.call_count == 2

    retried = manager.retry_notifications(flowcell=RUN_ID, kind="finalized", retry_uncertain=True)
    assert sender.call_count == 3
    assert retried["delivery_notifications"][-1]["status"] == "sent"
    assert processing_state(retried) == before_processing
    assert products(output) == before
    assert calls == ["analysis", "reporting", "finalization"]


def test_restart_recovers_completed_run_notification_without_reentering_pipeline(
    tmp_path, monkeypatch
):
    cfg, store, output, calls = prepare_pipeline(tmp_path, monkeypatch, start_stage="finalization")
    deliver = notification_delivery.deliver_pending
    monkeypatch.setattr(
        notification_delivery, "deliver_pending", Mock(side_effect=KeyboardInterrupt("restart"))
    )
    with pytest.raises(KeyboardInterrupt, match="restart"):
        run_pipeline(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert state["delivery_notifications"][-1]["status"] == "pending"
    before = products(output)
    before_processing = processing_state(state)
    monkeypatch.setattr(notification_delivery, "deliver_pending", deliver)
    sender = Mock()
    monkeypatch.setattr(notifications, "send_notification", sender)
    cfg.run.reset()

    notification_delivery.recover_pending(cfg, FlowcellStateStore(cfg.static.paths.manager_dir))
    notification_delivery.recover_pending(cfg, FlowcellStateStore(cfg.static.paths.manager_dir))

    assert sender.call_count == 1
    assert cli.candidate_flowcells(cfg, store) == []
    assert store.read(RUN_ID)["delivery_notifications"][-1]["status"] == "sent"
    assert processing_state(store.read(RUN_ID)) == before_processing
    assert calls == ["finalization"]
    assert products(output) == before


def completed_with_failed_notifications(tmp_path, monkeypatch):
    cfg, store, output, calls = prepare_pipeline(tmp_path, monkeypatch)
    monkeypatch.setattr(
        notifications, "send_notification", Mock(side_effect=OSError("SMTP offline"))
    )
    run_pipeline(cfg, store)
    return cfg, store, output, calls


def test_email_callbacks_observe_durably_completed_work(tmp_path, monkeypatch):
    cfg, store, _output, _calls = prepare_pipeline(tmp_path, monkeypatch)
    observed = []

    def inspect_state(_cfg, entry):
        # A separately constructed reader sees completion and the claimed intent,
        # even though SMTP delivery has not returned yet.
        observed.append((entry["kind"], FlowcellStateStore(store.manager_dir).read(RUN_ID)))

    monkeypatch.setattr(notifications, "send_notification", inspect_state)
    run_pipeline(cfg, store)

    (processed_kind, processed_state), (finalized_kind, finalized_state) = observed
    assert processed_kind == "processed"
    assert processed_state["stages"]["reporting"]["status"] == "completed"
    assert processed_state["stages"]["finalization"]["status"] == "queued"
    assert processed_state["delivery_notifications"][0]["status"] == "sending"
    assert finalized_kind == "finalized"
    assert finalized_state["status"] == "completed"
    assert finalized_state["stages"]["finalization"]["status"] == "completed"
    assert finalized_state["delivery_notifications"][-1]["status"] == "sending"
    assert {entry["status"] for entry in store.read(RUN_ID)["delivery_notifications"]} == {"sent"}


def test_finalization_rerun_preserves_processed_intent_and_replaces_finalized_intent(
    tmp_path, monkeypatch
):
    cfg, store, output, calls = completed_with_failed_notifications(tmp_path, monkeypatch)
    before = store.read(RUN_ID)["delivery_notifications"]
    state = manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", force=True)
    processed, finalized = state["delivery_notifications"]
    assert processed == before[0]
    assert finalized["status"] == "superseded"
    delivered = []
    monkeypatch.setattr(
        notifications, "send_notification", lambda _cfg, entry: delivered.append(entry["id"])
    )
    manager.retry_notifications(flowcell=RUN_ID, kind="processed")
    assert delivered == ["processed:1"]

    cfg.run.begin(cfg.static.paths.nova_base_dir / RUN_ID, cfg.static.paths)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep", "User": "test"},
        output / "SampleSheet.csv",
        output / "Sample-Submission-Form.xlsx",
    )
    run_pipeline(cfg, store)

    assert delivered == ["processed:1", "finalized:2"]
    assert calls == ["analysis", "reporting", "finalization", "finalization"]
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert [(entry["id"], entry["status"]) for entry in state["delivery_notifications"]] == [
        ("processed:1", "sent"),
        ("finalized:1", "superseded"),
        ("finalized:2", "sent"),
    ]


def test_active_retry_excludes_other_retries_and_forced_cleanup(tmp_path, monkeypatch):
    cfg, store, output, _calls = completed_with_failed_notifications(tmp_path, monkeypatch)
    started, release = Event(), Event()
    accepted = []

    def slow_send(_cfg, entry):
        started.set()
        assert release.wait(10), "test did not release blocked SMTP delivery"
        accepted.append(entry["id"])

    monkeypatch.setattr(notifications, "send_notification", slow_send)
    before = products(output)
    with ThreadPoolExecutor(max_workers=1) as executor:
        retry = executor.submit(manager.retry_notifications, flowcell=RUN_ID)
        try:
            assert started.wait(10), "retry did not start SMTP delivery"
            during = store.read(RUN_ID)
            with pytest.raises(StateConflictError, match="active"):
                manager.retry_notifications(flowcell=RUN_ID)
            for force in (False, True):
                with pytest.raises(StateConflictError, match="active|running"):
                    manager.rerun_flowcell(flowcell=RUN_ID, from_stage="analysis", force=force)
                with pytest.raises(StateConflictError, match="active|running"):
                    manager.archive_flowcell(flowcell=RUN_ID, force=force)
            notification_delivery.recover_pending(cfg, store)
            assert store.read(RUN_ID) == during
            assert products(output) == before
        finally:
            release.set()
        retried = retry.result(timeout=10)
    assert accepted == ["processed:1", "finalized:1"]
    assert retried["status"] == "completed"
    assert products(output) == before


@pytest.mark.parametrize("operation", ["rerun", "archive"])
def test_cleanup_invalidates_notifications_and_excludes_concurrent_retry(
    tmp_path, monkeypatch, operation
):
    _cfg, store, output, _calls = completed_with_failed_notifications(tmp_path, monkeypatch)
    started, release = Event(), Event()
    actual_cleanup = manager.apply_cleanup
    sender = Mock()
    monkeypatch.setattr(notifications, "send_notification", sender)

    def slow_cleanup(paths):
        started.set()
        assert release.wait(10), "test did not release blocked cleanup"
        return actual_cleanup(paths)

    monkeypatch.setattr(manager, "apply_cleanup", slow_cleanup)
    action = manager.rerun_flowcell if operation == "rerun" else manager.archive_flowcell
    kwargs = {"flowcell": RUN_ID, "force": True}
    if operation == "rerun":
        kwargs["from_stage"] = "reporting"
    with ThreadPoolExecutor(max_workers=1) as executor:
        cleanup = executor.submit(action, **kwargs)
        try:
            assert started.wait(10), "cleanup did not start"
            assert {entry["status"] for entry in store.read(RUN_ID)["delivery_notifications"]} == {
                "superseded"
            }
            with pytest.raises(StateConflictError, match="active"):
                manager.retry_notifications(flowcell=RUN_ID)
        finally:
            release.set()
        cleanup.result(timeout=10)
    with pytest.raises(StateConflictError, match="superseded outputs"):
        manager.retry_notifications(flowcell=RUN_ID)
    assert not sender.called
    assert not (output / f"{PROJECT}_260918.7za").exists()


def test_legacy_completed_run_does_not_acquire_new_notifications(tmp_path, monkeypatch):
    cfg, store, _output, _calls = prepare_pipeline(tmp_path, monkeypatch)
    monkeypatch.setattr(notifications, "send_notification", Mock())
    run_pipeline(cfg, store)
    legacy = store.read(RUN_ID)
    del legacy["delivery_notifications"]
    store.write(legacy)
    sender = Mock()
    monkeypatch.setattr(notifications, "send_notification", sender)
    before = store.state_path(RUN_ID).read_bytes()

    notification_delivery.recover_pending(cfg, store)
    with pytest.raises(StateConflictError, match="Legacy notifications are not reconstructed"):
        manager.retry_notifications(flowcell=RUN_ID)

    assert not sender.called
    assert store.state_path(RUN_ID).read_bytes() == before


@pytest.mark.parametrize("version", [1, 2])
def test_oversize_rejection_then_summary_only_retry_preserves_processing(
    tmp_path, monkeypatch, version
):
    cfg, store, output, calls = prepare_pipeline(tmp_path, monkeypatch)
    cfg.static.email["max_message_bytes"] = "200000"
    cfg.run.libprep = "Parse Biosciences"
    multiqc = output / f"multiqc_{PROJECT}_260918.html"
    supplementary = output / f"all_samples_web_summary_{PROJECT}_260918.html"

    def report(*_args):
        calls.append("reporting")
        multiqc.write_bytes(b"x" * 30_000)
        supplementary.write_bytes(b"y" * 20_000)
        (output / "Stats").mkdir(exist_ok=True)
        (output / "Stats" / f"sequencer_stats_{PROJECT}.html").write_text("Sequencing QC")
        return [PROJECT]

    monkeypatch.setattr(cli, "_run_reporting", report)
    monkeypatch.setattr(misc, "parseSampleSheetMetrics", lambda *_args, **_kwargs: "Sample groups")
    monkeypatch.setattr(misc, "analysisSampleMetrics", lambda *_args: "Analysis summary")
    monkeypatch.setattr(misc, "getFCmetricsImproved", lambda *_args: "Flowcell metrics")
    monkeypatch.setattr(afterFastq, "get_read_geometry", lambda *_args: "2x150")
    monkeypatch.setattr(afterFastq, "_disk_usage_message", lambda *_args: "Disk summary")
    make_payload = notifications.make_payload

    def versioned_payload(*args, **kwargs):
        payload = make_payload(*args, **kwargs)
        payload["email_version"] = version
        return payload

    monkeypatch.setattr(notifications, "make_payload", versioned_payload)
    smtp = Mock()
    smtp.esmtp_features = {}  # A downstream restriction need not appear in EHLO.
    smtp.send_message.side_effect = [smtplib.SMTPDataError(552, b"message too large"), {}]
    monkeypatch.setattr(notifications.smtplib, "SMTP", Mock(return_value=smtp))
    run_pipeline(cfg, store)

    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert [entry["status"] for entry in state["delivery_notifications"]] == ["failed", "sent"]
    assert "552" in state["delivery_notifications"][0]["last_error"]
    assert smtp.send_message.call_count == 2  # No implicit resend after rejection.
    before_products = products(output)
    before_processing = processing_state(state)
    fresh_cfg = PipelineConfig(
        static=replace(cfg.static, email={**cfg.static.email, "max_message_bytes": "15000"}),
        run=RunContext(),
    )
    PipelineConfig._instance = fresh_cfg
    smtp.send_message.side_effect = None
    smtp.send_message.return_value = {}
    retried = manager.retry_notifications(flowcell=RUN_ID, kind="processed")

    sent = smtp.send_message.call_args.args[0]
    assert email_attachments.message_size(sent) <= 15_000
    assert not any(
        part.get_filename() in {multiqc.name, supplementary.name}
        for part in sent.iter_attachments()
    )
    for kind in ("plain", "html"):
        body = sent.get_body(preferencelist=(kind,)).get_content()
        assert str(multiqc) in body
        assert str(supplementary) in body
        assert "not attached" in body
    assert retried["delivery_notifications"][0]["status"] == "sent"
    assert len(retried["delivery_notifications"][0]["attempts"]) == 2
    assert products(output) == before_products
    assert processing_state(retried) == before_processing
    assert calls == ["analysis", "reporting", "finalization"]
    assert not fresh_cfg.run.run_id
    manager.retry_notifications(flowcell=RUN_ID, kind="processed")
    notification_delivery.recover_pending(fresh_cfg, store)
    assert smtp.send_message.call_count == 3
