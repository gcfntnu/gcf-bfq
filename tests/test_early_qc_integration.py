"""Early reporting order, deduplication and explicit recovery; all mail is mocked."""

import copy
import logging

from unittest.mock import Mock

import flowcell_manager.flowcell_manager as manager
import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore, apply_cleanup, cleanup_plan, new_state
from test_state_integration import RUN_ID, configured_bfq, write_inputs

from bcl2fastq_pipeline import (
    afterFastq,
    cli,
    findFlowCells,
    makeFastq,
    misc,
    notification_delivery,
    notifications,
    sequencing_delivery,
    sequencing_qc,
)


def prepared(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    write_inputs(output)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep", "User": "test"},
        output / "SampleSheet.csv",
        output / "Sample-Submission-Form.xlsx",
    )
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(
        new_state(RUN_ID, source, output, origin="new", start_stage="demultiplexing", cfg=cfg)
    )
    calls = []
    monkeypatch.setattr(misc, "enoughFreeSpace", lambda: True)
    monkeypatch.setattr(makeFastq, "bcl2fq", lambda: ("bcl-convert", "4.2.4"))
    monkeypatch.setattr(makeFastq, "rename_fastqs", lambda: calls.append("rename"))
    monkeypatch.setattr(afterFastq, "md5sum_worker", lambda *_a, **_kw: calls.append("hash"))
    monkeypatch.setattr(afterFastq, "analysis_steps", lambda: calls.append("analysis"))
    monkeypatch.setattr(cli, "_run_reporting", lambda *_a: ["GCF-2026-001"])
    monkeypatch.setattr(afterFastq, "finalize", lambda: calls.append("finalize"))
    monkeypatch.setattr(findFlowCells, "markFinished", lambda: ["GCF-2026-001"])
    monkeypatch.setattr(misc, "write_error_report", lambda *_a: output / "test.error")
    monkeypatch.setattr(misc, "send_error_report", Mock())

    def report(saved, tool=None):
        calls.append("qc")
        assert tool == "bcl-convert"
        report_path = output / "Stats/sequencing_qc/sequencing_qc.html"
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text("<html>99.9% undetermined; operator decides</html>")
        return {
            "projects": ["GCF-2026-001"],
            "report_path": str(report_path),
            "summary_text": "99.9% undetermined",
            "summary_html": "99.9% undetermined",
        }

    monkeypatch.setattr(sequencing_qc, "generate", report)
    monkeypatch.setattr(
        notifications,
        "send_notification",
        lambda _cfg, entry: calls.append(entry["kind"] + "-mail"),
    )
    return cfg, store, output, calls


def run(cfg, store):
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("test-early"))


def test_qc_and_notification_precede_hashing_and_analysis_even_with_poor_qc(tmp_path, monkeypatch):
    cfg, store, _output, calls = prepared(tmp_path, monkeypatch)
    run(cfg, store)
    assert calls == [
        "rename",
        "qc",
        "sequencing-mail",
        "hash",
        "analysis",
        "processed-mail",
        "finalize",
        "finalized-mail",
    ]
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert state["sequencing_qc"]["duration_seconds"] >= 0
    timing = state["stages"]["reporting"]["metadata"]
    assert timing["sequencing_qc_duration_seconds"] == state["sequencing_qc"]["duration_seconds"]
    assert (
        timing["duration_seconds"]
        == timing["analysis_reporting_duration_seconds"] + timing["sequencing_qc_duration_seconds"]
    )


def test_failed_analysis_and_rerun_preserve_report_timing_and_do_not_resend(tmp_path, monkeypatch):
    cfg, store, output, calls = prepared(tmp_path, monkeypatch)
    monkeypatch.setattr(
        afterFastq, "analysis_steps", Mock(side_effect=RuntimeError("workflow broke"))
    )
    run(cfg, store)
    first = store.read(RUN_ID)
    assert first["status"] == "failed"
    assert first["delivery_notifications"][0]["status"] == "sent"
    qc = copy.deepcopy(first["sequencing_qc"])
    report = output / "Stats/sequencing_qc/sequencing_qc.html"
    before = report.stat().st_mtime_ns
    apply_cleanup(cleanup_plan(output, "analysis"))
    store.queue(RUN_ID, "analysis")
    cfg.run.begin(first["source_path"], cfg.static.paths)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep"},
        output / "SampleSheet.csv",
        output / "Sample-Submission-Form.xlsx",
    )
    monkeypatch.setattr(afterFastq, "analysis_steps", lambda: calls.append("analysis"))
    run(cfg, store)
    notification_delivery.recover_pending(cfg, store)
    assert calls.count("qc") == calls.count("sequencing-mail") == 1
    assert report.stat().st_mtime_ns == before
    assert store.read(RUN_ID)["sequencing_qc"] == qc


def test_report_failure_continues_and_explicit_recovery_uses_saved_sheet(tmp_path, monkeypatch):
    cfg, store, output, calls = prepared(tmp_path, monkeypatch)
    generate = sequencing_qc.generate
    monkeypatch.setattr(
        sequencing_qc, "generate", Mock(side_effect=RuntimeError("MultiQC unavailable"))
    )
    run(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert state["sequencing_qc"]["status"] == "failed"
    assert "MultiQC unavailable" in state["sequencing_qc"]["last_error"]
    assert "sequencing-mail" not in calls
    assert "duration_seconds" not in state["sequencing_qc"]
    (output / "SampleSheet.csv").write_text("edited sheet from another analysis")
    (output / "Sample-Submission-Form.xlsx").write_bytes(b"invalid workbook")

    def recovered(saved, tool=None):
        assert "sample,GCF-2026-001,ACGT" in saved.run.sample_sheet.read_text()
        return generate(saved, tool)

    monkeypatch.setattr(sequencing_qc, "generate", recovered)
    manager.retry_sequencing_qc(flowcell=RUN_ID)
    manager.retry_sequencing_qc(flowcell=RUN_ID)
    after = store.read(RUN_ID)
    assert after["status"] == "completed"
    assert after["attempt"] == state["attempt"]
    assert after["sequencing_qc"]["status"] == "completed"
    assert calls.count("sequencing-mail") == 1
    assert calls.count("hash") == calls.count("analysis") == 1


def test_smtp_failure_does_not_fail_pipeline_or_regenerate_qc(tmp_path, monkeypatch):
    cfg, store, _output, calls = prepared(tmp_path, monkeypatch)
    sender = notifications.send_notification

    def unavailable(_cfg, entry):
        if entry["kind"] == "sequencing":
            raise OSError("relay unavailable")
        sender(_cfg, entry)

    monkeypatch.setattr(notifications, "send_notification", unavailable)
    run(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert state["delivery_notifications"][0]["status"] == "failed"
    notification_delivery.recover_pending(cfg, store)
    monkeypatch.setattr(notifications, "send_notification", sender)
    manager.retry_notifications(flowcell=RUN_ID, kind="sequencing")
    assert calls.count("qc") == calls.count("sequencing-mail") == 1


def test_hash_failure_keeps_early_success_and_retry_eligible(tmp_path, monkeypatch):
    cfg, store, _output, calls = prepared(tmp_path, monkeypatch)
    monkeypatch.setattr(afterFastq, "md5sum_worker", Mock(side_effect=OSError("hash read failure")))
    monkeypatch.setattr(notifications, "send_notification", Mock(side_effect=OSError("relay down")))
    run(cfg, store)
    state = store.read(RUN_ID)
    assert state["stages"]["demultiplexing"]["status"] == "failed"
    assert state["sequencing_qc"]["status"] == "completed"
    assert state["delivery_notifications"][0]["status"] == "failed"
    monkeypatch.setattr(
        notifications,
        "send_notification",
        lambda _cfg, entry: calls.append(entry["kind"] + "-mail"),
    )
    manager.retry_notifications(flowcell=RUN_ID, kind="sequencing")
    assert store.read(RUN_ID)["delivery_notifications"][0]["status"] == "sent"
    assert "analysis" not in calls


def test_new_conversion_invalidates_qc_and_permits_new_notification(tmp_path, monkeypatch):
    cfg, store, output, calls = prepared(tmp_path, monkeypatch)
    run(cfg, store)
    state = store.read(RUN_ID)
    store.queue(RUN_ID, "demultiplexing")
    assert "sequencing_qc" not in store.read(RUN_ID)
    apply_cleanup(cleanup_plan(output, "demultiplexing"))
    cfg.run.begin(state["source_path"], cfg.static.paths)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep"},
        output / "SampleSheet.csv",
        output / "Sample-Submission-Form.xlsx",
    )
    run(cfg, store)
    after = store.read(RUN_ID)
    assert after["status"] == "completed"
    early = [e for e in after["delivery_notifications"] if e["kind"] == "sequencing"]
    assert [e["id"] for e in early] == ["sequencing:1", "sequencing:2"]
    assert [e["status"] for e in early] == ["superseded", "sent"]
    assert calls.count("qc") == calls.count("sequencing-mail") == 2
    assert len(after["sequencing_qc_history"]) == 1


def test_report_recovery_keeps_demultiplexing_identity_after_analysis_attempt(
    tmp_path, monkeypatch
):
    cfg, store, output, _calls = prepared(tmp_path, monkeypatch)
    store.begin_attempt(RUN_ID)
    sequencing_delivery.record_conversion(cfg, store, "bcl-convert", "4")
    store.complete_stage(RUN_ID, "demultiplexing")
    store.queue(RUN_ID, "analysis")
    store.begin_attempt(RUN_ID)
    with store.execution_lease(RUN_ID):
        sequencing_delivery.ensure_report(cfg, store, RUN_ID, retry=True)
    state = store.read(RUN_ID)
    assert state["attempt"] == 2
    assert state["delivery_notifications"][0]["id"] == "sequencing:1"


@pytest.mark.parametrize("stage", ["analysis", "reporting", "finalization"])
def test_all_downstream_cleanup_preserves_early_artifacts(tmp_path, stage):
    artifacts = [
        tmp_path / "Stats/sequencing_qc/sequencing_qc.html",
        tmp_path / "Stats/sequencing_qc_samplesheet.csv",
    ]
    for artifact in artifacts:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text("early result")
    apply_cleanup(cleanup_plan(tmp_path, stage))
    assert all(p.exists() for p in artifacts)


def test_recovery_after_reporting_refreshes_only_early_timing(tmp_path, monkeypatch):
    cfg, store, _output, _calls = prepared(tmp_path, monkeypatch)
    generate = sequencing_qc.generate
    monkeypatch.setattr(sequencing_qc, "generate", Mock(side_effect=RuntimeError("report failed")))
    run(cfg, store)
    before = store.read(RUN_ID)["stages"]["reporting"]["metadata"]
    assert before["sequencing_qc_duration_seconds"] is None
    monkeypatch.setattr(sequencing_qc, "generate", generate)
    manager.retry_sequencing_qc(flowcell=RUN_ID)
    state = store.read(RUN_ID)
    after = state["stages"]["reporting"]["metadata"]
    assert (
        after["analysis_reporting_duration_seconds"]
        == before["analysis_reporting_duration_seconds"]
    )
    assert after["sequencing_qc_duration_seconds"] == state["sequencing_qc"]["duration_seconds"]
    assert (
        after["duration_seconds"]
        == before["duration_seconds"] + state["sequencing_qc"]["duration_seconds"]
    )


def test_failed_report_rebuild_preserves_downstream_mail_and_remains_recoverable(
    tmp_path, monkeypatch
):
    cfg, store, output, calls = prepared(tmp_path, monkeypatch)
    sender = notifications.send_notification
    generate = sequencing_qc.generate

    def failed_early_mail(config, entry):
        if entry["kind"] == "sequencing":
            raise OSError("relay unavailable")
        sender(config, entry)

    monkeypatch.setattr(notifications, "send_notification", failed_early_mail)
    run(cfg, store)
    before = store.read(RUN_ID)
    early_attempts = before["delivery_notifications"][0]["attempts"]
    (output / "Stats/sequencing_qc/sequencing_qc.html").unlink()
    monkeypatch.setattr(sequencing_qc, "generate", Mock(side_effect=RuntimeError("report failed")))
    with pytest.raises(RuntimeError, match="Sequencing QC remains unavailable"):
        manager.retry_sequencing_qc(flowcell=RUN_ID)
    with pytest.raises(RuntimeError, match="Notification delivery remains incomplete"):
        manager.retry_notifications(flowcell=RUN_ID, kind="sequencing")
    failed = store.read(RUN_ID)
    assert failed["delivery_notifications"][0]["status"] == "failed"
    assert failed["delivery_notifications"][0]["attempts"] == early_attempts
    assert failed["delivery_notifications"][1:] == before["delivery_notifications"][1:]
    timing = failed["stages"]["reporting"]["metadata"]
    assert timing["sequencing_qc_duration_seconds"] is None
    assert timing["duration_seconds"] == timing["analysis_reporting_duration_seconds"]

    def regenerated(saved, tool=None):
        result = generate(saved, tool)
        result["summary_text"] = "Recovered sequencing summary"
        return result

    monkeypatch.setattr(sequencing_qc, "generate", regenerated)
    manager.retry_sequencing_qc(flowcell=RUN_ID)
    refreshed = store.read(RUN_ID)
    assert (
        refreshed["delivery_notifications"][0]["payload"]["sequencing_qc"]["summary_text"]
        == "Recovered sequencing summary"
    )
    assert refreshed["delivery_notifications"][0]["status"] == "failed"
    monkeypatch.setattr(notifications, "send_notification", sender)
    manager.retry_notifications(flowcell=RUN_ID, kind="sequencing")
    assert calls.count("sequencing-mail") == 1
    assert store.read(RUN_ID)["delivery_notifications"][1:] == before["delivery_notifications"][1:]


@pytest.mark.parametrize("status", ["sent", "sending", "uncertain"])
def test_report_rebuild_preserves_sent_or_uncertain_smtp_identity(tmp_path, monkeypatch, status):
    cfg, store, output, calls = prepared(tmp_path, monkeypatch)
    generate = sequencing_qc.generate
    run(cfg, store)

    def set_delivery_status(state):
        early = state["delivery_notifications"][0]
        early["status"] = status
        early["attempts"][-1]["status"] = status
        return state

    store.mutate(RUN_ID, set_delivery_status)
    before = store.read(RUN_ID)["delivery_notifications"][0]
    (output / "Stats/sequencing_qc/sequencing_qc.html").unlink()
    monkeypatch.setattr(sequencing_qc, "generate", Mock(side_effect=RuntimeError("report failed")))
    with pytest.raises(RuntimeError, match="Sequencing QC remains unavailable"):
        manager.retry_sequencing_qc(flowcell=RUN_ID)
    with store.execution_lease(RUN_ID):
        notification_delivery.deliver_pending(cfg, store, RUN_ID, kind="sequencing", retry=True)
    assert store.read(RUN_ID)["delivery_notifications"][0]["status"] == status
    monkeypatch.setattr(sequencing_qc, "generate", generate)
    manager.retry_sequencing_qc(flowcell=RUN_ID)
    after = store.read(RUN_ID)["delivery_notifications"][0]
    assert after["status"] == ("uncertain" if status == "sending" else status)
    assert after["payload"] == before["payload"]
    assert len(after["attempts"]) == 1
    assert calls.count("sequencing-mail") == 1


def test_report_duration_and_later_reporting_exclude_smtp_time(tmp_path, monkeypatch):
    cfg, store, _output, _calls = prepared(tmp_path, monkeypatch)
    now = [0.0]
    generate = sequencing_qc.generate
    monkeypatch.setattr(sequencing_delivery.time, "monotonic", lambda: now[0])

    def timed_report(*args, **kwargs):
        result = generate(*args, **kwargs)
        now[0] += 7.0
        return result

    def slow_smtp(*_args):
        now[0] += 100.0

    monkeypatch.setattr(sequencing_qc, "generate", timed_report)
    monkeypatch.setattr(notifications, "send_notification", slow_smtp)
    run(cfg, store)
    state = store.read(RUN_ID)
    assert state["sequencing_qc"]["duration_seconds"] == 7.0
    assert state["stages"]["reporting"]["metadata"]["duration_seconds"] == 7.0
    processed = next(e for e in state["delivery_notifications"] if e["kind"] == "processed")
    assert processed["payload"]["run_time"] == "0:00:07"


def test_report_recovery_after_interruption_and_missing_sheet_uses_stats(tmp_path, monkeypatch):
    cfg, store, output, calls = prepared(tmp_path, monkeypatch)
    store.begin_attempt(RUN_ID)
    (output / "SampleSheet.csv").unlink()
    sequencing_delivery.record_conversion(cfg, store, "bcl-convert", "4.2.4", run_time="0:03:00")
    store.start_sequencing_qc(RUN_ID)
    generate = sequencing_qc.generate

    def stats_only(saved, tool=None):
        assert saved.run.sample_sheet is None
        return generate(saved, tool)

    monkeypatch.setattr(sequencing_qc, "generate", stats_only)
    manager.retry_sequencing_qc(flowcell=RUN_ID)
    state = store.read(RUN_ID)
    assert state["sequencing_qc"]["attempts"][0]["status"] == "interrupted"
    assert state["sequencing_qc"]["attempts"][1]["status"] == "completed"
    assert state["delivery_notifications"][0]["payload"]["run_time"] == "0:03:00"
    assert calls == ["qc", "sequencing-mail"]
