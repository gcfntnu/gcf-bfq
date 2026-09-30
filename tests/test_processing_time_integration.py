"""Exercise real timing boundaries and email composition with controlled work/SMTP."""

import copy

import flowcell_manager.flowcell_manager as manager
import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore
from test_early_qc_integration import prepared, run
from test_state_integration import RUN_ID, write_fastq

from bcl2fastq_pipeline import (
    afterFastq,
    analysis_qc,
    cli,
    makeFastq,
    misc,
    notification_delivery,
    notifications,
    processing_times,
    sequencing_delivery,
    sequencing_qc,
)

REAL_FINALIZE = afterFastq.finalize
REAL_MD5_WORKER = afterFastq.md5sum_worker
PROJECT = "GCF-2026-001"


@pytest.fixture
def pipeline(tmp_path, monkeypatch):
    cfg, store, output, _calls = prepared(tmp_path, monkeypatch)
    cfg.static.email.update(
        host="smtp.invalid",
        from_address="bfq@example.test",
        finished_to="finished@example.test",
        error_to="operator@example.test",
    )
    now = [100.0]
    messages = []
    monkeypatch.setattr(processing_times.time, "monotonic", lambda: now[0])

    def advance(seconds):
        now[0] += seconds

    def conversion():
        advance(10)
        write_fastq(output / PROJECT / "sample_R1.fastq.gz")
        return "bcl-convert", "4.2.4"

    def report(*_args):
        advance(5)
        (output / f"multiqc_{PROJECT}_260918.html").write_text("<html>Project QC</html>")
        return [PROJECT]

    early_report = sequencing_qc.generate
    write_manifest = afterFastq._write_fastq_manifest

    def early(*args, **kwargs):
        advance(6)
        return early_report(*args, **kwargs)

    def hashing(*args):
        advance(3)
        write_manifest(*args)

    def archive(_cfg):
        advance(7)
        (output / f"{PROJECT}_260918.7za").write_text("test archive")

    def hash_archives(_cfg):
        assert store.read(RUN_ID)["processing_timings"]["archiving"][-1]["outcome"] == "completed"
        advance(8)
        (output / f"md5sum_{PROJECT}_260918_archive.txt").write_text("test archive checksum")

    def send(saved_cfg, entry):
        # Real MIME construction runs only after successful durations are durable.
        state = store.read(RUN_ID)
        if entry["kind"] != "sequencing":
            assert entry["payload"]["processing_timing"] == processing_times.snapshot(state)
        messages.append((copy.deepcopy(entry), notifications.build_message(saved_cfg, entry)))
        advance(1000)  # SMTP must never inflate any category.

    monkeypatch.setattr(makeFastq, "bcl2fq", conversion)
    monkeypatch.setattr(makeFastq, "rename_fastqs", lambda: advance(2))
    monkeypatch.setattr(afterFastq, "md5sum_worker", REAL_MD5_WORKER)
    monkeypatch.setattr(afterFastq, "_write_fastq_manifest", hashing)
    monkeypatch.setattr(afterFastq, "analysis_steps", lambda: advance(4))
    monkeypatch.setattr(cli, "_run_reporting", report)
    monkeypatch.setattr(sequencing_qc, "generate", early)
    monkeypatch.setattr(afterFastq, "finalize", REAL_FINALIZE)
    monkeypatch.setattr(afterFastq, "archive_worker", archive)
    monkeypatch.setattr(afterFastq, "md5sum_archive_worker", hash_archives)
    monkeypatch.setattr(analysis_qc, "collect", lambda *_a: {"summary_html": "Sample QC"})
    monkeypatch.setattr(misc, "parseSampleSheetMetrics", lambda *_a, **_kw: "Sample metadata")
    monkeypatch.setattr(afterFastq, "_disk_usage_message", lambda *_a: "Disk space")
    monkeypatch.setattr(notifications, "send_notification", send)
    return cfg, store, output, now, messages


def restart(cfg, store, boundary):
    manager.rerun_flowcell(flowcell=RUN_ID, from_stage=boundary, force=True)
    saved = store.read(RUN_ID)
    cfg.run.begin(saved["source_path"], cfg.static.paths)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep"},
        cfg.output_path / "SampleSheet.csv",
        cfg.output_path / "Sample-Submission-Form.xlsx",
    )


def test_fresh_run_and_both_email_formats_have_exact_successful_totals(pipeline):
    cfg, store, _output, _now, messages = pipeline
    run(cfg, store)
    assert store.read(RUN_ID)["status"] == "completed"
    assert [entry["kind"] for entry, _message in messages] == [
        "sequencing",
        "processed",
        "finalized",
    ]
    early, processed, finalized = messages
    assert "Total processing time" not in early[1].get_body(preferencelist=("plain",)).get_content()
    for mime in ("plain", "html"):
        text = processed[1].get_body(preferencelist=(mime,)).get_content()
        assert "Reporting: 0:00:11" in text
        assert "Archiving: Not completed" in text
        assert "Archive MD5 checksums: Not completed" in text
        assert "Total processing time (subtotal of completed steps): 0:00:30" in text
        assert "elapsed time" not in text
    final_text = finalized[1].get_content()
    for label, seconds in zip(
        (s[0] for s in processing_times.STEPS.values()), [12, 3, 4, 11, 7, 8], strict=True
    ):
        assert f"{label}: 0:00:{seconds:02}" in final_text
    assert "Total processing time: 0:00:45" in final_text
    assert "md5sum and 7zip runtime" not in final_text
    assert "Total runtime" not in final_text


def test_analysis_rerun_across_daemon_restart_preserves_inputs_and_replaces_success(pipeline):
    cfg, store, output, now, messages = pipeline
    run(cfg, store)
    before = store.read(RUN_ID)
    manifest = output / f"md5sum_{PROJECT}_fastq.txt"
    unchanged = (manifest.read_bytes(), manifest.stat().st_mtime_ns)
    now[0] += 90000
    store = FlowcellStateStore(store.manager_dir)
    restart(cfg, store, "analysis")
    run(cfg, store)
    state = store.read(RUN_ID)
    for step in ("demultiplexing", "fastq_checksums"):
        assert state["processing_timings"][step] == before["processing_timings"][step]
    for step in ("analysis", "reporting", "archiving", "archive_checksums"):
        records = state["processing_timings"][step]
        assert len(records) == 2 and records[0]["invalidated_at"]
        assert records[-1]["attempt"] == 2
    assert (manifest.read_bytes(), manifest.stat().st_mtime_ns) == unchanged
    assert messages[-1][0]["payload"]["processing_timing"]["total_seconds"] == 45


@pytest.mark.parametrize("boundary", ["demultiplexing", "reporting", "finalization"])
def test_other_reruns_match_actual_output_cleanup(pipeline, boundary):
    cfg, store, _output, _now, messages = pipeline
    run(cfg, store)
    before = store.read(RUN_ID)
    restart(cfg, store, boundary)
    run(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    for step, (_, stage) in processing_times.STEPS.items():
        records = state["processing_timings"][step]
        if processing_times.STAGES.index(stage) >= processing_times.STAGES.index(boundary):
            assert len(records) == 2 and records[0]["invalidated_at"]
        else:
            assert records == before["processing_timings"][step]
    assert messages[-1][0]["payload"]["processing_timing"]["total_seconds"] == 45


def test_failed_analysis_replacement_then_success_excludes_failure_and_idle_time(
    pipeline, monkeypatch
):
    cfg, store, _output, now, messages = pipeline
    run(cfg, store)
    successful = afterFastq.analysis_steps
    restart(cfg, store, "analysis")

    def failure():
        now[0] += 700
        raise RuntimeError("controlled workflow failure")

    monkeypatch.setattr(afterFastq, "analysis_steps", failure)
    run(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "failed"
    assert processing_times.snapshot(state)["total_seconds"] == 21  # conversion/hash/early QC
    assert processing_times.snapshot(state)["steps"]["analysis"]["status"] == "not_completed"
    assert state["processing_timings"]["analysis"][-1]["duration_seconds"] == 700
    now[0] += 50000
    monkeypatch.setattr(afterFastq, "analysis_steps", successful)
    restarted = FlowcellStateStore(store.manager_dir)
    restart(cfg, restarted, "analysis")
    run(cfg, restarted)
    assert messages[-1][0]["payload"]["processing_timing"]["total_seconds"] == 45
    assert [r["outcome"] for r in restarted.read(RUN_ID)["processing_timings"]["analysis"]] == [
        "completed",
        "failed",
        "completed",
    ]


def test_archive_success_persists_before_checksum_failure_and_is_invalidated_on_rerun(
    pipeline, monkeypatch
):
    cfg, store, _output, now, messages = pipeline
    successful = afterFastq.md5sum_archive_worker

    def failure(_cfg):
        now[0] += 300
        raise OSError("checksum failure")

    monkeypatch.setattr(afterFastq, "md5sum_archive_worker", failure)
    run(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "failed"
    assert state["processing_timings"]["archiving"][-1]["outcome"] == "completed"
    assert state["processing_timings"]["archive_checksums"][-1]["duration_seconds"] == 300
    assert processing_times.snapshot(state)["total_seconds"] == 37
    restart(cfg, store, "finalization")
    assert processing_times.snapshot(store.read(RUN_ID))["total_seconds"] == 30
    monkeypatch.setattr(afterFastq, "md5sum_archive_worker", successful)
    run(cfg, store)
    assert messages[-1][0]["payload"]["processing_timing"]["total_seconds"] == 45


def test_reused_hashes_keep_timing_and_repair_records_actual_work(pipeline):
    cfg, store, output, _now, _messages = pipeline
    run(cfg, store)
    before = store.read(RUN_ID)["processing_timings"]["fastq_checksums"]
    (output / f"md5sum_{PROJECT}_fastq.txt").unlink()
    restart(cfg, store, "analysis")
    run(cfg, store)
    records = store.read(RUN_ID)["processing_timings"]["fastq_checksums"]
    assert len(records) == len(before) + 1
    assert records[0]["invalidated_at"]
    assert records[-1]["duration_seconds"] == 3


def test_early_failure_and_recovery_contribute_only_success(pipeline, monkeypatch):
    cfg, store, _output, now, _messages = pipeline
    successful = sequencing_qc.generate

    def failure(*_a, **_kw):
        now[0] += 400
        raise OSError("QC tool failed")

    monkeypatch.setattr(sequencing_qc, "generate", failure)
    run(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "completed"  # Existing non-blocking early QC behavior.
    assert state["sequencing_qc"]["attempts"][-1]["duration_seconds"] == 400
    assert state["sequencing_qc"]["attempts"][-1]["attempt"] == 1
    partial = processing_times.snapshot(state)
    assert partial["total_seconds"] == 39
    assert partial["steps"]["reporting"] == {"status": "not_completed", "duration_seconds": 5}
    monkeypatch.setattr(sequencing_qc, "generate", successful)
    sequencing_delivery.ensure_report(cfg, store, RUN_ID, retry=True)
    assert processing_times.snapshot(store.read(RUN_ID))["total_seconds"] == 45


def test_notification_retry_reuses_saved_timing_snapshot(pipeline, monkeypatch):
    cfg, store, _output, now, _messages = pipeline
    monkeypatch.setattr(
        notifications,
        "send_notification",
        lambda *_a: (_ for _ in ()).throw(OSError("SMTP offline")),
    )
    run(cfg, store)
    before = store.read(RUN_ID)
    now[0] += 90000
    delivered = []
    monkeypatch.setattr(
        notifications, "send_notification", lambda _cfg, entry: delivered.append(entry["payload"])
    )
    notification_delivery.deliver_pending(cfg, store, RUN_ID, retry=True)
    assert delivered == [e["payload"] for e in before["delivery_notifications"]]
    assert store.read(RUN_ID)["processing_timings"] == before["processing_timings"]
