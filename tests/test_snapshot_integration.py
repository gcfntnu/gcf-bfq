"""Retained snapshots through the real daemon and operator state transitions."""

import shutil

from unittest.mock import Mock

import flowcell_manager.flowcell_manager as manager
import pytest

from bcl2fastq_pipeline.state import cleanup_plan
from test_analysis_snapshots import PROJECT, archive_path, put, read_member, ready
from test_notification_integration import prepare_pipeline, run_pipeline
from test_state_integration import RUN_ID, write_fastq, write_inputs

from bcl2fastq_pipeline import afterFastq, analysis_snapshots as snapshots, cli, notifications


def pipeline(tmp_path, monkeypatch):
    cfg, store, output, calls = prepare_pipeline(tmp_path, monkeypatch)
    cfg.run.pipeline = "test-workflow"
    monkeypatch.setenv("TMPDIR", str(tmp_path / "work"))
    monkeypatch.setattr(notifications, "send_notification", lambda *_args: None)

    def analyze():
        calls.append("analysis")
        work = snapshots.workdir_path(RUN_ID, PROJECT)
        put(work, "config.yaml", "delivered config")
        put(work, "Snakefile", "delivered workflow")
        put(work, ".snakemake/log/run.log", "completed log")
        put(work, "data/fastq", "excluded")
        return {PROJECT: snapshots.identify_workdir(work, RUN_ID, PROJECT)}

    monkeypatch.setattr(afterFastq, "analysis_steps", analyze)
    return cfg, store, output, calls


def restart(cfg, stage):
    manager.rerun_flowcell(flowcell=RUN_ID, from_stage=stage, force=True)
    cfg.run.begin(cfg.static.paths.nova_base_dir / RUN_ID, cfg.static.paths)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep", "User": "test"},
        cfg.output_path / "SampleSheet.csv",
        cfg.output_path / "Sample-Submission-Form.xlsx",
    )
    cfg.run.pipeline = "test-workflow"


def test_daemon_publishes_success_despite_email_failure(tmp_path, monkeypatch):
    cfg, store, output, calls = pipeline(tmp_path, monkeypatch)
    monkeypatch.setattr(notifications, "send_notification", Mock(side_effect=OSError("SMTP down")))
    run_pipeline(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert calls == ["analysis", "reporting", "finalization"]
    assert read_member(output) == "delivered config"
    assert {entry["status"] for entry in state["delivery_notifications"]} == {"failed"}
    assert "pending" not in state["analysis_snapshots"][PROJECT]
    assert (
        "pending" not in state["stages"]["finalization"]["metadata"]["analysis_snapshots"][PROJECT]
    )


@pytest.mark.parametrize("failure", ["analysis", "reporting", "finalization", "snapshot", "commit"])
def test_failed_rerun_preserves_previous_delivery_snapshot(tmp_path, monkeypatch, failure):
    cfg, store, output, _calls = pipeline(tmp_path, monkeypatch)
    run_pipeline(cfg, store)
    before = archive_path(output).read_bytes()
    old_identity = store.read(RUN_ID)["analysis_snapshots"][PROJECT]
    restart(cfg, "analysis")
    error = Mock(side_effect=OSError("injected failure"))
    if failure == "analysis":
        monkeypatch.setattr(afterFastq, "analysis_steps", error)
    elif failure == "reporting":
        monkeypatch.setattr(cli, "_run_reporting", error)
    elif failure == "finalization":
        monkeypatch.setattr(afterFastq, "finalize", error)
    elif failure == "snapshot":
        monkeypatch.setattr(snapshots, "_write_archive", error)
    else:
        monkeypatch.setattr(store, "complete_run", error)
    run_pipeline(cfg, store)
    assert store.read(RUN_ID)["status"] == "failed"
    assert store.read(RUN_ID)["analysis_snapshots"][PROJECT] == old_identity
    assert archive_path(output).read_bytes() == before
    assert not list((output / "provenance/.staging").glob("*"))


@pytest.mark.parametrize("stage", ["reporting", "finalization"])
def test_downstream_retry_reuses_matching_snapshot_without_workdir(tmp_path, monkeypatch, stage):
    cfg, store, output, calls = pipeline(tmp_path, monkeypatch)
    run_pipeline(cfg, store)
    before = archive_path(output).read_bytes(), archive_path(output).stat().st_mtime_ns
    old_analysis = store.read(RUN_ID)["stages"]["analysis"]
    shutil.rmtree(snapshots.workdir_path(RUN_ID, PROJECT))
    # A match must bypass source lookup entirely, irrespective of installed code.
    monkeypatch.setattr(snapshots, "_source", Mock(side_effect=AssertionError("source lookup")))
    restart(cfg, stage)
    run_pipeline(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert state["stages"]["analysis"] == old_analysis
    assert calls.count("analysis") == 1
    assert (archive_path(output).read_bytes(), archive_path(output).stat().st_mtime_ns) == before


@pytest.mark.parametrize("stage", ["demultiplexing", "analysis", "reporting", "finalization"])
def test_operator_cleanup_preserves_provenance_at_every_boundary(tmp_path, monkeypatch, stage):
    cfg, store, output, _calls = pipeline(tmp_path, monkeypatch)
    run_pipeline(cfg, store)
    before = archive_path(output).read_bytes()
    assert all(output / "provenance" != path for path in cleanup_plan(output, stage))
    restart(cfg, stage)
    assert archive_path(output).read_bytes() == before
    assert store.read(RUN_ID)["analysis_snapshots"]


def test_archive_recovers_committed_snapshot_before_removing_delivery_data(tmp_path, monkeypatch):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    write_inputs(output)
    write_fastq(output / PROJECT / "sample_R1.fastq.gz")
    put(output, f"{PROJECT}_260918.7za")
    results = snapshots.prepare(store.read(RUN_ID), [PROJECT])
    store.complete_run(RUN_ID, [PROJECT], snapshots=results)
    protected = put(output, "provenance/other.fastq.gz", "protected namespace")
    state = manager.archive_flowcell(flowcell=RUN_ID, force=True)
    assert state["status"] == "archived"
    assert read_member(output) == "original"
    assert protected.exists()
    assert not (output / PROJECT).exists()
    assert not list(output.glob("*.7za"))
    assert "pending" not in state["analysis_snapshots"][PROJECT]


def test_legacy_run_without_workdir_reports_unavailable(tmp_path, monkeypatch, caplog):
    cfg, store, output, _calls = prepare_pipeline(tmp_path, monkeypatch, start_stage="finalization")
    monkeypatch.setenv("TMPDIR", str(tmp_path / "absent"))
    monkeypatch.setattr(notifications, "send_notification", lambda *_args: None)
    run_pipeline(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert (
        state["stages"]["finalization"]["metadata"]["analysis_snapshots"][PROJECT]["status"]
        == "unavailable"
    )
    assert "Analysis snapshot unavailable" in caplog.text
    assert not archive_path(output).exists()


def test_postcommit_state_error_does_not_fail_committed_delivery(tmp_path, monkeypatch, caplog):
    cfg, store, output, _calls = pipeline(tmp_path, monkeypatch)
    complete = store.complete_run

    def committed_then_error(*args, **kwargs):
        complete(*args, **kwargs)
        raise OSError("directory fsync failed after rename")

    monkeypatch.setattr(store, "complete_run", committed_then_error)
    run_pipeline(cfg, store)
    assert store.read(RUN_ID)["status"] == "completed"
    assert read_member(output) == "delivered config"
    assert "Finalization was committed" in caplog.text


def test_publication_failure_keeps_completion_and_recovers_next_scan(tmp_path, monkeypatch):
    cfg, store, output, _calls = pipeline(tmp_path, monkeypatch)
    replace = snapshots.os.replace

    def fail_publication(src, dst):
        if str(dst).endswith("_analysis.tar.gz"):
            raise OSError("publication unavailable")
        return replace(src, dst)

    with monkeypatch.context() as patch:
        patch.setattr(snapshots.os, "replace", fail_publication)
        run_pipeline(cfg, store)
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert (output / state["analysis_snapshots"][PROJECT]["pending"]).exists()
    assert not archive_path(output).exists()
    snapshots.recover_pending(store)
    assert read_member(output) == "delivered config"


def test_active_execution_prevents_recovery_from_discarding_uncommitted_archive(
    tmp_path, monkeypatch
):
    _cfg, store, output = ready(tmp_path, monkeypatch)
    with store.execution_lease(RUN_ID):
        result = snapshots.prepare(store.read(RUN_ID), [PROJECT])
        snapshots.recover_pending(store)
        assert (output / result[PROJECT]["pending"]).exists()
    snapshots.recover_pending(store)
    assert not list((output / "provenance/.staging").glob("*"))


def test_new_success_replaces_one_snapshot_after_failed_attempt(tmp_path, monkeypatch):
    cfg, store, output, _calls = pipeline(tmp_path, monkeypatch)
    run_pipeline(cfg, store)
    first = store.read(RUN_ID)["analysis_snapshots"][PROJECT]["analysis_id"]
    restart(cfg, "analysis")
    with monkeypatch.context() as patch:
        patch.setattr(afterFastq, "finalize", Mock(side_effect=OSError("delivery failure")))
        run_pipeline(cfg, store)
    assert store.read(RUN_ID)["status"] == "failed"
    assert store.read(RUN_ID)["analysis_snapshots"][PROJECT]["analysis_id"] == first
    restart(cfg, "finalization")
    run_pipeline(cfg, store)
    assert store.read(RUN_ID)["analysis_snapshots"][PROJECT]["analysis_id"] != first
    assert len(list((output / "provenance").glob("*.tar.gz"))) == 1
    assert not list((output / "provenance/.staging").glob("*"))


def test_status_reports_retained_unavailable_and_pending_snapshots(tmp_path, monkeypatch, capsys):
    cfg, store, output, _calls = pipeline(tmp_path, monkeypatch)
    run_pipeline(cfg, store)
    manager.status_flowcell(flowcell=RUN_ID)
    assert f"Analysis snapshot {PROJECT}: retained (provenance/" in capsys.readouterr().out
    restart(cfg, "analysis")

    # Create a valid new analysis, then lose its workdir before finalization.
    def missing_workdir():
        shutil.rmtree(snapshots.workdir_path(RUN_ID, PROJECT))

    monkeypatch.setattr(afterFastq, "finalize", missing_workdir)
    run_pipeline(cfg, store)
    manager.status_flowcell(flowcell=RUN_ID)
    status = capsys.readouterr().out
    assert f"Analysis snapshot {PROJECT}: unavailable" in status
    assert "retained from previous successful finalization" in status
    assert archive_path(output).exists()
    state = store.read(RUN_ID)
    state["analysis_snapshots"][PROJECT]["pending"] = "provenance/.staging/pending.tar.gz"
    store.write(state)
    manager.status_flowcell(flowcell=RUN_ID)
    assert "publication pending; BFQ will retry" in capsys.readouterr().out
