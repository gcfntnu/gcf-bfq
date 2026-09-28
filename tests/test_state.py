import gzip
import json
import multiprocessing as mp
from pathlib import Path

import pytest

from bcl2fastq_pipeline.state import (
    ExecutionLeaseError,
    FlowcellStateStore,
    StateConflictError,
    StateValidationError,
    apply_cleanup,
    cleanup_plan,
    new_state,
    validate_restored_fastqs,
    validate_state,
)


RUN_ID = "260918_MN00686_0026_A000HCMFHF"


def make_store(tmp_path):
    return FlowcellStateStore(tmp_path / "manager")


def make_state(tmp_path, **kwargs):
    return new_state(
        RUN_ID,
        tmp_path / "instrument" / RUN_ID,
        tmp_path / "output" / RUN_ID,
        origin=kwargs.pop("origin", "new"),
        start_stage=kwargs.pop("start_stage", "demultiplexing"),
        **kwargs,
    )


def write_restored_inputs(output):
    output.mkdir(parents=True)
    (output / "SampleSheet.csv").write_text(
        "[CustomOptions]\nLibprep,Illumina DNA Prep\nUser,test\n",
        encoding="utf-8",
    )
    (output / "Sample-Submission-Form.xlsx").write_bytes(b"placeholder")


def write_fastq(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write("@read\nACGT\n+\n!!!!\n")


def _hold_execution_lease(manager_dir, started, release):
    store = FlowcellStateStore(manager_dir)
    with store.execution_lease(RUN_ID):
        started.set()
        release.wait(10)


def test_state_round_trip_and_schema(tmp_path):
    store = make_store(tmp_path)
    state = make_state(tmp_path)

    store.create(state)

    loaded = store.read(RUN_ID)
    assert loaded["schema_version"] == 1
    assert loaded["run_id"] == RUN_ID
    assert loaded["origin"] == "new"
    assert loaded["status"] == "queued"
    assert loaded["current_stage"] == "demultiplexing"
    assert loaded["stages"]["demultiplexing"]["status"] == "queued"
    assert loaded["notification"]["notified"] is False


def test_schema_rejects_missing_and_invalid_fields(tmp_path):
    state = make_state(tmp_path)
    del state["notification"]
    with pytest.raises(StateValidationError, match="notification"):
        validate_state(state)

    state = make_state(tmp_path)
    state["schema_version"] = 99
    with pytest.raises(StateValidationError, match="Unsupported state schema"):
        validate_state(state)


def test_create_refuses_overwrite_and_atomic_write_leaves_no_temp(tmp_path):
    store = make_store(tmp_path)
    state = make_state(tmp_path)
    store.create(state)

    with pytest.raises(StateConflictError, match="already exists"):
        store.create(state)

    saved = json.loads(store.state_path(RUN_ID).read_text())
    assert saved["run_id"] == RUN_ID
    assert not list(store.states_dir.glob("*.tmp"))


def test_execution_lease_excludes_another_process(tmp_path):
    store = make_store(tmp_path)
    started = mp.Event()
    release = mp.Event()
    process = mp.Process(
        target=_hold_execution_lease,
        args=(store.manager_dir, started, release),
    )
    process.start()
    try:
        assert started.wait(5)
        assert store.execution_active(RUN_ID)
        with pytest.raises(ExecutionLeaseError):
            with store.execution_lease(RUN_ID):
                pass
    finally:
        release.set()
        process.join(5)
        if process.is_alive():
            process.kill()
    assert not store.execution_active(RUN_ID)


def test_restored_state_infers_demultiplexing_completed(tmp_path):
    state = make_state(tmp_path, origin="restored_legacy_fastq", start_stage="analysis")

    assert state["stages"]["demultiplexing"]["status"] == "completed"
    assert state["stages"]["demultiplexing"]["inferred"] is True
    assert state["stages"]["analysis"]["status"] == "queued"


def test_attempt_and_stage_transitions_preserve_history(tmp_path):
    store = make_store(tmp_path)
    store.create(make_state(tmp_path))

    state = store.begin_attempt(RUN_ID)
    assert state["attempt"] == 1
    assert state["status"] == "running"
    assert state["stages"]["demultiplexing"]["status"] == "running"

    state = store.complete_stage(
        RUN_ID,
        "demultiplexing",
        {"tool": "bcl-convert", "version": "4.2.4"},
    )
    assert state["demultiplexing"] == {"tool": "bcl-convert", "version": "4.2.4"}
    assert state["current_stage"] == "analysis"
    assert state["stages"]["analysis"]["status"] == "queued"

    store.start_stage(RUN_ID, "analysis")
    store.complete_stage(RUN_ID, "analysis")
    store.start_stage(RUN_ID, "reporting")
    store.complete_stage(RUN_ID, "reporting")
    store.start_stage(RUN_ID, "finalization")
    state = store.complete_run(RUN_ID, ["GCF-2026-002", "GCF-2026-001"])

    assert state["status"] == "completed"
    assert state["projects"] == ["GCF-2026-001", "GCF-2026-002"]
    assert state["attempts"][0]["outcome"] == "completed"
    assert state["attempts"][0]["completed_at"]


def test_illegal_transitions_are_refused(tmp_path):
    store = make_store(tmp_path)
    store.create(make_state(tmp_path))

    with pytest.raises(StateConflictError, match="not running"):
        store.complete_stage(RUN_ID, "demultiplexing")

    store.begin_attempt(RUN_ID)
    with pytest.raises(StateConflictError, match="not queued"):
        store.start_stage(RUN_ID, "analysis")


def test_failure_records_signature_report_and_attempt(tmp_path):
    store = make_store(tmp_path)
    store.create(make_state(tmp_path))
    store.begin_attempt(RUN_ID)

    state = store.fail_stage(
        RUN_ID,
        "demultiplexing",
        summary="converter failed",
        report_path=tmp_path / "reports" / f"{RUN_ID}.error",
    )

    assert state["status"] == "failed"
    assert state["last_error"]["summary"] == "converter failed"
    assert state["last_error"]["report_path"].endswith(f"{RUN_ID}.error")
    assert state["last_error"]["failure_signature"]
    assert state["notification"]["failure_signature"] == state["last_error"]["failure_signature"]
    assert state["notification"]["notified"] is False
    assert state["attempts"][0]["outcome"] == "failed"


def test_stale_running_state_is_marked_interrupted(tmp_path):
    store = make_store(tmp_path)
    store.create(make_state(tmp_path))
    store.begin_attempt(RUN_ID)

    state = store.recover_interrupted(RUN_ID)

    assert state["status"] == "interrupted"
    assert state["stages"]["demultiplexing"]["status"] == "failed"
    assert "explicit rerun required" in state["last_error"]["summary"]
    assert state["attempts"][0]["outcome"] == "interrupted"


def test_active_running_state_is_not_recovered(tmp_path):
    store = make_store(tmp_path)
    store.create(make_state(tmp_path))
    store.begin_attempt(RUN_ID)

    with store.execution_lease(RUN_ID):
        state = store.recover_interrupted(RUN_ID)

    assert state["status"] == "running"


def test_rerun_preparation_and_queue_records_operator_request(tmp_path):
    store = make_store(tmp_path)
    store.create(make_state(tmp_path))
    store.begin_attempt(RUN_ID)
    store.complete_stage(RUN_ID, "demultiplexing")
    store.start_stage(RUN_ID, "analysis")
    store.fail_stage(RUN_ID, "analysis", summary="bad", report_path=None)

    state = store.set_preparing(
        RUN_ID,
        "analysis",
        reason="corrected sample sheet",
        refresh_inputs=True,
        hostname="operator-host",
    )
    assert state["status"] == "preparing"
    assert state["restart_request"]["from"] == "analysis"
    assert state["restart_request"]["hostname"] == "operator-host"
    assert state["restart_request"]["reason"] == "corrected sample sheet"
    assert state["restart_request"]["refresh_inputs"] is True

    state = store.queue(RUN_ID, "analysis")
    assert state["status"] == "queued"
    assert state["stages"]["analysis"]["status"] == "queued"
    assert state["stages"]["reporting"]["status"] == "pending"


def test_queue_rejects_restart_past_failed_upstream_stage(tmp_path):
    store = make_store(tmp_path)
    store.create(make_state(tmp_path))
    store.begin_attempt(RUN_ID)
    store.fail_stage(RUN_ID, "demultiplexing", summary="bad", report_path=None)

    with pytest.raises(StateConflictError, match="upstream stages are not complete"):
        store.queue(RUN_ID, "analysis")


def test_validate_restored_fastqs_accepts_recognized_tree(tmp_path):
    output = tmp_path / RUN_ID
    write_restored_inputs(output)
    write_fastq(output / "GCF-2026-001" / "sample_R1.fastq.gz")

    recognized, detail = validate_restored_fastqs(output)

    assert recognized is True
    assert detail == "recognized restored FASTQs"


@pytest.mark.parametrize("problem", ["no-inputs", "no-project-fastq", "empty-fastq", "bad-gzip"])
def test_validate_restored_fastqs_rejects_incomplete_trees(tmp_path, problem):
    output = tmp_path / RUN_ID
    if problem != "no-inputs":
        write_restored_inputs(output)
    else:
        output.mkdir()

    if problem == "no-project-fastq":
        write_fastq(output / "Undetermined_R1.fastq.gz")
    elif problem == "empty-fastq":
        path = output / "GCF-2026-001" / "sample_R1.fastq.gz"
        path.parent.mkdir(parents=True)
        path.touch()
    elif problem == "bad-gzip":
        path = output / "GCF-2026-001" / "sample_R1.fastq.gz"
        path.parent.mkdir(parents=True)
        path.write_text("not gzip")

    recognized, _detail = validate_restored_fastqs(output)
    assert recognized is False


def test_cleanup_from_demultiplexing_preserves_curated_inputs(tmp_path):
    output = tmp_path / RUN_ID
    write_restored_inputs(output)
    write_fastq(output / "GCF-2026-001" / "sample_R1.fastq.gz")
    (output / "QC_GCF-2026-001").mkdir()
    (output / "GCF-2026-001_260918.7za").touch()

    paths = cleanup_plan(output, "demultiplexing")
    apply_cleanup(paths)

    assert (output / "SampleSheet.csv").exists()
    assert (output / "Sample-Submission-Form.xlsx").exists()
    assert not (output / "GCF-2026-001").exists()
    assert not (output / "QC_GCF-2026-001").exists()
    assert not (output / "GCF-2026-001_260918.7za").exists()


def test_cleanup_boundaries_preserve_expected_products(tmp_path):
    output = tmp_path / RUN_ID
    write_restored_inputs(output)
    fastq = output / "GCF-2026-001" / "sample_R1.fastq.gz"
    write_fastq(fastq)
    qc = output / "QC_GCF-2026-001"
    qc.mkdir()
    report = output / "multiqc_GCF-2026-001_260918.html"
    report.touch()
    archive = output / "GCF-2026-001_260918.7za"
    archive.touch()
    checksum = output / "md5sum_GCF-2026-001_fastq.txt"
    checksum.touch()

    analysis_paths = cleanup_plan(output, "analysis")
    assert qc in analysis_paths
    assert report in analysis_paths
    assert archive in analysis_paths
    assert checksum in analysis_paths
    assert fastq not in analysis_paths

    reporting_paths = cleanup_plan(output, "reporting")
    assert qc not in reporting_paths
    assert report in reporting_paths
    assert archive in reporting_paths
    assert fastq not in reporting_paths

    final_paths = cleanup_plan(output, "finalization")
    assert qc not in final_paths
    assert report not in final_paths
    assert archive in final_paths
    assert checksum in final_paths


def test_state_survives_output_cleanup(tmp_path):
    store = make_store(tmp_path)
    state = make_state(tmp_path)
    store.create(state)
    output = Path(state["output_path"])
    output.mkdir(parents=True)
    (output / "throwaway").write_text("data")

    apply_cleanup(cleanup_plan(output, "demultiplexing"))

    assert store.read(RUN_ID)["run_id"] == RUN_ID
    assert store.state_path(RUN_ID).exists()


def test_invalid_run_id_cannot_escape_state_directory(tmp_path):
    with pytest.raises(StateValidationError):
        new_state(
            "../escape",
            tmp_path / "source",
            tmp_path / "output",
            origin="new",
            start_stage="demultiplexing",
        )
