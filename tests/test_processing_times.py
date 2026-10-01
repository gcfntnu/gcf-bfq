"""Controlled-clock timing semantics; no worker processes or SMTP."""

import copy

import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore, StateValidationError, validate_state
from test_state import RUN_ID, make_state, make_store

from bcl2fastq_pipeline import processing_times as timing


@pytest.fixture
def clock(monkeypatch):
    now = [100.0]
    monkeypatch.setattr(timing.time, "monotonic", lambda: now[0])
    return now


@pytest.fixture
def active(tmp_path):
    store = make_store(tmp_path)
    store.create(make_state(tmp_path))
    store.begin_attempt(RUN_ID)
    return store


def record(store, clock, step, seconds):
    with timing.measure(store, RUN_ID, step):
        clock[0] += seconds


def test_six_successes_survive_restart_and_exclude_idle_time(active, clock):
    for seconds, step in enumerate(timing.STEPS, 1):
        clock[0] += 1000
        record(active, clock, step, seconds)
    restarted = FlowcellStateStore(active.manager_dir)
    summary = timing.snapshot(restarted.read(RUN_ID))
    assert summary["total_seconds"] == 21
    assert not summary["incomplete"] and not summary["unavailable"]
    assert [item["duration_seconds"] for item in summary["steps"].values()] == [1, 2, 3, 4, 5, 6]
    text = timing.format_summary(summary)
    assert len(text.splitlines()) == 7
    assert text.endswith("Total processing time: 0:00:21")
    for records in restarted.read(RUN_ID)["processing_timings"].values():
        assert records[0]["attempt"] == 1
        assert records[0]["started_at"].endswith("+00:00")
        assert records[0]["completed_at"].endswith("+00:00")


@pytest.mark.parametrize("boundary", ["demultiplexing", "analysis", "reporting", "finalization"])
def test_restart_invalidation_preserves_upstream_and_history(active, clock, boundary):
    for step in timing.STEPS:
        record(active, clock, step, 10)
    active.set_preparing(RUN_ID, boundary, reason=None, refresh_inputs=False)
    state = active.read(RUN_ID)
    for step, (_, stage) in timing.STEPS.items():
        invalidated = timing.STAGES.index(stage) >= timing.STAGES.index(boundary)
        assert bool(state["processing_timings"][step][0].get("invalidated_at")) == invalidated
        assert timing.snapshot(state)["steps"][step]["status"] == (
            "not_completed" if invalidated else "completed"
        )


def test_failed_replacement_never_falls_back_to_old_success(active, clock):
    record(active, clock, "demultiplexing", 30)
    active.complete_stage(RUN_ID, "demultiplexing")
    active.start_stage(RUN_ID, "analysis")
    record(active, clock, "analysis", 80)
    active.fail_stage(RUN_ID, "analysis", summary="later failure", report_path=None)
    active.queue(RUN_ID, "analysis")
    active.begin_attempt(RUN_ID)
    with pytest.raises(RuntimeError, match="failed work"):
        with timing.measure(active, RUN_ID, "analysis"):
            clock[0] += 900
            raise RuntimeError("failed work")
    snapshot = timing.snapshot(active.read(RUN_ID))
    assert snapshot["total_seconds"] == 30
    assert snapshot["steps"]["analysis"]["status"] == "not_completed"
    record(active, clock, "analysis", 5)
    state = active.read(RUN_ID)
    assert timing.snapshot(state)["total_seconds"] == 35
    assert [r["duration_seconds"] for r in state["processing_timings"]["analysis"]] == [80, 900, 5]
    assert [r["outcome"] for r in state["processing_timings"]["analysis"]] == [
        "completed",
        "failed",
        "completed",
    ]
    assert state["processing_timings"]["analysis"][-1]["attempt"] == 2


def test_success_is_durable_before_another_step_in_same_stage_fails(active, clock):
    record(active, clock, "demultiplexing", 7)
    with pytest.raises(OSError):
        with timing.measure(active, RUN_ID, "fastq_checksums"):
            clock[0] += 11
            raise OSError("hash failed")
    active.fail_stage(RUN_ID, "demultiplexing", summary="hash failed", report_path=None)
    state = active.read(RUN_ID)
    assert timing.snapshot(state)["total_seconds"] == 7
    assert state["processing_timings"]["fastq_checksums"][0]["duration_seconds"] == 11
    active.queue(RUN_ID, "demultiplexing")
    assert timing.snapshot(active.read(RUN_ID))["total_seconds"] == 0


def test_hard_crash_has_no_recovery_duration(active, clock):
    record(active, clock, "demultiplexing", 3)
    active.mutate(RUN_ID, lambda state: timing.start(state, "fastq_checksums"))
    clock[0] += 90000
    state = active.recover_interrupted(RUN_ID)
    entry = state["processing_timings"]["fastq_checksums"][-1]
    assert entry["outcome"] == "interrupted"
    assert entry["duration_seconds"] is None
    assert entry["completed_at"] is None
    assert entry["recovered_at"]
    assert timing.snapshot(state)["total_seconds"] == 3


def test_controlled_interruption_keeps_diagnostic_duration(active, clock):
    with pytest.raises(KeyboardInterrupt):
        with timing.measure(active, RUN_ID, "demultiplexing"):
            clock[0] += 13
            raise KeyboardInterrupt
    entry = active.read(RUN_ID)["processing_timings"]["demultiplexing"][-1]
    assert entry["outcome"] == "interrupted"
    assert entry["duration_seconds"] == 13
    assert timing.snapshot(active.read(RUN_ID))["total_seconds"] == 0


def test_legacy_durations_are_unavailable_not_zero(tmp_path):
    state = make_state(tmp_path, origin="restored_legacy_fastq", start_stage="analysis")
    state.pop("processing_timings")
    validate_state(state)
    snapshot = timing.snapshot(state)
    assert snapshot["steps"]["demultiplexing"]["status"] == "unavailable"
    assert snapshot["steps"]["fastq_checksums"]["duration_seconds"] is None
    text = timing.format_summary(snapshot)
    assert "Demultiplexing: Timing unavailable" in text
    assert "Analysis: Not completed" in text
    assert "subtotal of completed steps; partial; some timings unavailable" in text
    timing.invalidate(state, "demultiplexing")
    assert timing.snapshot(state)["steps"]["demultiplexing"]["status"] == "not_completed"


def test_reporting_adds_early_success_once_and_marks_missing_parts(active, clock):
    active.mutate(
        RUN_ID,
        lambda state: {
            **state,
            "sequencing_qc": {
                "execution": 1,
                "status": "completed",
                "context": {},
                "attempts": [],
                "duration_seconds": 7,
                "result": {"report_path": "report.html"},
            },
        },
    )
    first = timing.snapshot(active.read(RUN_ID))
    assert first["steps"]["reporting"] == {"status": "not_completed", "duration_seconds": 7}
    record(active, clock, "reporting", 4)
    assert timing.snapshot(active.read(RUN_ID))["steps"]["reporting"] == {
        "status": "completed",
        "duration_seconds": 11,
    }
    active.mutate(RUN_ID, lambda state: timing.invalidate(state, "reporting"))
    assert timing.snapshot(active.read(RUN_ID))["total_seconds"] == 7
    record(active, clock, "reporting", 2)
    assert timing.snapshot(active.read(RUN_ID))["total_seconds"] == 9
    state = active.read(RUN_ID)
    state["sequencing_qc"]["status"] = "failed"
    assert timing.snapshot(state)["steps"]["reporting"] == {
        "status": "not_completed",
        "duration_seconds": 2,
    }


@pytest.mark.parametrize("bad", [True, -1, float("inf"), float("nan"), "12", None])
def test_validation_rejects_invalid_success_duration(active, clock, bad):
    record(active, clock, "analysis", 1)
    state = copy.deepcopy(active.read(RUN_ID))
    state["processing_timings"]["analysis"][0]["duration_seconds"] = bad
    with pytest.raises(StateValidationError, match="timing"):
        validate_state(state)
