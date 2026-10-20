"""Current project work must remain distinct from retained successful delivery."""

import copy

import pytest

from bcl2fastq_pipeline.state import new_state
from flowcell_manager.project_listing import (
    COLUMNS,
    UNAVAILABLE,
    UNIDENTIFIED,
    legacy_project_row,
    project_rows,
)

RUN = "990101_SYNTHETIC_0001_TEST"
A, B, C = (f"GCF-2099-{number:03}" for number in range(1, 4))
COMPLETED = "2099-01-02T12:00:00+00:00"


@pytest.fixture
def state(tmp_path):
    result = new_state(
        RUN,
        tmp_path / "source" / RUN,
        tmp_path / "output" / RUN,
        origin="new",
        start_stage="analysis",
    )
    result.update(status="running", attempt=2, projects=[A, B, C])
    result["stages"]["analysis"]["status"] = "running"
    return result


def analysis(status, generation="new", **extra):
    return {"status": status, "generation": generation, "attempt": 2, **extra}


def completion(generation="old", **extra):
    return {"attempt": 1, "generation": generation, "completed_at": COMPLETED, **extra}


def by_project(state):
    return {row["project"]: row for row in project_rows(state)}


def test_sibling_failure_does_not_invalidate_completed_delivery_or_start_waiting_work(state):
    state.update(status="failed")
    state["project_progress"] = {
        A: {
            "analysis": analysis("completed"),
            "finalization": {"status": "completed", "generation": "new", "completed_at": COMPLETED},
        },
        B: {"analysis": analysis("failed")},
    }
    rows = by_project(state)
    assert (rows[A]["status"], rows[A]["analysis"], rows[A]["delivery"]) == (
        "completed",
        "completed",
        "completed",
    )
    assert (rows[B]["status"], rows[B]["stage"]) == ("failed", "analysis")
    assert (rows[C]["status"], rows[C]["analysis"], rows[C]["delivery"]) == (
        "pending",
        "pending",
        "pending",
    )


@pytest.mark.parametrize("parent", ["running", "failed", "interrupted"])
def test_analysis_commits_wait_for_their_own_archive_start(state, parent):
    state.update(status=parent, current_stage="finalization")
    state["stages"]["reporting"]["status"] = "completed"
    state["project_progress"] = {
        A: {"analysis": analysis("completed")},
        B: {"analysis": analysis("completed"), "finalization": {"status": "failed"}},
    }
    rows = by_project(state)
    assert (rows[A]["status"], rows[A]["stage"], rows[A]["delivery"]) == (
        "pending",
        "finalization",
        "pending",
    )
    assert rows[B]["status"] == ("interrupted" if parent == "interrupted" else "failed")


def test_interruption_maps_only_the_affected_attempt(state):
    state["status"] = "interrupted"
    state["project_progress"] = {
        A: {"analysis": analysis("completed")},
        B: {"analysis": analysis("failed")},
        C: {"analysis": analysis("failed", attempt=1)},
    }
    rows = by_project(state)
    assert rows[A]["analysis"] == "completed"
    assert rows[B]["status"] == rows[B]["analysis"] == "interrupted"
    assert rows[C]["status"] == rows[C]["analysis"] == "failed"


@pytest.mark.parametrize("parent", ["running", "failed", "interrupted"])
def test_scoped_reporting_outcome_does_not_apply_to_unselected_analyzed_projects(state, parent):
    state.update(
        status=parent,
        current_stage="reporting",
        restart_request={"from": "analysis", "projects": [B]},
    )
    state["stages"]["analysis"]["status"] = "completed"
    state["stages"]["reporting"]["status"] = "running" if parent == "running" else "failed"
    state["project_progress"] = {
        A: {"analysis": analysis("completed", generation="retained", attempt=1)},
        B: {"analysis": analysis("completed")},
    }
    rows = by_project(state)
    assert (rows[A]["status"], rows[A]["stage"], rows[A]["analysis"], rows[A]["delivery"]) == (
        "pending",
        "reporting",
        "completed",
        "pending",
    )
    assert (rows[B]["status"], rows[B]["stage"], rows[B]["analysis"], rows[B]["delivery"]) == (
        parent,
        "reporting",
        "completed",
        "pending",
    )


@pytest.mark.parametrize("parent", ["preparing", "queued", "running", "failed"])
def test_scoped_reanalysis_keeps_old_delivery_separate_and_retained_project_complete(state, parent):
    state.update(status=parent, restart_request={"from": "analysis", "projects": [B]})
    state["project_completion"] = {A: completion(), B: completion()}
    state["project_progress"] = {B: {}}
    if parent in {"running", "failed"}:
        state["project_progress"][B]["analysis"] = analysis(parent)
    rows = by_project(state)
    assert (rows[A]["status"], rows[A]["analysis_generation"]) == ("completed", "old")
    assert rows[B]["status"] == parent
    assert rows[B]["delivery"] == "pending"
    assert rows[B]["last_delivery"] == COMPLETED
    assert rows[B]["last_delivery_generation"] == "old"
    assert rows[B]["analysis_generation"] == (
        "new" if parent in {"running", "failed"} else UNAVAILABLE
    )


def test_prior_delivery_generation_cannot_complete_new_analysis_even_with_aggregate_success(state):
    state["stages"]["analysis"]["status"] = "completed"
    state["stages"]["finalization"]["status"] = "completed"
    state["project_progress"] = {
        B: {
            "analysis": analysis("completed"),
            "finalization": {
                "status": "completed",
                "generation": "old",
                "completed_at": COMPLETED,
            },
        },
    }
    row = by_project(state)[B]
    assert row["analysis"] == "completed"
    assert row["delivery"] == "pending"
    assert row["status"] == "pending"
    assert row["analysis_generation"] == "new"
    assert row["last_delivery_generation"] == "old"


def test_frozen_priority_uses_producing_project_record_and_current_plan_for_unstarted_projects(
    state,
):
    explicit = {"priority": "0.12345678901234567890", "priority_source": "samplesheet"}
    default = {"priority": "1.0", "priority_source": "default"}
    state["project_progress"] = {A: {"analysis": analysis("completed", scheduling=explicit)}}
    state["attempts"] = [
        {
            "attempt": 2,
            "project_schedule": {
                "projects": {
                    A: {"priority": "99", "priority_source": "samplesheet"},
                    B: default,
                }
            },
        }
    ]
    rows = by_project(state)
    assert (rows[A]["priority"], rows[A]["priority_source"]) == (
        explicit["priority"],
        "samplesheet",
    )
    assert (rows[B]["priority"], rows[B]["priority_source"]) == ("1.0", "default")
    assert rows[C]["priority"] == rows[C]["priority_source"] == UNAVAILABLE


def test_new_queue_does_not_reuse_previous_attempt_priority(state):
    state.update(status="queued", restart_request={"from": "analysis"})
    state["attempts"] = [
        {
            "attempt": 2,
            "project_schedule": {
                "projects": {
                    A: {"priority": "9", "priority_source": "samplesheet"},
                }
            },
        }
    ]
    assert by_project(state)[A]["priority"] == UNAVAILABLE
    state["restart_request"]["analysis_resume"] = {"workflow": "synthetic"}
    assert by_project(state)[A]["priority"] == "9"


def test_retained_priority_does_not_use_new_correction_plan(state):
    state["restart_request"] = {"from": "analysis", "projects": [B]}
    state["project_completion"] = {
        A: completion(
            stages={
                "analysis": {
                    "metadata": {
                        "project_schedule": {
                            "projects": {A: {"priority": "0.5", "priority_source": "default"}}
                        },
                    }
                }
            }
        )
    }
    state["attempts"] = [
        {
            "attempt": 2,
            "project_schedule": {
                "projects": {
                    A: {"priority": "9", "priority_source": "samplesheet"},
                }
            },
        }
    ]
    assert by_project(state)[A]["priority"] == "0.5"
    del state["project_completion"][A]["stages"]
    assert by_project(state)[A]["priority"] == UNAVAILABLE


@pytest.mark.parametrize("parent", ["queued", "running", "failed", "interrupted"])
def test_shared_demultiplexing_affects_only_selected_projects(state, parent):
    state.update(status=parent, current_stage="demultiplexing")
    state["restart_request"] = {"from": "demultiplexing", "projects": [B]}
    state["project_completion"] = {A: completion()}
    rows = by_project(state)
    assert rows[A]["status"] == "completed"
    assert (rows[B]["status"], rows[B]["stage"], rows[B]["analysis"]) == (
        parent,
        "demultiplexing",
        "pending",
    )


def test_old_failure_without_project_evidence_remains_unavailable(state):
    state["status"] = "failed"
    assert by_project(state)[A]["status"] == UNAVAILABLE
    assert by_project(state)[A]["analysis"] == UNAVAILABLE
    state["stages"]["analysis"]["metadata"]["input_preflight"] = {"status": "failed"}
    row = by_project(state)[A]
    assert (row["status"], row["analysis"]) == ("failed", "pending")


def test_old_completed_stage_evidence_survives_missing_output_and_archive(state):
    state.update(status="archived", current_stage="finalization", completed_at=COMPLETED)
    state["archive"] = {"status": "archived", "archived_at": COMPLETED}
    for detail in state["stages"].values():
        detail["status"] = "completed"
    row = by_project(state)[A]
    assert (row["status"], row["analysis"], row["delivery"], row["last_delivery"]) == (
        "archived",
        "completed",
        "completed",
        COMPLETED,
    )
    assert row["priority"] == row["analysis_generation"] == UNAVAILABLE


def test_queued_unidentified_run_and_legacy_project_have_honest_detail(state):
    state.update(status="queued", projects=[], current_stage="demultiplexing")
    (row,) = project_rows(state)
    assert row["project"] == UNIDENTIFIED
    assert row["status"] == "queued"
    legacy = legacy_project_row(A, RUN, RUN, timestamp=COMPLETED)
    assert list(legacy) == COLUMNS
    assert (legacy["status"], legacy["stage"], legacy["last_delivery"]) == (
        "completed",
        "legacy",
        COMPLETED,
    )
    assert legacy["analysis"] == legacy["delivery"] == UNAVAILABLE


def test_projection_is_pure_and_does_not_inspect_paths(state, monkeypatch):
    state["project_progress"] = {B: {"analysis": analysis("failed")}}
    before = copy.deepcopy(state)

    def forbid_path_access(*args, **kwargs):
        raise AssertionError("Listing must only project recorded evidence")

    monkeypatch.setattr("pathlib.Path.read_text", forbid_path_access)
    monkeypatch.setattr("pathlib.Path.exists", forbid_path_access)
    rows = project_rows(state)
    assert state == before
    assert [row["project"] for row in rows] == [A, B, C]
    assert all(list(row) == COLUMNS for row in rows)
