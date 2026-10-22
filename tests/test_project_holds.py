"""Holds use real state/lifecycle with synthetic scientific work and guarded mail."""

import copy
import hashlib
import subprocess
import sys

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Event

import pytest

from bcl2fastq_pipeline.state import (
    FlowcellStateStore,
    StateConflictError,
    StateValidationError,
    validate_state,
)
from flowcell_manager import flowcell_manager as manager
from flowcell_manager.project_listing import project_rows
from test_project_lifecycle import restart_context, run
from test_project_scheduling import PROJECTS, RUN, configure_execution, set_priority

from bcl2fastq_pipeline import (
    afterFastq,
    cli,
    findFlowCells,
    initial_holds,
    notifications,
    preflight,
    processing_times,
    project_holds,
    workflow_config,
)

A, B = PROJECTS[:2]


@pytest.fixture
def execution(tmp_path, monkeypatch):
    return configure_execution(tmp_path, monkeypatch, ["2", "0"])


def control(command, projects, **options):
    return manager.control_projects(
        command=command,
        flowcell=RUN,
        projects=projects,
        reason="Synthetic operator decision",
        **options,
    )


def fingerprint(paths):
    return {
        path: (
            hashlib.sha256(path.read_bytes()).hexdigest(),
            path.stat().st_mtime_ns,
            path.stat().st_ino,
        )
        for path in paths
    }


def test_held_high_priority_project_waits_while_sibling_delivers_and_release_preserves_it(
    execution, monkeypatch
):
    cfg, store, events = execution
    before = store.read(RUN)
    control("hold", [B], dry_run=True)
    assert store.read(RUN) == before
    control("hold", [B])
    history = copy.deepcopy(store.read(RUN)["project_control_history"])
    control("hold", [B])
    assert store.read(RUN)["project_control_history"] == history
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting", waiting.get("last_error")
    assert waiting["completed_at"] is None
    assert events == [("analysis", A), ("archiving", A)]
    assert waiting["project_controls"][B]["applied"]
    rows = {row["project"]: row for row in project_rows(waiting)}
    assert rows[A]["status"] == "completed"
    assert rows[B]["status"] == "held"
    assert rows[B]["stage"] == "analysis"
    assert rows[B]["analysis"] == "pending"
    assert rows[B]["priority"] == "0"
    assert manager.combined_list(status="held")["project"].tolist() == [B]
    assert manager.combined_list(query=A).iloc[0]["status"] == "completed"
    output = Path(waiting["output_path"])
    paths = [
        p for p in output.rglob("*") if p.is_file() and (A in str(p) or p.name.startswith("md5sum"))
    ]
    preserved = fingerprint(paths)
    progress = copy.deepcopy(waiting["project_progress"][A])
    completion = copy.deepcopy(waiting["project_completion"][A])
    snapshot = copy.deepcopy(waiting["analysis_snapshots"][A])
    delivered = copy.deepcopy(waiting["delivery_notifications"])
    assert [(n["kind"], n["project"], n["status"]) for n in delivered] == [
        ("processed", A, "sent"),
        ("finalized", A, "sent"),
    ]
    restart_context(cfg, waiting)
    assert findFlowCells.flowCellProcessed()
    # Live curated inputs and the installed libprep configuration are not the
    # authority for a release of the same immutable attempt.
    workflow_config.AUTHORITATIVE_CONFIG.write_text("broken: [\n")
    (output / f"SampleSheet_{B}.csv").write_text("broken metadata\n")
    control("release", [B])
    released = store.read(RUN)
    assert released["status"] == "queued"
    assert Path(released["source_path"]) in cli.candidate_flowcells(cfg, store)
    events.clear()
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert events == [("analysis", B), ("archiving", B)]
    assert completed["attempt"] == waiting["attempt"]
    assert completed["projects"] == [A, B]
    assert completed["project_progress"][A] == progress
    assert completed["project_completion"][A] == completion
    assert completed["analysis_snapshots"][A] == snapshot
    assert completed["delivery_notifications"][:2] == delivered
    assert fingerprint(paths) == preserved


def test_all_held_waits_without_external_work_and_release_is_idempotent(execution):
    cfg, store, events = execution
    control("hold", [A, B, A])
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting", waiting.get("last_error")
    assert events == []
    assert "held" in waiting["waiting_reason"]
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    assert store.read(RUN)["status"] == "waiting"
    restart_context(cfg, waiting)
    control("release", [A, B])
    history = copy.deepcopy(store.read(RUN)["project_control_history"])
    control("release", [A, B])
    assert store.read(RUN)["project_control_history"] == history
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert events == [("analysis", B), ("analysis", A), ("archiving", B), ("archiving", A)]


def test_batch_timings_accumulate_without_held_or_smtp_waits(execution, monkeypatch):
    cfg, store, _events = execution
    clock = [0]
    monkeypatch.setattr(processing_times.time, "monotonic", lambda: clock[0])
    external = afterFastq.run_logged_command

    def work(command, cwd, log_path):
        clock[0] += 3
        return external(command, cwd=cwd, log_path=log_path)

    def smtp(_cfg, _entry):
        clock[0] += 1000

    monkeypatch.setattr(afterFastq, "run_logged_command", work)
    monkeypatch.setattr(notifications, "send_notification", smtp)
    control("hold", [B])
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting"
    assert waiting["processing_timings"]["analysis"][-1]["duration_seconds"] == 3
    clock[0] += 100000  # Time held is never active processing time.
    restart_context(cfg, waiting)
    control("release", [B])
    completed = run(cfg, store)
    assert completed["status"] == "completed"
    assert completed["processing_timings"]["analysis"][-1]["duration_seconds"] == 6
    assert completed["processing_timings"]["archiving"][-1]["duration_seconds"] == 12


def test_release_during_archiving_joins_subsequent_analysis_batch(execution, monkeypatch):
    cfg, store, events = execution
    control("hold", [B])
    external = afterFastq.run_logged_command
    released = False

    def release_in_archive(command, cwd, log_path):
        nonlocal released
        if command[0] == "7za" and not released:
            released = True
            control("release", [B])
            assert store.execution_active(RUN)
            assert store.read(RUN)["project_batch"]["projects"] == [A]
            row = next(row for row in project_rows(store.read(RUN)) if row["project"] == B)
            assert (row["status"], row["stage"], row["analysis"], row["eligibility"]) == (
                "queued",
                "analysis",
                "pending",
                "released",
            )
            assert manager.combined_list(query=B, stage="analysis")["project"].tolist() == [B]
        return external(command, cwd=cwd, log_path=log_path)

    monkeypatch.setattr(afterFastq, "run_logged_command", release_in_archive)
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert events == [("analysis", A), ("archiving", A), ("analysis", B), ("archiving", B)]


@pytest.mark.parametrize("kind", ["processed", "finalized"])
def test_notification_retry_remains_independent_while_sibling_is_held(execution, monkeypatch, kind):
    cfg, store, events = execution
    control("hold", [B])

    def fail_mail(_cfg, entry):
        if entry["kind"] == kind:
            raise OSError("Synthetic relay outage")

    monkeypatch.setattr(notifications, "send_notification", fail_mail)
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting"
    assert (
        next(n for n in waiting["delivery_notifications"] if n["kind"] == kind)["status"]
        == "failed"
    )
    observed = list(events)
    monkeypatch.setattr(notifications, "send_notification", lambda _cfg, _entry: None)
    manager.retry_notifications(flowcell=RUN, kind=kind)
    retried = store.read(RUN)
    assert retried["status"] == "waiting"
    assert retried["project_controls"][B]["held"]
    assert (
        next(n for n in retried["delivery_notifications"] if n["kind"] == kind)["status"] == "sent"
    )
    assert events == observed


@pytest.mark.parametrize("stage", ["analysis", "archiving"])
def test_hold_during_claimed_work_finishes_transaction_and_blocks_next_stage(
    execution, monkeypatch, stage
):
    cfg, store, events = execution
    external = afterFastq.run_logged_command
    requested = False

    def hold_running(command, cwd, log_path):
        nonlocal requested
        if not requested and ((command[0] == "7za") == (stage == "archiving")):
            requested = True
            control("hold", [B])  # priority zero: B starts first
            state = store.read(RUN)
            assert state["project_controls"][B]["held"]
            assert not state["project_controls"][B]["applied"]
            row = next(row for row in project_rows(state) if row["project"] == B)
            assert row["status"] == "running"
            assert row["eligibility"] == "hold_requested"
        return external(command, cwd=cwd, log_path=log_path)

    monkeypatch.setattr(afterFastq, "run_logged_command", hold_running)
    state = run(cfg, store)
    assert state["project_controls"][B]["applied"]
    if stage == "archiving":
        assert state["status"] == "completed", state.get("last_error")
        assert state["project_progress"][B]["finalization"]["status"] == "completed"
    else:
        assert state["status"] == "waiting", state.get("last_error")
        assert events == [("analysis", B), ("analysis", A), ("archiving", A)]
        restart_context(cfg, state)
        control("release", [B])
        events.clear()
        completed = run(cfg, store)
        assert completed["status"] == "completed", completed.get("last_error")
        assert events == [("archiving", B)]


def test_failed_analysis_stays_failed_under_controls_and_new_attempt_retains_hold(
    execution, monkeypatch
):
    cfg, store, _events = execution
    external = afterFastq.run_logged_command

    def fail(command, cwd, log_path):
        if command[0] != "7za":
            raise RuntimeError("Synthetic scientific failure")
        return external(command, cwd=cwd, log_path=log_path)

    monkeypatch.setattr(afterFastq, "run_logged_command", fail)
    failed = run(cfg, store)
    assert failed["status"] == "failed"
    control("hold", [B])
    assert store.read(RUN)["last_error"] == failed["last_error"]
    assert next(r for r in project_rows(store.read(RUN)) if r["project"] == B)["status"] == "failed"
    control("release", [B])
    assert store.read(RUN)["status"] == "failed"
    control("hold", [B])
    restart_context(cfg, failed)
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", force=True)
    assert store.read(RUN)["project_controls"][B]["held"]
    monkeypatch.setattr(afterFastq, "run_logged_command", external)
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting", waiting.get("last_error")
    assert waiting["attempt"] == failed["attempt"] + 1


def test_claim_and_hold_are_serialized_by_state_lock(execution):
    cfg, store, _events = execution
    state = store.begin_attempt(RUN, cfg=cfg)
    initial_holds.initialize(
        store,
        RUN,
        preflight.validate_selection(preflight.select_run_inputs(cfg.output_path, cfg.output_path)),
    )
    control("hold", [B])
    reached = Event()

    def claim():
        reached.set()
        return project_holds.claim(store, RUN, B, "analysis")

    with ThreadPoolExecutor(max_workers=1) as pool:
        with store.lock(RUN):
            pending = pool.submit(claim)
            assert reached.wait(5)
            assert not pending.done()
        assert not pending.result(timeout=5)
    control("release", [B])
    assert project_holds.claim(store, RUN, B, "analysis")
    control("hold", [B])
    assert not store.read(RUN)["project_controls"][B]["applied"]
    assert state["attempt"] == store.read(RUN)["attempt"]


def test_scoped_correction_waits_and_releases_without_replacing_staging(execution):
    cfg, store, events = execution
    original = run(cfg, store)
    assert original["status"] == "completed"
    first = copy.deepcopy(original["project_completion"][A])
    restart_context(cfg, original)
    control("hold", [B])
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[B], force=True)
    events.clear()
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting", waiting.get("last_error")
    assert events == []
    assert waiting["project_completion"][A] == first
    restart_context(cfg, waiting)
    control("release", [B])
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert completed["attempt"] == waiting["attempt"]
    assert completed["project_completion"][A] == first
    assert events == [("analysis", B), ("archiving", B)]


def test_scoped_correction_of_delivered_sibling_retains_held_pending_work(execution):
    cfg, store, events = execution
    control("hold", [B])
    first = run(cfg, store)
    assert first["status"] == "waiting"
    original_inputs = copy.deepcopy(first["attempts"][-1]["execution_inputs"])
    # An unrelated scoped attempt may validate newer curated priorities, but
    # cannot change a suspended project's original producing plan.
    set_priority(Path(first["output_path"]) / f"SampleSheet_{B}.csv", "3")
    set_priority(Path(first["output_path"]) / f"SampleSheet_{A}.csv", "2")
    restart_context(cfg, first)
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[A], force=True)
    events.clear()
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting", waiting.get("last_error")
    assert events == [("analysis", A), ("archiving", A)]
    assert waiting["continuation_context"]["execution_inputs"] == original_inputs
    assert waiting["attempts"][-1]["project_schedule"]["order"] == [B, A]
    assert waiting["attempts"][-1]["project_schedule"]["projects"][B]["priority"] == "0"
    completed_a = copy.deepcopy(waiting["project_completion"][A])
    restart_context(cfg, waiting)
    control("release", [B])
    events.clear()
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert completed["project_completion"][A] == completed_a
    assert events == [("analysis", B), ("archiving", B)]


def test_held_pair_still_participates_in_full_preflight(execution):
    cfg, store, events = execution
    control("hold", [B])
    sheet = Path(store.read(RUN)["output_path"]) / f"SampleSheet_{B}.csv"
    sheet.write_text(sheet.read_text().replace("Priority,0", "Priority,NaN"))
    failed = run(cfg, store)
    assert failed["status"] == "failed"
    assert failed["project_controls"][B]["held"]
    assert events == []
    assert "scheduling.priority" in str(failed["stages"]["analysis"]["metadata"])


def test_eligible_decimal_order_and_equivalent_ties_with_held_zero(tmp_path, monkeypatch):
    cfg, store, events = configure_execution(tmp_path, monkeypatch, ["1.00", None, "0.5", "0", "2"])
    control("hold", [PROJECTS[3]])
    waiting = run(cfg, store)
    expected = [PROJECTS[i] for i in (2, 0, 1, 4)]
    assert waiting["status"] == "waiting", waiting.get("last_error")
    assert events == [(s, p) for s in ("analysis", "archiving") for p in expected]


def test_archive_failure_then_held_recovery_retains_first_delivery(execution, monkeypatch):
    cfg, store, events = execution
    external = afterFastq.run_logged_command

    def fail_archive(command, cwd, log_path):
        if command[0] == "7za" and Path(command[2]).name.startswith(A):
            raise RuntimeError("Synthetic archive failure")
        return external(command, cwd=cwd, log_path=log_path)

    monkeypatch.setattr(afterFastq, "run_logged_command", fail_archive)
    failed = run(cfg, store)
    assert failed["status"] == "failed"
    first = copy.deepcopy(failed["project_completion"][B])
    control("hold", [A])
    control("release", [A])
    assert store.read(RUN)["status"] == "failed"
    control("hold", [A])
    restart_context(cfg, failed)
    manager.rerun_flowcell(flowcell=RUN, from_stage="finalization", force=True)
    monkeypatch.setattr(afterFastq, "run_logged_command", external)
    events.clear()
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting", waiting.get("last_error")
    assert events == []
    assert waiting["project_completion"][B] == first
    restart_context(cfg, waiting)
    control("release", [A])
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert completed["project_completion"][B] == first
    assert events == [("archiving", A)]


def test_downstream_retry_after_last_commit_failure_finishes_without_rework(execution, monkeypatch):
    cfg, store, events = execution
    control("release", [A, B])
    write = store._atomic_write_unlocked
    raised = False

    def fail_after_commit(state):
        nonlocal raised
        write(state)
        if (
            not raised
            and state["status"] == "running"
            and all(
                state.get("project_progress", {}).get(p, {}).get("finalization", {}).get("status")
                == "completed"
                for p in (A, B)
            )
        ):
            raised = True
            raise OSError("Synthetic durability failure after the last project commit")

    monkeypatch.setattr(store, "_atomic_write_unlocked", fail_after_commit)
    first = run(cfg, store)
    assert first["status"] == "failed"
    monkeypatch.setattr(store, "_atomic_write_unlocked", write)
    completion = copy.deepcopy(first["project_completion"])
    control("hold", [A, B])
    control("release", [A, B])
    assert store.read(RUN)["status"] == "failed"
    restart_context(cfg, first)
    manager.rerun_flowcell(flowcell=RUN, from_stage="finalization", force=True)
    events.clear()
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert completed["project_completion"] == completion
    assert events == []


def test_scoped_success_does_not_queue_a_held_siblings_failed_archive(execution, monkeypatch):
    cfg, store, events = execution
    external = afterFastq.run_logged_command

    def fail_archive(command, cwd, log_path):
        if command[0] == "7za" and Path(command[2]).name.startswith(A):
            raise RuntimeError("Synthetic sibling archive failure")
        return external(command, cwd=cwd, log_path=log_path)

    monkeypatch.setattr(afterFastq, "run_logged_command", fail_archive)
    first = run(cfg, store)
    assert first["status"] == "failed"
    control("hold", [A])
    restart_context(cfg, first)
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[B], force=True)
    monkeypatch.setattr(afterFastq, "run_logged_command", external)
    events.clear()
    state = run(cfg, store)
    assert state["status"] == "failed"
    assert state["project_progress"][A]["finalization"]["status"] == "failed"
    assert state["project_progress"][B]["finalization"]["status"] == "completed"
    assert "explicit recovery" in state["last_error"]["summary"]
    assert events == [("analysis", B), ("archiving", B)]


def test_controls_reject_unknown_projects_and_malformed_state_without_mutation(execution):
    _cfg, store, _events = execution
    before = store.read(RUN)
    with pytest.raises(StateConflictError, match="known projects"):
        control("hold", ["GCF-2099-999"])
    assert store.read(RUN) == before
    control("hold", [B])
    malformed = store.read(RUN)
    malformed["project_controls"][B]["held"] = "true"
    with pytest.raises(StateValidationError, match="project hold"):
        validate_state(malformed)


@pytest.mark.parametrize("name", ["fm", "flowcell-manager"])
def test_manager_command_parser_controls_and_preview_agree(execution, monkeypatch, capsys, name):
    _cfg, store, _events = execution
    monkeypatch.setattr(
        sys, "argv", [name, "hold", RUN, "--project", B, "--reason", "Review", "--dry-run"]
    )
    before = store.read(RUN)
    manager.main()
    assert store.read(RUN) == before
    assert "Preview only" in capsys.readouterr().out
    monkeypatch.setattr(sys, "argv", [name, "hold", RUN, "--project", B, "--reason", "Review"])
    manager.main()
    assert store.read(RUN)["project_controls"][B]["held"]


@pytest.mark.parametrize("name", ["fm", "flowcell-manager"])
def test_installed_controls_persist_and_project_filters_agree(execution, tmp_path, name):
    cfg, store, _events = execution
    config = tmp_path / "manager-test.ini"
    config.write_text(
        "[Paths]\n"
        + "\n".join(f"{name}={path}" for name, path in vars(cfg.static.paths).items())
        + "\n"
    )
    driver = """
import runpy, smtplib, sys
from bcl2fastq_pipeline.config import PipelineConfig
def forbidden_smtp(*args, **kwargs):
    raise AssertionError("Scheduling must never send mail")
smtplib.SMTP = forbidden_smtp
smtplib.SMTP_SSL = forbidden_smtp
PipelineConfig.load(sys.argv[1])
sys.argv = sys.argv[2:]
runpy.run_path(sys.argv[0], run_name="__main__")
"""

    def invoke(*args):
        return subprocess.run(
            [
                sys.executable,
                "-c",
                driver,
                str(config),
                str(Path(sys.executable).parent / name),
                *args,
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=False,
        )

    before = store.read(RUN)
    preview = invoke("hold", RUN, "--project", B, "--reason", "Review", "--dry-run")
    assert preview.returncode == 0, preview.stderr
    assert store.read(RUN) == before
    changed = invoke("hold", RUN, "--project", B, "--reason", "Review")
    assert changed.returncode == 0, changed.stderr
    assert store.read(RUN)["project_controls"][B]["held"]
    listed = invoke("list", "--status", "held")
    searched = invoke("search", B, "--status", "held")
    assert listed.returncode == searched.returncode == 0
    assert listed.stdout == searched.stdout
    assert B in listed.stdout and A not in listed.stdout
    assert invoke("release", RUN, "--project", B, "--reason", "Ready").returncode == 0
    assert not store.read(RUN)["project_controls"][B]["held"]
