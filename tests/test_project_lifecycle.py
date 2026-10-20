"""Project boundaries use real state, reports, publication and snapshot recovery."""

import copy
import hashlib
import logging

from pathlib import Path
from unittest.mock import Mock

import pytest

from bcl2fastq_pipeline.state import DeliveryUncertainError
from flowcell_manager import flowcell_manager as manager
from test_project_scheduling import PROJECTS, RUN, configure_execution

from bcl2fastq_pipeline import (
    afterFastq,
    cli,
    notification_delivery,
    notifications,
    project_corrections,
)

A, B = PROJECTS[:2]


@pytest.fixture
def execution(tmp_path, monkeypatch, request):
    return configure_execution(tmp_path, monkeypatch, request.param)


def run(cfg, store):
    cli._run_state_backed_flowcell(
        cfg, store, logging.getLogger("project-boundary-test"), prepare=True
    )
    return store.read(RUN)


def restart_context(cfg, state):
    cfg.run.begin(Path(state["source_path"]), cfg.static.paths)


def mail_double(cfg, sent):
    cfg.static.email.update(
        host="relay.example.org",
        from_address="bfq@example.org",
        finished_to="analysis@example.org",
        error_to="delivery@example.org",
    )

    def send(current_cfg, entry):
        message = notifications.build_message(current_cfg, entry)
        if entry["kind"] == "processed":
            attachments = list(message.iter_attachments())
            assert any(entry["project"] in item.get_filename() for item in attachments)
            assert len(entry["payload"]["project_contexts"]) == 1
            assert entry["project"] in Path(entry["payload"]["sample_sheet"]).read_text()
        sent.append(copy.deepcopy(entry))

    return send


@pytest.mark.parametrize("execution", [["0", "1"]], indirect=True)
def test_analysis_and_delivery_notifications_precede_next_project_work(execution, monkeypatch):
    cfg, store, _events = execution
    sent = []
    monkeypatch.setattr(notifications, "send_notification", mail_double(cfg, sent))
    external = afterFastq.run_logged_command

    def ordered_work(command, cwd, log_path):
        if command[0] != "7za" and cwd.name.startswith(B):
            assert [(entry["kind"], entry["project"]) for entry in sent] == [("processed", A)]
            assert store.read(RUN)["project_progress"][A]["analysis"]["status"] == "completed"
        if command[0] == "7za":
            assert all(
                store.read(RUN)["project_progress"][p]["analysis"]["status"] == "completed"
                for p in (A, B)
            )
            if Path(command[2]).name.startswith(B):
                state = store.read(RUN)
                assert state["status"] == "running"
                assert (
                    state["project_completion"][A]["stages"]["finalization"]["status"]
                    == "completed"
                )
                snapshot = state["analysis_snapshots"][A]
                assert "pending" not in snapshot
                assert (Path(state["output_path"]) / snapshot["archive"]).is_file()
                assert any(entry["kind"] == "finalized" and entry["project"] == A for entry in sent)
        return external(command, cwd=cwd, log_path=log_path)

    monkeypatch.setattr(afterFastq, "run_logged_command", ordered_work)
    state = run(cfg, store)
    assert state["status"] == "completed", state.get("last_error")
    assert [(entry["kind"], entry["project"]) for entry in sent] == [
        ("processed", A),
        ("processed", B),
        ("finalized", A),
        ("finalized", B),
    ]
    assert len({entry["id"] for entry in state["delivery_notifications"]}) == 4
    assert all(entry["status"] == "sent" for entry in state["delivery_notifications"])
    for step in ("archiving", "archive_checksums"):
        assert state["processing_timings"][step][-1]["duration_seconds"] == sum(
            state["project_progress"][project]["timings"][step] for project in (A, B)
        )


@pytest.mark.parametrize("execution", [["0", "1"]], indirect=True)
@pytest.mark.parametrize("kind", ["processed", "finalized"])
@pytest.mark.parametrize("uncertain", [False, True])
def test_mail_failure_never_stops_projects_and_requires_explicit_retry(
    execution, monkeypatch, kind, uncertain
):
    cfg, store, events = execution
    calls = []

    def send(_cfg, entry):
        calls.append(entry["id"])
        if entry["kind"] == kind and entry["project"] == A:
            raise (
                DeliveryUncertainError("synthetic interrupted SMTP")
                if uncertain
                else OSError("synthetic SMTP failure")
            )

    monkeypatch.setattr(notifications, "send_notification", send)
    state = run(cfg, store)
    assert state["status"] == "completed", state.get("last_error")
    failed = next(
        entry
        for entry in state["delivery_notifications"]
        if entry["kind"] == kind and entry["project"] == A
    )
    assert failed["status"] == ("uncertain" if uncertain else "failed")
    notification_delivery.recover_pending(cfg, store)
    assert calls.count(failed["id"]) == 1
    observed = copy.deepcopy(events)
    monkeypatch.setattr(notifications, "send_notification", Mock())
    result = manager.retry_notifications(flowcell=RUN, kind=kind, retry_uncertain=uncertain)
    assert result is not None
    assert events == observed
    assert (
        next(
            entry
            for entry in store.read(RUN)["delivery_notifications"]
            if entry["id"] == failed["id"]
        )["status"]
        == "sent"
    )


@pytest.mark.parametrize("execution", [["0", "1"]], indirect=True)
def test_later_analysis_failure_and_scoped_retry_preserve_first_analysis_and_mail(
    execution, monkeypatch
):
    cfg, store, events = execution
    sent = []
    monkeypatch.setattr(
        notifications, "send_notification", lambda _cfg, entry: sent.append(copy.deepcopy(entry))
    )
    external = afterFastq.run_logged_command

    def fail_second(command, cwd, log_path):
        if command[0] != "7za" and cwd.name.startswith(B):
            raise RuntimeError("synthetic second analysis failure")
        return external(command, cwd=cwd, log_path=log_path)

    monkeypatch.setattr(afterFastq, "run_logged_command", fail_second)
    state = run(cfg, store)
    assert state["status"] == "failed"
    first = copy.deepcopy(state["project_progress"][A]["analysis"])
    first_mail = copy.deepcopy(state["delivery_notifications"][0])
    assert first_mail["status"] == "sent"
    retained = {
        p: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
        for p in Path(first["reports"]).iterdir()
    }
    restart_context(cfg, state)
    monkeypatch.setattr(afterFastq, "run_logged_command", external)
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[B], force=True)
    events.clear()
    retried = run(cfg, store)
    assert retried["status"] == "queued", retried.get("last_error")
    assert retried["current_stage"] == "finalization"
    assert retried["project_progress"][A]["analysis"] == first
    assert retried["delivery_notifications"][0] == first_mail
    assert events == [("analysis", B), ("archiving", B)]
    restart_context(cfg, retried)
    events.clear()
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert events == [("archiving", A)]
    assert (
        len([entry for entry in sent if entry["kind"] == "processed" and entry["project"] == A])
        == 1
    )
    assert {
        p: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns) for p in retained
    } == retained


@pytest.mark.parametrize("execution", [["0", "1"]], indirect=True)
def test_partial_delivery_retry_skips_committed_archive_snapshot_and_notifications(
    execution, monkeypatch
):
    cfg, store, events = execution
    sent = []
    monkeypatch.setattr(
        notifications, "send_notification", lambda _cfg, entry: sent.append(copy.deepcopy(entry))
    )
    external = afterFastq.run_logged_command

    def fail_archive(command, cwd, log_path):
        if command[0] == "7za" and Path(command[2]).name.startswith(B):
            raise RuntimeError("synthetic second archive failure")
        return external(command, cwd=cwd, log_path=log_path)

    monkeypatch.setattr(afterFastq, "run_logged_command", fail_archive)
    state = run(cfg, store)
    assert state["status"] == "failed"
    output = Path(state["output_path"])
    paths = [
        output / f"{A}_990101.7za",
        output / f"QC_{A}_990101.7za",
        output / state["analysis_snapshots"][A]["archive"],
    ]
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ino) for p in paths}
    completion = copy.deepcopy(state["project_completion"][A])
    first_mail = next(
        entry
        for entry in state["delivery_notifications"]
        if entry["kind"] == "finalized" and entry["project"] == A
    )
    restart_context(cfg, state)
    manager.rerun_flowcell(flowcell=RUN, from_stage="finalization", force=True)
    monkeypatch.setattr(afterFastq, "run_logged_command", external)
    events.clear()
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert events == [("archiving", B)]
    assert completed["project_completion"][A] == completion
    assert {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ino) for p in paths} == before
    assert (
        next(
            entry
            for entry in completed["delivery_notifications"]
            if entry["id"] == first_mail["id"]
        )
        == first_mail
    )


@pytest.mark.parametrize("execution", [["0", "1"]], indirect=True)
def test_interrupted_smtp_claim_is_not_resent_automatically(execution, monkeypatch):
    cfg, store, _events = execution
    sent = []
    monkeypatch.setattr(
        notifications, "send_notification", lambda _cfg, entry: sent.append(entry["id"])
    )
    finish = store._finish_delivery_unlocked
    raised = False

    def interrupt(current, entry, status, error):
        nonlocal raised
        if entry.get("project") == A and entry["kind"] == "processed" and not raised:
            raised = True
            raise OSError("synthetic persistence failure after SMTP acceptance")
        return finish(current, entry, status, error)

    monkeypatch.setattr(store, "_finish_delivery_unlocked", interrupt)
    state = run(cfg, store)
    assert state["status"] == "completed", state.get("last_error")
    first = next(
        entry
        for entry in state["delivery_notifications"]
        if entry.get("project") == A and entry["kind"] == "processed"
    )
    assert first["status"] == "uncertain"
    notification_delivery.recover_pending(cfg, store)
    assert sent.count(first["id"]) == 1


@pytest.mark.parametrize("execution", [["0", "1"]], indirect=True)
def test_staged_publication_failure_preserves_old_delivery_and_analysis_generation(
    execution, monkeypatch
):
    cfg, store, events = execution
    initial = run(cfg, store)
    assert initial["status"] == "completed"
    restart_context(cfg, initial)
    output = Path(initial["output_path"])
    archive = output / f"{B}_990101.7za"
    before = archive.read_bytes(), archive.stat().st_ino
    old_snapshot = (output / initial["analysis_snapshots"][B]["archive"]).read_bytes()
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[B], force=True)
    commit = store.complete_project_finalization
    monkeypatch.setattr(
        store,
        "complete_project_finalization",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("synthetic precommit failure")),
    )
    failed = run(cfg, store)
    assert failed["status"] == "failed"
    assert (archive.read_bytes(), archive.stat().st_ino) == before
    assert (output / failed["analysis_snapshots"][B]["archive"]).read_bytes() == old_snapshot
    analysis = copy.deepcopy(failed["project_progress"][B]["analysis"])
    assert analysis["status"] == "completed"
    assert Path(analysis["reports"]).is_dir()
    processed = next(
        entry
        for entry in failed["delivery_notifications"]
        if entry.get("generation") == analysis["generation"] and entry["kind"] == "processed"
    )
    assert processed["status"] == "sent"
    assert project_corrections.current(failed)["project_publications"][B]["status"] == "rolled_back"
    monkeypatch.setattr(store, "complete_project_finalization", commit)
    restart_context(cfg, failed)
    manager.rerun_flowcell(flowcell=RUN, from_stage="finalization", projects=[B], force=True)
    events.clear()
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert events == [("archiving", B)]
    assert completed["project_progress"][B]["analysis"] == analysis
    assert (
        next(
            entry for entry in completed["delivery_notifications"] if entry["id"] == processed["id"]
        )
        == processed
    )


@pytest.mark.parametrize("execution", [["0", "1"]], indirect=True)
def test_partial_correction_publication_keeps_first_commit_on_second_failure(
    execution, monkeypatch
):
    cfg, store, events = execution
    original = run(cfg, store)
    restart_context(cfg, original)
    output = Path(original["output_path"])
    retained_b = (output / f"{B}_990101.7za").stat().st_ino
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[A, B], force=True)
    commit = store.complete_project_finalization

    def fail_second(run_id, project, payload, snapshot):
        if project == B:
            raise OSError("synthetic second project precommit failure")
        return commit(run_id, project, payload, snapshot)

    monkeypatch.setattr(store, "complete_project_finalization", fail_second)
    failed = run(cfg, store)
    assert failed["status"] == "failed"
    journals = project_corrections.current(failed)["project_publications"]
    assert journals[A]["status"] == "committed"
    assert journals[B]["status"] == "rolled_back"
    assert (output / f"{B}_990101.7za").stat().st_ino == retained_b
    first = copy.deepcopy(failed["project_completion"][A])
    snapshot_a = (output / failed["analysis_snapshots"][A]["archive"]).read_bytes()
    project_corrections.recover_pending(store)
    assert store.read(RUN)["project_completion"][A] == first
    restart_context(cfg, failed)
    monkeypatch.setattr(store, "complete_project_finalization", commit)
    manager.rerun_flowcell(flowcell=RUN, from_stage="finalization", force=True)
    events.clear()
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert events == [("archiving", B)]
    assert completed["project_completion"][A] == first
    assert (output / completed["analysis_snapshots"][A]["archive"]).read_bytes() == snapshot_a


@pytest.mark.parametrize("execution", [["0", "1"]], indirect=True)
def test_interruption_keeps_committed_analysis_and_notification_across_recovery(
    execution, monkeypatch
):
    cfg, store, _events = execution
    external = afterFastq.run_logged_command

    def interrupt(command, cwd, log_path):
        if command[0] != "7za" and cwd.name.startswith(B):
            raise KeyboardInterrupt("synthetic worker interruption")
        return external(command, cwd=cwd, log_path=log_path)

    monkeypatch.setattr(afterFastq, "run_logged_command", interrupt)
    with pytest.raises(KeyboardInterrupt):
        run(cfg, store)
    before = store.read(RUN)
    recovered = store.recover_interrupted(RUN)
    assert recovered["status"] == "interrupted"
    assert recovered["project_progress"][A]["analysis"] == before["project_progress"][A]["analysis"]
    assert recovered["project_progress"][B]["analysis"]["status"] == "failed"
    assert recovered["delivery_notifications"][0] == before["delivery_notifications"][0]
    notification_delivery.recover_pending(cfg, store)
    assert store.read(RUN)["delivery_notifications"][0] == before["delivery_notifications"][0]


@pytest.mark.parametrize("execution", [["0", "1"]], indirect=True)
def test_state_error_after_project_commit_preserves_publication_and_skips_rework(
    execution, monkeypatch
):
    cfg, store, events = execution
    initial = run(cfg, store)
    restart_context(cfg, initial)
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[B], force=True)
    write = store._atomic_write_unlocked
    raised = False

    def fail_after_commit(state):
        nonlocal raised
        write(state)
        progress = state.get("project_progress", {}).get(B, {}).get("finalization", {})
        if (
            not raised
            and state["status"] == "running"
            and progress.get("status") == "completed"
            and progress.get("attempt") == state["attempt"]
        ):
            raised = True
            raise OSError("synthetic directory fsync failure after project commit")

    monkeypatch.setattr(store, "_atomic_write_unlocked", fail_after_commit)
    failed = run(cfg, store)
    assert failed["status"] == "failed"
    assert failed["project_progress"][B]["finalization"]["status"] == "completed"
    assert project_corrections.current(failed)["project_publications"][B]["status"] == "committed"
    assert "pending" not in failed["analysis_snapshots"][B]
    notification_delivery.recover_pending(cfg, store)
    completed_mail = next(
        entry
        for entry in store.read(RUN)["delivery_notifications"]
        if entry["kind"] == "finalized"
        and entry.get("project") == B
        and entry["attempt"] == failed["attempt"]
    )
    assert completed_mail["status"] == "sent"
    restart_context(cfg, failed)
    manager.rerun_flowcell(flowcell=RUN, from_stage="finalization", force=True)
    events.clear()
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert events == []
