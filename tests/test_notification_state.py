"""Failure injection for persisted completion mail, without any SMTP access."""

import copy

from concurrent.futures import ThreadPoolExecutor
from threading import Event

import pytest

from bcl2fastq_pipeline.state import (
    DeliveryUncertainError,
    FlowcellStateStore,
    StateConflictError,
    StateValidationError,
    new_state,
    validate_state,
)

RUN_ID = "260918_MN00686_0026_A000HCMFHF"


def reporting_store(tmp_path):
    store = FlowcellStateStore(tmp_path / "manager")
    state = new_state(
        RUN_ID,
        tmp_path / "source",
        tmp_path / "output",
        origin="restored_legacy_fastq",
        start_stage="reporting",
    )
    store.create(state)
    store.begin_attempt(RUN_ID)
    return store


def processed(store):
    return store.complete_stage(
        RUN_ID,
        "reporting",
        notification={"kind": "processed", "payload": {"elapsed": "0:01:00"}},
    )


def completed_store(tmp_path):
    store = reporting_store(tmp_path)
    processed(store)
    store.start_stage(RUN_ID, "finalization")
    store.complete_run(
        RUN_ID,
        ["GCF-2026-001"],
        notification={"kind": "finalized", "payload": {"elapsed": "0:02:00"}},
    )
    return store


def entry(store, identifier="processed:1"):
    return next(
        item for item in store.read(RUN_ID)["delivery_notifications"] if item["id"] == identifier
    )


def processing_state(store):
    state = store.read(RUN_ID)
    state.pop("delivery_notifications", None)
    state.pop("updated_at")
    return state


def test_completion_and_notification_are_committed_together(tmp_path):
    store = reporting_store(tmp_path)
    state = processed(store)
    assert state["stages"]["reporting"]["status"] == "completed"
    assert state["stages"]["finalization"]["status"] == "queued"
    notification = entry(store)
    assert notification["status"] == "pending"
    assert notification["attempt"] == 1
    assert notification["attempts"] == []
    assert notification["payload"] == {"elapsed": "0:01:00"}
    assert notification["stage"] == "reporting"


@pytest.mark.parametrize("finalization", [False, True])
def test_interrupted_completion_write_keeps_previous_state(tmp_path, monkeypatch, finalization):
    store = reporting_store(tmp_path)
    if finalization:
        processed(store)
        store.start_stage(RUN_ID, "finalization")
    before = store.state_path(RUN_ID).read_bytes()

    def fail_replace(*args):
        raise OSError("interrupted atomic rename")

    monkeypatch.setattr("bcl2fastq_pipeline.state.os.replace", fail_replace)
    with pytest.raises(OSError, match="interrupted atomic rename"):
        if finalization:
            store.complete_run(RUN_ID, [], notification={"kind": "finalized", "payload": {}})
        else:
            processed(store)
    assert store.state_path(RUN_ID).read_bytes() == before
    assert not list(store.states_dir.glob("*.tmp"))


def test_send_observes_durable_completion_and_claim_and_cannot_mutate_payload(tmp_path):
    store = completed_store(tmp_path)
    before = processing_state(store)
    observed = []

    def send(notification):
        snapshot = FlowcellStateStore(store.manager_dir).read(RUN_ID)
        observed.append(snapshot)
        assert snapshot["status"] == "completed"
        assert snapshot["delivery_notifications"][1]["status"] == "sending"
        notification["payload"]["elapsed"] = "mutated by sender"

    assert store.deliver_notification(RUN_ID, "finalized:1", send)
    assert len(observed) == 1
    notification = entry(store, "finalized:1")
    assert notification["status"] == "sent"
    assert notification["sent_at"]
    assert notification["payload"]["elapsed"] == "0:02:00"
    assert notification["attempts"][0]["status"] == "sent"
    assert notification["attempts"][0]["completed_at"]
    assert processing_state(store) == before
    assert not store.deliver_notification(RUN_ID, "finalized:1", send, retry=True)
    assert len(observed) == 1


@pytest.mark.parametrize("exception", [ValueError("bad settings"), OSError("connection refused")])
def test_failed_delivery_requires_explicit_retry_and_preserves_processing(tmp_path, exception):
    store = completed_store(tmp_path)
    before = processing_state(store)

    def fail(notification):
        raise exception

    assert not store.deliver_notification(RUN_ID, "processed:1", fail)
    failed = entry(store)
    assert failed["status"] == "failed"
    assert str(exception) in failed["last_error"]
    assert failed["attempts"][0]["status"] == "failed"
    calls = []
    restarted = FlowcellStateStore(store.manager_dir)
    assert not restarted.deliver_notification(RUN_ID, "processed:1", calls.append)
    assert calls == []
    assert restarted.deliver_notification(RUN_ID, "processed:1", calls.append, retry=True)
    assert len(calls) == 1
    notification = entry(store)
    assert notification["status"] == "sent"
    assert notification["last_error"] is None
    assert [attempt["status"] for attempt in notification["attempts"]] == ["failed", "sent"]
    assert processing_state(store) == before


def test_uncertain_delivery_requires_explicit_duplicate_acknowledgement(tmp_path):
    store = completed_store(tmp_path)

    def timeout(notification):
        raise DeliveryUncertainError("Disconnected after DATA; relay may have accepted it")

    assert not store.deliver_notification(RUN_ID, "processed:1", timeout)
    assert entry(store)["status"] == "uncertain"
    calls = []
    assert not store.deliver_notification(RUN_ID, "processed:1", calls.append, retry=True)
    assert not store.deliver_notification(RUN_ID, "processed:1", calls.append, retry_uncertain=True)
    assert calls == []
    assert store.deliver_notification(
        RUN_ID, "processed:1", calls.append, retry=True, retry_uncertain=True
    )
    assert len(calls) == 1
    assert [attempt["status"] for attempt in entry(store)["attempts"]] == ["uncertain", "sent"]


def test_crash_during_delivery_is_not_automatically_retried(tmp_path):
    store = completed_store(tmp_path)
    before = processing_state(store)
    accepted = []

    def crash(notification):
        accepted.append(notification["id"])
        raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        store.deliver_notification(RUN_ID, "processed:1", crash)
    assert entry(store)["status"] == "sending"
    restarted = FlowcellStateStore(store.manager_dir)
    assert not restarted.deliver_notification(RUN_ID, "processed:1", accepted.append)
    notification = entry(store)
    assert notification["status"] == "uncertain"
    assert notification["attempts"][0]["status"] == "uncertain"
    assert "interrupted" in notification["last_error"]
    assert accepted == ["processed:1"]
    assert processing_state(store) == before


def test_claim_write_failure_prevents_sending(tmp_path, monkeypatch):
    store = completed_store(tmp_path)
    before = store.state_path(RUN_ID).read_bytes()
    calls = []

    def fail_replace(*args):
        raise OSError("full disk")

    monkeypatch.setattr("bcl2fastq_pipeline.state.os.replace", fail_replace)
    with pytest.raises(OSError, match="full disk"):
        store.deliver_notification(RUN_ID, "processed:1", calls.append)
    assert calls == []
    assert store.state_path(RUN_ID).read_bytes() == before


def test_acceptance_before_success_write_failure_remains_uncertain(tmp_path, monkeypatch):
    store = completed_store(tmp_path)
    original_write = store._atomic_write_unlocked
    calls = []

    def fail_success_write(state):
        if state["delivery_notifications"][0]["status"] == "sent":
            raise OSError("cannot persist accepted mail")
        original_write(state)

    monkeypatch.setattr(store, "_atomic_write_unlocked", fail_success_write)
    with pytest.raises(OSError, match="cannot persist accepted mail"):
        store.deliver_notification(RUN_ID, "processed:1", calls.append)
    assert len(calls) == 1
    assert entry(store)["status"] == "sending"
    restarted = FlowcellStateStore(store.manager_dir)
    assert not restarted.deliver_notification(RUN_ID, "processed:1", calls.append)
    assert entry(store)["status"] == "uncertain"
    assert len(calls) == 1
    assert restarted.deliver_notification(
        RUN_ID, "processed:1", calls.append, retry=True, retry_uncertain=True
    )
    assert len(calls) == 2


@pytest.mark.parametrize("from_stage", ["demultiplexing", "analysis", "reporting", "finalization"])
def test_preparing_invalidates_only_notifications_dependent_on_restart(tmp_path, from_stage):
    store = completed_store(tmp_path)
    state = store.set_preparing(RUN_ID, from_stage, reason="rerun", refresh_inputs=False)
    statuses = [item["status"] for item in state["delivery_notifications"]]
    if from_stage == "finalization":
        assert statuses == ["pending", "superseded"]
    else:
        assert statuses == ["superseded", "superseded"]
    calls = []
    assert not store.deliver_notification(
        RUN_ID, "finalized:1", calls.append, retry=True, retry_uncertain=True
    )
    assert calls == []


def test_queue_defensively_invalidates_and_new_attempt_gets_new_identity(tmp_path):
    store = completed_store(tmp_path)
    assert store.deliver_notification(RUN_ID, "processed:1", lambda notification: None)
    sent = entry(store)
    store.queue(RUN_ID, "reporting")
    superseded = entry(store)
    assert superseded["status"] == "superseded"
    assert superseded["sent_at"] == sent["sent_at"]
    assert superseded["attempts"] == sent["attempts"]
    assert superseded["invalidated_at"]
    store.begin_attempt(RUN_ID)
    processed(store)
    assert entry(store, "processed:2")["status"] == "pending"
    assert entry(store)["status"] == "superseded"


def test_archive_invalidation_precedes_cleanup_and_preserves_sent_history(tmp_path):
    store = completed_store(tmp_path)
    assert store.deliver_notification(RUN_ID, "processed:1", lambda notification: None)
    state = store.invalidate_notifications(RUN_ID, reason="Archive cleanup", unsent_only=True)
    assert state["archive"]["status"] == "active"
    assert entry(store)["status"] == "sent"
    assert entry(store, "finalized:1")["status"] == "superseded"
    state = store.mark_archived(RUN_ID)
    assert state["status"] == "archived"
    assert entry(store)["status"] == "sent"


def test_mark_archived_defensively_invalidates_pending_entries(tmp_path):
    store = completed_store(tmp_path)
    state = store.mark_archived(RUN_ID)
    assert {item["status"] for item in state["delivery_notifications"]} == {"superseded"}


def test_concurrent_delivery_only_calls_sender_once(tmp_path):
    store = completed_store(tmp_path)
    sending = Event()
    release = Event()
    competing = Event()
    calls = []

    def send(notification):
        calls.append(notification["id"])
        sending.set()
        assert release.wait(5)

    def competing_delivery():
        competing.set()
        return FlowcellStateStore(store.manager_dir).deliver_notification(
            RUN_ID, "processed:1", calls.append, retry=True
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        first = executor.submit(store.deliver_notification, RUN_ID, "processed:1", send)
        try:
            assert sending.wait(5)
            second = executor.submit(competing_delivery)
            assert competing.wait(5)
        finally:
            release.set()
        assert first.result(timeout=5)
        assert not second.result(timeout=5)
    assert calls == ["processed:1"]
    assert len(entry(store)["attempts"]) == 1


def test_cleanup_transition_waits_for_inflight_delivery(tmp_path):
    store = completed_store(tmp_path)
    sending = Event()
    release = Event()
    preparing = Event()

    def send(notification):
        sending.set()
        assert release.wait(5)

    def prepare():
        preparing.set()
        return FlowcellStateStore(store.manager_dir).set_preparing(
            RUN_ID, "analysis", reason="rerun", refresh_inputs=False
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        first = executor.submit(store.deliver_notification, RUN_ID, "processed:1", send)
        try:
            assert sending.wait(5)
            second = executor.submit(prepare)
            assert preparing.wait(5)
            assert store.read(RUN_ID)["status"] == "completed"
            assert entry(store)["status"] == "sending"
        finally:
            release.set()
        assert first.result(timeout=5)
        assert second.result(timeout=5)["status"] == "preparing"
    notification = entry(store)
    assert notification["status"] == "superseded"
    assert notification["sent_at"]


def test_legacy_state_load_does_not_fabricate_delivery_intents(tmp_path):
    store = completed_store(tmp_path)
    state = store.read(RUN_ID)
    del state["delivery_notifications"]
    store.write(state)
    before = store.state_path(RUN_ID).read_bytes()
    assert "delivery_notifications" not in store.read(RUN_ID)
    assert store.state_path(RUN_ID).read_bytes() == before
    store.queue(RUN_ID, "reporting")
    store.begin_attempt(RUN_ID)
    processed(store)
    notifications = store.read(RUN_ID)["delivery_notifications"]
    assert [item["id"] for item in notifications] == ["processed:2"]


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("status", "invented"),
        ("status", []),
        ("stage", "email"),
        ("payload", []),
        ("attempt", True),
        ("attempt", 0),
        ("attempts", {}),
        ("attempts", [{"status": "sent"}]),
        ("kind", ""),
        ("sent_at", 12),
        ("last_error", {}),
    ],
)
def test_invalid_notification_schema_is_rejected(tmp_path, field, value):
    store = completed_store(tmp_path)
    state = store.read(RUN_ID)
    state["delivery_notifications"][0][field] = value
    with pytest.raises(StateValidationError):
        validate_state(state)


def test_duplicate_notification_identity_is_rejected(tmp_path):
    store = completed_store(tmp_path)
    state = store.read(RUN_ID)
    state["delivery_notifications"].append(copy.deepcopy(state["delivery_notifications"][0]))
    with pytest.raises(StateValidationError, match="Duplicate"):
        validate_state(state)


def test_unknown_notification_does_not_invent_delivery(tmp_path):
    store = completed_store(tmp_path)
    calls = []
    with pytest.raises(StateConflictError, match="No notification"):
        store.deliver_notification(RUN_ID, "processed:999", calls.append, retry=True)
    assert calls == []


def test_invalid_notification_prevents_completion_write(tmp_path):
    store = reporting_store(tmp_path)
    before = store.state_path(RUN_ID).read_bytes()
    with pytest.raises(StateValidationError):
        store.complete_stage(RUN_ID, "reporting", notification={"kind": "processed"})
    assert store.state_path(RUN_ID).read_bytes() == before
