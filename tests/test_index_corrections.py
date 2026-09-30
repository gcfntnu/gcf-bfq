"""Exercise durable input edits and crash boundaries without a demultiplexer."""

import hashlib

from pathlib import Path

import pytest

from bcl2fastq_pipeline.preflight import InputSelection
from bcl2fastq_pipeline.state import FlowcellStateStore, StateConflictError, new_state
from test_state_integration import RUN_ID, configured_bfq, write_inputs

from bcl2fastq_pipeline import index_corrections as corrections

SHEET = b"[Data]\r\nSample_ID,Sample_Project,index,index2\r\nsample,GCF-2026-001,AAGC,AGTC\r\n"


@pytest.fixture
def prepared(tmp_path):
    cfg, source, output = configured_bfq(tmp_path)
    for directory in (source, output):
        write_inputs(directory)
        (directory / "SampleSheet.csv").write_bytes(SHEET)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(
        new_state(
            RUN_ID, source, output, origin="new", start_stage="demultiplexing", preparing=True
        )
    )
    return store, source, output


def make_plan(prepared, indexes=("index2",)):
    store, source, output = prepared
    selection = InputSelection(output / "SampleSheet.csv", output / "Sample-Submission-Form.xlsx")
    return corrections.plan(store, RUN_ID, selection, source, indexes, reason="operator correction")


def test_plan_is_readonly_and_apply_retains_backup_history_and_attempt(prepared):
    store, source, output = prepared
    before = store.read(RUN_ID)
    correction = make_plan(prepared, ("index1", "index2"))
    assert not correction.backup_path.exists()
    assert store.read(RUN_ID) == before
    assert (output / "SampleSheet.csv").read_bytes() == SHEET

    with store.execution_lease(RUN_ID):
        applied = corrections.apply(store, RUN_ID, correction)
        queued = store.queue(RUN_ID, "demultiplexing")
        attempt = store.begin_attempt(RUN_ID)

    assert applied["status"] == "preparing"
    assert queued["status"] == "queued"
    entry = applied["index_corrections"][0]
    assert entry["status"] == "applied"
    assert entry["indexes"] == ["index1", "index2"]
    assert entry["reason"] == "operator correction"
    assert entry["changed_rows"] == {"index1": 1, "index2": 1}
    assert entry["before_sha256"] == hashlib.sha256(SHEET).hexdigest()
    assert entry["after_sha256"] == hashlib.sha256(correction.after).hexdigest()
    assert entry["orientation"]["source_sha256"] == entry["before_sha256"]
    assert correction.backup_path.read_bytes() == SHEET
    assert not correction.backup_path.is_relative_to(output)
    assert (source / "SampleSheet.csv").read_bytes() == SHEET
    assert (output / "SampleSheet.csv").read_bytes() == SHEET.replace(b"AAGC,AGTC", b"GCTT,GACT")
    assert attempt["attempts"][-1]["index_correction_id"] == correction.operation_id


def test_new_request_always_toggles_but_same_operation_is_idempotent(prepared):
    store, _, output = prepared
    first = make_plan(prepared)
    with store.execution_lease(RUN_ID):
        corrections.apply(store, RUN_ID, first)
        corrections.apply(store, RUN_ID, first)
    assert len(store.read(RUN_ID)["index_corrections"]) == 1
    second = make_plan(prepared)
    assert first.operation_id != second.operation_id
    with store.execution_lease(RUN_ID):
        corrections.apply(store, RUN_ID, second)
    assert (output / "SampleSheet.csv").read_bytes() == SHEET
    assert first.backup_path.read_bytes() == SHEET
    assert second.backup_path.read_bytes() == first.after
    assert len(store.read(RUN_ID)["index_corrections"]) == 2


def test_backup_failure_leaves_effective_sheet_and_history_untouched(prepared, monkeypatch):
    store, _, output = prepared
    correction = make_plan(prepared)
    before = store.read(RUN_ID)

    def fail(*_):
        raise OSError("backup storage unavailable")

    monkeypatch.setattr(corrections, "_write_backup", fail)
    with pytest.raises(OSError, match="backup storage"):
        corrections.apply(store, RUN_ID, correction)
    assert (output / "SampleSheet.csv").read_bytes() == SHEET
    assert store.read(RUN_ID) == before


def test_interruption_before_replace_reconciles_without_toggling(prepared, monkeypatch):
    store, _, output = prepared
    correction = make_plan(prepared)

    def fail(*_):
        raise OSError("replace unavailable")

    monkeypatch.setattr(corrections, "_replace_sheet", fail)
    with pytest.raises(OSError, match="replace unavailable"):
        corrections.apply(store, RUN_ID, correction)
    assert store.read(RUN_ID)["index_corrections"][0]["status"] == "pending"
    with pytest.raises(StateConflictError, match="pending index"):
        store.queue(RUN_ID, "demultiplexing")
    recovered = corrections.recover(store, RUN_ID)
    assert recovered["status"] == "preparing"
    assert recovered["index_corrections"][0]["status"] == "not_applied"
    assert (output / "SampleSheet.csv").read_bytes() == SHEET
    assert corrections.recover(store, RUN_ID) == recovered


def test_interruption_after_replace_reconciles_committed_sheet(prepared, monkeypatch):
    store, _, output = prepared
    correction = make_plan(prepared)

    def fail(*_, **__):
        raise OSError("history completion failed")

    with monkeypatch.context() as patch:
        patch.setattr(corrections, "_finish", fail)
        with pytest.raises(OSError, match="history completion"):
            corrections.apply(store, RUN_ID, correction)
    assert (output / "SampleSheet.csv").read_bytes() == correction.after
    recovered = corrections.recover(store, RUN_ID)
    assert recovered["index_corrections"][0]["status"] == "applied"
    assert recovered["restart_request"]["index_correction_id"] == correction.operation_id
    assert (output / "SampleSheet.csv").read_bytes() == correction.after
    assert corrections.recover(store, RUN_ID) == recovered
    # A new explicit command intentionally toggles back after recovery.
    corrections.apply(store, RUN_ID, make_plan(prepared))
    assert (output / "SampleSheet.csv").read_bytes() == SHEET


def test_interrupted_history_does_not_overwrite_later_manual_edits(prepared, monkeypatch):
    store, _, output = prepared

    def fail(*_):
        raise OSError("replace unavailable")

    with monkeypatch.context() as patch:
        patch.setattr(corrections, "_replace_sheet", fail)
        with pytest.raises(OSError):
            corrections.apply(store, RUN_ID, make_plan(prepared))
    edited = SHEET.replace(b"AGTC", b"CCAA")
    (output / "SampleSheet.csv").write_bytes(edited)
    result = corrections.recover(store, RUN_ID)
    assert result["index_corrections"][0]["status"] == "interrupted_unknown"
    assert (output / "SampleSheet.csv").read_bytes() == edited


def test_changed_effective_sheet_rejected_before_backup(prepared):
    store, _, output = prepared
    correction = make_plan(prepared)
    (output / "SampleSheet.csv").write_bytes(SHEET.replace(b"AGTC", b"CCCC"))
    with pytest.raises(StateConflictError, match="changed after"):
        corrections.apply(store, RUN_ID, correction)
    assert not correction.backup_path.exists()
    assert "index_corrections" not in store.read(RUN_ID)


def test_source_change_between_preview_and_confirmation_rejected(prepared):
    _, source, _ = prepared
    preview = make_plan(prepared)
    (source / "SampleSheet.csv").write_bytes(SHEET.replace(b"AGTC", b"CCCC"))
    current = make_plan(prepared)
    with pytest.raises(StateConflictError, match="changed after preview"):
        corrections.verify_plan(preview, current)


def test_source_unavailable_does_not_block_toggle(prepared):
    store, source, output = prepared
    (source / "SampleSheet.csv").unlink()
    correction = make_plan(prepared)
    assert correction.orientation["per_index"]["index2"]["status"] == "unavailable"
    corrections.apply(store, RUN_ID, correction)
    assert (output / "SampleSheet.csv").read_bytes() == correction.after


@pytest.mark.parametrize("link_type", ["symlink", "hardlink"])
def test_atomic_replace_detaches_output_link_and_preserves_source(prepared, link_type):
    store, source, output = prepared
    target = output / "SampleSheet.csv"
    target.unlink()
    if link_type == "symlink":
        target.symlink_to(source / "SampleSheet.csv")
    else:
        target.hardlink_to(source / "SampleSheet.csv")
    correction = make_plan(prepared)
    corrections.apply(store, RUN_ID, correction)
    assert not target.is_symlink()
    assert not target.samefile(source / "SampleSheet.csv")
    assert (source / "SampleSheet.csv").read_bytes() == SHEET


def test_live_orientation_uses_current_source_and_effective_files(prepared):
    store, source, _ = prepared
    original = corrections.current_orientation(store.read(RUN_ID))
    assert original["per_index"]["index2"]["status"] == "original"
    corrections.apply(store, RUN_ID, make_plan(prepared))
    reversed_report = corrections.current_orientation(store.read(RUN_ID))
    assert reversed_report["per_index"]["index2"]["status"] == "reversed"
    (source / "SampleSheet.csv").write_bytes(SHEET.replace(b"AGTC", b"GACT"))
    latest = corrections.current_orientation(store.read(RUN_ID))
    assert latest["per_index"]["index2"]["status"] == "original"
    assert latest["source_sha256"] != original["source_sha256"]
    historical = store.read(RUN_ID)["index_corrections"][0]["orientation"]
    assert historical["per_index"]["index2"]["status"] == "reversed"
    assert historical["source_sha256"] == original["source_sha256"]


def test_history_symlink_rejected(prepared, tmp_path):
    store, _, output = prepared
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (store.manager_dir / "input-history").symlink_to(elsewhere, target_is_directory=True)
    with pytest.raises(StateConflictError, match="symlink"):
        corrections.apply(store, RUN_ID, make_plan(prepared))
    assert (output / "SampleSheet.csv").read_bytes() == SHEET
    assert not list(elsewhere.iterdir())


def test_active_probe_is_readonly_even_without_state_directories(tmp_path):
    store = FlowcellStateStore(tmp_path / "absent")
    assert not store.execution_active(RUN_ID)
    assert not store.manager_dir.exists()
    with store.execution_lease(RUN_ID):
        lock = store._execution_path(RUN_ID)
        before = lock.read_bytes(), lock.stat().st_mtime_ns
        assert store.execution_active(RUN_ID)
        assert (lock.read_bytes(), lock.stat().st_mtime_ns) == before
    assert not store.execution_active(RUN_ID)
    assert (lock.read_bytes(), lock.stat().st_mtime_ns) == before


def test_durable_backup_survives_output_cleanup(prepared):
    store, _, output = prepared
    correction = make_plan(prepared)
    corrections.apply(store, RUN_ID, correction)
    for path in output.iterdir():
        path.unlink()
    output.rmdir()
    assert Path(store.read(RUN_ID)["index_corrections"][0]["backup_path"]).read_bytes() == SHEET
