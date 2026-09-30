"""Cleanup safety and restoration checks, using isolated output/state trees."""

import flowcell_manager.flowcell_manager as manager
import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore, StateConflictError
from test_state_integration import (
    RUN_ID,
    completed_state,
    configured_bfq,
    write_fastq,
    write_inputs,
)

from bcl2fastq_pipeline import fastq_cleanup

PROJECT = "GCF-2026-001"


@pytest.fixture
def delivered(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    write_inputs(output)
    write_inputs(source)
    write_fastq(output / PROJECT / "nested" / "sample_R1.fastq.gz")
    write_fastq(output / "Undetermined_S0_R1.fastq.gz")
    for name in (
        f"{PROJECT}_260918.7za",
        f"md5sum_{PROJECT}_260918_archive.txt",
        f"QC_{PROJECT}_260918.7za",
        f"md5sum_{PROJECT}_fastq.txt",
        f"encryption.{PROJECT}",
        "config.yaml",
        "report.html",
    ):
        (output / name).write_text("retained product")
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(completed_state(cfg, source, output))

    def no_external_commands(*args, **kwargs):
        raise AssertionError("Cleanup must not run checksum/archive/subprocess commands")

    monkeypatch.setattr("subprocess.check_call", no_external_commands)
    monkeypatch.setattr("subprocess.check_output", no_external_commands)
    monkeypatch.setattr("subprocess.run", no_external_commands)
    return cfg, store, output


def state_bytes(store):
    return store.state_path(RUN_ID).read_bytes()


def retained_bytes(output):
    return {
        str(p.relative_to(output)): p.read_bytes()
        for p in output.rglob("*")
        if p.is_file() and not p.name.endswith(fastq_cleanup.FASTQ_SUFFIXES)
    }


@pytest.mark.parametrize("absolute", [False, True])
@pytest.mark.parametrize("mode", ["dry-run", "decline", "confirm"])
def test_preview_confirmation_and_repeat(delivered, monkeypatch, capsys, absolute, mode):
    _cfg, store, output = delivered
    before = state_bytes(store)
    retained = retained_bytes(output)
    files = fastq_cleanup.plan(output)
    size = sum(info["size"] for info in files.values())
    monkeypatch.setattr("builtins.input", lambda _: "yes" if mode == "confirm" else "no")
    result = manager.clean_fastqs(
        flowcell=str(output) + "/" if absolute else RUN_ID, dry_run=mode == "dry-run"
    )
    preview = capsys.readouterr().out
    assert str(output) in preview
    assert f"{size:,} bytes" in preview
    assert all(str(output / name) in preview for name in files)
    assert retained_bytes(output) == retained
    if mode != "confirm":
        assert state_bytes(store) == before
        assert fastq_cleanup.plan(output) == files
        return
    assert fastq_cleanup.plan(output) == {}
    assert result["status"] == "completed"
    assert result["fastq_cleanup"]["status"] == "completed"
    assert set(result["fastq_cleanup"]["required_files"]) == set(files)
    original = __import__("json").loads(before)
    for key in ("stages", "completed_at", "attempts", "delivery_notifications"):
        assert result[key] == original[key]
    repeated = manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert repeated["fastq_cleanup"]["required_files"] == result["fastq_cleanup"]["required_files"]
    assert repeated["fastq_cleanup"]["removed_files"] == []


def test_symlink_targets_and_directory_links_never_deleted(delivered, tmp_path):
    _cfg, store, output = delivered
    external = tmp_path / "external"
    external.mkdir()
    target = external / "outside.fastq.gz"
    target.write_bytes(b"a large external target" * 100)
    directory_link = output / "linked-directory"
    directory_link.symlink_to(external, target_is_directory=True)
    file_link = output / PROJECT / "linked.fastq.gz"
    file_link.symlink_to(target)
    dangling = output / "dangling.fastq.gz"
    dangling.symlink_to(external / "absent")
    # Even a directory whose name ends in .fastq.gz must not be removed.
    (output / "directory.fastq.gz").symlink_to(external, target_is_directory=True)
    expected = sum(i["size"] for i in fastq_cleanup.plan(output).values() if not i["symlink"])
    result = manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert target.is_file()
    assert directory_link.is_symlink()
    assert (output / "directory.fastq.gz").is_symlink()
    assert not file_link.is_symlink()
    assert not dangling.is_symlink()
    assert result["fastq_cleanup"]["estimated_removed_bytes"] == expected


@pytest.mark.parametrize("suffix", fastq_cleanup.FASTQ_SUFFIXES)
def test_recursive_recognized_suffixes(delivered, suffix):
    _cfg, _store, output = delivered
    fastq = output / PROJECT / "a" / "b" / "c" / f"reads{suffix}"
    fastq.parent.mkdir(parents=True)
    fastq.write_bytes(b"reads")
    manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert not fastq.exists()


@pytest.mark.parametrize(
    "missing", [f"{PROJECT}_260918.7za", f"md5sum_{PROJECT}_260918_archive.txt"]
)
def test_missing_required_delivery_product_prevents_deletion(delivered, missing):
    _cfg, store, output = delivered
    (output / missing).unlink()
    before = state_bytes(store)
    files = fastq_cleanup.plan(output)
    with pytest.raises(StateConflictError, match="Missing required delivery"):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert state_bytes(store) == before
    assert fastq_cleanup.plan(output) == files


def test_all_projects_require_archives_even_with_no_remaining_reads(delivered):
    _cfg, _store, output = delivered
    (output / "GCF-2026-002").mkdir()
    with pytest.raises(StateConflictError, match="GCF-2026-002_260918.7za"):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)


@pytest.mark.parametrize("status", ["queued", "running", "failed", "archived"])
def test_ineligible_processing_status(delivered, status):
    _cfg, store, output = delivered
    store.mutate(RUN_ID, lambda s: {**s, "status": status})
    with pytest.raises(StateConflictError, match="completed finalization"):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert fastq_cleanup.plan(output)


def test_incomplete_finalization(delivered):
    _cfg, store, output = delivered

    def update(state):
        state["stages"]["finalization"]["status"] = "pending"
        return state

    store.mutate(RUN_ID, update)
    with pytest.raises(StateConflictError, match="completed finalization"):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert fastq_cleanup.plan(output)


def test_conflicting_argument_and_state_locations(delivered, tmp_path):
    _cfg, store, output = delivered
    other = tmp_path / RUN_ID
    other.mkdir()
    with pytest.raises(StateConflictError, match="Output path mismatch"):
        manager.clean_fastqs(flowcell=str(other), force=True)
    store.mutate(RUN_ID, lambda s: {**s, "output_path": str(other)})
    with pytest.raises(StateConflictError, match="Output path mismatch"):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert fastq_cleanup.plan(output)


def test_unavailable_output(delivered):
    _cfg, _store, output = delivered
    output.rename(output.with_name("unmounted"))
    with pytest.raises(StateConflictError, match="Cannot inspect FASTQs"):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)


def test_active_execution_blocks_even_force(delivered):
    _cfg, store, output = delivered
    before = state_bytes(store)
    with store.execution_lease(RUN_ID):
        with pytest.raises(StateConflictError, match="active execution lease"):
            manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert state_bytes(store) == before
    assert fastq_cleanup.plan(output)


@pytest.mark.parametrize("change", ["archive", "files", "state"])
def test_revalidation_after_confirmation(delivered, monkeypatch, change):
    _cfg, store, output = delivered

    def confirm(_):
        if change == "archive":
            (output / f"{PROJECT}_260918.7za").unlink()
        elif change == "files":
            (output / "new.fastq.gz").write_bytes(b"new")
        else:
            store.mutate(RUN_ID, lambda s: s)
        return "yes"

    monkeypatch.setattr("builtins.input", confirm)
    with pytest.raises(StateConflictError):
        manager.clean_fastqs(flowcell=RUN_ID)
    assert (output / "Undetermined_S0_R1.fastq.gz").exists()
    assert "fastq_cleanup" not in store.read(RUN_ID)


def test_partial_failure_is_recorded_and_retryable(delivered, monkeypatch):
    _cfg, store, output = delivered
    original = fastq_cleanup.unlink_fastq

    def failing(root_fd, name, info):
        if name.startswith(PROJECT):
            raise PermissionError("simulated failure")
        return original(root_fd, name, info)

    monkeypatch.setattr(fastq_cleanup, "unlink_fastq", failing)
    with pytest.raises(StateConflictError, match="partially failed"):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)
    partial = store.read(RUN_ID)
    assert partial["status"] == "completed"
    assert partial["fastq_cleanup"]["status"] == "partial"
    assert partial["fastq_cleanup"]["removed_files"] == ["Undetermined_S0_R1.fastq.gz"]
    monkeypatch.setattr(fastq_cleanup, "unlink_fastq", original)
    result = manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert result["fastq_cleanup"]["status"] == "completed"
    assert not fastq_cleanup.plan(output)
    assert len(result["fastq_cleanup"]["required_files"]) == 2


def test_interruption_has_durable_manifest_and_holds_lease(delivered, monkeypatch):
    _cfg, store, output = delivered

    def interrupt(root_fd, name, info):
        assert store.execution_active(RUN_ID)
        current = store.read(RUN_ID)
        assert current["fastq_cleanup"]["status"] == "in_progress"
        assert name in current["fastq_cleanup"]["required_files"]
        raise KeyboardInterrupt()

    monkeypatch.setattr(fastq_cleanup, "unlink_fastq", interrupt)
    with pytest.raises(KeyboardInterrupt):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert store.read(RUN_ID)["fastq_cleanup"]["status"] == "interrupted"
    assert not store.execution_active(RUN_ID)
    assert fastq_cleanup.plan(output)


@pytest.mark.parametrize("stage", ["analysis", "reporting", "finalization"])
def test_rerun_requires_all_restored_reads_before_archive_invalidation(delivered, stage):
    _cfg, store, output = delivered
    saved = {name: (output / name).read_bytes() for name in fastq_cleanup.plan(output)}
    manager.clean_fastqs(flowcell=RUN_ID, force=True)
    before = state_bytes(store)
    retained = retained_bytes(output)
    with pytest.raises(StateConflictError, match="Restore all cleaned FASTQs"):
        manager.rerun_flowcell(flowcell=RUN_ID, from_stage=stage, force=True)
    assert state_bytes(store) == before
    assert retained_bytes(output) == retained
    first = next(iter(saved))
    (output / first).write_bytes(saved[first])
    with pytest.raises(StateConflictError, match="Restore all cleaned FASTQs"):
        manager.rerun_flowcell(flowcell=RUN_ID, from_stage=stage, force=True)
    for name, data in saved.items():
        (output / name).write_bytes(data)
    # No changes or checksum commands needed to recognize restoration.
    fastq_cleanup.require_restored(store.read(RUN_ID), output, stage)
    if stage == "finalization":
        queued = manager.rerun_flowcell(flowcell=RUN_ID, from_stage=stage, force=True)
        assert queued["status"] == "queued"


def test_demultiplexing_can_regenerate_and_clear_old_restoration_requirement(
    delivered, monkeypatch
):
    _cfg, store, output = delivered
    manager.clean_fastqs(flowcell=RUN_ID, force=True)
    # Preflight is tested elsewhere; this case checks actual restart/state lifecycle.
    monkeypatch.setattr(manager.preflight, "require_valid_inputs", lambda *a, **k: (None, None))
    queued = manager.rerun_flowcell(flowcell=RUN_ID, from_stage="demultiplexing", force=True)
    assert queued["status"] == "queued"
    store.begin_attempt(RUN_ID)
    store.complete_stage(RUN_ID, "demultiplexing")
    state = store.read(RUN_ID)
    assert state["fastq_cleanup"]["status"] == "regenerated"
    assert state["fastq_cleanup"]["required_files"] == {}
    fastq_cleanup.require_restored(state, output, "analysis")


def test_unlink_refuses_parent_replaced_by_external_symlink(delivered, tmp_path, monkeypatch):
    _cfg, store, output = delivered
    external = tmp_path / "external"
    external.mkdir()
    (external / "nested").mkdir()
    target = external / "nested" / "sample_R1.fastq.gz"
    target.write_bytes(b"external data")
    original = fastq_cleanup.unlink_fastq

    def swap(root_fd, name, info):
        if name.startswith(PROJECT):
            (output / PROJECT).rename(output / "saved-project")
            (output / PROJECT).symlink_to(external, target_is_directory=True)
        return original(root_fd, name, info)

    monkeypatch.setattr(fastq_cleanup, "unlink_fastq", swap)
    with pytest.raises(StateConflictError, match="partially failed"):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert target.read_bytes() == b"external data"
    assert store.read(RUN_ID)["fastq_cleanup"]["errors"]


def test_file_replaced_after_preview_is_not_removed(delivered, monkeypatch):
    _cfg, _store, output = delivered
    original = fastq_cleanup.unlink_fastq

    def replace(root_fd, name, info):
        (output / name).write_bytes(b"changed file")
        return original(root_fd, name, info)

    monkeypatch.setattr(fastq_cleanup, "unlink_fastq", replace)
    with pytest.raises(StateConflictError, match="FASTQ changed after preview"):
        manager.clean_fastqs(flowcell=RUN_ID, force=True)
    assert all(
        (output / name).read_bytes() == b"changed file" for name in fastq_cleanup.plan(output)
    )
