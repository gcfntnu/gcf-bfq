"""Operator cleanup must target the same output directory as the daemon."""

import logging

from dataclasses import replace
from pathlib import Path

import flowcell_manager.flowcell_manager as manager
import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore, StateConflictError, new_state
from test_notification_integration import prepare_pipeline
from test_state_integration import RUN_ID, configured_bfq, write_fastq, write_inputs

from bcl2fastq_pipeline import afterFastq, cli, findFlowCells, misc, notifications

PROJECT = "GCF-2026-001"


def prepared_outputs(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    write_inputs(output)
    write_fastq(output / PROJECT / "sample_R1.fastq.gz")
    afterFastq.md5sum_worker(cfg)
    (output / f"{PROJECT}_260918_fastq.7za").write_bytes(b"original archive")
    (output / f"md5sum_{PROJECT}_archive.txt").write_text("original archive checksum\n")
    (output / f"multiqc_{PROJECT}_260918.html").write_text("<html>report</html>")
    # Reproduce a manager invoked from /opt rather than the output root.
    unrelated = tmp_path / "opt"
    unrelated.mkdir()
    decoy = unrelated / RUN_ID
    decoy.mkdir()
    (decoy / "unrelated.7za").write_bytes(b"must not remove")
    monkeypatch.chdir(unrelated)
    return cfg, source, output


def completed_bare_path_state(tmp_path, monkeypatch):
    cfg, source, output = prepared_outputs(tmp_path, monkeypatch)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(
        new_state(
            RUN_ID, source, RUN_ID, origin="restored_legacy_fastq", start_stage="reporting", cfg=cfg
        )
    )
    store.begin_attempt(RUN_ID)
    store.complete_stage(
        RUN_ID, "reporting", notification={"kind": "processed", "payload": {"run_id": RUN_ID}}
    )
    store.start_stage(RUN_ID, "finalization")
    store.complete_run(
        RUN_ID, [PROJECT], notification={"kind": "finalized", "payload": {"run_id": RUN_ID}}
    )
    return cfg, store, output


def contents(output):
    return {
        str(path.relative_to(output)): (path.read_bytes(), path.stat().st_mtime_ns)
        for path in output.rglob("*")
        if path.is_file()
    }


@pytest.mark.parametrize("argument", ["bare", "absolute"])
def test_finalization_rerun_normalizes_bare_state_and_cleans_real_products(
    tmp_path, monkeypatch, capsys, argument
):
    _cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    before = store.read(RUN_ID)
    old_files = contents(output)

    queued = manager.rerun_flowcell(
        flowcell=RUN_ID if argument == "bare" else str(output),
        from_stage="finalization",
        force=True,
    )

    assert queued["output_path"] == str(output)
    assert store.read(RUN_ID)["output_path"] == str(output)
    assert queued["status"] == "queued"
    assert queued["current_stage"] == "finalization"
    for stage in ("demultiplexing", "analysis", "reporting"):
        assert queued["stages"][stage] == before["stages"][stage]
    assert queued["delivery_notifications"][0] == before["delivery_notifications"][0]
    assert queued["delivery_notifications"][1]["status"] == "superseded"
    assert contents(output) == {
        name: value
        for name, value in old_files.items()
        if not name.endswith(".7za") and not name.endswith("_archive.txt")
    }
    preview = capsys.readouterr().out
    assert str(output) in preview
    assert str(output / f"{PROJECT}_260918_fastq.7za") in preview
    assert "(none)" not in preview
    assert (Path.cwd() / RUN_ID / "unrelated.7za").read_bytes() == b"must not remove"


@pytest.mark.parametrize("argument", ["bare", "absolute"])
@pytest.mark.parametrize("decision", ["dry-run", "decline"])
def test_preview_or_decline_does_not_normalize_or_delete_anything(
    tmp_path, monkeypatch, capsys, argument, decision
):
    _cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    state_before = store.state_path(RUN_ID).read_bytes()
    files_before = contents(output)
    monkeypatch.setattr("builtins.input", lambda _prompt: "no")

    manager.rerun_flowcell(
        flowcell=RUN_ID if argument == "bare" else str(output),
        from_stage="finalization",
        dry_run=decision == "dry-run",
    )

    assert store.state_path(RUN_ID).read_bytes() == state_before
    assert contents(output) == files_before
    assert str(output) in capsys.readouterr().out


@pytest.mark.parametrize("decision", ["confirm", "dry-run", "decline"])
def test_legacy_bare_inventory_path_is_resolved_before_migration(
    tmp_path, monkeypatch, capsys, decision
):
    cfg, _source, output = prepared_outputs(tmp_path, monkeypatch)
    manager.add_flowcell(project=PROJECT, path=RUN_ID)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    files_before = contents(output)
    inventory_path = cfg.static.paths.manager_dir / "flowcells.processed"
    inventory_before = inventory_path.read_bytes()
    monkeypatch.setattr("builtins.input", lambda _prompt: "no")

    result = manager.rerun_flowcell(
        flowcell=RUN_ID,
        from_stage="finalization",
        force=decision == "confirm",
        dry_run=decision == "dry-run",
    )

    assert str(output) in capsys.readouterr().out
    if decision == "confirm":
        assert result["origin"] == "legacy_rerun"
        assert result["output_path"] == str(output)
        assert store.read(RUN_ID)["output_path"] == str(output)
        assert result["status"] == "queued"
        assert not (output / f"{PROJECT}_260918_fastq.7za").exists()
        assert not (output / f"md5sum_{PROJECT}_archive.txt").exists()
        assert (output / f"md5sum_{PROJECT}_fastq.txt").is_file()
    else:
        assert not store.exists(RUN_ID)
        assert contents(output) == files_before
        assert inventory_path.read_bytes() == inventory_before


@pytest.mark.parametrize("path_kind", ["absolute-mismatch", "relative-components"])
def test_ambiguous_output_path_is_refused_without_state_or_file_changes(
    tmp_path, monkeypatch, path_kind
):
    _cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    foreign_root = tmp_path if path_kind == "absolute-mismatch" else Path.cwd()
    foreign = foreign_root / "different-output" / RUN_ID
    foreign.mkdir(parents=True)
    sentinel = foreign / "unrelated.7za"
    sentinel.write_bytes(b"must not remove")
    invalid = str(foreign) if path_kind == "absolute-mismatch" else f"different-output/{RUN_ID}"
    state = store.read(RUN_ID)
    state["output_path"] = invalid
    store.write(state)
    state_before = store.state_path(RUN_ID).read_bytes()
    files_before = contents(output)

    with pytest.raises(StateConflictError):
        manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", force=True)

    assert store.state_path(RUN_ID).read_bytes() == state_before
    assert contents(output) == files_before
    assert sentinel.read_bytes() == b"must not remove"


@pytest.mark.parametrize("path_kind", ["absolute-mismatch", "relative-components"])
def test_ambiguous_legacy_inventory_is_refused_before_state_creation(
    tmp_path, monkeypatch, path_kind
):
    cfg, _source, output = prepared_outputs(tmp_path, monkeypatch)
    foreign_root = tmp_path if path_kind == "absolute-mismatch" else Path.cwd()
    foreign = foreign_root / "different-output" / RUN_ID
    foreign.mkdir(parents=True)
    sentinel = foreign / "unrelated.7za"
    sentinel.write_bytes(b"must not remove")
    invalid = str(foreign) if path_kind == "absolute-mismatch" else f"different-output/{RUN_ID}"
    manager.add_flowcell(project=PROJECT, path=invalid)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    inventory_path = cfg.static.paths.manager_dir / "flowcells.processed"
    inventory_before = inventory_path.read_bytes()
    files_before = contents(output)

    with pytest.raises(StateConflictError):
        manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", force=True)

    assert not store.exists(RUN_ID)
    assert inventory_path.read_bytes() == inventory_before
    assert contents(output) == files_before
    assert sentinel.read_bytes() == b"must not remove"


@pytest.mark.parametrize("stage", ["analysis", "reporting", "finalization"])
def test_missing_output_directory_blocks_downstream_rerun_without_mutation(
    tmp_path, monkeypatch, stage
):
    _cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    retained = output.with_name(f"{RUN_ID}.moved")
    output.rename(retained)
    before = store.state_path(RUN_ID).read_bytes()

    with pytest.raises(StateConflictError):
        manager.rerun_flowcell(flowcell=RUN_ID, from_stage=stage, force=True)

    assert store.state_path(RUN_ID).read_bytes() == before
    assert not output.exists()
    assert (retained / f"{PROJECT}_260918_fastq.7za").exists()


@pytest.mark.parametrize("location", ["state", "argument"])
def test_absolute_alias_to_the_same_output_directory_is_accepted(tmp_path, monkeypatch, location):
    _cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    alias_root = tmp_path / "output-alias"
    alias_root.symlink_to(output.parent, target_is_directory=True)
    alias = alias_root / RUN_ID
    argument = RUN_ID
    if location == "state":
        state = store.read(RUN_ID)
        state["output_path"] = str(alias)
        store.write(state)
    else:
        argument = str(alias)

    queued = manager.rerun_flowcell(flowcell=argument, from_stage="finalization", force=True)

    assert Path(queued["output_path"]).is_absolute()
    assert Path(queued["output_path"]).samefile(output)
    assert not (output / f"{PROJECT}_260918_fastq.7za").exists()
    assert not (output / f"md5sum_{PROJECT}_archive.txt").exists()
    assert (output / f"md5sum_{PROJECT}_fastq.txt").is_file()


def test_preview_identifies_output_directory_when_no_products_need_cleanup(
    tmp_path, monkeypatch, capsys
):
    _cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    (output / f"{PROJECT}_260918_fastq.7za").unlink()
    (output / f"md5sum_{PROJECT}_archive.txt").unlink()
    before = store.state_path(RUN_ID).read_bytes()

    manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", dry_run=True)

    preview = capsys.readouterr().out
    assert str(output) in preview
    assert "(none)" in preview
    assert store.state_path(RUN_ID).read_bytes() == before


def test_instrument_path_argument_identifies_run_without_overriding_output(tmp_path, monkeypatch):
    cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    source = cfg.static.paths.nova_base_dir / RUN_ID
    source_product = source / "unrelated.7za"
    source_product.write_bytes(b"instrument data")

    queued = manager.rerun_flowcell(flowcell=str(source), from_stage="finalization", force=True)

    assert queued["output_path"] == str(output)
    assert store.read(RUN_ID)["output_path"] == str(output)
    assert source_product.read_bytes() == b"instrument data"
    assert not (output / f"{PROJECT}_260918_fastq.7za").exists()


@pytest.mark.parametrize("conflict", [False, True])
def test_multiple_inventory_locations_are_validated_before_cleanup(tmp_path, monkeypatch, conflict):
    cfg, _source, output = prepared_outputs(tmp_path, monkeypatch)
    manager.add_flowcell(project=PROJECT, path=RUN_ID)
    other = tmp_path / "foreign" / RUN_ID if conflict else output
    manager.add_flowcell(project=PROJECT, path=str(other))
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    if conflict:
        before = contents(output)
        with pytest.raises(StateConflictError, match="Output path mismatch"):
            manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", force=True)
        assert not store.exists(RUN_ID)
        assert contents(output) == before
    else:
        queued = manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", force=True)
        assert queued["output_path"] == str(output)
        assert not (output / f"{PROJECT}_260918_fastq.7za").exists()


def test_relative_configured_output_root_is_refused(tmp_path, monkeypatch):
    cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    cfg.static = replace(cfg.static, paths=replace(cfg.static.paths, output_dir=Path("output")))
    before = store.state_path(RUN_ID).read_bytes()
    files = contents(output)
    with pytest.raises(StateConflictError, match="outputDir must be absolute"):
        manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", force=True)
    assert store.state_path(RUN_ID).read_bytes() == before
    assert contents(output) == files


@pytest.mark.parametrize("state_backed", [False, True])
def test_archive_normalizes_bare_paths_and_updates_legacy_inventory(
    tmp_path, monkeypatch, state_backed
):
    cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    if not state_backed:
        store.state_path(RUN_ID).unlink()
    manager.add_flowcell(project=PROJECT, path=RUN_ID)
    manager.archive_flowcell(flowcell=RUN_ID, force=True)
    assert not (output / PROJECT).exists()
    assert not (output / f"{PROJECT}_260918_fastq.7za").exists()
    assert (output / "SampleSheet.csv").is_file()
    inventory = manager._read_inventory(cfg)
    assert list(inventory["archived"]) != ["0"]
    if state_backed:
        state = store.read(RUN_ID)
        assert state["output_path"] == str(output)
        assert state["status"] == "archived"
        assert all(n["status"] == "superseded" for n in state["delivery_notifications"])
    assert (Path.cwd() / RUN_ID / "unrelated.7za").read_bytes() == b"must not remove"


def test_daemon_normalizes_before_preparation_with_one_execution_lease(tmp_path, monkeypatch):
    cfg, store, output, calls = prepare_pipeline(tmp_path, monkeypatch, start_stage="finalization")
    state = store.read(RUN_ID)
    state["output_path"] = RUN_ID
    store.write(state)
    prepare = findFlowCells.newFlowCell

    def check_preparation():
        assert store.execution_active(RUN_ID)
        assert store.read(RUN_ID)["output_path"] == str(output)
        with pytest.raises(StateConflictError, match="active"):
            manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", force=True)
        prepare()

    monkeypatch.setattr(findFlowCells, "newFlowCell", check_preparation)
    monkeypatch.setattr(cli.workflow_config, "prepare_execution", lambda *_args: None)
    monkeypatch.setattr(notifications, "send_notification", lambda *_args: None)
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("test"), prepare=True)
    assert store.read(RUN_ID)["status"] == "completed"
    assert calls == ["finalization"]
    assert (output / f"{PROJECT}_260918.7za").exists()


@pytest.mark.parametrize("problem", ["mismatch", "missing"])
def test_daemon_refuses_bad_output_before_preparing_inputs(tmp_path, monkeypatch, problem):
    cfg, store, output, calls = prepare_pipeline(tmp_path, monkeypatch, start_stage="finalization")
    state = store.read(RUN_ID)
    state["output_path"] = str(tmp_path / "foreign" / RUN_ID) if problem == "mismatch" else RUN_ID
    store.write(state)
    if problem == "missing":
        output.rename(output.with_suffix(".retained"))
    state_before = store.state_path(RUN_ID).read_bytes()

    def unexpected_prepare():
        pytest.fail("Inputs must not be prepared for a conflicting output path")

    monkeypatch.setattr(findFlowCells, "newFlowCell", unexpected_prepare)
    with pytest.raises(StateConflictError):
        cli._run_state_backed_flowcell(cfg, store, logging.getLogger("test"), prepare=True)
    assert store.state_path(RUN_ID).read_bytes() == state_before
    assert not calls


@pytest.mark.parametrize("after_replace", [False, True])
def test_normalization_write_failure_never_deletes_outputs(tmp_path, monkeypatch, after_replace):
    _cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    before = store.state_path(RUN_ID).read_bytes()
    files = contents(output)
    write = FlowcellStateStore._atomic_write_unlocked

    def fail_write(self, state):
        if after_replace:
            write(self, state)
        raise OSError("injected persistence failure")

    monkeypatch.setattr(FlowcellStateStore, "_atomic_write_unlocked", fail_write)
    with pytest.raises(OSError, match="injected persistence"):
        manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", force=True)
    assert contents(output) == files
    if after_replace:
        state = store.read(RUN_ID)
        assert state["output_path"] == str(output)
        assert state["status"] == "preparing"
        assert state["delivery_notifications"][1]["status"] == "superseded"
    else:
        assert store.state_path(RUN_ID).read_bytes() == before


def test_preparation_failure_preserves_previous_successful_attempt(tmp_path, monkeypatch):
    cfg, store, output = completed_bare_path_state(tmp_path, monkeypatch)
    attempts = store.read(RUN_ID)["attempts"]
    manager.rerun_flowcell(flowcell=RUN_ID, from_stage="finalization", force=True)
    files = contents(output)

    def preparation_error():
        raise OSError("injected input preparation failure")

    monkeypatch.setattr(findFlowCells, "newFlowCell", preparation_error)
    monkeypatch.setattr(misc, "write_error_report", lambda *_args: None)
    monkeypatch.setattr(misc, "send_error_report", lambda *_args, **_kwargs: None)
    log = logging.getLogger("test")
    try:
        cli._run_state_backed_flowcell(cfg, store, log, prepare=True)
    except OSError:
        cli.report_run_error(cfg, log, "Preparation failed", store=store, stage="finalization")
    else:
        pytest.fail("Preparation should fail")
    state = store.read(RUN_ID)
    assert state["status"] == "failed"
    assert state["attempts"] == attempts
    assert state["attempt"] == 1
    assert contents(output) == files
