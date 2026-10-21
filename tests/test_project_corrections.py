"""Synthetic correction publication and retention; no external processing or mail."""

import hashlib
import os
import shutil

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore, StateConflictError, cleanup_plan
from flowcell_manager import flowcell_manager as manager

from bcl2fastq_pipeline import preflight, project_corrections
from tests.support.fixtures import RUN_ID, completed_state, configured_bfq

A = "GCF-2026-001"
B = "GCF-2026-002"


@dataclass(frozen=True)
class Pair:
    project_id: str
    sample_sheet: Path
    submission_form: Path


@dataclass(frozen=True)
class Selection:
    pairs: tuple

    @property
    def sample_sheet(self):
        return self.pairs[0].sample_sheet

    @property
    def submission_form(self):
        return self.pairs[0].submission_form


def _select(directory, _output=None, **_kwargs):
    return Selection(
        tuple(
            Pair(
                project,
                directory / f"SampleSheet_{project}.csv",
                directory / f"{project}_Sample-Submission-Form.xlsx",
            )
            for project in (A, B)
        )
    )


def _copy_selection(selection, output):
    output.mkdir(parents=True, exist_ok=True)
    for pair in selection.pairs:
        for path in (pair.sample_sheet, pair.submission_form):
            shutil.copy2(path, output / path.name)
    return _select(output)


def identity(path):
    stat = path.stat()
    return hashlib.sha256(path.read_bytes()).hexdigest(), stat.st_mtime_ns, stat.st_ino


@pytest.fixture
def correction(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    output.mkdir()
    for directory in (source, output):
        for project in (A, B):
            (directory / f"SampleSheet_{project}.csv").write_text(
                f"[CustomOptions]\nLibprep,Prep {project}\nUser,test\n"
                "[Data]\nSample_ID,Sample_Project,index\n"
                f"sample-{project},{project},ACGT\n"
            )
            (directory / f"{project}_Sample-Submission-Form.xlsx").write_bytes(b"synthetic form")
    for project in (A, B):
        (output / project).mkdir()
        (output / project / "reads.fastq.gz").write_bytes(f"original {project}".encode())
        (output / f"md5sum_{project}_fastq.txt").write_text("valid fastq checksum\n")
        (output / f"{project}_260918.7za").write_bytes(f"archive {project}".encode())
        (output / f"multiqc_{project}_260918.html").write_text(f"report {project}")
        (output / f"QC_{project}").mkdir()
        (output / f"QC_{project}" / "result.txt").write_text(f"analysis {project}")
    (output / "Stats").mkdir()
    (output / "Stats" / "conversion.json").write_text("original conversion")
    (output / "provenance").mkdir()
    (output / "provenance" / "demultiplexing").mkdir()
    (output / "provenance" / "demultiplexing" / "manifest.json").write_text("original plan")
    (output / "provenance" / f"{A}_analysis.tar.gz").write_bytes(b"prior snapshot A")
    (output / "provenance" / f"{B}_analysis.tar.gz").write_bytes(b"prior snapshot B")
    state = completed_state(cfg, source, output)
    state["attempt"] = 1
    state["projects"] = [A, B]
    state["stages"]["analysis"]["metadata"] = {
        "project_workflows": {A: "workflow-a", B: "workflow-b"}
    }
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(state)
    monkeypatch.setattr(preflight, "copy_run_inputs", _copy_selection)
    monkeypatch.setattr(preflight, "select_run_inputs", _select)
    monkeypatch.setattr(
        preflight, "validate_selection", lambda _selection: SimpleNamespace(ok=True)
    )
    return cfg, store, source, output


def queue(correction, *, stage="demultiplexing"):
    cfg, store, source, output = correction
    # This represents explicit selected-project refresh: A stays curated, B comes
    # from the instrument. Both effective pairs are frozen together.
    source_selection = _select(source)
    output_selection = _select(output)
    selection = Selection((output_selection.pairs[0], source_selection.pairs[1]))
    store.set_preparing(RUN_ID, stage, reason="correct B", refresh_inputs=True, projects=[B])
    project_corrections.prepare(store, RUN_ID, selection, [B], stage, refresh=True)
    store.queue(RUN_ID, stage)
    store.begin_attempt(RUN_ID, cfg=cfg)
    return project_corrections.prepare_execution(store, RUN_ID, cfg)


def finish_state(store):
    def update(state):
        for detail in state["stages"].values():
            detail["status"] = "completed"
        state["stages"]["analysis"]["metadata"] = {"project_workflows": {B: "corrected-workflow"}}
        state["stages"]["finalization"]["status"] = "running"
        return state

    store.mutate(RUN_ID, update)
    return store.complete_run(RUN_ID, [B])


def test_success_publishes_only_selected_and_retains_prior_success(correction):
    _cfg, store, source, output = correction
    shared_report = output / f"sequencer_stats_{A}_{B}_260918.html"
    shared_report.write_text("original shared conversion report")
    protected = [path for path in output.rglob("*") if path.is_file() and A in str(path)]
    before = {path: identity(path) for path in protected}
    (source / f"SampleSheet_{B}.csv").write_text(
        (source / f"SampleSheet_{B}.csv").read_text().replace("User,test", "User,corrected")
    )
    working = queue(correction)
    (working / B).mkdir()
    (working / B / "reads.fastq.gz").write_bytes(b"corrected B")
    (working / f"{B}_260918.7za").write_bytes(b"corrected archive B")
    (working / f"multiqc_{B}_260918.html").write_text("corrected report B")
    # A was regenerated by a combined conversion, but is never installed.
    (working / A).mkdir()
    (working / A / "reads.fastq.gz").write_bytes(b"unused regenerated A")
    extra = working / "GCF-2099-003"
    extra.mkdir()
    (extra / "reads.fastq.gz").write_bytes(b"unused new project")
    (working / "Stats").mkdir()
    (working / "Stats" / "conversion.json").write_text("correction conversion")
    (working / shared_report.name).write_text("corrected conversion report")
    project_corrections.publish(store, RUN_ID)
    finish_state(store)
    project_corrections.commit_publication(store, RUN_ID)

    assert {path: identity(path) for path in protected} == before
    assert (output / B / "reads.fastq.gz").read_bytes() == b"corrected B"
    assert (working / B).is_symlink()
    assert (working / B / "reads.fastq.gz").read_bytes() == b"corrected B"
    assert not (working / A).exists()
    assert not extra.exists()
    assert shared_report.read_text() == "original shared conversion report"
    assert (working / shared_report.name).read_text() == "corrected conversion report"
    report = working / f"multiqc_{B}_260918.html"
    assert report.is_symlink()
    assert report.read_text() == "corrected report B"
    assert (output / "Stats" / "conversion.json").read_text() == "original conversion"

    assert (output / "provenance" / f"{B}_analysis.tar.gz").read_bytes() == b"prior snapshot B"
    state = store.read(RUN_ID)
    record = project_corrections.current(state)
    assert record["status"] == "committed"
    assert (
        Path(record["root"]) / "previous" / B / "reads.fastq.gz"
    ).read_bytes() == f"original {B}".encode()
    assert state["projects"] == [A, B]
    assert state["project_completion"][A]["attempt"] == 1
    assert state["project_completion"][A]["workflow"] == "workflow-a"
    assert state["project_completion"][B]["attempt"] == 2
    assert state["project_completion"][B]["workflow"] == "corrected-workflow"


def test_startup_recovers_uncommitted_publication_before_discovery(correction):
    _cfg, store, _source, output = correction
    original = identity(output / B / "reads.fastq.gz")
    working = queue(correction)
    (working / B).mkdir()
    (working / B / "reads.fastq.gz").write_bytes(b"uncommitted replacement")
    project_corrections.publish(store, RUN_ID)
    assert (output / B / "reads.fastq.gz").read_bytes() == b"uncommitted replacement"
    project_corrections.recover_pending(store)
    assert identity(output / B / "reads.fastq.gz") == original
    assert project_corrections.current(store.read(RUN_ID))["status"] == "rolled_back"
    assert store.recover_interrupted(RUN_ID)["status"] == "interrupted"


def test_cleanup_failure_retries_after_a_later_restart_request(correction, monkeypatch):
    _cfg, store, _source, _output = correction
    working = queue(correction)
    (working / B).mkdir()
    (working / B / "reads.fastq.gz").write_bytes(b"replacement")
    (working / A).mkdir()
    (working / A / "reads.fastq.gz").write_bytes(b"unused")
    project_corrections.publish(store, RUN_ID)
    finish_state(store)
    remove = project_corrections.shutil.rmtree

    def fail_cleanup(_path):
        raise OSError("busy")

    monkeypatch.setattr(project_corrections.shutil, "rmtree", fail_cleanup)
    project_corrections.commit_publication(store, RUN_ID)
    assert store.read(RUN_ID)["status"] == "completed"
    assert (
        project_corrections.current(store.read(RUN_ID))["unused_output_cleanup"]["status"]
        == "failed"
    )
    monkeypatch.setattr(project_corrections.shutil, "rmtree", remove)
    store.set_preparing(RUN_ID, "analysis", reason="later rerun", refresh_inputs=False)
    assert project_corrections.current(store.read(RUN_ID)) is None
    project_corrections.recover_pending(store)
    assert not (working / A).exists()
    assert (
        store.read(RUN_ID)["project_corrections"][-1]["unused_output_cleanup"]["status"]
        == "completed"
    )


def test_failure_and_restart_restore_previous_selected_outputs(correction, monkeypatch):
    _cfg, store, _source, output = correction
    before = {path: identity(path) for path in output.rglob("*") if path.is_file()}
    working = queue(correction)
    (working / B).mkdir()
    (working / B / "reads.fastq.gz").write_bytes(b"failed replacement")
    real_replace = os.replace
    moved = []

    def interrupted(source, target):
        if str(target).startswith(str(output)) and "/states/" not in str(target):
            moved.append(target)
            if len(moved) == 3:
                raise OSError("synthetic interrupted publication")
        return real_replace(source, target)

    with monkeypatch.context() as patch:
        patch.setattr(project_corrections.os, "replace", interrupted)
        with pytest.raises(OSError, match="interrupted publication"):
            project_corrections.publish(store, RUN_ID)
    project_corrections.recover(store, RUN_ID)
    project_corrections.recover(store, RUN_ID)
    assert {path: identity(path) for path in before} == before
    assert project_corrections.current(store.read(RUN_ID))["status"] == "rolled_back"


def test_post_publication_failure_rolls_back_until_completion(correction):
    _cfg, store, _source, output = correction
    before = identity(output / B / "reads.fastq.gz")
    working = queue(correction)
    (working / B).mkdir()
    (working / B / "reads.fastq.gz").write_bytes(b"replacement")
    project_corrections.publish(store, RUN_ID)
    project_corrections.recover(store, RUN_ID)
    assert identity(output / B / "reads.fastq.gz") == before
    assert (working / B / "reads.fastq.gz").read_bytes() == b"replacement"
    assert not (working / B).is_symlink()


def test_analysis_correction_preserves_fastq_and_manifest(correction):
    _cfg, store, _source, output = correction
    fastq = output / B / "reads.fastq.gz"
    manifest = output / f"md5sum_{B}_fastq.txt"
    before = (
        fastq.read_bytes(),
        fastq.stat().st_mtime_ns,
        manifest.read_bytes(),
        manifest.stat().st_mtime_ns,
    )
    working = queue(correction, stage="analysis")
    assert not (working / A).exists()
    assert (working / B / "reads.fastq.gz").read_bytes() == before[0]
    assert not (working / f"{B}_260918.7za").exists()
    assert not (working / f"QC_{B}").exists()
    assert (
        working / "provenance" / "demultiplexing" / "manifest.json"
    ).read_text() == "original plan"
    (working / f"QC_{B}").mkdir()
    (working / f"QC_{B}" / "result.txt").write_text("corrected analysis")
    project_corrections.publish(store, RUN_ID)
    finish_state(store)
    project_corrections.recover(store, RUN_ID)
    assert (
        fastq.read_bytes(),
        fastq.stat().st_mtime_ns,
        manifest.read_bytes(),
        manifest.stat().st_mtime_ns,
    ) == before


def test_frozen_inputs_ignore_later_source_edit(correction):
    cfg, store, source, _output = correction
    working = queue(correction)
    frozen = (working / f"SampleSheet_{B}.csv").read_bytes()
    (source / f"SampleSheet_{B}.csv").write_text("later instrument edit")
    assert cfg.run.input_selection.pairs[1].sample_sheet.read_bytes() == frozen
    assert cfg.run.custom["User"] == "test"
    assert (
        project_corrections.current(store.read(RUN_ID))["inputs"][1]["sample_sheet_sha256"]
        == hashlib.sha256(frozen).hexdigest()
    )


def test_cleanup_scope_excludes_curated_inputs_and_other_projects(correction):
    _cfg, _store, _source, output = correction
    paths = cleanup_plan(output, "demultiplexing", projects=[B])
    assert output / B in paths
    assert output / f"{B}_260918.7za" in paths
    assert output / f"SampleSheet_{B}.csv" not in paths
    assert output / f"{B}_Sample-Submission-Form.xlsx" not in paths
    assert not any(A in path.name or path.name == "Stats" for path in paths)


def test_manager_preview_is_read_only_and_previews_selected_scope(correction, monkeypatch, capsys):
    _cfg, store, source, output = correction
    before = {path: identity(path) for path in output.rglob("*") if path.is_file()}
    state_before = store.read(RUN_ID)
    seen = []

    def valid(source_path, output_path, **kwargs):
        seen.append(kwargs)
        return _select(source), SimpleNamespace(ok=True)

    monkeypatch.setattr(preflight, "require_valid_inputs", valid)
    manager.rerun_flowcell(
        flowcell=RUN_ID, from_stage="analysis", projects=[B], refresh_inputs=True, dry_run=True
    )
    text = capsys.readouterr().out
    assert f"Selected projects: {B}" in text
    assert f"Retained projects: {A}" in text
    assert "Restart boundary: analysis" in text
    assert seen == [{"refresh": True, "projects": [B], "state": state_before}]
    assert store.read(RUN_ID) == state_before
    assert {path: identity(path) for path in before} == before
    assert not (output / "provenance" / "project-corrections").exists()


def test_manager_refuses_unknown_projects_before_mutation(correction):
    _cfg, store, _source, output = correction
    with pytest.raises(StateConflictError, match="Unknown selected"):
        project_corrections.validate_scope(store.read(RUN_ID), _select(output), ["GCF-2026-999"])


def test_analysis_completion_merges_retained_workflow_metadata(correction):
    _cfg, store, _source, _output = correction
    queue(correction, stage="analysis")
    store.complete_stage(
        RUN_ID,
        "analysis",
        {"project_workflows": {B: "corrected-workflow"}, "workdirs": {B: {"path": "staged"}}},
    )
    metadata = store.read(RUN_ID)["stages"]["analysis"]["metadata"]
    assert metadata["project_workflows"] == {A: "workflow-a", B: "corrected-workflow"}
    assert store.read(RUN_ID)["project_completion"][A]["workflow"] == "workflow-a"


def test_manager_queues_frozen_correction_without_changing_success(correction, monkeypatch):
    _cfg, store, source, output = correction
    before = {path: identity(path) for path in output.rglob("*") if path.is_file()}
    monkeypatch.setattr(
        preflight,
        "require_valid_inputs",
        lambda *_args, **_kwargs: (_select(source), SimpleNamespace(ok=True)),
    )
    result = manager.rerun_flowcell(
        flowcell=RUN_ID, from_stage="analysis", projects=[B], refresh_inputs=True, force=True
    )
    assert result["status"] == "queued"
    assert result["restart_request"]["projects"] == [B]
    assert result["project_completion"][A]["workflow"] == "workflow-a"
    assert project_corrections.current(store.read(RUN_ID))["status"] == "prepared"
    assert {path: identity(path) for path in before} == before


def test_multi_pair_index_toggle_rejected_before_mutation(correction, monkeypatch):
    _cfg, store, source, output = correction
    state_before = store.read(RUN_ID)
    before = {path: identity(path) for path in output.rglob("*") if path.is_file()}
    monkeypatch.setattr(
        preflight,
        "require_valid_inputs",
        lambda *_args, **_kwargs: (_select(source), SimpleNamespace(ok=True)),
    )
    with pytest.raises(
        StateConflictError, match="Index-toggle options require a single input pair"
    ):
        manager.rerun_flowcell(
            flowcell=RUN_ID, from_stage="demultiplexing", force=True, reverse_complement_index1=True
        )
    assert store.read(RUN_ID) == state_before
    assert {path: identity(path) for path in before} == before


def test_legacy_selection_with_empty_pairs_still_freezes_inputs(correction, monkeypatch):
    cfg, store, source, _output = correction
    legacy = SimpleNamespace(
        pairs=(),
        sample_sheet=source / f"SampleSheet_{A}.csv",
        submission_form=source / f"{A}_Sample-Submission-Form.xlsx",
    )

    def copy_legacy(selection, output):
        shutil.copy2(selection.sample_sheet, output / "SampleSheet.csv")
        shutil.copy2(selection.submission_form, output / "Sample-Submission-Form.xlsx")
        return select_legacy(output, output)

    def select_legacy(directory, _output):
        return SimpleNamespace(
            pairs=(),
            sample_sheet=directory / "SampleSheet.csv",
            submission_form=directory / "Sample-Submission-Form.xlsx",
        )

    monkeypatch.setattr(preflight, "copy_run_inputs", copy_legacy)
    monkeypatch.setattr(preflight, "select_run_inputs", select_legacy)
    store.set_preparing(RUN_ID, "analysis", reason="legacy B", refresh_inputs=False, projects=[B])
    project_corrections.prepare(store, RUN_ID, legacy, [B], "analysis")
    store.queue(RUN_ID, "analysis")
    store.begin_attempt(RUN_ID, cfg=cfg)
    working = project_corrections.prepare_execution(store, RUN_ID, cfg)
    assert (working / "SampleSheet.csv").read_bytes() == legacy.sample_sheet.read_bytes()
    assert len(project_corrections.current(store.read(RUN_ID))["inputs"]) == 1
