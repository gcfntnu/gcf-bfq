"""Scoped reruns through manager, preflight, orchestration and publication.

Scientific processing and SMTP are explicit doubles. Pairing, shared metadata
validation, libprep selection, FASTQ checksums, persistent state and snapshots
are real and use only generated synthetic inputs.
"""

import copy
import gzip
import hashlib
import json
import logging
import tarfile

from pathlib import Path
from unittest.mock import Mock

import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore, new_state
from flowcell_manager import flowcell_manager as manager
from test_project_preflight import pair

from bcl2fastq_pipeline import (
    afterFastq,
    analysis_snapshots,
    cli,
    demultiplexing,
    makeFastq,
    misc,
    notification_delivery,
    notifications,
    preflight,
    project_corrections,
    sequencing_qc,
    workflow_config,
)
from tests.support.fixtures import configured_bfq, write_fastq

RUN = "990101_TEST_0001_SYNTHETIC"
A, B = "GCF-2099-901", "GCF-2099-902"
DATE = "990101"


def identity(path):
    stat = path.stat()
    return hashlib.sha256(path.read_bytes()).hexdigest(), stat.st_mtime_ns, stat.st_ino


def tree_identity(directory):
    return {path: identity(path) for path in directory.rglob("*") if path.is_file()}


@pytest.fixture
def execution(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path, run_id=RUN)
    output.mkdir()
    geometry = (
        '<RunInfo><Run><Reads><Read Number="1" NumCycles="75" IsIndexedRead="N"/>'
        '<Read Number="2" NumCycles="8" IsIndexedRead="Y"/>'
        '<Read Number="3" NumCycles="75" IsIndexedRead="N"/>'
        '</Reads><FlowcellLayout LaneCount="1"/></Run></RunInfo>'
    )
    for directory in (source, output):
        (directory / "RunInfo.xml").write_text(geometry)
        for project, sample, index, kit in (
            (A, "alpha", "AAAACCCC", "Kit A"),
            (B, "beta", "TTTTGGGG", "Kit B"),
        ):
            sheet, _form = pair(directory, project, sample, index)
            content = sheet.read_text().replace("Illumina DNA Prep", kit)
            if project == B:
                content = content.replace("[Data]", "[Settings]\nOverrideCycles,Y75;I8;N75\n[Data]")
            sheet.write_text(content)

    authoritative = tmp_path / "workflows" / "libprep.config"
    authoritative.parent.mkdir()
    authoritative.write_text(
        "Kit A PE:\n  workflow: synthetic-alpha\nKit B SE:\n  workflow: synthetic-beta\n"
    )
    monkeypatch.setattr(workflow_config, "AUTHORITATIVE_CONFIG", authoritative)
    monkeypatch.setenv("TMPDIR", str(cfg.static.paths.analysis_dir))
    monkeypatch.setattr(misc, "enoughFreeSpace", lambda: True)
    monkeypatch.setattr(misc, "write_error_report", lambda *_args: tmp_path / "synthetic.error")
    monkeypatch.setattr(misc, "send_error_report", Mock())

    workdirs = {}
    for project, sample, workflow in (
        (A, "alpha", "synthetic-alpha"),
        (B, "beta", "synthetic-beta"),
    ):
        write_fastq(output / project / f"{sample}_R1.fastq.gz")
        report = output / f"QC_{project}" / "result.txt"
        report.parent.mkdir()
        report.write_text(f"original analysis {project}")
        (output / f"multiqc_{project}_{DATE}.html").write_text(f"original report {project}")
        (output / f"{project}_{DATE}.7za").write_bytes(f"original archive {project}".encode())
        (output / f"md5sum_{project}_{DATE}_archive.txt").write_text(
            f"original archive manifest {project}"
        )
        work = cfg.static.paths.analysis_dir / f"{project}_{DATE}"
        work.mkdir()
        (work / "config.yaml").write_text(f"workflow: {workflow}\noriginal: true\n")
        (work / "Snakefile").write_text("# original synthetic workflow\n")
        (work / "data").mkdir()
        (work / "data" / "fastq").symlink_to(output / project, target_is_directory=True)
        workdirs[project] = analysis_snapshots.identify_workdir(work, RUN, project)
    afterFastq.md5sum_worker(cfg)
    stats = output / "Reports" / "legacy" / "Stats"
    stats.mkdir(parents=True)
    (stats / "Stats.json").write_text('{"original_flowcell_total": 12}')
    (output / "Stats").symlink_to(stats, target_is_directory=True)
    cfg.run.input_selection = preflight.select_run_inputs(source, output)
    demultiplexing.prepare(cfg)
    demultiplexing.record(cfg, version="original", status="completed")
    for project in (A, B):
        demultiplexing.delivery_provenance(cfg, project)

    state = new_state(RUN, source, output, origin="new", start_stage="finalization", cfg=cfg)
    state.update(status="running", current_stage="finalization", attempt=1, projects=[A, B])
    for detail in state["stages"].values():
        detail.update(status="completed", completed_at="2099-01-01T00:00:00+00:00")
    state["stages"]["finalization"]["status"] = "running"
    state["stages"]["analysis"]["metadata"] = {
        "project_workflows": {A: "synthetic-alpha", B: "synthetic-beta"},
        "workdirs": workdirs,
    }
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(state)
    snapshots = analysis_snapshots.prepare(state, [A, B])
    store.complete_run(RUN, [A, B], snapshots=snapshots)
    analysis_snapshots.recover(store, RUN)
    for project in (A, B):
        manager.add_flowcell(project=project, path=str(output), timestamp="2099-01-01T00:00:00")

    calls = []
    sent = Mock()
    monkeypatch.setattr(notifications, "send_notification", sent)

    def convert():
        calls.append(("conversion", cfg.output_path))
        assert cfg.output_path != output
        demultiplexing.prepare(cfg)
        assert cfg.run.demultiplexing_plan["passes"] == 1
        for project, sample in ((A, "alpha"), (B, "beta")):
            write_fastq(cfg.output_path / project / f"{sample}_R1.fastq.gz")
            with gzip.open(cfg.output_path / project / f"{sample}_R1.fastq.gz", "wt") as reads:
                reads.write("@regenerated\nTTTT\n+\n!!!!\n")
        converted = cfg.output_path / "Reports" / "legacy" / "Stats"
        converted.mkdir(parents=True)
        (converted / "Stats.json").write_text('{"correction_flowcell_total": 12}')
        (cfg.output_path / "Stats").symlink_to(converted, target_is_directory=True)
        demultiplexing.record(cfg, version="synthetic", status="completed")
        for project in (A, B):
            demultiplexing.delivery_provenance(cfg, project, refresh=True)
        return "bcl-convert", "synthetic"

    def analyze():
        projects = afterFastq.get_project_names(afterFastq.get_project_dirs(cfg))
        assert projects == {B}
        calls.append(("analysis", B))
        selected = workflow_config.select_workflow(cfg, B)
        assert selected.workflow == "synthetic-beta"
        assert cfg.run.project_configs[B].custom["User"] == "corrected"
        assert cfg.run.project_configs[A].custom["User"] == "test"
        work = workflow_config.project_workdir(cfg, B)
        assert work != Path(workdirs[B]["path"])
        work.mkdir(parents=True)
        (work / "config.yaml").write_text("workflow: synthetic-beta\noriginal: false\n")
        (work / "Snakefile").write_text("# corrected synthetic workflow\n")
        (work / "data").mkdir()
        (work / "data" / "fastq").symlink_to(cfg.output_path / B, target_is_directory=True)
        report = cfg.output_path / f"QC_{B}" / "result.txt"
        report.parent.mkdir(exist_ok=True)
        report.write_text("corrected analysis")
        (cfg.output_path / f"multiqc_{B}_{DATE}.html").write_text("corrected report")
        cfg.run.project_workdirs[B] = work
        return {B: analysis_snapshots.identify_workdir(work, RUN, B)}

    def report(*_args):
        calls.append(("reporting", B))
        return [B]

    def finalize(**_kwargs):
        calls.append(("finalization", B))
        (cfg.output_path / f"{B}_{DATE}.7za").write_bytes(b"corrected archive")

    def early_report(saved, *, tool):
        assert tool == "bcl-convert"
        report = saved.output_path / "Stats" / "sequencing_qc" / "report.html"
        report.parent.mkdir()
        report.write_text("synthetic conversion report")
        return {
            "projects": [B],
            "report_path": str(report),
            "summary_text": "synthetic",
            "summary_html": "synthetic",
        }

    monkeypatch.setattr(makeFastq, "bcl2fq", convert)
    monkeypatch.setattr(makeFastq, "rename_fastqs", lambda: None)
    monkeypatch.setattr(afterFastq, "analysis_steps", analyze)
    monkeypatch.setattr(cli, "_run_reporting", report)
    monkeypatch.setattr(afterFastq, "finalize", finalize)
    monkeypatch.setattr(sequencing_qc, "generate", early_report)
    return cfg, source, output, store, workdirs, calls, sent


def queue(execution, stage):
    cfg, source, output, store, _workdirs, _calls, _sent = execution
    sheet = source / f"SampleSheet_{B}.csv"
    sheet.write_text(sheet.read_text().replace("User,test", "User,corrected"))
    before = tree_identity(output)
    manager.rerun_flowcell(
        flowcell=RUN, from_stage=stage, projects=[B], refresh_inputs=True, force=True
    )
    assert {path: identity(path) for path in before} == before
    # Neither this selected instrument edit nor the unselected instrument's
    # invalid form can change inputs frozen by the successful queue operation.
    sheet.write_text("later unvalidated instrument edit")
    (source / f"{A}_Sample-Submission-Form.xlsx").write_bytes(b"unselected instrument edit")
    cfg.run.begin(source, cfg.static.paths)
    return before


@pytest.mark.parametrize("stage", ["analysis", "demultiplexing"])
def test_scoped_execution_preserves_completed_project_and_notification_independence(
    execution, monkeypatch, stage
):
    cfg, _source, output, store, workdirs, calls, sent = execution
    retained_workdir = tree_identity(Path(workdirs[A]["path"]))
    old_completion = copy.deepcopy(store.read(RUN)["project_completion"][A])
    before = queue(execution, stage)
    retained = {
        path: record for path, record in before.items() if A in path.name or A in str(path.parent)
    }
    old_stats = identity(output / "Stats" / "Stats.json")
    old_b_fastq = identity(output / B / "beta_R1.fastq.gz")
    sent.side_effect = OSError("synthetic SMTP failure")
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger(__name__), prepare=True)

    state = store.read(RUN)
    assert state["status"] == "completed", state["last_error"]
    assert {path: identity(path) for path in retained} == retained
    assert tree_identity(Path(workdirs[A]["path"])) == retained_workdir
    assert identity(output / "Stats" / "Stats.json") == old_stats
    assert state["project_completion"][A] == old_completion
    assert state["stages"]["analysis"]["metadata"]["project_workflows"] == {
        A: "synthetic-alpha",
        B: "synthetic-beta",
    }
    assert (output / f"multiqc_{B}_{DATE}.html").read_text() == "corrected report"
    assert (output / f"{B}_{DATE}.7za").read_bytes() == b"corrected archive"
    assert "User,corrected" in (output / f"SampleSheet_{B}.csv").read_text()
    assert state["stages"]["analysis"]["metadata"]["input_preflight"]["status"] == "passed"
    checked = state["stages"]["analysis"]["metadata"]["input_preflight"]["inputs"]
    assert {f"SampleSheet_{project}.csv" for project in (A, B)} <= {
        Path(item["path"]).name for item in checked
    }
    with tarfile.open(output / "provenance" / f"{B}_analysis.tar.gz") as snapshot:
        assert b"original: false" in snapshot.extractfile("config.yaml").read()
    correction = project_corrections.current(state)
    assert correction["status"] == "committed"
    assert Path(correction["working"]).name == RUN
    if stage == "analysis":
        assert identity(output / B / "beta_R1.fastq.gz")[:2] == old_b_fastq[:2]
        assert not any(call[0] == "conversion" for call in calls)
    else:
        assert any(call[0] == "conversion" for call in calls)
        assert not (Path(correction["working"]) / A).exists()
        assert correction["unused_output_cleanup"]["status"] == "completed"
        assert (
            Path(correction["working"]) / "provenance" / "demultiplexing" / "plan.json"
        ).is_file()
    own_provenance = output / "provenance" / "project-demultiplexing" / B / "plan.json"
    assert json.loads(own_provenance.read_text())["tool_version"] == (
        "original" if stage == "analysis" else "synthetic"
    )
    assert sent.call_count >= 2
    assert {item["status"] for item in state["delivery_notifications"] if item["attempt"] == 2} == {
        "failed"
    }


@pytest.mark.parametrize("failure", ["analysis", "finalization", "commit"])
def test_failed_scoped_execution_retains_both_prior_successes(execution, monkeypatch, failure):
    cfg, _source, output, store, workdirs, _calls, sent = execution
    old_success = copy.deepcopy(store.read(RUN)["project_completion"])
    old_workdirs = {
        project: tree_identity(Path(record["path"])) for project, record in workdirs.items()
    }
    before = queue(execution, "demultiplexing")
    inventory_before = manager._read_inventory(cfg).copy(deep=True)
    error = Mock(side_effect=OSError("synthetic correction failure"))
    if failure == "analysis":
        monkeypatch.setattr(afterFastq, "analysis_steps", error)
    elif failure == "finalization":
        monkeypatch.setattr(afterFastq, "finalize", error)
    else:
        monkeypatch.setattr(store, "complete_run", error)
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger(__name__), prepare=True)
    state = store.read(RUN)
    assert state["status"] == "failed"
    assert {path: identity(path) for path in before} == before
    assert state["project_completion"] == old_success
    assert manager._read_inventory(cfg).equals(inventory_before)
    assert {
        project: tree_identity(Path(record["path"])) for project, record in workdirs.items()
    } == old_workdirs
    assert not list((output / "provenance" / ".staging").glob("analysis-*.tar.gz"))
    assert sent.call_count == 0
    notification_delivery.recover_pending(cfg, store)
    assert sent.call_count == 0
