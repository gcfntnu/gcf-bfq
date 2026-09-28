import gzip
import logging

from unittest.mock import Mock
import pytest

import flowcell_manager.flowcell_manager as manager

from bcl2fastq_pipeline import afterFastq, cli, findFlowCells, makeFastq, misc
from bcl2fastq_pipeline.config import Paths, PipelineConfig, RunContext, StaticConfig
from bcl2fastq_pipeline.state import FlowcellStateStore, StateConflictError, new_state


RUN_ID = "260918_MN00686_0026_A000HCMFHF"


def configured_bfq(tmp_path):
    paths = Paths(
        ekista_base_dir=tmp_path / "ekista",
        nova_base_dir=tmp_path / "nova",
        output_dir=tmp_path / "output",
        log_dir=tmp_path / "logs",
        manager_dir=tmp_path / "manager",
        report_dir=tmp_path / "reports",
        analysis_dir=tmp_path / "analysis",
    )
    for path in (
        paths.ekista_base_dir,
        paths.nova_base_dir,
        paths.output_dir,
        paths.log_dir,
        paths.manager_dir,
        paths.report_dir,
        paths.analysis_dir,
    ):
        path.mkdir(parents=True, exist_ok=True)
    cfg = PipelineConfig(
        static=StaticConfig(
            paths=paths,
            system={"minspace": "1", "sleeptime": "1"},
            version={"pipeline": "0.3.1"},
        ),
        run=RunContext(),
    )
    PipelineConfig._instance = cfg
    source = paths.nova_base_dir / RUN_ID
    source.mkdir()
    cfg.run.begin(source, paths)
    return cfg, source, paths.output_dir / RUN_ID


def write_inputs(directory, suffix=""):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "SampleSheet.csv").write_text(
        "[CustomOptions]\n"
        "Libprep,Illumina DNA Prep\n"
        f"User,test{suffix}\n",
        encoding="utf-8",
    )
    (directory / "Sample-Submission-Form.xlsx").write_bytes(f"form{suffix}".encode())


def write_fastq(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write("@read\nACGT\n+\n!!!!\n")


def completed_state(cfg, source, output):
    state = new_state(
        RUN_ID,
        source,
        output,
        origin="new",
        start_stage="demultiplexing",
        cfg=cfg,
    )
    state["status"] = "completed"
    state["current_stage"] = "finalization"
    state["completed_at"] = "2026-09-27T12:00:00+00:00"
    state["projects"] = ["GCF-2026-001"]
    for detail in state["stages"].values():
        detail["status"] = "completed"
        detail["completed_at"] = "2026-09-27T12:00:00+00:00"
    return state


def test_discovery_json_overrides_inventory_and_legacy_markers(tmp_path):
    cfg, source, output = configured_bfq(tmp_path)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(
        new_state(
            RUN_ID,
            source,
            output,
            origin="legacy_rerun",
            start_stage="analysis",
            cfg=cfg,
        )
    )
    output.mkdir()
    (output / "bcl.done").touch()
    (output / "files.renamed").touch()
    (output / "analysis.made").touch()
    (output / "fastq.made").touch()
    manager.add_flowcell(project="GCF-2026-001", path=str(output))

    assert findFlowCells.flowCellProcessed() is False
    assert store.read(RUN_ID)["current_stage"] == "analysis"


def test_discovery_inventory_only_flowcell_is_protected(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    manager.add_flowcell(project="GCF-2026-001", path=str(output))

    assert findFlowCells.flowCellProcessed() is True
    assert not FlowcellStateStore(cfg.static.paths.manager_dir).exists(RUN_ID)


def test_discovery_restored_fastqs_bootstrap_analysis(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    write_inputs(output)
    write_fastq(output / "GCF-2026-001" / "sample_R1.fastq.gz")
    (output / "fastq.made").touch()

    assert findFlowCells.flowCellProcessed() is False

    state = FlowcellStateStore(cfg.static.paths.manager_dir).read(RUN_ID)
    assert state["origin"] == "restored_legacy_fastq"
    assert state["current_stage"] == "analysis"
    assert state["stages"]["demultiplexing"]["status"] == "completed"


def test_discovery_new_run_creates_state_before_output(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    assert not output.exists()

    assert findFlowCells.flowCellProcessed() is False

    assert FlowcellStateStore(cfg.static.paths.manager_dir).exists(RUN_ID)
    assert not output.exists()


def test_discovery_ambiguous_output_is_refused_with_actionable_report(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    output.mkdir()
    (output / "partial.txt").write_text("partial")

    assert findFlowCells.flowCellProcessed() is True
    assert not FlowcellStateStore(cfg.static.paths.manager_dir).exists(RUN_ID)
    report = cfg.static.paths.report_dir / f"{RUN_ID}.error"
    assert report.exists()
    assert "flowcell-manager initialize" in report.read_text()


def test_inventory_completion_update_is_duplicate_free(tmp_path):
    _cfg, _source, output = configured_bfq(tmp_path)

    manager.add_flowcell(project="GCF-2026-001", path=str(output), timestamp="first")
    manager.add_flowcell(project="GCF-2026-001", path=str(output), timestamp="second")

    inventory = manager.list_all()
    assert len(inventory) == 1
    assert inventory.iloc[0]["timestamp"] == "second"


def test_rerun_analysis_preserves_fastqs_and_inputs_but_invalidates_downstream(
    tmp_path, monkeypatch
):
    cfg, source, output = configured_bfq(tmp_path)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(completed_state(cfg, source, output))
    write_inputs(output)
    fastq = output / "GCF-2026-001" / "sample_R1.fastq.gz"
    write_fastq(fastq)
    qc = output / "QC_GCF-2026-001"
    qc.mkdir()
    report = output / "multiqc_GCF-2026-001_260918.html"
    report.touch()
    archive = output / "GCF-2026-001_260918.7za"
    archive.touch()
    work = tmp_path / "work" / "GCF-2026-001_260918"
    work.mkdir(parents=True)
    monkeypatch.setenv("TMPDIR", str(tmp_path / "work"))

    state = manager.rerun_flowcell(
        flowcell=RUN_ID,
        from_stage="analysis",
        force=True,
        dry_run=False,
        refresh_inputs=False,
        reason="repeat workflow",
    )

    assert state["status"] == "queued"
    assert state["current_stage"] == "analysis"
    assert state["stages"]["demultiplexing"]["status"] == "completed"
    assert fastq.exists()
    assert (output / "SampleSheet.csv").exists()
    assert not qc.exists()
    assert not report.exists()
    assert not archive.exists()
    assert not work.exists()


def test_rerun_dry_run_does_not_mutate_state_or_files(tmp_path):
    cfg, source, output = configured_bfq(tmp_path)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(completed_state(cfg, source, output))
    write_inputs(output)
    archive = output / "GCF-2026-001_260918.7za"
    archive.touch()
    before = store.read(RUN_ID)

    manager.rerun_flowcell(
        flowcell=RUN_ID,
        from_stage="finalization",
        force=True,
        dry_run=True,
        refresh_inputs=False,
        reason=None,
    )

    assert archive.exists()
    assert store.read(RUN_ID) == before


def test_rerun_refresh_inputs_is_explicit(tmp_path):
    cfg, source, output = configured_bfq(tmp_path)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(completed_state(cfg, source, output))
    write_inputs(source, "-instrument")
    write_inputs(output, "-corrected")

    manager.rerun_flowcell(
        flowcell=RUN_ID,
        from_stage="demultiplexing",
        force=True,
        dry_run=False,
        refresh_inputs=True,
        reason="take instrument inputs",
    )

    assert "test-instrument" in (output / "SampleSheet.csv").read_text()
    assert (output / "Sample-Submission-Form.xlsx").read_bytes() == b"form-instrument"


def test_initialize_analysis_requires_recognized_restored_fastqs(tmp_path):
    _cfg, _source, output = configured_bfq(tmp_path)
    output.mkdir()

    with pytest.raises(StateConflictError, match="Cannot initialize from analysis"):
        manager.initialize_flowcell(
            flowcell=RUN_ID,
            from_stage="analysis",
            force=True,
            dry_run=False,
            refresh_inputs=False,
            reason=None,
        )


def test_active_flowcell_rerun_is_refused_without_force(tmp_path):
    cfg, source, output = configured_bfq(tmp_path)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(
        new_state(
            RUN_ID,
            source,
            output,
            origin="new",
            start_stage="demultiplexing",
            cfg=cfg,
        )
    )
    store.begin_attempt(RUN_ID)

    with store.execution_lease(RUN_ID):
        with pytest.raises(StateConflictError, match="currently running"):
            manager.rerun_flowcell(
                flowcell=RUN_ID,
                from_stage="demultiplexing",
                force=False,
                dry_run=False,
                refresh_inputs=False,
                reason=None,
            )


def test_daemon_executes_only_from_queued_analysis_boundary(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    write_inputs(output)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep", "User": "test"},
        output / "SampleSheet.csv",
        output / "Sample-Submission-Form.xlsx",
    )
    cfg.run.pipeline = "rnaseq"
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(
        new_state(
            RUN_ID,
            source,
            output,
            origin="restored_legacy_fastq",
            start_stage="analysis",
            cfg=cfg,
        )
    )
    calls = []
    monkeypatch.setattr(misc, "enoughFreeSpace", lambda: True)
    monkeypatch.setattr(makeFastq, "bcl2fq", lambda: calls.append("demultiplexing"))
    monkeypatch.setattr(makeFastq, "rename_fastqs", lambda: calls.append("rename"))
    monkeypatch.setattr(afterFastq, "analysis_steps", lambda: calls.append("analysis"))
    monkeypatch.setattr(cli, "_run_reporting", lambda *_args: calls.append("reporting"))
    monkeypatch.setattr(afterFastq, "finalize", lambda: calls.append("finalization"))
    monkeypatch.setattr(misc, "finalizedEmail", lambda *_args: calls.append("final-email"))
    monkeypatch.setattr(findFlowCells, "markFinished", lambda: ["GCF-2026-001"])

    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("test"))

    assert "demultiplexing" not in calls
    assert "rename" not in calls
    assert calls == ["analysis", "reporting", "finalization", "final-email"]
    state = store.read(RUN_ID)
    assert state["status"] == "completed"
    assert state["attempt"] == 1
    assert state["projects"] == ["GCF-2026-001"]


def test_daemon_failure_is_recorded_in_state(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    write_inputs(output)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep", "User": "test"},
        output / "SampleSheet.csv",
        output / "Sample-Submission-Form.xlsx",
    )
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(
        new_state(
            RUN_ID,
            source,
            output,
            origin="restored_legacy_fastq",
            start_stage="analysis",
            cfg=cfg,
        )
    )
    monkeypatch.setattr(misc, "enoughFreeSpace", lambda: True)
    monkeypatch.setattr(afterFastq, "analysis_steps", Mock(side_effect=RuntimeError("workflow bad")))
    report = cfg.static.paths.report_dir / f"{RUN_ID}.error"
    monkeypatch.setattr(misc, "errorEmail", lambda *_args: report)

    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("test"))

    state = store.read(RUN_ID)
    assert state["status"] == "failed"
    assert state["current_stage"] == "analysis"
    assert "workflow bad" in state["last_error"]["summary"]
    assert state["last_error"]["report_path"] == str(report)
