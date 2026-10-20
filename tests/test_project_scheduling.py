"""Scheduling policy with real metadata/state and synthetic external processing."""

import copy
import json
import logging
import re

from decimal import localcontext
from pathlib import Path
from unittest.mock import Mock

import pytest
import yaml

from bcl2fastq_pipeline.state import FlowcellStateStore, StateConflictError, new_state
from flowcell_manager import flowcell_manager as manager
from openpyxl import load_workbook
from test_project_preflight import pair

from bcl2fastq_pipeline import (
    afterFastq,
    cli,
    misc,
    notifications,
    preflight,
    project_scheduling,
    workflow_config,
)
from tests.support.fixtures import configured_bfq, write_fastq, write_inputs

RUN = "990101_SYNTHETIC_0001_TEST"
PROJECTS = tuple(f"GCF-2099-{index:03d}" for index in range(1, 7))
INDEXES = ("AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT", "ACACACAC", "TGTGTGTG")


def set_priority(sheet, value):
    content = re.sub(r"^Priority,.*\n", "", sheet.read_text(), flags=re.MULTILINE)
    sheet.write_text(content.replace("[Data]", f"Priority,{value}\n[Data]"))


def paired_inputs(directory, priorities):
    directory.mkdir(parents=True, exist_ok=True)
    for project, index, priority in zip(PROJECTS, INDEXES, priorities):
        sheet, _form = pair(directory, project, f"sample-{project}", index)
        if priority is not None:
            set_priority(sheet, priority)
    return preflight.select_run_inputs(directory, directory)


@pytest.mark.parametrize("custom", [{}, {"Priority": ""}, {"Priority": "  "}])
def test_missing_and_blank_priority_use_default(custom):
    assert project_scheduling.resolve_priority(custom) == {
        "priority": "1.0",
        "priority_source": "default",
    }


@pytest.mark.parametrize(
    "value", ["-1", "-0.00001", "NaN", "sNaN", "Infinity", "-inf", "true", "1_0", "a", "1,5", "١"]
)
def test_invalid_explicit_priority_is_rejected(value):
    with pytest.raises(ValueError, match="finite, non-negative decimal"):
        project_scheduling.resolve_priority({"Priority": value})


def test_decimal_order_is_exact_even_with_low_decimal_context_precision(tmp_path):
    selection = paired_inputs(
        tmp_path / "inputs", ["1.00", "1", "1.00000000000000000001", "0.5", "0", "2"]
    )
    with localcontext() as context:
        context.prec = 2
        result = preflight.validate_selection(selection)
        assert result.ok, result.render_text()
        plan = result.summary["project_schedule"]
        assert plan["order"] == [PROJECTS[i] for i in (4, 3, 0, 1, 2, 5)]
        assert project_scheduling.ordered(json.loads(json.dumps(plan))) == plan["order"]
    assert plan["projects"][PROJECTS[2]]["priority"] == "1.00000000000000000001"
    assert all(len(item["sample_sheet_sha256"]) == 64 for item in plan["projects"].values())


def test_defaulted_priority_plan_is_read_only_and_preview_matches(tmp_path, capsys):
    source = tmp_path / "source"
    selection = paired_inputs(source, ["2", None, "", "0.5", "0", "1.00"])
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in source.iterdir()}
    _selection, result = preflight.require_valid_inputs(source, tmp_path / "output")
    plan = result.summary["project_schedule"]
    assert plan["order"] == [PROJECTS[i] for i in (4, 3, 1, 2, 5, 0)]
    assert plan["projects"][PROJECTS[1]]["priority_source"] == "default"
    assert plan["projects"][PROJECTS[2]]["priority_source"] == "default"
    assert plan["projects"][PROJECTS[5]]["priority_source"] == "samplesheet"
    expected = " -> ".join(plan["order"])
    assert f"Analysis order: {expected}" in capsys.readouterr().out
    assert {
        path: (path.read_bytes(), path.stat().st_mtime_ns) for path in source.iterdir()
    } == before
    assert not (tmp_path / "output").exists()
    assert preflight.validate_selection(selection).summary["project_schedule"] == plan


def test_one_legacy_pair_applies_its_priority_to_all_projects(tmp_path):
    write_inputs(tmp_path)
    sheet = tmp_path / "SampleSheet.csv"
    sheet.write_text(sheet.read_text() + "sample-2,GCF-2026-002,TGCA\n")
    form = tmp_path / "Sample-Submission-Form.xlsx"
    workbook = load_workbook(form)
    workbook["Sample-Submission-Form"].append(["sample-2", "group"])
    workbook["INFO (GCF-lab only)"].append(["sample-2"])
    workbook.save(form)
    set_priority(sheet, "0.5")
    selection = preflight.select_run_inputs(tmp_path, tmp_path)
    assert not selection.multi_pair
    result = preflight.validate_selection(selection)
    assert result.ok, result.render_text()
    plan = result.summary["project_schedule"]
    assert plan["order"] == ["GCF-2026-001", "GCF-2026-002"]
    assert all(item["priority"] == "0.5" for item in plan["projects"].values())


@pytest.mark.parametrize("value", ["NaN", "-1", "true", "broken"])
def test_invalid_curated_priority_blocks_scoped_refresh_before_state_or_cleanup(tmp_path, value):
    cfg, source, output = configured_bfq(tmp_path, run_id=RUN)
    paired_inputs(source, ["0", "1"])
    selection = paired_inputs(output, ["0", "1"])
    set_priority(selection.pairs[0].sample_sheet, value)
    for project in PROJECTS[:2]:
        write_fastq(output / project / "synthetic_R1.fastq.gz")
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    state = new_state(RUN, source, output, origin="new", start_stage="analysis", cfg=cfg)
    state["projects"] = list(PROJECTS[:2])
    store.create(state)
    before = store.read(RUN)
    retained = {
        path: (path.read_bytes(), path.stat().st_mtime_ns)
        for path in output.rglob("*")
        if path.is_file()
    }
    with pytest.raises(preflight.PreflightValidationError, match="scheduling.priority"):
        manager.rerun_flowcell(
            flowcell=RUN,
            from_stage="analysis",
            projects=[PROJECTS[1]],
            refresh_inputs=True,
            force=True,
        )
    assert store.read(RUN) == before
    assert {
        path: (path.read_bytes(), path.stat().st_mtime_ns)
        for path in output.rglob("*")
        if path.is_file()
    } == retained


@pytest.fixture
def execution(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path, run_id=RUN)
    paired_inputs(source, ["2", None, "0.5", "0"])
    geometry = '<RunInfo><Run><Reads><Read Number="1" NumCycles="75" IsIndexedRead="N"/><Read Number="2" NumCycles="8" IsIndexedRead="Y"/><Read Number="3" NumCycles="75" IsIndexedRead="N"/></Reads><FlowcellLayout LaneCount="1"/></Run></RunInfo>'
    (source / "RunInfo.xml").write_text(geometry)
    for project in PROJECTS[:4]:
        write_fastq(output / project / "synthetic_R1.fastq.gz")
    authoritative = tmp_path / "workflows/libprep.config"
    authoritative.parent.mkdir()
    authoritative.write_text("Illumina DNA Prep PE:\n  workflow: default\n")
    monkeypatch.setattr(workflow_config, "AUTHORITATIVE_CONFIG", authoritative)
    monkeypatch.setenv("TMPDIR", str(cfg.static.paths.analysis_dir))
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(tmp_path / "cache"))
    monkeypatch.setattr(misc, "enoughFreeSpace", lambda: True)
    monkeypatch.setattr(misc, "send_error_report", Mock())
    monkeypatch.setattr(afterFastq, "multiqc_stats", lambda _cfg: None)

    real_check_call = afterFastq.subprocess.check_call

    def configmaker(command, cwd, **kwargs):
        if command[0] == "md5sum":
            return real_check_call(command, cwd=cwd, **kwargs)
        assert command[1] == "/opt/conda/bin/configmaker.py"
        (cwd / "config.yaml").write_text("workflow: default\n")
        (cwd / "Snakefile").write_text("# Synthetic, never executed\n")

    monkeypatch.setattr(afterFastq.subprocess, "check_call", configmaker)
    monkeypatch.setattr(notifications, "send_notification", Mock())
    events = []
    edited = False

    def external_work(command, cwd, log_path):
        nonlocal edited
        if command[0] == "7za":
            archive = Path(command[2])
            if archive.name.startswith("QC_"):
                archive.write_bytes(b"synthetic QC archive")
            else:
                project = archive.name.rsplit("_", 1)[0]
                events.append(("archiving", project))
                archive.write_bytes(b"synthetic FASTQ archive")
            return
        project = cwd.name.rsplit("_", 1)[0]
        events.append(("analysis", project))
        # Edits to the original curated sheet must not reprioritize the attempt.
        if not edited:
            original = output / f"SampleSheet_{PROJECTS[0]}.csv"
            original.write_text(original.read_text().replace("Priority,2", "Priority,0"))
            edited = True
        report_dir = cwd / "data/tmp/default/bfq"
        report_dir.mkdir(parents=True)
        (report_dir / f"multiqc_{project}.html").write_text("synthetic report")
        (report_dir / ".multiqc_config.yaml").write_text("{}")
        (cwd / "data/tmp/sample_info.tsv").write_text("synthetic sample info")

    monkeypatch.setattr(afterFastq, "run_logged_command", external_work)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    preflight.copy_run_inputs(preflight.select_run_inputs(source, output), output)
    store.create(new_state(RUN, source, output, origin="new", start_stage="analysis", cfg=cfg))
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", force=True)
    return cfg, store, events


def test_real_stage_loops_use_frozen_priority_and_keep_analysis_phase_first(execution, capsys):
    cfg, store, events = execution
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("scheduling-test"), prepare=True)
    state = store.read(RUN)
    assert state["status"] == "completed", state.get("last_error")
    cfg.run.begin(Path(state["source_path"]), cfg.static.paths)
    expected = [PROJECTS[i] for i in (3, 2, 1, 0)]
    assert events == [
        (stage, project) for stage in ("analysis", "archiving") for project in expected
    ]
    schedule = state["attempts"][-1]["project_schedule"]
    assert schedule["order"] == expected
    for stage in ("analysis", "archiving"):
        assert [item["project"] for item in schedule["actual_order"][stage]] == expected
    assert schedule["projects"][PROJECTS[0]]["priority"] == "2"
    cfg_file = yaml.safe_load((cfg.output_path / "bcl2fastq.ini").read_text())
    assert cfg_file["run"]["project_schedule"]["order"] == expected
    for project in expected:
        project_json = workflow_config.project_workdir(cfg, project) / "inputs/project.json"
        assert json.loads(project_json.read_text())["scheduling"] == schedule["projects"][project]
    manager.status_flowcell(flowcell=RUN)
    assert "Actual archiving starts: " + " -> ".join(expected) in capsys.readouterr().out

    # A downstream retry uses the completed analysis plan rather than the edited sheet.
    manager.rerun_flowcell(flowcell=RUN, from_stage="finalization", force=True)
    events.clear()
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("retry-test"), prepare=True)
    retried = store.read(RUN)
    assert retried["status"] == "completed", retried.get("last_error")
    assert events == [("archiving", project) for project in expected]
    retry_plan = retried["attempts"][-1]["project_schedule"]
    assert retry_plan["resolved_in_attempt"] == schedule["attempt"]
    assert retry_plan["actual_order"]["analysis"] == []
    assert retry_plan["projects"] == schedule["projects"]


def test_stage_start_rejects_projects_outside_frozen_plan_and_stale_attempt(execution):
    cfg, store, _events = execution
    state = store.begin_attempt(RUN, cfg=cfg)
    cfg.run.input_selection = preflight.select_run_inputs(cfg.output_path, cfg.output_path)
    result = preflight.run_preflight(cfg, store, "analysis")
    project_scheduling.freeze(cfg, store, result)
    with pytest.raises(StateConflictError, match="absent from the frozen"):
        project_scheduling.project_order(cfg, [PROJECTS[-1]])
    previous = copy.deepcopy(cfg.run.project_schedule)
    store.fail_stage(RUN, "analysis", summary="synthetic failure", report_path=None)
    with pytest.raises(StateConflictError, match="running attempt"):
        project_scheduling.record_start(cfg, "analysis", PROJECTS[0])
    assert store.read(RUN)["attempts"][-1]["project_schedule"] == previous
    assert state["attempt"] == previous["attempt"]


def test_changed_effective_inputs_cannot_replace_an_active_plan(execution):
    cfg, store, _events = execution
    store.begin_attempt(RUN, cfg=cfg)
    cfg.run.input_selection = preflight.select_run_inputs(cfg.output_path, cfg.output_path)
    result = preflight.run_preflight(cfg, store, "analysis")
    project_scheduling.freeze(cfg, store, result)
    original_plan = copy.deepcopy(cfg.run.project_schedule)
    set_priority(cfg.run.input_selection.pairs[0].sample_sheet, "0")
    changed = preflight.run_preflight(cfg, store, "analysis")
    assert changed.ok
    with pytest.raises(StateConflictError, match="Effective inputs differ"):
        project_scheduling.freeze(cfg, store, changed)
    assert store.read(RUN)["attempts"][-1]["project_schedule"] == original_plan


def test_scoped_correction_preserves_unselected_producing_priorities(execution):
    cfg, store, events = execution
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("initial-test"), prepare=True)
    original = store.read(RUN)
    assert original["status"] == "completed", original.get("last_error")
    cfg.run.begin(Path(original["source_path"]), cfg.static.paths)
    selected = PROJECTS[1]
    retained = {
        project: copy.deepcopy(original["project_completion"][project])
        for project in PROJECTS[:4]
        if project != selected
    }
    # The external double edited the first project's canonical priority to zero.
    # Correcting a different project must retain the first's producing value (2).
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[selected], force=True)
    events.clear()
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("correction-test"), prepare=True)
    corrected = store.read(RUN)
    assert corrected["status"] == "completed", corrected.get("last_error")
    assert events == [("analysis", selected), ("archiving", selected)]
    assert {project: corrected["project_completion"][project] for project in retained} == retained
    saved = corrected["stages"]["analysis"]["metadata"]["project_schedule"]
    assert saved["projects"][PROJECTS[0]]["priority"] == "2"
    assert set(saved["retained_projects"]) == set(retained)

    cfg.run.begin(Path(original["source_path"]), cfg.static.paths)
    manager.rerun_flowcell(flowcell=RUN, from_stage="finalization", force=True)
    events.clear()
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("retry-test"), prepare=True)
    final = store.read(RUN)
    assert final["status"] == "completed", final.get("last_error")
    assert events == [("archiving", PROJECTS[i]) for i in (3, 2, 1, 0)]
