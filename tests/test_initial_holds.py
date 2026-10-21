"""Optional holds use shared parsing, atomic state and the existing project lifecycle."""

import copy
import json
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
from openpyxl import load_workbook
from test_project_holds import A, B, control, fingerprint
from test_project_lifecycle import restart_context, run
from test_project_scheduling import RUN, configure_execution, paired_inputs

from bcl2fastq_pipeline import (
    initial_holds,
    makeFastq,
    preflight,
    project_holds,
    sequencing_delivery,
)
from tests.support.fixtures import write_inputs


def set_initial(sheet, value):
    sheet.write_text(sheet.read_text().replace("[Data]", f"InitialHold,{value}\n[Data]"))


@pytest.mark.parametrize(
    "value,held,source",
    [
        (None, False, "default"),
        ("", False, "default"),
        ("  ", False, "default"),
        ("false", False, "samplesheet"),
        (" fAlSe ", False, "samplesheet"),
        ("true", True, "samplesheet"),
        (" TRUE ", True, "samplesheet"),
    ],
)
def test_effective_optional_value_and_read_only_preview(tmp_path, capsys, value, held, source):
    directory = tmp_path / "inputs"
    selection = paired_inputs(directory, ["2", "0"])
    if value is not None:
        set_initial(selection.pairs[1].sample_sheet, value)
    before = fingerprint(directory.iterdir())
    _selection, result = preflight.require_valid_inputs(directory, tmp_path / "output")
    entries = result.summary["project_initial_holds"]
    assert entries[B]["held"] is held
    assert entries[B]["source"] == source
    assert len(entries[B]["sample_sheet_sha256"]) == 64
    assert result.summary["project_schedule"]["order"] == [B, A]
    assert f"InitialHold={str(held).lower()} ({source})" in capsys.readouterr().out
    assert fingerprint(directory.iterdir()) == before
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("value", ["1", "0", "yes", "no", "null", "truthy"])
def test_invalid_curated_value_blocks_selected_refresh_before_cleanup(tmp_path, monkeypatch, value):
    cfg, store, _events = configure_execution(tmp_path, monkeypatch, ["2", "0"])
    set_initial(cfg.output_path / f"SampleSheet_{A}.csv", value)
    before = store.read(RUN)
    files = fingerprint(p for p in cfg.output_path.rglob("*") if p.is_file())
    with pytest.raises(preflight.PreflightValidationError, match="scheduling.initial_hold"):
        manager.rerun_flowcell(
            flowcell=RUN, from_stage="analysis", projects=[B], refresh_inputs=True, force=True
        )
    assert store.read(RUN) == before
    assert fingerprint(files) == files


def test_legacy_shared_sheet_seeds_all_represented_projects(tmp_path):
    write_inputs(tmp_path)
    sheet = tmp_path / "SampleSheet.csv"
    sheet.write_text(sheet.read_text() + "sample-2,GCF-2099-002,TGCA\n")
    workbook = load_workbook(tmp_path / "Sample-Submission-Form.xlsx")
    workbook["Sample-Submission-Form"].append(["sample-2", "group"])
    workbook["INFO (GCF-lab only)"].append(["sample-2"])
    workbook.save(tmp_path / "Sample-Submission-Form.xlsx")
    set_initial(sheet, "true")
    selection = preflight.select_run_inputs(tmp_path, tmp_path)
    assert not selection.multi_pair
    result = preflight.validate_selection(selection)
    assert result.ok, result.render_text()
    entries = result.summary["project_initial_holds"]
    assert len(entries) == 2
    assert all(entry["held"] for entry in entries.values())
    assert len({entry["sample_sheet_sha256"] for entry in entries.values()}) == 1


def test_selected_refresh_uses_each_effective_hold(tmp_path):
    source = tmp_path / "source"
    output = tmp_path / "output"
    paired_inputs(source, ["2", "0"])
    paired_inputs(output, ["2", "0"])
    set_initial(source / f"SampleSheet_{A}.csv", "true")
    set_initial(output / f"SampleSheet_{B}.csv", "true")
    result = preflight.validate_selection(
        preflight.select_run_inputs(source, output, refresh=True, projects=[B])
    )
    assert result.ok, result.render_text()
    assert not any(entry["held"] for entry in result.summary["project_initial_holds"].values())


def test_initial_hold_shared_conversion_wait_release_and_retention(tmp_path, monkeypatch):
    cfg, store, events = configure_execution(tmp_path, monkeypatch, ["2", "0"])
    set_initial(cfg.output_path / f"SampleSheet_{B}.csv", "true")
    store.queue(RUN, "demultiplexing")
    converted = []

    def conversion():
        state = store.read(RUN)
        assert state["project_initial_holds"][B]["held"]
        assert project_holds.held(state, B)
        assert not project_holds.held(state, A)
        converted.extend(pair.project_id for pair in cfg.run.input_selection.pairs)
        return "synthetic-converter", "test"

    monkeypatch.setattr(makeFastq, "bcl2fq", conversion)
    monkeypatch.setattr(makeFastq, "rename_fastqs", lambda: None)
    monkeypatch.setattr(sequencing_delivery, "record_conversion", lambda *a, **k: None)
    monkeypatch.setattr(sequencing_delivery, "ensure_report", lambda *a, **k: None)
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting", waiting.get("last_error")
    assert converted == [A, B]
    assert events == [("analysis", A), ("archiving", A)]
    assert waiting["project_initial_holds"][A]["held"] is False
    assert waiting["project_initial_holds"][A]["source"] == "default"
    assert waiting["project_controls"][B]["source"] == "samplesheet"
    assert {n["project"] for n in waiting["delivery_notifications"] if n["status"] == "sent"} == {A}
    snapshot = {
        k: copy.deepcopy(waiting[k])
        for k in (
            "project_completion",
            "analysis_snapshots",
            "project_initial_holds",
            "delivery_notifications",
        )
    }
    paths = [
        p
        for p in Path(waiting["output_path"]).rglob("*")
        if p.is_file() and (A in str(p) or p.name.startswith("md5sum"))
    ]
    preserved = fingerprint(paths)
    restart_context(cfg, waiting)
    # Original input edits cannot change a frozen continuation or its seed.
    (Path(waiting["output_path"]) / f"SampleSheet_{B}.csv").write_text("broken metadata\n")
    control("release", [B])
    events.clear()
    completed = run(cfg, FlowcellStateStore(cfg.static.paths.manager_dir))
    assert completed["status"] == "completed", completed.get("last_error")
    assert completed["attempt"] == waiting["attempt"]
    assert events == [("analysis", B), ("archiving", B)]
    for key in ("project_completion", "analysis_snapshots"):
        assert completed[key][A] == snapshot[key][A]
    assert completed["delivery_notifications"][:2] == snapshot["delivery_notifications"]
    assert completed["project_initial_holds"] == snapshot["project_initial_holds"]
    assert fingerprint(paths) == preserved


def test_all_initially_held_and_manual_precedence_across_scoped_rerun(tmp_path, monkeypatch):
    cfg, store, events = configure_execution(tmp_path, monkeypatch, ["2", "0"])
    for project in (A, B):
        set_initial(cfg.output_path / f"SampleSheet_{project}.csv", "true")
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting", waiting.get("last_error")
    assert events == []
    initial = copy.deepcopy(waiting["project_initial_holds"])
    restart_context(cfg, waiting)
    control("release", [A, B])
    completed = run(cfg, store)
    assert completed["status"] == "completed"
    restart_context(cfg, completed)
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[B], force=True)
    events.clear()
    rerun = run(cfg, store)
    assert rerun["status"] == "completed", rerun.get("last_error")
    assert events == [("analysis", B), ("archiving", B)]
    assert rerun["project_initial_holds"] == initial
    assert not project_holds.held(rerun, B)


@pytest.mark.parametrize("manual_hold", [True, False])
def test_manager_before_seed_wins_and_new_sources_do_not_reseed(tmp_path, monkeypatch, manual_hold):
    cfg, store, _events = configure_execution(tmp_path, monkeypatch, ["2", "0"])
    set_initial(cfg.output_path / f"SampleSheet_{B}.csv", str(not manual_hold))
    control("hold" if manual_hold else "release", [B])
    store.begin_attempt(RUN, cfg=cfg)
    selection = preflight.select_run_inputs(cfg.output_path, cfg.output_path)
    result = preflight.validate_selection(selection)
    initial_holds.initialize(store, RUN, result)
    state = store.read(RUN)
    assert project_holds.held(state, B) is manual_hold
    assert state["project_initial_holds"][B]["decision"] == "existing_control"
    initial = copy.deepcopy(state["project_initial_holds"])
    control("release" if manual_hold else "hold", [B])
    initial_holds.initialize(store, RUN, result)
    assert store.read(RUN)["project_initial_holds"] == initial
    assert project_holds.held(store.read(RUN), B) is not manual_hold


@pytest.mark.parametrize("origin,started", [("new", True), ("legacy_rerun", False)])
def test_older_active_and_legacy_work_is_not_retroactively_held(
    tmp_path, monkeypatch, origin, started
):
    cfg, store, events = configure_execution(tmp_path, monkeypatch, ["2", "0"])
    for project in (A, B):
        set_initial(cfg.output_path / f"SampleSheet_{project}.csv", "true")
    state = store.read(RUN)
    state.pop("project_initial_holds", None)  # State written before this optional extension.
    state["origin"] = origin
    if started:
        state["started_at"] = state["created_at"]
    store.write(state)
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert events == [(s, p) for s in ("analysis", "archiving") for p in (B, A)]
    assert all(
        record["decision"] == "historical" for record in completed["project_initial_holds"].values()
    )


def test_default_false_is_durable_through_source_edit_and_whole_rerun(tmp_path, monkeypatch):
    cfg, store, events = configure_execution(tmp_path, monkeypatch, ["2", "0"])
    first = run(cfg, store)
    assert first["status"] == "completed"
    assert events == [(s, p) for s in ("analysis", "archiving") for p in (B, A)]
    initial = copy.deepcopy(first["project_initial_holds"])
    assert all(not record["held"] for record in initial.values())
    restart_context(cfg, first)
    for project in (A, B):
        set_initial(Path(first["output_path"]) / f"SampleSheet_{project}.csv", "true")
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", force=True)
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert completed["project_initial_holds"] == initial
    assert not any(project_holds.held(completed, p) for p in (A, B))


def test_old_completed_state_preserves_eligibility_across_explicit_rerun(tmp_path, monkeypatch):
    cfg, store, _events = configure_execution(tmp_path, monkeypatch, ["2", "0"])
    first = run(cfg, store)
    first.pop("project_initial_holds")
    store.write(first)
    restart_context(cfg, first)
    set_initial(Path(first["output_path"]) / f"SampleSheet_{B}.csv", "true")
    manager.rerun_flowcell(flowcell=RUN, from_stage="analysis", projects=[B], force=True)
    completed = run(cfg, store)
    assert completed["status"] == "completed", completed.get("last_error")
    assert completed["project_initial_holds"][B]["decision"] == "historical"


def test_seed_serializes_with_claim_and_manager_requests(tmp_path, monkeypatch):
    cfg, store, _events = configure_execution(tmp_path, monkeypatch, ["2", "0"])
    set_initial(cfg.output_path / f"SampleSheet_{B}.csv", "true")
    store.begin_attempt(RUN, cfg=cfg)
    result = preflight.validate_selection(
        preflight.select_run_inputs(cfg.output_path, cfg.output_path)
    )
    for action in (
        lambda: project_holds.admit(store, RUN),
        lambda: project_holds.claim(store, RUN, B, "analysis"),
    ):
        with pytest.raises(StateConflictError, match="initialized before"):
            action()
    reached = Event()

    def seed():
        reached.set()
        initial_holds.initialize(store, RUN, result)

    with ThreadPoolExecutor(max_workers=1) as pool:
        with store.lock(RUN):
            pending = pool.submit(seed)
            assert reached.wait(5)
            assert not pending.done()
        pending.result(timeout=5)
    assert not project_holds.claim(store, RUN, B, "analysis")
    control("release", [B])
    assert project_holds.claim(store, RUN, B, "analysis")
    initial_holds.initialize(store, RUN, result)
    assert not project_holds.held(store.read(RUN), B)
    assert (
        len(
            [
                e
                for e in store.read(RUN)["project_control_history"]
                if e["project"] == B and e["action"] == "initial_eligibility"
            ]
        )
        == 1
    )


@pytest.mark.parametrize("name", ["fm", "flowcell-manager"])
def test_installed_previews_status_and_history_agree(tmp_path, monkeypatch, name):
    cfg, store, _events = configure_execution(tmp_path, monkeypatch, ["2", "0"])
    set_initial(cfg.output_path / f"SampleSheet_{B}.csv", "true")
    waiting = run(cfg, store)
    assert waiting["status"] == "waiting"
    config = tmp_path / "manager.ini"
    config.write_text(
        "[Paths]\n"
        + "\n".join(
            f"{k}={v}"
            for k, v in {
                "nova_baseDir": cfg.static.paths.nova_base_dir,
                "ekista_baseDir": cfg.static.paths.ekista_base_dir,
                "outputDir": cfg.static.paths.output_dir,
                "logDir": cfg.static.paths.log_dir,
                "reportDir": cfg.static.paths.report_dir,
                "analysisDir": cfg.static.paths.analysis_dir,
                "manager_dir": cfg.static.paths.manager_dir,
            }.items()
        )
        + "\n"
    )
    driver = """
import runpy, sys
from bcl2fastq_pipeline.config import PipelineConfig
PipelineConfig.load(sys.argv[1])
sys.argv = sys.argv[2:]
runpy.run_path(sys.argv[0], run_name="__main__")
"""

    def invoke(*args):
        completed = subprocess.run(
            [
                sys.executable,
                "-c",
                driver,
                str(config),
                str(Path(sys.executable).parent / name),
                *args,
            ],
            capture_output=True,
            text=True,
            check=False,
            cwd=tmp_path,
        )
        assert completed.returncode == 0, completed.stderr
        return completed.stdout

    before = store.read(RUN)
    assert "InitialHold=true (samplesheet); held (persisted)" in invoke("validate", RUN)
    assert "InitialHold=true (samplesheet)" in invoke("status", RUN)
    assert "initial_eligibility" in invoke("show", RUN)
    assert "InitialHold=true (samplesheet); held (persisted)" in invoke(
        "rerun", RUN, "--from", "analysis", "--project", B, "--dry-run"
    )
    assert store.read(RUN) == before
    invoke("release", RUN, "--project", B, "--reason", "Reviewed")
    assert "InitialHold=true (samplesheet); released (persisted)" in invoke("validate", RUN)
    assert "manager" in invoke("show", RUN)
    malformed = store.read(RUN)
    malformed["project_initial_holds"][B]["held"] = "true"
    with pytest.raises(StateValidationError, match="initial project hold"):
        validate_state(json.loads(json.dumps(malformed)))
