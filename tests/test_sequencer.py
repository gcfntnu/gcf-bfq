"""Known sequencers are mandatory at BFQ's processing boundaries."""

import hashlib
import logging
import sys

from types import SimpleNamespace
from unittest.mock import Mock

import flowcell_manager.flowcell_manager as manager
import pytest

from bcl2fastq_pipeline.sequencer import SequencerError, require_sequencer
from bcl2fastq_pipeline.state import FlowcellStateStore, new_state
from configmaker.configmaker import SEQUENCERS, find_machine
from support.fixtures import completed_state, configured_bfq, write_fastq, write_inputs
from test_input_preflight import simulated_execution, snapshot

from bcl2fastq_pipeline import afterFastq, cli, preflight

KNOWN_RUN = "261101_TESTSEQ_0001_TESTFLOWCELL"
UNKNOWN_RUN = "261101_UNREGISTERED_0001_TESTFLOWCELL"
INVALID_RUNS = [UNKNOWN_RUN, "malformed", "261101__0001_TESTFLOWCELL"]


@pytest.fixture
def registry(monkeypatch):
    monkeypatch.setitem(SEQUENCERS, "TESTSEQ", "NovaSeq 6000")


@pytest.mark.parametrize("run_id", [KNOWN_RUN, UNKNOWN_RUN])
def test_lookup_matches_configmaker(registry, run_id, caplog):
    assert afterFastq.get_sequencer(run_id) == find_machine([run_id])
    if run_id == UNKNOWN_RUN:
        assert "UNREGISTERED" in caplog.text
        assert "SEQUENCERS" in caplog.text
        assert run_id in caplog.text
    else:
        assert require_sequencer(run_id) == "NovaSeq 6000"
        assert not caplog.records


@pytest.mark.parametrize("run_id", ["malformed", "261101__0001_TESTFLOWCELL", "", None])
def test_malformed_identity_has_clear_error_and_empty_lookup(run_id, caplog):
    with pytest.raises(SequencerError, match="malformed run ID"):
        require_sequencer(run_id)
    assert afterFastq.get_sequencer(run_id) == ""
    assert "second underscore-separated field" in caplog.text


@pytest.mark.parametrize("run_id", INVALID_RUNS)
@pytest.mark.parametrize("command", ["fm", "flowcell-manager"])
def test_manual_validation_fails_read_only(tmp_path, monkeypatch, capsys, run_id, command):
    cfg, source, _output = configured_bfq(tmp_path, run_id=run_id)
    write_inputs(source)
    cfg.static.paths.manager_dir.rmdir()
    before = snapshot(tmp_path)
    monkeypatch.setattr(sys, "argv", [command, "validate", run_id])
    with pytest.raises(SystemExit) as error:
        manager.main()
    assert error.value.code == 1
    assert "Cannot determine sequencer" in capsys.readouterr().err
    assert snapshot(tmp_path) == before
    assert not cfg.static.paths.manager_dir.exists()


@pytest.mark.parametrize("boundary", ["demultiplexing", "analysis"])
@pytest.mark.parametrize("run_id", INVALID_RUNS)
def test_execution_records_failure_before_processing(tmp_path, monkeypatch, boundary, run_id):
    cfg, source, output = configured_bfq(tmp_path, run_id=run_id)
    write_inputs(output)
    cfg.run.sample_sheet = output / "SampleSheet.csv"
    cfg.run.sample_submission_form = output / "Sample-Submission-Form.xlsx"
    fastq = output / "GCF-2026-001/sample_R1.fastq.gz"
    write_fastq(fastq)
    checksum = hashlib.md5(fastq.read_bytes()).hexdigest()
    (output / "md5sum_GCF-2026-001_fastq.txt").write_text(
        f"{checksum}  {fastq.relative_to(output)}\n"
    )
    (output / "provenance").mkdir()
    (output / "provenance/successful.tar").write_bytes(b"prior successful snapshot")
    before = snapshot(output)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(new_state(run_id, source, output, origin="new", start_stage=boundary, cfg=cfg))
    calls = simulated_execution(cfg, monkeypatch)
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("sequencer-test"))
    state = store.read(run_id)
    assert calls == []
    assert state["status"] == "failed"
    assert state["current_stage"] == boundary
    assert run_id in state["last_error"]["summary"]
    report = state["stages"][boundary]["metadata"]["input_preflight"]
    assert any(error["code"] == "run.sequencer" for error in report["errors"])
    assert snapshot(output) == before


@pytest.mark.parametrize("boundary", ["demultiplexing", "analysis"])
@pytest.mark.parametrize("refresh", [False, True])
@pytest.mark.parametrize("projects", [None, ["GCF-2026-001"]])
def test_unknown_instrument_blocks_rerun_before_cleanup(  # noqa: PLR0913
    tmp_path, monkeypatch, boundary, refresh, projects
):
    cfg, source, output = configured_bfq(tmp_path)
    monkeypatch.delitem(SEQUENCERS, cfg.run.run_id.split("_")[1])
    write_inputs(source)
    write_inputs(output)
    write_fastq(output / "GCF-2026-001/sample_R1.fastq.gz")
    (output / "multiqc_GCF-2026-001.html").write_text("prior analysis")
    (output / "GCF-2026-001.7za").write_bytes(b"prior delivery")
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(completed_state(cfg, source, output))
    before = snapshot(output), snapshot(source), store.state_path(cfg.run.run_id).read_bytes()
    with pytest.raises(preflight.PreflightValidationError, match="SEQUENCERS"):
        manager.rerun_flowcell(
            flowcell=cfg.run.run_id,
            from_stage=boundary,
            force=True,
            refresh_inputs=refresh,
            projects=projects,
        )
    assert (
        snapshot(output),
        snapshot(source),
        store.state_path(cfg.run.run_id).read_bytes(),
    ) == before


def test_unknown_instrument_blocks_initialization(tmp_path):
    cfg, source, output = configured_bfq(tmp_path, run_id=UNKNOWN_RUN)
    write_inputs(source)
    write_inputs(output)
    (output / "existing.7za").write_bytes(b"prior delivery")
    before = snapshot(tmp_path)
    with pytest.raises(preflight.PreflightValidationError, match="SEQUENCERS"):
        manager.initialize_flowcell(flowcell=str(source), from_stage="demultiplexing", force=True)
    assert snapshot(tmp_path) == before
    assert not FlowcellStateStore(cfg.static.paths.manager_dir).exists(UNKNOWN_RUN)


def test_preflight_uses_state_identity_for_staged_inputs(tmp_path, registry):
    inputs = tmp_path / "inputs"
    write_inputs(inputs)
    state = new_state(
        KNOWN_RUN,
        tmp_path / KNOWN_RUN,
        tmp_path / "output" / KNOWN_RUN,
        origin="new",
        start_stage="analysis",
    )
    _selection, result = preflight.require_valid_inputs(inputs, inputs, state=state)
    assert result.ok


@pytest.mark.parametrize("boundary", ["reporting", "finalization"])
@pytest.mark.parametrize("refresh", [False, True])
def test_delivery_only_retry_does_not_require_registered_instrument(
    tmp_path, monkeypatch, boundary, refresh
):
    cfg, source, output = configured_bfq(tmp_path)
    write_inputs(source)
    write_inputs(output)
    monkeypatch.delitem(SEQUENCERS, cfg.run.run_id.split("_")[1])
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(completed_state(cfg, source, output))
    state = manager.rerun_flowcell(
        flowcell=cfg.run.run_id, from_stage=boundary, force=True, refresh_inputs=refresh
    )
    assert state["status"] == "queued"
    assert state["stages"]["analysis"]["status"] == "completed"


@pytest.mark.parametrize("resume", [None, {"projects": {}}])
def test_direct_analysis_refuses_unknown_machine_before_side_effects(monkeypatch, resume):
    cfg = SimpleNamespace(run=SimpleNamespace(run_id=UNKNOWN_RUN))
    prepare = Mock()
    command = Mock()
    monkeypatch.setattr(afterFastq.workflow_config, "prepare_analysis_projects", prepare)
    monkeypatch.setattr(afterFastq.subprocess, "check_call", command)
    monkeypatch.setattr(afterFastq, "run_logged_command", command)
    with pytest.raises(SequencerError, match="UNREGISTERED"):
        afterFastq.full_align(cfg, resume=resume)
    prepare.assert_not_called()
    command.assert_not_called()
