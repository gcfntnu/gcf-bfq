"""Index corrections across real CLI selection, preflight and daemon preparation."""

import logging
import shutil
import sys

import flowcell_manager.flowcell_manager as manager
import pytest

from bcl2fastq_pipeline.state import (
    ExecutionLeaseError,
    FlowcellStateStore,
    StateConflictError,
)
from test_input_preflight import simulated_execution, snapshot
from test_state_integration import RUN_ID, completed_state, configured_bfq, write_inputs

from bcl2fastq_pipeline import cli, index_corrections, makeFastq, preflight

ORIGINAL_ROW = b"sample,GCF-2026-001,ACGA,TAGC"


def indexed_inputs(directory, suffix=""):
    write_inputs(directory, suffix)
    path = directory / "SampleSheet.csv"
    prefix = path.read_bytes().split(b"[Data]", 1)[0]
    path.write_bytes(
        prefix + b"[Settings]\nReverseComplementIndexP5,0\nReverseComplementIndexP7,1\n"
        b"[Data]\nSample_ID,Sample_Project,index,index2\n" + ORIGINAL_ROW + b"\n"
    )
    return path.read_bytes()


def prepared_run(tmp_path):
    cfg, source, output = configured_bfq(tmp_path)
    original = indexed_inputs(source)
    indexed_inputs(output)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(completed_state(cfg, source, output))
    return cfg, store, source, output, original


def invoke(monkeypatch, *arguments):
    monkeypatch.setattr(sys, "argv", ["flowcell-manager", *map(str, arguments)])
    manager.main()


@pytest.mark.parametrize(
    ("flags", "row"),
    [
        (["--reverse-complement-index1"], b"sample,GCF-2026-001,TCGT,TAGC"),
        (["--reverse-complement-index2"], b"sample,GCF-2026-001,ACGA,GCTA"),
        (["--tom-mode"], b"sample,GCF-2026-001,ACGA,GCTA"),
        (["--tom-mode", "--reverse-complement-index2"], b"sample,GCF-2026-001,ACGA,GCTA"),
        (
            ["--reverse-complement-index1", "--reverse-complement-index2"],
            b"sample,GCF-2026-001,TCGT,GCTA",
        ),
        (["--reverse-complement-index1", "--tom-mode"], b"sample,GCF-2026-001,TCGT,GCTA"),
    ],
)
def test_public_cli_options_apply_once_and_preserve_instrument_and_unselected_bytes(
    tmp_path, monkeypatch, flags, row
):
    _cfg, store, source, output, original = prepared_run(tmp_path)
    source_before = snapshot(source)
    (output / "obsolete.7za").write_bytes(b"old archive")

    invoke(monkeypatch, "rerun", RUN_ID, "--force", *flags)

    assert (output / "SampleSheet.csv").read_bytes() == original.replace(ORIGINAL_ROW, row)
    assert snapshot(source) == source_before
    assert store.read(RUN_ID)["status"] == "queued"
    assert not (output / "obsolete.7za").exists()


def test_each_explicit_toggle_uses_current_sheet_and_plain_rerun_preserves_it(tmp_path):
    _cfg, store, source, output, original = prepared_run(tmp_path)
    manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    reversed_sheet = (output / "SampleSheet.csv").read_bytes()
    manager.rerun_flowcell(flowcell=RUN_ID, force=True)
    assert (output / "SampleSheet.csv").read_bytes() == reversed_sheet
    manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    assert (output / "SampleSheet.csv").read_bytes() == original
    manager.rerun_flowcell(flowcell=RUN_ID, force=True, reverse_complement_index1=True)
    manager.rerun_flowcell(flowcell=RUN_ID, force=True, reverse_complement_index2=True)
    assert (output / "SampleSheet.csv").read_bytes() == original.replace(
        ORIGINAL_ROW, b"sample,GCF-2026-001,TCGT,GCTA"
    )
    assert (source / "SampleSheet.csv").read_bytes() == original
    assert store.read(RUN_ID)["status"] == "queued"


def test_explicit_refresh_then_toggle_uses_source_pair(tmp_path):
    _cfg, _store, source, output, _original = prepared_run(tmp_path)
    source_original = indexed_inputs(source, "-instrument")
    indexed_inputs(output, "-curated")
    manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    manager.rerun_flowcell(flowcell=RUN_ID, force=True, refresh_inputs=True, tom_mode=True)
    assert (output / "SampleSheet.csv").read_bytes() == source_original.replace(b"TAGC", b"GCTA")
    assert (output / "Sample-Submission-Form.xlsx").read_bytes() == (
        source / "Sample-Submission-Form.xlsx"
    ).read_bytes()


@pytest.mark.parametrize("decision", ["dry-run", "decline"])
@pytest.mark.parametrize("running", [False, True])
def test_correction_preview_never_mutates_files_state_or_lease(
    tmp_path, monkeypatch, decision, running
):
    _cfg, store, _source, output, _original = prepared_run(tmp_path)
    if running:
        store.queue(RUN_ID, "demultiplexing")
        store.begin_attempt(RUN_ID)
    (output / "obsolete.7za").write_bytes(b"old archive")
    monkeypatch.setattr("builtins.input", lambda _prompt: "no")
    before = snapshot(tmp_path)
    state_before = store.read(RUN_ID)
    manager.rerun_flowcell(flowcell=RUN_ID, tom_mode=True, dry_run=decision == "dry-run")
    assert store.read(RUN_ID) == state_before
    assert snapshot(tmp_path) == before


@pytest.mark.parametrize("command", ["rerun", "initialize"])
@pytest.mark.parametrize("boundary", ["analysis", "reporting", "finalization"])
def test_invalid_boundary_rejected_before_any_recovery_or_input_access(
    tmp_path, monkeypatch, command, boundary
):
    _cfg, store, _source, _output, _original = prepared_run(tmp_path)
    store.queue(RUN_ID, "demultiplexing")
    store.begin_attempt(RUN_ID)
    monkeypatch.setattr(
        manager, "get_cfg", lambda: pytest.fail("configuration loaded before guard")
    )
    before = snapshot(tmp_path)
    with pytest.raises(StateConflictError, match="require --from demultiplexing"):
        getattr(manager, f"{command}_flowcell")(
            flowcell=RUN_ID, from_stage=boundary, force=True, tom_mode=True
        )
    assert snapshot(tmp_path) == before


def test_missing_source_does_not_prevent_toggle_of_curated_output(tmp_path, capsys):
    _cfg, store, source, output, original = prepared_run(tmp_path)
    shutil.rmtree(source)
    manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    assert (output / "SampleSheet.csv").read_bytes() == original.replace(b"TAGC", b"GCTA")
    assert store.read(RUN_ID)["status"] == "queued"
    assert "unavailable" in capsys.readouterr().out.lower()


def test_initialize_explicit_path_selects_exact_source_and_defaults_to_demultiplexing(
    tmp_path, monkeypatch
):
    cfg, source, output = configured_bfq(tmp_path)
    indexed_inputs(source, "-first-root")
    exact = cfg.static.paths.ekista_base_dir / RUN_ID
    original = indexed_inputs(exact, "-chosen-root")
    source_before = snapshot(exact)
    invoke(monkeypatch, "initialize", exact, "--tom-mode", "--force")
    state = FlowcellStateStore(cfg.static.paths.manager_dir).read(RUN_ID)
    assert state["source_path"] == str(exact)
    assert state["current_stage"] == "demultiplexing"
    assert state["status"] == "queued"
    assert cli.candidate_flowcells(cfg, FlowcellStateStore(cfg.static.paths.manager_dir)) == [exact]
    assert (output / "SampleSheet.csv").read_bytes() == original.replace(b"TAGC", b"GCTA")
    assert snapshot(exact) == source_before


def test_initialize_bare_id_errors_if_roots_are_ambiguous_without_mutation(tmp_path):
    cfg, source, _output = configured_bfq(tmp_path)
    indexed_inputs(source)
    indexed_inputs(cfg.static.paths.ekista_base_dir / RUN_ID)
    before = snapshot(tmp_path)
    with pytest.raises(StateConflictError, match="Ambiguous input run"):
        manager.initialize_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    assert snapshot(tmp_path) == before


@pytest.mark.parametrize("decision", ["dry-run", "decline"])
def test_initialize_preview_is_read_only_and_displays_correction(
    tmp_path, monkeypatch, capsys, decision
):
    cfg, source, output = configured_bfq(tmp_path)
    indexed_inputs(source)
    cfg.static.paths.manager_dir.rmdir()
    monkeypatch.setattr("builtins.input", lambda _prompt: "no")
    before = snapshot(tmp_path)
    manager.initialize_flowcell(flowcell=RUN_ID, tom_mode=True, dry_run=decision == "dry-run")
    assert snapshot(tmp_path) == before
    assert not cfg.static.paths.manager_dir.exists()
    assert not output.exists()
    assert "index2" in capsys.readouterr().out


def test_initialize_preserves_curated_output_precedence_and_rejects_existing_state(tmp_path):
    cfg, source, output = configured_bfq(tmp_path)
    indexed_inputs(source, "-instrument")
    original = indexed_inputs(output, "-curated")
    manager.initialize_flowcell(flowcell=RUN_ID, force=True, reverse_complement_index1=True)
    assert (output / "SampleSheet.csv").read_bytes() == original.replace(b"ACGA", b"TCGT")
    before = snapshot(tmp_path)
    with pytest.raises(StateConflictError, match="use rerun instead"):
        manager.initialize_flowcell(flowcell=source, force=True, tom_mode=True)
    assert snapshot(tmp_path) == before
    assert FlowcellStateStore(cfg.static.paths.manager_dir).exists(RUN_ID)


@pytest.mark.parametrize("source_kind", ["same-directory", "symlink", "missing", "file"])
def test_initialize_rejects_unusable_or_output_source_before_mutating(tmp_path, source_kind):
    _cfg, source, output = configured_bfq(tmp_path)
    indexed_inputs(output)
    if source_kind == "same-directory":
        selected = output
    elif source_kind == "symlink":
        source.rmdir()
        source.symlink_to(output, target_is_directory=True)
        selected = source
    elif source_kind == "missing":
        selected = tmp_path / "missing" / RUN_ID
    else:
        source.rmdir()
        source.write_bytes(b"not a directory")
        selected = source
    before = snapshot(tmp_path)
    with pytest.raises(StateConflictError):
        manager.initialize_flowcell(flowcell=selected, force=True, tom_mode=True)
    assert snapshot(tmp_path) == before


def test_daemon_demultiplexer_receives_corrected_output_and_never_reapplies_toggle(
    tmp_path, monkeypatch
):
    cfg, store, source, output, original = prepared_run(tmp_path)
    manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    calls = simulated_execution(cfg, monkeypatch)
    observed = []

    def demultiplex():
        observed.append(cfg.run.sample_sheet.read_bytes())
        return "test", "1"

    monkeypatch.setattr(makeFastq, "bcl2fq", demultiplex)
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("index-cli-test"), prepare=True)
    assert observed == [original.replace(b"TAGC", b"GCTA")]
    assert "finalize" in calls
    assert store.read(RUN_ID)["status"] == "completed"
    assert (source / "SampleSheet.csv").read_bytes() == original
    assert (output / "SampleSheet.csv").read_bytes() == observed[0]


def test_show_reports_live_orientation_without_persisting_derived_field(tmp_path, capsys):
    _cfg, store, _source, _output, _original = prepared_run(tmp_path)
    manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    result = manager.show_flowcell(flowcell=RUN_ID)
    assert result["index_orientation"] == index_corrections.current_orientation(store.read(RUN_ID))
    assert "index_orientation" not in store.read(RUN_ID)
    assert "index_orientation" in capsys.readouterr().out


@pytest.mark.parametrize("mutated", ["effective", "source"])
def test_changed_sheet_after_confirmation_blocks_cleanup_and_toggle(tmp_path, monkeypatch, mutated):
    _cfg, store, source, output, _original = prepared_run(tmp_path)
    (output / "obsolete.7za").write_bytes(b"old archive")
    before_state = store.read(RUN_ID)

    def confirm(_prompt, _force):
        sheet = (source if mutated == "source" else output) / "SampleSheet.csv"
        sheet.write_bytes(sheet.read_bytes().replace(b"ACGA", b"CCCC"))
        return True

    monkeypatch.setattr(manager, "_confirm", confirm)
    with pytest.raises(StateConflictError):
        manager.rerun_flowcell(flowcell=RUN_ID, tom_mode=True)
    assert store.read(RUN_ID) == before_state
    assert (output / "obsolete.7za").exists()
    assert b",TAGC\n" in (output / "SampleSheet.csv").read_bytes()


def test_invalid_metadata_preflight_blocks_correction_and_cleanup(tmp_path):
    _cfg, store, _source, output, _original = prepared_run(tmp_path)
    (output / "Sample-Submission-Form.xlsx").write_bytes(b"invalid workbook")
    before = snapshot(tmp_path)
    before_state = store.read(RUN_ID)
    with pytest.raises(preflight.PreflightValidationError):
        manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    assert store.read(RUN_ID) == before_state
    assert snapshot(tmp_path) == before


def test_active_lease_blocks_even_forced_index_correction(tmp_path):
    _cfg, store, _source, output, original = prepared_run(tmp_path)
    with store.execution_lease(RUN_ID):
        with pytest.raises(ExecutionLeaseError):
            manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    assert (output / "SampleSheet.csv").read_bytes() == original
    assert store.read(RUN_ID)["status"] == "completed"


@pytest.mark.parametrize("command", ["initialize", "rerun"])
def test_help_documents_canonical_options_but_hides_alias(monkeypatch, capsys, command):
    with pytest.raises(SystemExit) as error:
        invoke(monkeypatch, command, "--help")
    assert error.value.code == 0
    text = capsys.readouterr().out
    assert "--reverse-complement-index1" in text
    assert "--reverse-complement-index2" in text
    assert "--tom-mode" not in text


@pytest.mark.parametrize("replaced", [False, True])
def test_plain_rerun_recovers_interrupted_toggle_without_reapplying(
    tmp_path, monkeypatch, replaced
):
    _cfg, store, _source, output, original = prepared_run(tmp_path)
    real_replace = index_corrections._replace_sheet

    def interrupt_replacement(target, content):
        if replaced:
            real_replace(target, content)
        raise OSError("simulated interruption")

    monkeypatch.setattr(index_corrections, "_replace_sheet", interrupt_replacement)
    with pytest.raises(OSError, match="simulated interruption"):
        manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    pending = store.read(RUN_ID)
    assert pending["status"] == "preparing"
    assert pending["index_corrections"][-1]["status"] == "pending"
    current = (output / "SampleSheet.csv").read_bytes()
    assert current == (original.replace(b"TAGC", b"GCTA") if replaced else original)

    monkeypatch.setattr(index_corrections, "_replace_sheet", real_replace)
    manager.rerun_flowcell(flowcell=RUN_ID, force=True)
    assert (output / "SampleSheet.csv").read_bytes() == current
    resumed = store.read(RUN_ID)
    assert resumed["status"] == "queued"
    assert resumed["index_corrections"][-1]["status"] == ("applied" if replaced else "not_applied")

    manager.rerun_flowcell(flowcell=RUN_ID, force=True, tom_mode=True)
    assert (output / "SampleSheet.csv").read_bytes() == (
        original if replaced else original.replace(b"TAGC", b"GCTA")
    )


def test_initialize_exact_relative_input_and_normalized_dot_path(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    original = indexed_inputs(source)
    monkeypatch.chdir(tmp_path)
    manager.initialize_flowcell(flowcell=f"./nova/{RUN_ID}/../{RUN_ID}", force=True, tom_mode=True)
    assert FlowcellStateStore(cfg.static.paths.manager_dir).read(RUN_ID)["source_path"] == str(
        source
    )
    assert (output / "SampleSheet.csv").read_bytes() == original.replace(b"TAGC", b"GCTA")
