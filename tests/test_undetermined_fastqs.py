"""Synthetic converter trees; no converter, workflow or mail is executed."""

import logging
import shutil
import subprocess

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import flowcell_manager.flowcell_manager as manager
import pytest

from bcl2fastq_pipeline.state import StateConflictError, find_fastqs, validate_restored_fastqs
from support.fixtures import write_fastq
from test_early_qc_integration import RUN_ID, prepared, run

from bcl2fastq_pipeline import (
    afterFastq,
    fastq_cleanup,
    makeFastq,
    notifications,
    sequencing_qc,
    undetermined_fastqs,
)

ROOTS = [
    "",
    "nested/bcl-convert",
    "outs/fastq_path",
    "FLOWCELL/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ/fork0/chnk0/files/fastq_path",
]
NAMES = [
    f"Undetermined_S0{lane}_{read}_{chunk}.fastq.gz"
    for lane in ("", "_L001", "_L008")
    for read in ("R1", "R2", "R3", "I1", "I2")
    for chunk in ("001", "002")
]
SAMPLE_NAMES = [
    "Undetermined_S1_L001_R1_001.fastq.gz",
    "Undetermined_control_S2_L001_I1_001.fastq.gz",
    "Undetermined_S0_control_R1.fastq.gz",
    "Undetermined_R1.fastq.gz",
    "Undetermined_S0_S1_R1_001.fastq.gz",
    "sample_Undetermined_S0_R1_001.fastq.gz",
    "Undetermined_S01_R1_001.fastq.gz",
]


@pytest.mark.parametrize("root", ROOTS)
def test_supported_layouts_are_idempotent_and_preserve_other_outputs(tmp_path, caplog, root):
    output = tmp_path / "output"
    for name in NAMES:
        path = output / root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"unassigned reads")
    kept = {}
    for relative in [
        *[f"GCF-2026-001/{name}" for name in SAMPLE_NAMES],
        "GCF-2026-002/sample_R1.fastq.gz",
        "InterOp/IndexMetricsOut.bin",
        "Reports/Demultiplex_Stats.csv",
        "Reports/Top_Unknown_Barcodes.csv",
        "Stats/Stats.json",
        "Stats/DemultiplexingStats.xml",
        "sequencer_stats_synthetic.html",
        "md5sum_GCF-2026-001_fastq.txt",
        "conversion.log",
        "Undetermined_S0_R1_001.fastq.gz.md5",
        "Undetermined_S0_R1_001.fastq.gz.html",
        "completed.7za",
    ]:
        path = output / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"retained synthetic product")
        kept[path] = path.read_bytes(), path.stat().st_mtime_ns
    caplog.set_level(logging.INFO)
    undetermined_fastqs.cleanup(output)
    undetermined_fastqs.cleanup(output)
    assert not any((output / root / name).exists() for name in NAMES)
    assert {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in kept} == kept
    size = len(NAMES) * len(b"unassigned reads")
    assert (
        f"removed {len(NAMES)} paths, {size} regular-file bytes; 0 unlink failures" in caplog.text
    )
    assert "removed 0 paths, 0 regular-file bytes; 0 unlink failures" in caplog.text
    assert all(name in caplog.text for name in NAMES)


def test_symlinks_are_unlinked_without_visiting_external_trees(tmp_path, monkeypatch):
    output, external = tmp_path / "output", tmp_path / "external"
    output.mkdir()
    external.mkdir()
    target = external / NAMES[0]
    target.write_bytes(b"external reads")
    (output / NAMES[0]).symlink_to(target)
    (output / NAMES[1]).symlink_to(external, target_is_directory=True)
    (output / NAMES[2]).symlink_to(external / "absent.fastq.gz")
    (output / "external-directory").symlink_to(external, target_is_directory=True)
    original = undetermined_fastqs.os.scandir

    def scan(path):
        assert path != external and "external-directory" not in str(path)
        return original(path)

    monkeypatch.setattr(undetermined_fastqs.os, "scandir", scan)
    undetermined_fastqs.cleanup(output)
    assert target.read_bytes() == b"external reads"
    assert not any((output / name).is_symlink() for name in NAMES[:3])
    assert (output / "external-directory").is_symlink()


def test_cleanup_rejects_linked_root(tmp_path):
    external = tmp_path / "external"
    external.mkdir()
    target = external / NAMES[0]
    target.write_bytes(b"external reads")
    output = tmp_path / "output"
    output.symlink_to(external, target_is_directory=True)
    with pytest.raises(StateConflictError, match="must not be a symlink"):
        undetermined_fastqs.cleanup(output)
    assert target.read_bytes() == b"external reads"


def test_parent_swapped_for_symlink_after_plan_cannot_delete_external_file(tmp_path, monkeypatch):
    output, external = tmp_path / "output", tmp_path / "external"
    (output / "nested").mkdir(parents=True)
    external.mkdir()
    relative = Path("nested") / NAMES[0]
    (output / relative).write_bytes(b"managed reads")
    (external / NAMES[0]).write_bytes(b"external reads")
    original = fastq_cleanup.unlink_fastq

    def swap(root_fd, name, info):
        (output / "nested").rename(output / "saved")
        (output / "nested").symlink_to(external, target_is_directory=True)
        return original(root_fd, name, info)

    monkeypatch.setattr(fastq_cleanup, "unlink_fastq", swap)
    with pytest.raises(StateConflictError, match="retry from analysis"):
        undetermined_fastqs.cleanup(output)
    assert (external / NAMES[0]).read_bytes() == b"external reads"
    assert (output / "saved" / NAMES[0]).read_bytes() == b"managed reads"


def test_files_missing_between_plan_and_unlink_are_normal(tmp_path, monkeypatch, caplog):
    path = tmp_path / NAMES[0]
    path.write_bytes(b"reads")
    original = fastq_cleanup.unlink_fastq

    def vanished(root_fd, name, info):
        path.unlink()
        return original(root_fd, name, info)

    monkeypatch.setattr(fastq_cleanup, "unlink_fastq", vanished)
    caplog.set_level(logging.INFO)
    undetermined_fastqs.cleanup(tmp_path)
    assert "removed 0 paths, 0 regular-file bytes; 0 unlink failures" in caplog.text


def test_rename_preserves_reserved_identity_and_renames_real_samples(tmp_path, monkeypatch):
    output = tmp_path
    cfg = SimpleNamespace(output_path=output, run=SimpleNamespace(libprep="Illumina DNA Prep"))
    monkeypatch.setattr(makeFastq.PipelineConfig, "get", lambda: cfg)
    for root in ("nested", "GCF-2026-001/sample"):
        for name in [NAMES[0], "Undetermined_S1_L001_R1_001.fastq.gz"]:
            path = output / root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"reads")
    makeFastq.rename_fastqs()
    undetermined_fastqs.cleanup(output)
    for root in ("nested", "GCF-2026-001/sample"):
        assert not (output / root / NAMES[0]).exists()
        assert (output / root / "Undetermined_L001_R1.fastq.gz").read_bytes() == b"reads"


@pytest.mark.parametrize("failure", [PermissionError("unlink denied"), KeyboardInterrupt()])
def test_partial_failure_or_interruption_retries_analysis_without_conversion(
    tmp_path, monkeypatch, failure
):
    cfg, store, output, calls = prepared(tmp_path, monkeypatch)
    convert = makeFastq.bcl2fq
    converter = Mock(side_effect=convert)
    monkeypatch.setattr(makeFastq, "bcl2fq", converter)
    # Sorting ensures one unlink succeeds before the failing second file.
    paths = [output / name for name in sorted(NAMES[:2])]
    for path in paths:
        path.write_bytes(b"unassigned reads")
    original = fastq_cleanup.unlink_fastq

    def fail_second(root_fd, name, info):
        state = store.read(RUN_ID)
        assert store.execution_active(RUN_ID)
        assert state["stages"]["demultiplexing"]["status"] == "completed"
        assert state["sequencing_qc"]["status"] == "completed"
        assert "hash" in calls and "analysis" not in calls
        if name == paths[1].name:
            raise failure
        return original(root_fd, name, info)

    monkeypatch.setattr(fastq_cleanup, "unlink_fastq", fail_second)
    if isinstance(failure, KeyboardInterrupt):
        with pytest.raises(KeyboardInterrupt):
            run(cfg, store)
        first = store.recover_interrupted(RUN_ID)
        assert first["status"] == "interrupted"
    else:
        run(cfg, store)
        first = store.read(RUN_ID)
        assert first["status"] == "failed"
        assert "unlink denied" in first["last_error"]["summary"]
    assert first["current_stage"] == "analysis"
    assert not paths[0].exists() and paths[1].exists()
    qc = first["sequencing_qc"]
    # Exercise existing manager planning, state and execution-lease recovery.
    manager.rerun_flowcell(flowcell=RUN_ID, from_stage="analysis", force=True)
    cfg.run.begin(first["source_path"], cfg.static.paths)
    cfg.run.apply_custom(
        {"Libprep": "Illumina DNA Prep"},
        output / "SampleSheet.csv",
        output / "Sample-Submission-Form.xlsx",
    )
    monkeypatch.setattr(fastq_cleanup, "unlink_fastq", original)
    run(cfg, store)
    assert store.read(RUN_ID)["status"] == "completed"
    assert store.read(RUN_ID)["sequencing_qc"] == qc
    assert not any(path.exists() for path in paths)
    converter.assert_called_once()
    assert calls.count("qc") == calls.count("sequencing-mail") == calls.count("analysis") == 1


@pytest.mark.parametrize("qc_fails", [False, True])
def test_mail_or_report_failure_does_not_block_cleanup(tmp_path, monkeypatch, qc_fails):
    cfg, store, output, _calls = prepared(tmp_path, monkeypatch)
    path = output / NAMES[0]
    path.write_bytes(b"unassigned reads")
    (output / "Stats").mkdir()
    stats = output / "Stats/Stats.json"
    stats.write_text('{"undetermined": 25}')
    monkeypatch.setattr(notifications, "send_notification", Mock(side_effect=OSError("no relay")))
    if qc_fails:
        monkeypatch.setattr(sequencing_qc, "generate", Mock(side_effect=RuntimeError("no report")))

    def analysis():
        assert not path.exists()
        assert stats.read_text() == '{"undetermined": 25}'

    monkeypatch.setattr(afterFastq, "analysis_steps", analysis)
    run(cfg, store)
    assert store.read(RUN_ID)["status"] == "completed"


def test_failed_conversion_retains_diagnostic_reads_and_never_cleans(tmp_path, monkeypatch):
    cfg, store, output, calls = prepared(tmp_path, monkeypatch)
    path = output / NAMES[0]
    path.write_bytes(b"diagnostic reads")
    monkeypatch.setattr(makeFastq, "bcl2fq", Mock(side_effect=RuntimeError("conversion failed")))
    cleaner = Mock(side_effect=AssertionError("must not clean failed conversion"))
    monkeypatch.setattr(undetermined_fastqs, "cleanup", cleaner)
    run(cfg, store)
    cleaner.assert_not_called()
    assert path.read_bytes() == b"diagnostic reads"
    assert store.read(RUN_ID)["stages"]["demultiplexing"]["status"] == "failed"
    assert not calls


def test_restored_validation_and_old_cleanup_records_do_not_require_undetermined(
    tmp_path, monkeypatch
):
    _cfg, _store, output, _calls = prepared(tmp_path, monkeypatch)
    sample = output / "GCF-2026-001" / SAMPLE_NAMES[0]
    write_fastq(sample)
    (output / NAMES[0]).write_bytes(b"empty or damaged unassigned diagnostic")
    assert find_fastqs(output) == [sample]
    assert validate_restored_fastqs(output)[0]
    state = {
        "fastq_cleanup": {
            "required_files": {
                NAMES[0]: 123,
                str(sample.relative_to(output)): sample.stat().st_size,
            }
        }
    }
    (output / NAMES[0]).unlink()
    fastq_cleanup.require_restored(state, output, "analysis")
    sample.unlink()
    with pytest.raises(StateConflictError, match="Restore all cleaned FASTQs"):
        fastq_cleanup.require_restored(state, output, "analysis")


@pytest.mark.skipif(shutil.which("7za") is None, reason="optional real 7-Zip contract check")
def test_real_archiver_excludes_nested_reads_but_keeps_similar_samples(tmp_path):
    # This is a tiny local archive, never a scientific tool or operational input.
    root = tmp_path / "mkfastq"
    for parent in (root, root / "outs/fastq_path", root / "staging/deep/files"):
        parent.mkdir(parents=True, exist_ok=True)
        for name in [*NAMES, *SAMPLE_NAMES, "Stats.json", "unknown-index.csv"]:
            (parent / name).write_bytes(b"synthetic")
    archive = tmp_path / "delivery.7z"
    subprocess.run(
        ["7za", "a", str(archive), *undetermined_fastqs.archive_exclusions(), str(root)],
        check=True,
        capture_output=True,
    )
    listing = subprocess.check_output(["7za", "l", "-slt", str(archive)], text=True)
    members = [
        line.removeprefix("Path = ") for line in listing.splitlines() if line.startswith("Path = ")
    ]
    assert not any(undetermined_fastqs.is_undetermined(name) for name in members)
    for name in [*SAMPLE_NAMES, "Stats.json", "unknown-index.csv"]:
        assert sum(Path(member).name == name for member in members) == 3
