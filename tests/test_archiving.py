import json
import os
import shutil
import subprocess
import sys
import textwrap

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from bcl2fastq_pipeline.state import apply_cleanup, cleanup_plan

from bcl2fastq_pipeline import afterFastq, misc, undetermined_fastqs

PROJECT = "GCF-2026-001"
RUN_ID = "260923_A01990_0001_AFLOWCELL"
RUN_DATE = "260923"
REAL_ARCHIVER = shutil.which("7za")


@pytest.fixture
def archive_run(tmp_path, monkeypatch):
    output = tmp_path / "output with spaces;not-a-command"
    output.mkdir()
    inputs = {
        "project": output / PROJECT,
        "stats": output / "Stats",
        "project_samplesheet": output / f"{PROJECT}_samplesheet.tsv",
        "samplesheet": output / "SampleSheet.csv",
        "submission_form": output / "Sample-Submission-Form.xlsx",
        "fastq_checksum": output / f"md5sum_{PROJECT}_fastq.txt",
    }
    inputs["project"].mkdir()
    inputs["stats"].mkdir()
    for path in list(inputs.values())[2:]:
        path.write_text(f"synthetic input: {path.name}\n")
    fastq = inputs["project"] / "sample_R1.fastq.gz"
    fastq.write_bytes(b"synthetic FASTQ contents")
    report = output / f"sequencer_stats_{PROJECT}_{RUN_DATE}.html"
    report.write_text("<html>completed sequencing report</html>\n")
    analysis_report = output / f"multiqc_{PROJECT}_{RUN_DATE}.html"
    analysis_report.write_text("<html>completed analysis report</html>\n")
    work_root = tmp_path / "work with spaces"
    workdir = work_root / f"{PROJECT}_{RUN_DATE}"
    qc_dir = workdir / "data" / "tmp" / "singlecell" / "bfq"
    qc_dir.mkdir(parents=True)
    qc_result = qc_dir / "analysis.txt"
    qc_result.write_text("completed scientific output\n")
    workflow = workdir / "workflow.smk"
    workflow.write_text("# retained workflow identity\n")
    inputs["qc"] = qc_dir
    cfg = SimpleNamespace(
        output_path=output,
        run=SimpleNamespace(
            run_id=RUN_ID,
            sensitive=False,
            libprep="10X Genomics Chromium 3' Gene Expression",
            pipeline="singlecell",
        ),
        static=SimpleNamespace(
            paths=SimpleNamespace(log_dir=tmp_path / "logs", report_dir=tmp_path / "reports")
        ),
    )

    # Exercise BFQ's actual process/logging boundary without invoking an archiver.
    binary_dir = tmp_path / "bin"
    binary_dir.mkdir()
    executable = binary_dir / "7za"
    executable.write_text(
        f"#!{sys.executable}\n"
        + textwrap.dedent(
            """\
            import json
            import os
            import sys
            from pathlib import Path

            args = sys.argv[1:]
            assert args.pop(0) == "a"
            while args and args[0].startswith("-p"):
                args.pop(0)
            archive = Path(args.pop(0))
            # Exclusion switches are archive options, not filesystem inputs.
            args = [arg for arg in args if not arg.startswith("-xr!")]
            kind = "qc" if archive.name.startswith("QC_") else "fastq"
            with open(os.environ["FAKE_ARCHIVER_CALLS"], "a") as calls:
                calls.write(json.dumps({"argv": sys.argv[1:], "cwd": os.getcwd()}) + "\\n")
            print(f"{kind}: standard output diagnostic", flush=True)
            print(f"{kind}: standard error diagnostic", file=sys.stderr, flush=True)
            archive.write_text(json.dumps(args))
            missing = [path for path in args if not Path(path).exists()]
            if missing:
                for path in missing:
                    print(f"WARNING: missing input: {path}", file=sys.stderr)
                sys.exit(1)
            if os.environ.get("FAKE_ARCHIVER_FAIL_KIND") == kind:
                print(f"{kind}: deliberate archive failure", file=sys.stderr)
                sys.exit(int(os.environ["FAKE_ARCHIVER_EXIT_CODE"]))
            """
        )
    )
    executable.chmod(0o755)
    calls = tmp_path / "archiver-calls.jsonl"
    monkeypatch.setenv("PATH", f"{binary_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("TMPDIR", str(work_root))
    monkeypatch.setenv("FAKE_ARCHIVER_CALLS", str(calls))
    monkeypatch.delenv("FAKE_ARCHIVER_FAIL_KIND", raising=False)
    monkeypatch.setattr(afterFastq, "get_project_names", Mock(return_value={PROJECT}))
    monkeypatch.setattr(afterFastq, "get_project_dirs", Mock(return_value=set()))
    monkeypatch.setattr(misc.PipelineConfig, "get", Mock(return_value=cfg))
    return SimpleNamespace(
        cfg=cfg,
        inputs=inputs,
        extra=output / "FLOWCELL",
        fastq_archive=output / f"{PROJECT}_{RUN_DATE}.7za",
        qc_archive=output / f"QC_{PROJECT}_{RUN_DATE}.7za",
        calls=calls,
        retained=[fastq, inputs["fastq_checksum"], report, analysis_report, qc_result, workflow],
    )


@pytest.mark.parametrize("extra_layout", ["absent", "directory", "symlink"])
@pytest.mark.parametrize("restart", [False, True])
def test_tenx_archive_accepts_restored_fastqs_and_preserves_optional_output(
    archive_run, tmp_path, extra_layout, restart
):
    run = archive_run
    if extra_layout == "directory":
        run.extra.mkdir()
        (run.extra / "mkfastq-report.txt").write_text("completed mkfastq report\n")
    elif extra_layout == "symlink":
        linked_output = tmp_path / "mkfastq output"
        linked_output.mkdir()
        (linked_output / "mkfastq-report.txt").write_text("completed mkfastq report\n")
        run.extra.symlink_to(linked_output, target_is_directory=True)
    retained_contents = {path: path.read_bytes() for path in run.retained}
    if restart:
        run.fastq_archive.write_text("partial old archive\n")
        run.qc_archive.write_text("stale QC archive\n")
        for prefix in [PROJECT, f"QC_{PROJECT}"]:
            (run.cfg.output_path / f"md5sum_{prefix}_{RUN_DATE}_archive.txt").write_text(
                "stale archive checksum\n"
            )
        apply_cleanup(cleanup_plan(run.cfg.output_path, "finalization"))
        assert not run.fastq_archive.exists()
        assert not run.qc_archive.exists()
        assert not list(run.cfg.output_path.glob("md5sum_*_archive.txt"))

    afterFastq.archive_worker(run.cfg)

    archived_inputs = json.loads(run.fastq_archive.read_text())
    assert (str(run.extra) in archived_inputs) == (extra_layout != "absent")
    assert {str(path) for name, path in run.inputs.items() if name != "qc"}.issubset(
        archived_inputs
    )
    assert json.loads(run.qc_archive.read_text()) == [str(run.inputs["qc"])]
    assert {path: path.read_bytes() for path in run.retained} == retained_contents
    calls = [json.loads(line) for line in run.calls.read_text().splitlines()]
    assert len(calls) == 2
    assert set(undetermined_fastqs.archive_exclusions()).issubset(calls[0]["argv"])
    assert all(call["cwd"] == str(run.cfg.output_path) for call in calls)
    for kind in ["fastq", "qc"]:
        command_log = run.cfg.static.paths.log_dir / f"{RUN_ID}_{PROJECT}_archive_{kind}.log"
        assert f"{kind}: standard output diagnostic" in command_log.read_text()
        assert f"{kind}: standard error diagnostic" in command_log.read_text()


@pytest.mark.parametrize("layout", ["directory", "symlink"])
@pytest.mark.skipif(REAL_ARCHIVER is None, reason="optional real 7-Zip contract check")
def test_archive_worker_filters_nested_mkfastq_inputs_and_preserves_sample_reads(
    archive_run, tmp_path, monkeypatch, layout
):
    run = archive_run
    if layout == "symlink":
        extra = tmp_path / "external-mkfastq"
        run.extra.symlink_to(extra, target_is_directory=True)
    else:
        extra = run.extra
    nested = extra / "outs/fastq_path/staging"
    nested.mkdir(parents=True)
    discarded = [
        run.cfg.output_path / "Undetermined_S0_R1_001.fastq.gz",
        nested / "Undetermined_S0_L001_I2_001.fastq.gz",
        run.inputs["project"] / "Undetermined_S0_L002_R3_001.fastq.gz",
    ]
    for path in discarded:
        path.write_bytes(b"unassigned")
    sample = nested / "Undetermined_control_S1_L001_R1_001.fastq.gz"
    sample.write_bytes(b"sample reads")
    stats = nested / "unknown-index.csv"
    stats.write_bytes(b"unknown index counts")
    original = afterFastq._run_archive_command

    def real_command(cmd, cfg, project, kind):
        if kind == "fastq":
            assert not any(str(path) in cmd for path in discarded)
            assert set(undetermined_fastqs.archive_exclusions()).issubset(cmd)
        return original([REAL_ARCHIVER, *cmd[1:]], cfg, project, kind)

    monkeypatch.setattr(afterFastq, "_run_archive_command", real_command)
    afterFastq.archive_worker(run.cfg)
    listing = subprocess.check_output(
        [REAL_ARCHIVER, "l", "-slt", str(run.fastq_archive)], text=True
    )
    assert "Undetermined_S0" not in listing
    assert sample.name in listing and stats.name in listing
    assert "sample_R1.fastq.gz" in listing
    assert f"md5sum_{PROJECT}_fastq.txt" in listing
    assert all(path.read_bytes() == b"unassigned" for path in discarded)


@pytest.mark.parametrize(
    "missing_input",
    [
        "project",
        "stats",
        "project_samplesheet",
        "samplesheet",
        "submission_form",
        "fastq_checksum",
        "qc",
    ],
)
def test_missing_required_archive_input_remains_fatal(archive_run, missing_input):
    run = archive_run
    path = run.inputs[missing_input]
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()

    with pytest.raises(subprocess.CalledProcessError) as error:
        afterFastq.archive_worker(run.cfg)

    assert error.value.returncode == 1
    assert str(path) in error.value.cmd
    assert f"WARNING: missing input: {path}" in error.value.output


def test_dangling_optional_tenx_symlink_is_not_silently_omitted(archive_run, tmp_path):
    run = archive_run
    run.extra.symlink_to(tmp_path / "unavailable mkfastq output", target_is_directory=True)

    with pytest.raises(subprocess.CalledProcessError) as error:
        afterFastq.archive_worker(run.cfg)

    assert error.value.returncode == 1
    assert str(run.extra) in error.value.cmd
    assert f"WARNING: missing input: {run.extra}" in error.value.output


def test_unreadable_optional_tenx_path_is_not_treated_as_absent(archive_run, monkeypatch):
    run = archive_run
    original_lstat = Path.lstat

    def denied_extra(path):
        if path == run.extra:
            raise PermissionError("Cannot inspect existing mkfastq output")
        return original_lstat(path)

    monkeypatch.setattr(Path, "lstat", denied_extra)
    with pytest.raises(PermissionError, match="Cannot inspect existing mkfastq output"):
        afterFastq.archive_worker(run.cfg)
    assert not run.calls.exists()


@pytest.mark.parametrize("kind", ["fastq", "qc"])
@pytest.mark.parametrize("returncode", [1, 2])
def test_archiver_warning_or_failure_retains_diagnostics_in_error_report(
    archive_run, monkeypatch, kind, returncode
):
    run = archive_run
    monkeypatch.setenv("FAKE_ARCHIVER_FAIL_KIND", kind)
    monkeypatch.setenv("FAKE_ARCHIVER_EXIT_CODE", str(returncode))

    with pytest.raises(subprocess.CalledProcessError) as error:
        try:
            afterFastq.archive_worker(run.cfg)
        except subprocess.CalledProcessError:
            report_path = misc.write_error_report(
                sys.exc_info(), "Got an error during finalization"
            )
            raise

    assert error.value.returncode == returncode
    report = report_path.read_text()
    command_log = run.cfg.static.paths.log_dir / f"{RUN_ID}_{PROJECT}_archive_{kind}.log"
    assert error.value.output == command_log.read_text()
    for expected in [
        f"{kind}: standard output diagnostic",
        f"{kind}: standard error diagnostic",
        f"{kind}: deliberate archive failure",
    ]:
        assert expected in report
        assert expected in error.value.output
    calls = [json.loads(line) for line in run.calls.read_text().splitlines()]
    assert len(calls) == (1 if kind == "fastq" else 2)


@pytest.mark.parametrize("kind", ["fastq", "qc"])
def test_sensitive_archive_failure_redacts_password_without_changing_invocation(
    archive_run, monkeypatch, kind
):
    run = archive_run
    run.cfg.run.sensitive = True
    password = "synthetic-secret-password"
    monkeypatch.setattr(afterFastq, "generate_password", Mock(return_value=password))
    monkeypatch.setenv("FAKE_ARCHIVER_FAIL_KIND", kind)
    monkeypatch.setenv("FAKE_ARCHIVER_EXIT_CODE", "2")

    with pytest.raises(subprocess.CalledProcessError) as error:
        try:
            afterFastq.archive_worker(run.cfg)
        except subprocess.CalledProcessError:
            report_path = misc.write_error_report(
                sys.exc_info(), "Got an error during finalization"
            )
            raise

    calls = [json.loads(line) for line in run.calls.read_text().splitlines()]
    assert all(f"-p{password}" in call["argv"] for call in calls)
    assert password not in str(error.value)
    assert password not in report_path.read_text()
    assert f"{kind}: deliberate archive failure" in report_path.read_text()
