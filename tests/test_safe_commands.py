import ast
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pandas as pd
import pytest

from bcl2fastq_pipeline import afterFastq, makeFastq
from flowcell_manager import flowcell_manager


def test_configured_commands_are_split_without_losing_quoted_values():
    assert afterFastq.command_args(
        "/opt/tools/multiqc", "--title 'Project with spaces' --force"
    ) == ["/opt/tools/multiqc", "--title", "Project with spaces", "--force"]
    assert makeFastq.command_args(
        "cellranger mkfastq", "--jobmode=local --localmem 55"
    ) == ["cellranger", "mkfastq", "--jobmode=local", "--localmem", "55"]


def test_bcl_convert_keeps_dynamic_paths_as_single_arguments(tmp_path, monkeypatch):
    flowcell_path = tmp_path / "flowcell with spaces;not-a-command"
    (flowcell_path / "InterOp").mkdir(parents=True)
    output_path = tmp_path / "output with spaces;not-a-command"
    output_path.mkdir()
    log_dir = tmp_path / "logs"
    log_dir.mkdir()
    sample_sheet = flowcell_path / "SampleSheet with spaces.csv"
    sample_sheet.write_text("[Data]\n")

    cfg = SimpleNamespace(
        output_path=output_path,
        run=SimpleNamespace(
            flowcell_path=flowcell_path,
            sample_sheet=sample_sheet,
            libprep="Illumina DNA Prep",
            run_id=flowcell_path.name,
        ),
        static=SimpleNamespace(
            paths=SimpleNamespace(log_dir=log_dir),
            commands={},
        ),
    )
    check_call = Mock()
    monkeypatch.setattr(makeFastq.PipelineConfig, "get", Mock(return_value=cfg))
    monkeypatch.setattr(makeFastq.subprocess, "check_call", check_call)
    monkeypatch.delenv("FORCE_BCL2FASTQ", raising=False)
    monkeypatch.setenv("BCL_CONVERT_VERSION", "4.2.4")

    assert makeFastq.bcl2fq() == ["bcl-convert", "4.2.4"]

    command = check_call.call_args.args[0]
    assert command[0] == "bcl-convert"
    assert command[command.index("--bcl-input-directory") + 1] == str(flowcell_path)
    assert command[command.index("--output-directory") + 1] == str(output_path)
    assert command[command.index("--sample-sheet") + 1] == str(sample_sheet)
    assert "shell" not in check_call.call_args.kwargs


def test_archive_commands_expand_inputs_without_shell_globbing(tmp_path, monkeypatch):
    output_path = tmp_path / "output with spaces;not-a-command"
    project = "GCF-2026-001 project;not-a-command"
    (output_path / project).mkdir(parents=True)
    (output_path / "Stats").mkdir()
    (output_path / "Reports").mkdir()
    (output_path / "Undetermined lane_R1.fastq.gz").touch()
    (output_path / f"{project}_samplesheet.tsv").touch()
    (output_path / "SampleSheet.csv").touch()
    (output_path / "Sample-Submission-Form.xlsx").touch()
    (output_path / f"md5sum_{project}_fastq.txt").touch()

    work_root = tmp_path / "work root"
    qc_dir = work_root / f"{project}_260923" / "data" / "tmp" / "rnaseq" / "bfq"
    qc_dir.mkdir(parents=True)

    cfg = SimpleNamespace(
        output_path=output_path,
        run=SimpleNamespace(
            run_id="260923_A01990_0001_ABC",
            sensitive=False,
            libprep="RNA-seq",
            pipeline="rnaseq",
        ),
    )
    check_call = Mock()
    monkeypatch.setattr(afterFastq, "get_project_names", Mock(return_value={project}))
    monkeypatch.setattr(afterFastq, "get_project_dirs", Mock(return_value=set()))
    monkeypatch.setattr(afterFastq.subprocess, "check_call", check_call)
    monkeypatch.setenv("TMPDIR", str(work_root))

    afterFastq.archive_worker(cfg)

    fastq_command = check_call.call_args_list[0].args[0]
    qc_command = check_call.call_args_list[1].args[0]
    assert fastq_command[:2] == ["7za", "a"]
    assert str(output_path / "Undetermined lane_R1.fastq.gz") in fastq_command
    assert str(output_path / project) in fastq_command
    assert qc_command == [
        "7za",
        "a",
        str(output_path / f"QC_{project}_260923.7za"),
        str(qc_dir),
    ]
    assert all("shell" not in call.kwargs for call in check_call.call_args_list)


def test_archive_md5_uses_file_handles_instead_of_redirection(tmp_path, monkeypatch):
    archive = tmp_path / "archive with spaces;not-a-command.7za"
    archive.write_bytes(b"archive")

    def write_checksum(command, stdout, cwd):
        stdout.write(f"checksum  {archive.name}\n")

    check_call = Mock(side_effect=write_checksum)
    monkeypatch.setattr(afterFastq.subprocess, "check_call", check_call)

    afterFastq.md5sum_archive(archive)

    assert check_call.call_args.args[0] == ["md5sum", archive.name]
    assert check_call.call_args.kwargs["cwd"] == archive.parent
    assert "shell" not in check_call.call_args.kwargs
    assert (tmp_path / "md5sum_archive with spaces;not-a-command_archive.txt").exists()


def test_fastq_md5_generation_handles_paths_with_shell_metacharacters(tmp_path, monkeypatch):
    output_path = tmp_path / "output with spaces"
    project = "GCF-2026-001 project;not-a-command"
    project_path = output_path / project
    project_path.mkdir(parents=True)
    fastqs = [project_path / "sample one_R1.fastq.gz", project_path / "sample two_R2.fastq.gz"]
    for fastq in fastqs:
        fastq.touch()

    cfg = SimpleNamespace(output_path=output_path)
    monkeypatch.setattr(afterFastq, "get_project_names", Mock(return_value={project}))
    monkeypatch.setattr(afterFastq, "get_project_dirs", Mock(return_value=set()))

    afterFastq.md5sum_worker(cfg)

    md5_file = output_path / f"md5sum_{project}_fastq.txt"
    assert md5_file.read_text().splitlines() == [
        f"d41d8cd98f00b204e9800998ecf8427e  {fastq.relative_to(output_path)}"
        for fastq in fastqs
    ]


def test_workflow_commands_keep_config_values_as_single_arguments(tmp_path, monkeypatch):
    output_path = tmp_path / "output with spaces;not-a-command"
    output_path.mkdir()
    work_root = tmp_path / "workflow cache with spaces"
    project = "GCF-2026-001 project;not-a-command"
    cfg = SimpleNamespace(
        output_path=output_path,
        static=SimpleNamespace(paths=SimpleNamespace(log_dir=tmp_path / "logs")),
        run=SimpleNamespace(
            run_id="260923_A01990_0001_ABC",
            libprep="RNA prep;not-a-command",
            pipeline="rnaseq",
        ),
    )
    check_call = Mock()
    run_logged_command = Mock()
    monkeypatch.setattr(afterFastq, "get_project_names", Mock(return_value={project}))
    monkeypatch.setattr(afterFastq, "get_project_dirs", Mock(return_value=set()))
    monkeypatch.setattr(afterFastq, "get_sequencer", Mock(return_value="Nova Seq;not-a-command"))
    monkeypatch.setattr(afterFastq.shutil, "copytree", Mock())
    monkeypatch.setattr(afterFastq.shutil, "copy2", Mock())
    monkeypatch.setattr(afterFastq.subprocess, "check_call", check_call)
    monkeypatch.setattr(afterFastq, "run_logged_command", run_logged_command)
    monkeypatch.setenv("TMPDIR", str(work_root))
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(tmp_path / "cache with spaces"))

    afterFastq.full_align(cfg)

    configmaker_command = check_call.call_args_list[0].args[0]
    snakemake_command = run_logged_command.call_args.args[0]
    assert configmaker_command[configmaker_command.index("--libkit") + 1] == cfg.run.libprep
    assert configmaker_command[configmaker_command.index("--machine") + 1] == (
        "Nova Seq;not-a-command"
    )
    assert snakemake_command[snakemake_command.index("--singularity-prefix") + 1] == str(
        tmp_path / "cache with spaces"
    )
    assert all("shell" not in call.kwargs for call in check_call.call_args_list)
    assert "shell" not in run_logged_command.call_args.kwargs
    assert run_logged_command.call_args.kwargs["log_path"] == (
        tmp_path / "logs" / f"{cfg.run.run_id}_{project}_snakemake.log"
    )


def test_logged_command_preserves_output_and_raises_with_tail(tmp_path, capsys):
    command = [
        sys.executable,
        "-c",
        "import sys; print('standard output'); print('workflow error', file=sys.stderr); sys.exit(3)",
    ]
    log_path = tmp_path / "snakemake.log"

    with pytest.raises(subprocess.CalledProcessError) as error:
        afterFastq.run_logged_command(command, cwd=tmp_path, log_path=log_path)

    assert error.value.returncode == 3
    assert set(error.value.output.splitlines()) == {"standard output", "workflow error"}
    assert log_path.read_text() == error.value.output
    assert capsys.readouterr().out == error.value.output


def test_flowcell_rerun_deletes_only_the_inventory_path(tmp_path, monkeypatch):
    manager_dir = tmp_path / "manager"
    manager_dir.mkdir()
    flowcell = tmp_path / "run;touch escaped"
    flowcell.mkdir()
    (flowcell / "result.txt").touch()
    unrelated = tmp_path / "escaped"
    unrelated.touch()

    inventory = pd.DataFrame(
        [
            {
                "project": "GCF-2026-001",
                "flowcell_path": str(flowcell),
                "timestamp": "2026-09-23T12:00:00",
                "archived": 0,
            }
        ]
    )
    inventory.to_csv(manager_dir / "flowcells.processed", index=False)
    cfg = SimpleNamespace(static=SimpleNamespace(paths=SimpleNamespace(manager_dir=manager_dir)))
    monkeypatch.setattr(flowcell_manager, "get_cfg", Mock(return_value=cfg))

    flowcell_manager.rerun_flowcell(flowcell=str(flowcell), force=True)

    assert not flowcell.exists()
    assert unrelated.exists()
    assert pd.read_csv(manager_dir / "flowcells.processed").empty


def test_python_sources_do_not_enable_shell_execution():
    repository_root = Path(__file__).parents[1]
    source_paths = list((repository_root / "bcl2fastq_pipeline").rglob("*.py"))
    for source_path in source_paths:
        if "build" in source_path.relative_to(repository_root).parts:
            continue
        tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            for keyword in node.keywords:
                assert not (
                    keyword.arg == "shell"
                    and isinstance(keyword.value, ast.Constant)
                    and keyword.value.value is True
                ), f"shell=True remains in {source_path}:{node.lineno}"
