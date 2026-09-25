from pathlib import Path

import pytest

from bcl2fastq_pipeline import containers


def write_config(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


def test_build_command_uses_image_from_workflow_config(tmp_path, monkeypatch):
    config = write_config(
        tmp_path / "docker.config",
        "docker:\n  cellranger: gcfntnu/cellranger:10.0.0\n",
    )
    monkeypatch.setenv("BFQ_APPTAINER_COMMAND", "singularity")

    command = containers.build_command("cellranger", ["mkfastq", "--localcores=8"], config)

    assert command == [
        "singularity",
        "exec",
        "docker://gcfntnu/cellranger:10.0.0",
        "cellranger",
        "mkfastq",
        "--localcores=8",
    ]


def test_config_path_can_be_overridden_for_development(tmp_path, monkeypatch):
    config = write_config(
        tmp_path / "docker.config",
        "docker:\n  multiqc: docker://gcfntnu/multiqc:test\n",
    )
    monkeypatch.setenv("GCF_WORKFLOWS_DOCKER_CONFIG", str(config))

    assert containers.build_command("multiqc", ["--version"])[2:] == [
        "docker://gcfntnu/multiqc:test",
        "multiqc",
        "--version",
    ]


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("not_docker: {}\n", "must contain a 'docker' mapping"),
        ("docker:\n  cellranger: [invalid]\n", "non-string image references"),
    ],
)
def test_invalid_config_is_reported(tmp_path, content, message):
    config = write_config(tmp_path / "docker.config", content)

    with pytest.raises(containers.ContainerConfigError, match=message):
        containers.load_images(config)


def test_missing_tool_image_is_reported(tmp_path):
    config = write_config(tmp_path / "docker.config", "docker:\n  multiqc: example/multiqc:1\n")

    with pytest.raises(containers.ContainerConfigError, match="'bcl-convert'.*is missing"):
        containers.build_command("bcl-convert", config_path=config)


def test_run_tool_preserves_argument_boundaries(tmp_path, monkeypatch):
    config = write_config(
        tmp_path / "docker.config",
        "docker:\n  spaceranger: gcfntnu/spaceranger:1.3.0\n",
    )
    calls = []
    monkeypatch.setattr(containers.subprocess, "call", lambda command: calls.append(command) or 17)

    result = containers.run_tool("spaceranger", ["mkfastq", "--id=run with spaces"], config)

    assert result == 17
    assert calls == [
        [
            "apptainer",
            "exec",
            "docker://gcfntnu/spaceranger:1.3.0",
            "spaceranger",
            "mkfastq",
            "--id=run with spaces",
        ]
    ]
