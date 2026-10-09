"""GPU opt-in must cover fresh analysis and retained-workdir resume alike."""

import pytest

from bcl2fastq_pipeline.afterFastq import snakemake_command


@pytest.mark.parametrize("resume", [False, True])
@pytest.mark.parametrize("setting", ["1", "true", "YES", "on"])
def test_gpu_arguments_reach_both_analysis_modes(monkeypatch, tmp_path, resume, setting):
    monkeypatch.setenv("BFQ_GPU", setting)
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(tmp_path / "cache with spaces"))

    command = snakemake_command(resume=resume)

    assert command.count("--singularity-args=--nv") == 1
    assert "--use-singularity" in command
    assert command[command.index("--singularity-prefix") + 1] == str(tmp_path / "cache with spaces")
    assert command[command.index("--cores") + 1] == "32"
    assert command[-1] == "multiqc_report"
    assert ("--rerun-incomplete" in command) == resume
    assert "--keep-incomplete" not in command


@pytest.mark.parametrize("resume", [False, True])
@pytest.mark.parametrize("setting", [None, "0", "false", "NO", "off"])
def test_cpu_execution_does_not_require_nvidia(monkeypatch, tmp_path, resume, setting):
    monkeypatch.delenv("BFQ_GPU", raising=False)
    if setting is not None:
        monkeypatch.setenv("BFQ_GPU", setting)
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(tmp_path))

    command = snakemake_command(resume=resume)

    assert not any("--nv" in arg or "--singularity-args" in arg for arg in command)
    assert ("--rerun-incomplete" in command) == resume


def test_invalid_gpu_setting_fails_explicitly(monkeypatch, tmp_path):
    monkeypatch.setenv("BFQ_GPU", "nvidia")
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(tmp_path))
    with pytest.raises(ValueError, match="BFQ_GPU must be"):
        snakemake_command()


def test_smoke_uses_same_gpu_command_with_its_own_target(monkeypatch, tmp_path):
    monkeypatch.setenv("BFQ_GPU", "1")
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(tmp_path))
    command = snakemake_command(target="gpu_smoke", cores=1)
    assert command[-1] == "gpu_smoke"
    assert command[command.index("--cores") + 1] == "1"
    assert "--singularity-args=--nv" in command
