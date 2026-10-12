"""GPU operator smoke wiring, with all external execution and CUDA calls doubled."""

import json
import subprocess
import sys

from importlib import metadata
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from bcl2fastq_pipeline import afterFastq, gpu_smoke


def gpu_result():
    return {
        "device_count": 1,
        "device_name": "Test GPU",
        "cuda_driver_version": 12090,
        "cuda_runtime_version": 12090,
        "cupy_version": "test",
        "result": 357389824,
        "expected": 357389824,
        "synchronized": True,
    }


@pytest.fixture
def configured_smoke(tmp_path, monkeypatch):
    config = tmp_path / "docker.config"
    config.write_text("docker:\n  rapids-scanpy: example/rapids:configured\n")
    monkeypatch.setenv("BFQ_GPU", "1")
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(tmp_path / "images"))
    monkeypatch.setattr(gpu_smoke, "preflight", lambda: {"nvidia_smi": "mock GPU, mock driver"})
    return config


def test_smoke_uses_shared_normal_and_resume_commands(configured_smoke, tmp_path, monkeypatch):
    commands = []

    def run(command, cwd, log_path):
        commands.append(command)
        (cwd / command[-1]).write_text(json.dumps(gpu_result()))
        log_path.write_text("mock external execution\n")

    monkeypatch.setattr(afterFastq, "run_logged_command", run)
    workdir = gpu_smoke.run_smoke(configured_smoke, tmp_path)

    assert workdir.parent == tmp_path
    assert len(commands) == 2
    for index, phase in enumerate(("normal", "resume")):
        command = commands[index]
        assert command == afterFastq.snakemake_command(
            resume=phase == "resume", target=f"{phase}.json", cores=1
        )
        assert "--singularity-args=--nv" in command
        assert ("--rerun-incomplete" in command) == (phase == "resume")
        assert json.loads((workdir / f"{phase}.command.json").read_text()) == command
        assert (workdir / f"{phase}.json").is_file()
    snakefile = (workdir / "Snakefile").read_text()
    assert "docker://example/rapids:configured" in snakefile
    assert "python {input:q} {output:q}" in snakefile
    evidence = json.loads((workdir / "preflight.json").read_text())
    assert evidence["image"] == "docker://example/rapids:configured"
    assert evidence["docker_config"] == str(configured_smoke)
    compile((workdir / "gpu_check.py").read_text(), "gpu_check.py", "exec")


def test_smoke_requires_explicit_gpu_enablement(configured_smoke, monkeypatch):
    monkeypatch.delenv("BFQ_GPU")
    with pytest.raises(gpu_smoke.GPUSmokeError, match="BFQ_GPU=1"):
        gpu_smoke.run_smoke(configured_smoke)


def test_smoke_rejects_missing_argument_propagation(configured_smoke, tmp_path, monkeypatch):
    monkeypatch.setattr(afterFastq, "snakemake_command", lambda **kwargs: ["snakemake"])
    with pytest.raises(gpu_smoke.GPUSmokeError, match="did not propagate --nv"):
        gpu_smoke.run_smoke(configured_smoke, tmp_path)


def test_smoke_requires_existing_cache_configuration(configured_smoke, monkeypatch):
    monkeypatch.delenv("SINGULARITY_CACHEDIR")
    with pytest.raises(gpu_smoke.GPUSmokeError, match="SINGULARITY_CACHEDIR"):
        gpu_smoke.run_smoke(configured_smoke)


def test_smoke_reports_missing_workflow_image(configured_smoke):
    configured_smoke.write_text("docker: {}\n")
    with pytest.raises(gpu_smoke.GPUSmokeError, match="Missing 'rapids-scanpy'"):
        gpu_smoke.run_smoke(configured_smoke)


@pytest.mark.parametrize("phase", ["normal", "resume"])
def test_smoke_reports_driver_failure_and_retains_log(
    configured_smoke, tmp_path, monkeypatch, phase
):
    commands = []

    def run(command, cwd, log_path):
        commands.append(command)
        if command[-1] == f"{phase}.json":
            output = "cudaErrorInsufficientDriver: CUDA driver version is insufficient\n"
            log_path.write_text(output)
            raise subprocess.CalledProcessError(1, command, output=output)
        (cwd / command[-1]).write_text(json.dumps(gpu_result()))

    monkeypatch.setattr(afterFastq, "run_logged_command", run)
    with pytest.raises(gpu_smoke.GPUSmokeError, match="runtime rejected the exposed driver"):
        gpu_smoke.run_smoke(configured_smoke, tmp_path)
    workdir = next(tmp_path.glob("bfq-gpu-smoke-*"))
    assert "cudaErrorInsufficientDriver" in (workdir / f"{phase}.log").read_text()
    assert len(commands) == (1 if phase == "normal" else 2)


@pytest.mark.parametrize("result", [None, {}, {**gpu_result(), "result": 0}])
def test_success_exit_without_verified_gpu_result_is_failure(
    configured_smoke, tmp_path, monkeypatch, result
):
    def run(command, cwd, log_path):
        if result is not None:
            (cwd / command[-1]).write_text(json.dumps(result))

    monkeypatch.setattr(afterFastq, "run_logged_command", run)
    with pytest.raises(gpu_smoke.GPUSmokeError, match="GPU result"):
        gpu_smoke.run_smoke(configured_smoke, tmp_path)


@pytest.fixture
def outer_gpu(monkeypatch):
    monkeypatch.setattr(Path, "glob", lambda self, pattern: [Path("/dev/nvidia0")])
    monkeypatch.setattr(gpu_smoke.ctypes, "CDLL", lambda name: object())
    monkeypatch.setattr(
        gpu_smoke.shutil,
        "which",
        lambda name, path=None: (
            None if name == "ldconfig.real" else f"/sbin/{name}" if path else f"/usr/bin/{name}"
        ),
    )
    outputs = {
        ("singularity", "exec", "--help"): "--nv enable NVIDIA",
        ("singularity", "config", "global", "--get", "binary path"): "/sbin:/usr/sbin:$PATH:/bin",
        ("singularity", "buildcfg"): "PREFIX=/opt/conda\nLIBEXECDIR=/opt/conda/libexec",
        ("/sbin/ldconfig", "-p"): "libcuda.so.1 => /usr/lib/nvidia/libcuda.so.1",
        ("singularity", "--version"): "apptainer mock",
        (
            "nvidia-smi",
            "--query-gpu=name,driver_version,uuid",
            "--format=csv,noheader",
        ): "GPU, driver, UUID",
    }
    monkeypatch.setattr(gpu_smoke, "command_output", lambda command: outputs[tuple(command)])
    return outputs


def test_preflight_records_driver_and_actual_apptainer_library_lookup(outer_gpu):
    evidence = gpu_smoke.preflight()
    assert evidence["devices"] == ["/dev/nvidia0"]
    assert evidence["ldconfig"] == "/sbin/ldconfig"
    assert evidence["nvidia_smi"] == "GPU, driver, UUID"


def test_apptainer_library_lookup_prefers_system_path_and_real_binary(outer_gpu, monkeypatch):
    calls = []
    monkeypatch.setenv("PATH", "/opt/conda/bin:/usr/bin")

    def which(name, path=None):
        calls.append((name, path))
        return "/sbin/ldconfig.real" if name == "ldconfig.real" else None

    monkeypatch.setattr(gpu_smoke.shutil, "which", which)
    executable, path = gpu_smoke.apptainer_ldconfig()
    assert executable == "/sbin/ldconfig.real"
    assert path == "/sbin:/usr/sbin:/opt/conda/libexec/apptainer/bin:/opt/conda/bin:/usr/bin:/bin"
    assert calls == [("ldconfig.real", path)]


def test_apptainer_library_lookup_uses_default_when_config_key_unset(outer_gpu, monkeypatch):
    outer_gpu[("singularity", "config", "global", "--get", "binary path")] = ""
    monkeypatch.setenv("PATH", "/opt/conda/bin")
    _executable, path = gpu_smoke.apptainer_ldconfig()
    assert path.startswith("/opt/conda/libexec/apptainer/bin:/opt/conda/bin:")


def test_preflight_distinguishes_missing_docker_devices(outer_gpu, monkeypatch):
    monkeypatch.setattr(Path, "glob", lambda self, pattern: [])
    with pytest.raises(gpu_smoke.GPUSmokeError, match="Docker GPU passthrough"):
        gpu_smoke.preflight()


def test_preflight_distinguishes_unloadable_outer_driver(outer_gpu, monkeypatch):
    def missing_library(name):
        raise OSError("missing shared object")

    monkeypatch.setattr(gpu_smoke.ctypes, "CDLL", missing_library)
    with pytest.raises(gpu_smoke.GPUSmokeError, match="BFQ cannot load libcuda.so.1"):
        gpu_smoke.preflight()


def test_preflight_distinguishes_missing_apptainer_gpu_support(outer_gpu):
    outer_gpu[("singularity", "exec", "--help")] = "ordinary exec help"
    with pytest.raises(gpu_smoke.GPUSmokeError, match="does not advertise --nv"):
        gpu_smoke.preflight()


def test_preflight_distinguishes_wrong_library_lookup(outer_gpu):
    outer_gpu[("/sbin/ldconfig", "-p")] = "libc.so.6"
    with pytest.raises(gpu_smoke.GPUSmokeError, match="cannot discover libcuda.so.1"):
        gpu_smoke.preflight()


def test_nested_diagnostics_distinguish_missing_driver_binding():
    assert "device/driver" in gpu_smoke.nested_failure_hint("libcuda.so.1: cannot open")
    assert "image download" in gpu_smoke.nested_failure_hint("connection refused")


@pytest.mark.parametrize("result", [357389824, 0])
def test_gpu_program_requires_compute_result_and_synchronization(tmp_path, monkeypatch, result):
    cupy = MagicMock()
    cupy.cuda.runtime.getDeviceCount.return_value = 1
    cupy.cuda.runtime.getDeviceProperties.return_value = {"name": b"Mock GPU"}
    cupy.cuda.runtime.driverGetVersion.return_value = 12090
    cupy.cuda.runtime.runtimeGetVersion.return_value = 12090
    cupy.sum.return_value.get.return_value = result
    cupy.__version__ = "mock"
    monkeypatch.setitem(sys.modules, "cupy", cupy)
    output = tmp_path / "result.json"
    monkeypatch.setattr(sys, "argv", ["gpu_check.py", str(output)])

    if result:
        exec(compile(gpu_smoke.GPU_PROGRAM, "gpu_check.py", "exec"), {})
        evidence = json.loads(output.read_text())
        assert evidence["device_name"] == "Mock GPU"
        assert evidence["result"] == evidence["expected"] == result
        assert evidence["synchronized"] is True
    else:
        with pytest.raises(RuntimeError, match="computation mismatch"):
            exec(compile(gpu_smoke.GPU_PROGRAM, "gpu_check.py", "exec"), {})
        assert not output.exists()
    cupy.arange.assert_called_once_with(1024, dtype=cupy.int64)
    cupy.cuda.runtime.deviceSynchronize.assert_called_once_with()


def test_cli_help_never_starts_gpu_work(monkeypatch, capsys):
    monkeypatch.setattr(gpu_smoke, "run_smoke", MagicMock(side_effect=AssertionError))
    with pytest.raises(SystemExit) as exit_info:
        gpu_smoke.main(["--help"])
    assert exit_info.value.code == 0
    assert "Snakemake" in capsys.readouterr().out


def test_cli_reports_preflight_failure_without_traceback(monkeypatch, capsys):
    monkeypatch.setattr(
        gpu_smoke, "run_smoke", MagicMock(side_effect=gpu_smoke.GPUSmokeError("missing GPU"))
    )
    assert gpu_smoke.main([]) == 1
    assert capsys.readouterr().err.strip() == "bfq-gpu-smoke: missing GPU"


def test_gpu_smoke_console_script_is_installed():
    scripts = metadata.entry_points(group="console_scripts", name="bfq-gpu-smoke")
    assert [script.value for script in scripts] == ["bcl2fastq_pipeline.gpu_smoke:main"]
