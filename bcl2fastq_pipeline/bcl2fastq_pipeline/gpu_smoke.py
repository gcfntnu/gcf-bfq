"""Operator-invoked GPU verification through BFQ's normal Snakemake command."""

from __future__ import annotations

import argparse
import ctypes
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

from pathlib import Path

from bcl2fastq_pipeline import containers
from bcl2fastq_pipeline.entrypoint import _remove_executable_directory_from_import_path

GPU_PROGRAM = """import json
import sys
from pathlib import Path
import cupy as cp

count = cp.cuda.runtime.getDeviceCount()
if count < 1:
    raise RuntimeError("CUDA reported no GPU devices")
with cp.cuda.Device(0):
    properties = cp.cuda.runtime.getDeviceProperties(0)
    values = cp.arange(1024, dtype=cp.int64)
    result = int(cp.sum(values * values).get())
    cp.cuda.runtime.deviceSynchronize()
expected = sum(value * value for value in range(1024))
if result != expected:
    raise RuntimeError(f"GPU computation mismatch: {result} != {expected}")
name = properties["name"]
if isinstance(name, bytes):
    name = name.decode()
evidence = {
    "device_count": count,
    "device_name": name,
    "cuda_driver_version": cp.cuda.runtime.driverGetVersion(),
    "cuda_runtime_version": cp.cuda.runtime.runtimeGetVersion(),
    "cupy_version": cp.__version__,
    "result": result,
    "expected": expected,
    "synchronized": True,
}
Path(sys.argv[1]).write_text(json.dumps(evidence, indent=2) + "\\n")
print(json.dumps(evidence, sort_keys=True))
"""


class GPUSmokeError(RuntimeError):
    """An actionable failure in the nested GPU execution path."""


def rapids_image(config_path: Path | None) -> str:
    """Use the scientific workflow's selected image, without a BFQ version pin."""
    images = containers.load_images(config_path)
    image = images.get("rapids-scanpy", "").strip()
    if not image:
        raise GPUSmokeError(
            f"Missing 'rapids-scanpy' image in {config_path or containers.docker_config_path()}"
        )
    return containers.apptainer_image(image)


def command_output(command: list[str]) -> str:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    output = result.stdout + result.stderr
    if result.returncode:
        raise GPUSmokeError(f"{shlex.join(command)} failed:\n{output.strip()}")
    return output.strip()


def apptainer_ldconfig() -> tuple[str, str]:
    """Resolve ldconfig using the pinned Apptainer 1.5 binary-search contract."""
    binary_path = command_output(["singularity", "config", "global", "--get", "binary path"])
    binary_path = (
        binary_path or "$PATH:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
    )
    buildcfg = command_output(["singularity", "buildcfg"])
    match = re.search(r"^LIBEXECDIR=(.+)$", buildcfg, re.M)
    if not match:
        raise GPUSmokeError("Cannot determine Apptainer LIBEXECDIR from singularity buildcfg")
    libexec = match[1].strip().strip("\"'")
    if not Path(libexec).is_absolute():
        raise GPUSmokeError(f"Apptainer LIBEXECDIR is not an absolute directory: {libexec!r}")
    internal = str(Path(libexec) / "apptainer" / "bin") + ":"
    user_path = os.environ.get("PATH", "")
    if user_path:
        user_path += ":"
    search_path = (
        binary_path.replace("$PATH:", internal + user_path, 1)
        if "$PATH:" in binary_path
        else internal + binary_path
    )
    # Apptainer prefers Ubuntu's real executable over an ldconfig wrapper.
    executable = shutil.which("ldconfig.real", path=search_path) or shutil.which(
        "ldconfig", path=search_path
    )
    if not executable:
        raise GPUSmokeError(f"Apptainer cannot find ldconfig using binary path {binary_path!r}")
    return executable, search_path


def preflight() -> dict:
    """Inspect the outer BFQ container before downloading or launching an image."""
    devices = sorted(str(path) for path in Path("/dev").glob("nvidia[0-9]*"))
    if not devices:
        raise GPUSmokeError(
            "No /dev/nvidia GPU devices in BFQ. Enable Docker GPU passthrough with "
            "--gpus all and NVIDIA_DRIVER_CAPABILITIES=compute,utility; verify the "
            "host NVIDIA driver and NVIDIA Container Toolkit."
        )
    try:
        ctypes.CDLL("libcuda.so.1")
    except OSError as error:
        raise GPUSmokeError(
            "BFQ cannot load libcuda.so.1. Check NVIDIA Container Toolkit injection, "
            "NVIDIA_DRIVER_CAPABILITIES=compute,utility and the outer library paths. "
            f"Loader error: {error}"
        ) from error
    for executable in ("singularity", "snakemake", "nvidia-smi"):
        if shutil.which(executable) is None:
            raise GPUSmokeError(f"Required executable {executable!r} is missing from BFQ PATH")
    if "--nv" not in command_output(["singularity", "exec", "--help"]):
        raise GPUSmokeError("The installed singularity/Apptainer does not advertise --nv support")
    ldconfig, binary_path = apptainer_ldconfig()
    if "libcuda.so.1" not in command_output([ldconfig, "-p"]):
        raise GPUSmokeError(
            f"Apptainer's ldconfig ({ldconfig}) cannot discover libcuda.so.1. "
            "Check injected driver library directories, /etc/ld.so.conf.d and the "
            "BFQ image's Apptainer 'binary path' setting."
        )
    return {
        "devices": devices,
        "apptainer": command_output(["singularity", "--version"]),
        "ldconfig": ldconfig,
        "apptainer_binary_path": binary_path,
        "nvidia_smi": command_output(
            ["nvidia-smi", "--query-gpu=name,driver_version,uuid", "--format=csv,noheader"]
        ),
    }


def nested_failure_hint(output: str) -> str:
    """Keep the original log and attach a diagnosis specific to common GPU failures."""
    if re.search(r"insufficient.?driver|unsupported.?ptx|system.?driver.?mismatch", output, re.I):
        return (
            "The nested CUDA runtime rejected the exposed driver. Compare the recorded "
            "host driver with the configured RAPIDS CUDA requirements; also verify that "
            "Apptainer bound the injected host libcuda.so.1 rather than a stub library."
        )
    if re.search(r"libcuda|no.?cuda.?capable|no.?device|no nvidia", output, re.I):
        return (
            "The nested image could not access the CUDA device/driver. Inspect the logged "
            "singularity exec --nv command, device bindings and Apptainer library discovery."
        )
    return (
        "Inspect the Snakemake log below for image download, Apptainer GPU binding, "
        "CUDA initialization or computation errors."
    )


def run_smoke(config_path: Path | None = None, workdir_parent: Path | None = None) -> Path:
    # Import only after CLI path sanitization, like the main BFQ entrypoint.
    from bcl2fastq_pipeline import afterFastq  # noqa: PLC0415

    if not afterFastq.gpu_enabled():
        raise GPUSmokeError("GPU smoke requires BFQ_GPU=1, matching the GPU-enabled BFQ launcher")
    if not os.environ.get("SINGULARITY_CACHEDIR"):
        raise GPUSmokeError("Set SINGULARITY_CACHEDIR to the BFQ launcher's writable image cache")
    image = rapids_image(config_path)
    evidence = preflight()
    workdir = Path(tempfile.mkdtemp(prefix="bfq-gpu-smoke-", dir=workdir_parent)).resolve()
    print(f"GPU smoke artifacts: {workdir}", flush=True)
    print(f"Workflow image: {image}", flush=True)
    (workdir / "gpu_check.py").write_text(GPU_PROGRAM)
    (workdir / "Snakefile").write_text(
        "rule gpu_smoke:\n"
        "    input: 'gpu_check.py'\n"
        "    output: '{phase}.json'\n"
        "    wildcard_constraints: phase='normal|resume'\n"
        "    resources: gpu=1\n"
        f"    container: {image!r}\n"
        "    shell: 'python {input:q} {output:q}'\n"
    )
    evidence.update(
        {"image": image, "docker_config": str(config_path or containers.docker_config_path())}
    )
    (workdir / "preflight.json").write_text(json.dumps(evidence, indent=2) + "\n")
    for phase in ("normal", "resume"):
        command = afterFastq.snakemake_command(
            resume=phase == "resume", target=f"{phase}.json", cores=1
        )
        if "--singularity-args=--nv" not in command:
            raise GPUSmokeError("BFQ's Snakemake command did not propagate --nv")
        (workdir / f"{phase}.command.json").write_text(json.dumps(command, indent=2) + "\n")
        log_path = workdir / f"{phase}.log"
        try:
            afterFastq.run_logged_command(command, cwd=workdir, log_path=log_path)
        except subprocess.CalledProcessError as error:
            output = error.output or ""
            raise GPUSmokeError(
                f"{phase} GPU check failed. {nested_failure_hint(output)}\n"
                f"Log: {log_path}\n{output}"
            ) from error
        try:
            result = json.loads((workdir / f"{phase}.json").read_text())
            valid = (
                result["device_count"] >= 1
                and result["result"] == result["expected"] == 357389824
                and result["synchronized"] is True
            )
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise GPUSmokeError(
                f"{phase} GPU result is missing or invalid; inspect {log_path}"
            ) from error
        if not valid:
            raise GPUSmokeError(f"{phase} GPU result failed validation; inspect {log_path}")
    print("GPU computation passed through BFQ/Snakemake/Apptainer for normal and resume commands.")
    return workdir


def main(argv: list[str] | None = None) -> int:
    _remove_executable_directory_from_import_path()
    parser = argparse.ArgumentParser(
        description="Verify NVIDIA GPU computation through BFQ, Snakemake and nested Apptainer."
    )
    parser.add_argument(
        "--docker-config", type=Path, help="workflow docker.config (default: installed workflow)"
    )
    parser.add_argument(
        "--workdir-parent",
        type=Path,
        help="parent for retained unique smoke artifacts (default: TMPDIR)",
    )
    args = parser.parse_args(argv)
    try:
        run_smoke(args.docker_config, args.workdir_parent)
    except (GPUSmokeError, containers.ContainerConfigError, OSError, ValueError) as error:
        print(f"bfq-gpu-smoke: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
