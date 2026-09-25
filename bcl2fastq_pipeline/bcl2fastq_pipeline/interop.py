import subprocess

from pathlib import Path


def prepare_index_metrics(output_path: Path):
    """Expose BCL Convert index metrics in the layout expected by InterOp."""
    index_metrics = output_path / "InterOp" / "IndexMetricsOut.bin"
    bcl_convert_index_metrics = output_path / "Reports" / "IndexMetricsOut.bin"

    if index_metrics.exists() or not bcl_convert_index_metrics.exists():
        return

    index_metrics.parent.mkdir(parents=True, exist_ok=True)
    index_metrics.symlink_to(Path("..") / "Reports" / "IndexMetricsOut.bin")


def run_interop_csv(command: str, run_path: Path, output_path: Path, cwd: Path):
    """Run an InterOp CSV command and publish only a successful, nonempty result."""
    temporary_output = output_path.with_suffix(f"{output_path.suffix}.tmp")
    temporary_error = output_path.with_suffix(f"{output_path.suffix}.stderr.tmp")
    try:
        with temporary_output.open("w") as output_fh, temporary_error.open("w") as error_fh:
            try:
                subprocess.check_call(
                    [command, str(run_path), "--csv=1"],
                    stdout=output_fh,
                    stderr=error_fh,
                    cwd=cwd,
                )
            except subprocess.CalledProcessError as error:
                error_fh.flush()
                details = temporary_error.read_text(errors="replace").strip()
                if not details:
                    details = "No stderr output was produced."
                raise RuntimeError(
                    f"{command} failed with exit code {error.returncode}:\n{details}"
                ) from error
        if temporary_output.stat().st_size == 0:
            raise RuntimeError(f"{command} produced an empty CSV")
        temporary_output.replace(output_path)
    finally:
        temporary_output.unlink(missing_ok=True)
        temporary_error.unlink(missing_ok=True)
