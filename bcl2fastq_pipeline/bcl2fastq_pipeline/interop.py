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
    try:
        with temporary_output.open("w") as output_fh:
            subprocess.check_call(
                [command, str(run_path), "--csv=1"],
                stdout=output_fh,
                cwd=cwd,
            )
        if temporary_output.stat().st_size == 0:
            raise RuntimeError(f"{command} produced an empty CSV")
        temporary_output.replace(output_path)
    finally:
        temporary_output.unlink(missing_ok=True)
