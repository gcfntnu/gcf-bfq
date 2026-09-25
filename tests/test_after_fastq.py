from pathlib import Path
from subprocess import CalledProcessError
from unittest.mock import ANY, Mock

import pytest

from bcl2fastq_pipeline.interop import prepare_index_metrics, run_interop_csv


def test_prepare_index_metrics_links_bcl_convert_report(tmp_path):
    source = tmp_path / "Reports" / "IndexMetricsOut.bin"
    source.parent.mkdir()
    source.write_bytes(b"index metrics")

    prepare_index_metrics(tmp_path)

    destination = tmp_path / "InterOp" / "IndexMetricsOut.bin"
    assert destination.is_symlink()
    assert destination.read_bytes() == b"index metrics"
    assert destination.readlink() == Path("..") / "Reports" / "IndexMetricsOut.bin"


def test_prepare_index_metrics_preserves_existing_interop_file(tmp_path):
    destination = tmp_path / "InterOp" / "IndexMetricsOut.bin"
    destination.parent.mkdir()
    destination.write_bytes(b"existing metrics")
    source = tmp_path / "Reports" / "IndexMetricsOut.bin"
    source.parent.mkdir()
    source.write_bytes(b"bcl convert metrics")

    prepare_index_metrics(tmp_path)

    assert not destination.is_symlink()
    assert destination.read_bytes() == b"existing metrics"


def test_run_interop_csv_publishes_nonempty_output(tmp_path, monkeypatch):
    output_path = tmp_path / "Stats" / "interop_index-summary.csv"
    output_path.parent.mkdir()

    def write_csv(command, stdout, stderr, cwd):
        stdout.write("Lane,Sample\n1,sample\n")

    check_call = Mock(side_effect=write_csv)
    monkeypatch.setattr("bcl2fastq_pipeline.interop.subprocess.check_call", check_call)

    run_interop_csv("interop_index-summary", tmp_path, output_path, output_path.parent)

    assert output_path.read_text() == "Lane,Sample\n1,sample\n"
    assert not output_path.with_suffix(".csv.tmp").exists()
    check_call.assert_called_once_with(
        ["interop_index-summary", str(tmp_path), "--csv=1"],
        stdout=ANY,
        stderr=ANY,
        cwd=output_path.parent,
    )


def test_run_interop_csv_rejects_empty_output(tmp_path, monkeypatch):
    output_path = tmp_path / "interop_index-summary.csv"
    output_path.write_text("previous report\n")
    monkeypatch.setattr("bcl2fastq_pipeline.interop.subprocess.check_call", Mock())

    with pytest.raises(RuntimeError, match="produced an empty CSV"):
        run_interop_csv("interop_index-summary", tmp_path, output_path, tmp_path)

    assert output_path.read_text() == "previous report\n"
    assert not output_path.with_suffix(".csv.tmp").exists()


def test_run_interop_csv_preserves_previous_output_on_command_failure(tmp_path, monkeypatch):
    output_path = tmp_path / "interop_index-summary.csv"
    output_path.write_text("previous report\n")

    def fail_with_diagnostics(command, stderr, **_kwargs):
        stderr.write("libgomp.so.1: cannot open shared object file\n")
        stderr.flush()
        raise CalledProcessError(127, command)

    monkeypatch.setattr(
        "bcl2fastq_pipeline.interop.subprocess.check_call",
        Mock(side_effect=fail_with_diagnostics),
    )

    with pytest.raises(RuntimeError, match="libgomp.so.1.*cannot open shared object file"):
        run_interop_csv("interop_index-summary", tmp_path, output_path, tmp_path)

    assert output_path.read_text() == "previous report\n"
    assert not output_path.with_suffix(".csv.tmp").exists()
    assert not output_path.with_suffix(".csv.stderr.tmp").exists()
