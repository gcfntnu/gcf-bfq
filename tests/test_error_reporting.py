import logging
import subprocess
import sys

from types import SimpleNamespace
from unittest.mock import Mock

from bcl2fastq_pipeline import cli, misc


def test_error_report_uses_flowcell_name_and_contains_traceback(tmp_path, monkeypatch):
    report_dir = tmp_path / "reports"
    cfg = SimpleNamespace(
        static=SimpleNamespace(paths=SimpleNamespace(report_dir=report_dir)),
        run=SimpleNamespace(run_id="260924_A01990_0221_TEST"),
    )
    monkeypatch.setattr(misc.PipelineConfig, "get", Mock(return_value=cfg))

    try:
        raise RuntimeError("workflow failed")
    except RuntimeError:
        report_path = misc.errorEmail(sys.exc_info(), "Got an error during postMakeSteps")

    assert report_path == report_dir / "260924_A01990_0221_TEST.error"
    report = report_path.read_text()
    assert "Got an error during postMakeSteps" in report
    assert "Traceback (most recent call last):" in report
    assert "RuntimeError: workflow failed" in report


def test_run_context_is_reset_after_error_report(monkeypatch):
    events = []
    cfg = SimpleNamespace(run=SimpleNamespace(reset=lambda: events.append("reset")))
    log = logging.getLogger("test-error-reporting")

    def write_report(error_info, message):
        assert error_info[0] is RuntimeError
        events.append("report")

    monkeypatch.setattr(misc, "errorEmail", write_report)

    try:
        raise RuntimeError("workflow failed")
    except RuntimeError:
        cli.report_run_error(cfg, log, "workflow failed")

    assert events == ["report", "reset"]


def test_error_report_includes_captured_command_output(tmp_path, monkeypatch):
    report_dir = tmp_path / "reports"
    cfg = SimpleNamespace(
        static=SimpleNamespace(paths=SimpleNamespace(report_dir=report_dir)),
        run=SimpleNamespace(run_id="260924_A01990_0221_TEST"),
    )
    monkeypatch.setattr(misc.PipelineConfig, "get", Mock(return_value=cfg))

    try:
        raise subprocess.CalledProcessError(
            1,
            ["snakemake", "multiqc_report"],
            output="/usr/bin/bash: wget: command not found\nError in rule ensembl_genome\n",
        )
    except subprocess.CalledProcessError:
        report_path = misc.errorEmail(sys.exc_info(), "Got an error during postMakeSteps")

    report = report_path.read_text()
    assert "Captured command output (last 400 lines):" in report
    assert "wget: command not found" in report
    assert "Error in rule ensembl_genome" in report
