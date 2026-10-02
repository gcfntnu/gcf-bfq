"""Durable diagnostics for failures in the documented development command."""

import importlib.util
import json
import os
import subprocess
import sys

from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "bfq_dev", Path(__file__).resolve().parents[1] / "scripts/dev.py"
)
dev = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dev)


def test_failed_check_retains_transcript_junit_and_summary(tmp_path):
    probe = tmp_path / "test_probe.py"
    probe.write_text('def test_probe():\n    assert False, "deliberate diagnostic probe"\n')
    # No repository conftest/config is involved in this disposable failure.
    with pytest.raises(subprocess.CalledProcessError):
        dev.pytest_check(Path(sys.executable), os.environ.copy(), tmp_path, [str(probe)])
    summary = json.loads((tmp_path / "tests-summary.json").read_text())
    assert summary["outcome"] == "failed"
    assert summary["elapsed_seconds"] > 0
    assert "deliberate diagnostic probe" in (tmp_path / "tests.log").read_text()
    assert 'failures="1"' in (tmp_path / "tests-junit.xml").read_text()
    assert probe.is_file()


def test_new_check_does_not_inherit_another_scenarios_config(tmp_path, monkeypatch):
    monkeypatch.setenv("BFQ_TEST_SCENARIOS", "1")
    monkeypatch.setenv("BFQ_TEST_SCENARIO_ROOT", "/another/invocation")
    env = dev.environment(tmp_path, guarded=True)
    assert not any(key.startswith("BFQ_TEST_") for key in env)
    assert env["PYTHONPATH"] == str(dev.GUARD)
