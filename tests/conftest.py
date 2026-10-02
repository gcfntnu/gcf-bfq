import os
import runpy

from pathlib import Path

import pytest

from bcl2fastq_pipeline.config import PipelineConfig

GUARD = Path(__file__).parent / "support"
# Cover direct pytest use as well as scripts/dev.py, including real SMTP on
# custom ports and subprocesses that inherit the test environment.
runpy.run_path(str(GUARD / "sitecustomize.py"))


@pytest.fixture(autouse=True)
def guard_child_processes(monkeypatch):
    previous = os.environ.get("PYTHONPATH")
    path = str(GUARD.resolve())
    monkeypatch.setenv("PYTHONPATH", path + (os.pathsep + previous if previous else ""))


@pytest.fixture(autouse=True)
def reset_pipeline_config_singleton():
    """Keep singleton state from leaking between tests."""
    PipelineConfig._instance = None
    yield
    PipelineConfig._instance = None


@pytest.fixture(autouse=True)
def forbid_real_smtp(monkeypatch):
    """Notification tests must install a mock; never reach a real relay."""

    def unexpected_smtp(*args, **kwargs):
        raise AssertionError("Tests must mock SMTP explicitly")

    monkeypatch.setattr("smtplib.SMTP", unexpected_smtp)
    monkeypatch.setattr("smtplib.SMTP_SSL", unexpected_smtp)
