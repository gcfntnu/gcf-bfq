import pytest

from bcl2fastq_pipeline.config import PipelineConfig


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
