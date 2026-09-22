import sys

from pathlib import Path

import pytest

PACKAGE_ROOT = Path(__file__).resolve().parents[1] / "bcl2fastq_pipeline"
sys.path.insert(0, str(PACKAGE_ROOT))

from bcl2fastq_pipeline.config import PipelineConfig  # noqa: E402


@pytest.fixture(autouse=True)
def reset_pipeline_config_singleton():
    """Keep singleton state from leaking between tests."""
    PipelineConfig._instance = None
    yield
    PipelineConfig._instance = None
