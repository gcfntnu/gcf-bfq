import pytest

from bcl2fastq_pipeline.config import PipelineConfig


@pytest.fixture(autouse=True)
def reset_pipeline_config_singleton():
    """Keep singleton state from leaking between tests."""
    PipelineConfig._instance = None
    yield
    PipelineConfig._instance = None
