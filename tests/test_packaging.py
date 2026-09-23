import importlib
import sys

from importlib import metadata
from unittest.mock import Mock

from bcl2fastq_pipeline.config import PipelineConfig


def test_distribution_contains_runtime_dependencies():
    requirements = metadata.requires("bcl2fastq-pipeline") or []
    requirement_names = {requirement.split()[0].split(">")[0] for requirement in requirements}

    assert {"gcf-tools", "pandas", "PyYAML", "urllib3"} <= requirement_names


def test_console_scripts_are_installed():
    scripts = {
        entry_point.name: entry_point.value
        for entry_point in metadata.entry_points(group="console_scripts")
    }

    expected = "bcl2fastq_pipeline.cli:main"
    assert scripts["bfq"] == expected
    assert scripts["bfq.py"] == expected

    flowcell_manager = "flowcell_manager.flowcell_manager:main"
    assert scripts["flowcell-manager"] == flowcell_manager


def test_importing_cli_does_not_start_pipeline(monkeypatch):
    load_config = Mock()
    monkeypatch.setattr(PipelineConfig, "load", load_config)
    sys.modules.pop("bcl2fastq_pipeline.cli", None)

    importlib.import_module("bcl2fastq_pipeline.cli")

    load_config.assert_not_called()
