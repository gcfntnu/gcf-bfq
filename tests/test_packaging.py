import importlib
import sys

from importlib import metadata
from unittest.mock import Mock

from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.entrypoint import _remove_executable_directory_from_import_path


def test_distribution_contains_runtime_dependencies():
    requirements = metadata.requires("bcl2fastq-pipeline") or []
    requirement_names = {requirement.split()[0].split(">")[0] for requirement in requirements}

    assert {"gcf-tools", "pandas", "PyYAML", "urllib3"} <= requirement_names


def test_console_scripts_are_installed():
    scripts = {
        entry_point.name: entry_point.value
        for entry_point in metadata.entry_points(group="console_scripts")
    }

    assert scripts["bfq"] == "bcl2fastq_pipeline.entrypoint:main"
    assert "bfq.py" not in scripts

    flowcell_manager = "flowcell_manager.flowcell_manager:main"
    assert scripts["flowcell-manager"] == flowcell_manager
    assert "flowcell_manager.py" not in scripts


def test_entrypoint_avoids_neighboring_script_shadowing(tmp_path, monkeypatch):
    executable_directory = tmp_path / "bin"
    executable_directory.mkdir()
    (executable_directory / "configmaker.py").write_text(
        'raise RuntimeError("neighboring script was imported")\n', encoding="utf-8"
    )
    monkeypatch.setattr(sys, "argv", [str(executable_directory / "bfq")])
    monkeypatch.setattr(sys, "path", [str(executable_directory), *sys.path])
    monkeypatch.delitem(sys.modules, "configmaker", raising=False)
    monkeypatch.delitem(sys.modules, "configmaker.configmaker", raising=False)

    _remove_executable_directory_from_import_path()

    configmaker = importlib.import_module("configmaker.configmaker")
    assert hasattr(configmaker, "SEQUENCERS")


def test_importing_cli_does_not_start_pipeline(monkeypatch):
    load_config = Mock()
    monkeypatch.setattr(PipelineConfig, "load", load_config)
    sys.modules.pop("bcl2fastq_pipeline.cli", None)

    importlib.import_module("bcl2fastq_pipeline.cli")

    load_config.assert_not_called()
