import sys

from importlib import metadata
from types import SimpleNamespace
from unittest.mock import Mock

import flowcell_manager.flowcell_manager as manager
import pytest
import yaml

from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.state import FlowcellStateStore, collect_versions, new_state

from bcl2fastq_pipeline import entrypoint, version


@pytest.mark.parametrize(
    ("command", "main"),
    [("bfq", entrypoint.main), ("fm", manager.main), ("flowcell-manager", manager.main)],
)
def test_version_exits_without_configuration_or_processing(command, main, monkeypatch, capsys):
    load = Mock(side_effect=AssertionError("Version must not load site configuration"))
    run = Mock(side_effect=AssertionError("Version must not start processing"))
    monkeypatch.setattr(PipelineConfig, "load", load)
    monkeypatch.setitem(sys.modules, "bcl2fastq_pipeline.cli", SimpleNamespace(main=run))
    monkeypatch.setattr(sys, "argv", [command, "--version"])
    monkeypatch.setattr(sys, "path", sys.path.copy())

    with pytest.raises(SystemExit) as result:
        main()

    assert result.value.code == 0
    assert capsys.readouterr().out == f"BFQ {metadata.version('bcl2fastq-pipeline')}\n"
    load.assert_not_called()
    run.assert_not_called()


def test_bfq_without_arguments_still_starts_pipeline(monkeypatch):
    run = Mock()
    monkeypatch.setitem(sys.modules, "bcl2fastq_pipeline.cli", SimpleNamespace(main=run))
    monkeypatch.setattr(sys, "argv", ["bfq"])
    monkeypatch.setattr(sys, "path", sys.path.copy())

    entrypoint.main()

    run.assert_called_once_with()


def test_missing_package_metadata_is_unknown(monkeypatch, capsys):
    monkeypatch.setattr(
        version.metadata, "version", Mock(side_effect=metadata.PackageNotFoundError)
    )
    monkeypatch.setattr(sys, "argv", ["bfq", "--version"])
    monkeypatch.setattr(sys, "path", sys.path.copy())

    assert collect_versions()["bfq"] is None
    with pytest.raises(SystemExit) as result:
        entrypoint.main()

    assert result.value.code == 0
    assert capsys.readouterr().out == "BFQ unknown (package not installed)\n"


def test_installed_version_overrides_legacy_label_and_preserves_history(tmp_path, monkeypatch):
    ini = tmp_path / "bfq.ini"
    ini.write_text(
        f"[Paths]\noutputDir = {tmp_path / 'output'}\n"
        "[Version]\npipeline = 0.3.1\ndeployment = site-label\n",
        encoding="utf-8",
    )
    cfg = PipelineConfig.load(ini)
    run_id = "260918_MN00686_0026_A000HCMFHF"
    output = cfg.static.paths.output_dir / run_id
    state = new_state(
        run_id, tmp_path / run_id, output, origin="new", start_stage="demultiplexing", cfg=cfg
    )
    installed_version = metadata.version("bcl2fastq-pipeline")
    assert state["versions"]["bfq"] == installed_version

    store = FlowcellStateStore(tmp_path / "manager")
    store.create(state)
    with monkeypatch.context() as previous_installation:
        previous_installation.setattr(version.metadata, "version", lambda name: "0.3.1")
        store.begin_attempt(run_id, cfg=cfg)
    store.fail_stage(run_id, "demultiplexing", summary="test failure", report_path=None)
    store.queue(run_id, "demultiplexing")

    updated = store.begin_attempt(run_id, cfg=cfg)

    assert updated["schema_version"] == 1
    assert updated["versions"]["bfq"] == installed_version
    assert updated["attempts"][0]["versions"]["bfq"] == "0.3.1"
    assert updated["attempts"][1]["versions"]["bfq"] == installed_version
    snapshot = tmp_path / "config.yaml"
    cfg.to_file(snapshot)
    assert yaml.safe_load(snapshot.read_text())["static"]["version"] == {
        "pipeline": "0.3.1",
        "deployment": "site-label",
    }
