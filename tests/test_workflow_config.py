import json
import logging

from pathlib import Path
from unittest.mock import Mock

import pytest
import yaml

from bcl2fastq_pipeline.state import FlowcellStateStore, new_state
from configmaker.configmaker import add_workflow
from configmaker.libprep import LibprepConfigError
from test_state_integration import RUN_ID, configured_bfq, write_fastq, write_inputs

from bcl2fastq_pipeline import (
    afterFastq,
    cli,
    findFlowCells,
    makeFastq,
    misc,
    notifications,
    workflow_config,
)

CONTENT = b"""# uncommitted local parameter edit
Illumina DNA Prep SE:
  workflow: rnaseq
  filter:
    trim:
      fastp:
        params: '-q 17'
Illumina DNA Prep PE:
  workflow: metagenome
  filter:
    trim:
      fastp:
        params: '-q 19'
Custom PE:
  workflow: default
"""


def write_stats(output, geometry):
    (output / "Stats").mkdir(exist_ok=True)
    (output / "Stats/Stats.json").write_text(
        json.dumps(
            {
                "ReadInfosForLanes": [
                    {
                        "ReadInfos": [
                            {"IsIndexedRead": False, "NumCycles": length} for length in geometry
                        ]
                    }
                ]
            }
        )
    )


@pytest.fixture
def workflow_run(tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    write_inputs(output)
    write_inputs(source)
    cfg.run.libprep = "Illumina DNA Prep"
    authoritative = tmp_path / "opt/gcf-workflows/libprep.config"
    authoritative.parent.mkdir(parents=True)
    authoritative.write_bytes(CONTENT)
    override = tmp_path / "ignored.config"
    override.write_bytes(CONTENT.replace(b"metagenome", b"wrong").replace(b"-q 19", b"-q 99"))
    monkeypatch.setenv("BFQ_LIBPREP_CONFIG", str(override))
    monkeypatch.setattr(workflow_config, "AUTHORITATIVE_CONFIG", authoritative)
    monkeypatch.setenv("TMPDIR", str(tmp_path / "analysis"))
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(tmp_path / "cache"))
    return cfg, source, output, authoritative


@pytest.mark.parametrize("stage", ["demultiplexing", "analysis"])
@pytest.mark.parametrize(
    "geometry,workflow,parameter", [([75], "rnaseq", "-q 17"), ([150, 150], "metagenome", "-q 19")]
)
def test_daemon_and_configmaker_share_snapshot_through_source_edits(  # noqa: PLR0913
    workflow_run, monkeypatch, caplog, stage, geometry, workflow, parameter
):
    cfg, source, output, authoritative = workflow_run
    caplog.set_level(logging.INFO)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(new_state(RUN_ID, source, output, origin="new", start_stage=stage, cfg=cfg))
    projects = ["GCF-2026-001", "GCF-2026-002"]
    for project in projects:
        write_fastq(output / project / "sample_R1.fastq.gz")
    write_stats(output, geometry)
    real_prepare = workflow_config.prepare_execution

    def prepare(cfg, stage):
        assert store.execution_active(RUN_ID)
        real_prepare(cfg, stage)
        # Edit after capture, before workflow copying/selection.
        authoritative.write_bytes(CONTENT.replace(b"-q 19", b"-q 88"))

    generated = []

    def configmaker(command, cwd):
        def value(flag):
            return command[command.index(flag) + 1]

        assert Path(value("--libprep-config")).read_bytes() == CONTENT
        with monkeypatch.context() as context:
            context.chdir(cwd)
            config = add_workflow(
                {
                    "libprepkit": value("--libkit"),
                    "read_geometry": geometry,
                    "filter": {"subsample_fastq": "skip"},
                },
                libprep_config=value("--libprep-config"),
                expected_sha256=value("--libprep-sha256"),
                expected_entry=value("--libprep-entry"),
                expected_read_geometry=[
                    int(n) for n in command[command.index("--expected-read-geometry") + 1 :]
                ],
            )
            Path("config.yaml").write_text(yaml.safe_dump(config))
        assert config["workflow"] == cfg.run.pipeline == workflow
        assert config["filter"]["trim"]["fastp"]["params"] == parameter
        assert config["libprep_selection"]["sha256"] == cfg.run.libprep_config.sha256
        generated.append(config)
        # A second project must still receive the original bytes.
        authoritative.write_bytes(b"broken: [")

    def snakemake(_command, cwd, log_path):
        project = cwd.name.rsplit("_", 1)[0]
        result = cwd / "data/tmp" / workflow / "bfq"
        result.mkdir(parents=True)
        (result / f"multiqc_{project}.html").write_text("report")
        (result / ".multiqc_config.yaml").write_text("{}")
        (cwd / "data/tmp/sample_info.tsv").write_text("sample")

    monkeypatch.setattr(workflow_config, "prepare_execution", prepare)
    monkeypatch.setattr(misc, "enoughFreeSpace", lambda: True)
    monkeypatch.setattr(makeFastq, "bcl2fq", lambda: ("test", "1"))
    monkeypatch.setattr(makeFastq, "rename_fastqs", lambda: None)
    monkeypatch.setattr(afterFastq.subprocess, "check_call", configmaker)
    monkeypatch.setattr(afterFastq, "run_logged_command", snakemake)
    monkeypatch.setattr(cli, "_run_reporting", lambda *_args: projects)
    monkeypatch.setattr(afterFastq, "finalize", lambda: None)
    monkeypatch.setattr(findFlowCells, "markFinished", lambda: projects)
    monkeypatch.setattr(notifications, "send_notification", lambda *_args: None)
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("test"), prepare=True)
    assert store.read(RUN_ID)["status"] == "completed"
    assert len(generated) == 2
    assert (output / workflow_config.SNAPSHOT_NAME).read_bytes() == CONTENT
    saved = json.loads((output / workflow_config.SELECTION_NAME).read_text())
    assert saved["source"] == str(authoritative)
    assert saved["entry"].endswith(" SE" if len(geometry) == 1 else " PE")
    assert "retired and ignored" in caplog.text
    assert saved["sha256"] in caplog.text


@pytest.mark.parametrize("problem", ["missing", "malformed", "unknown"])
def test_configuration_errors_prevent_analysis_launch(workflow_run, monkeypatch, problem):
    cfg, _source, output, authoritative = workflow_run
    write_stats(output, [150, 150])
    if problem == "missing":
        authoritative.unlink()
    elif problem == "malformed":
        authoritative.write_text("workflow: [")
    else:
        cfg.run.libprep = "Misspelled kit"
    launch = Mock()
    monkeypatch.setattr(afterFastq.subprocess, "check_call", launch)
    monkeypatch.setattr(afterFastq, "run_logged_command", launch)
    with pytest.raises(LibprepConfigError, match="configuration|Unsupported library kit"):
        afterFastq.analysis_steps()
    launch.assert_not_called()


def test_analysis_restart_refreshes_but_later_stages_restore(workflow_run):
    cfg, _source, output, authoritative = workflow_run
    write_stats(output, [150, 150])
    first = workflow_config.select_workflow(cfg)
    authoritative.write_bytes(CONTENT.replace(b"metagenome", b"rnaseq"))
    for stage in ("reporting", "finalization"):
        workflow_config.prepare_execution(cfg, stage)
        assert cfg.run.pipeline == "metagenome"
        assert cfg.run.libprep_config.sha256 == first.config.sha256
    workflow_config.prepare_execution(cfg, "analysis")
    second = workflow_config.select_workflow(cfg)
    assert second.workflow == "rnaseq"
    assert second.config.sha256 != first.config.sha256
    cfg.run.reset()
    assert cfg.run.libprep_config is None
    assert cfg.run.workflow_selection is None


def test_corrupt_retained_snapshot_fails_closed(workflow_run):
    cfg, _source, output, _authoritative = workflow_run
    write_stats(output, [150, 150])
    workflow_config.select_workflow(cfg)
    (output / workflow_config.SNAPSHOT_NAME).write_bytes(CONTENT.replace(b"-q 19", b"-q 20"))
    with pytest.raises(LibprepConfigError, match="snapshot files or restart from analysis"):
        workflow_config.prepare_execution(cfg, "finalization")


def test_legacy_later_stage_bootstrap_is_explicit(workflow_run, caplog):
    cfg, _source, output, _authoritative = workflow_run
    write_stats(output, [150, 150])
    workflow_config.prepare_execution(cfg, "reporting")
    assert cfg.run.pipeline == "metagenome"
    assert "resolving legacy results" in caplog.text
