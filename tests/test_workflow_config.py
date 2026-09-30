import json
import logging
import tarfile

from pathlib import Path
from unittest.mock import Mock

import pytest
import yaml

from bcl2fastq_pipeline.state import FlowcellStateStore, new_state
from configmaker.configmaker import add_workflow
from configmaker.libprep import LibprepConfigError
from test_state_integration import (
    RUN_ID,
    completed_state,
    configured_bfq,
    write_fastq,
    write_inputs,
)

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

    def prepare(cfg, stage, **kwargs):
        assert store.execution_active(RUN_ID)
        real_prepare(cfg, stage, **kwargs)
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
            Path("Snakefile").write_text("# generated workflow entry point\n")
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
    assert not (output / "bfq-libprep.config").exists()
    assert not (output / "bfq-libprep.json").exists()
    metadata = store.read(RUN_ID)["stages"]["analysis"]["metadata"]
    assert metadata["workflow"] == workflow
    assert metadata["input_preflight"]["status"] == "passed"
    selection = generated[0]["libprep_selection"]
    assert selection["entry"].endswith(" SE" if len(geometry) == 1 else " PE")
    assert str(authoritative) in caplog.text
    assert "retired and ignored" in caplog.text
    assert selection["sha256"] in caplog.text
    for project in projects:
        with tarfile.open(output / "provenance" / f"{project}_analysis.tar.gz") as archive:
            assert archive.extractfile("src/gcf-workflows/libprep.config").read() == CONTENT
            retained_config = yaml.safe_load(archive.extractfile("config.yaml"))
            assert retained_config["filter"]["trim"]["fastp"]["params"] == parameter
            assert archive.extractfile("Snakefile").read() == b"# generated workflow entry point\n"
            assert not any(
                name == "data" or name.startswith("data/") for name in archive.getnames()
            )


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


def test_analysis_restart_refreshes_but_later_stages_restore_only_workflow(workflow_run):
    cfg, _source, output, authoritative = workflow_run
    write_stats(output, [150, 150])
    first = workflow_config.select_workflow(cfg)
    authoritative.unlink()
    for stage in ("reporting", "finalization"):
        workflow_config.prepare_execution(cfg, stage, workflow=first.workflow)
        assert cfg.run.pipeline == "metagenome"
        assert cfg.run.libprep_config is None
        assert cfg.run.workflow_selection is None
    authoritative.write_bytes(CONTENT.replace(b"metagenome", b"rnaseq"))
    workflow_config.prepare_execution(cfg, "analysis")
    second = workflow_config.select_workflow(cfg)
    assert second.workflow == "rnaseq"
    assert second.config.sha256 != first.config.sha256
    cfg.run.reset()
    assert cfg.run.libprep_config is None
    assert cfg.run.workflow_selection is None


@pytest.mark.parametrize("stage", ["reporting", "finalization"])
@pytest.mark.parametrize("recorded", [False, True])
def test_downstream_daemon_recovers_only_workflow(workflow_run, monkeypatch, stage, recorded):
    cfg, source, output, authoritative = workflow_run
    authoritative.unlink()
    projects = ["GCF-2026-001", "GCF-2026-002"]
    for project in projects:
        write_fastq(output / project / "sample_R1.fastq.gz")
        if not recorded:
            work = cfg.static.paths.analysis_dir / f"{project}_260918"
            work.mkdir(parents=True)
            (work / "config.yaml").write_text("workflow: metagenome\n")
    # Artifacts from the first PR revision must neither be read nor be needed.
    for name in ("bfq-libprep.config", "bfq-libprep.json"):
        (output / name).write_text("obsolete: [")
    state = new_state(RUN_ID, source, output, origin="legacy_rerun", start_stage=stage, cfg=cfg)
    if recorded:
        state["stages"]["analysis"]["metadata"]["workflow"] = "metagenome"
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(state)
    calls = []

    def check_workflow(*_args):
        assert store.execution_active(RUN_ID)
        assert cfg.run.pipeline == "metagenome"
        assert cfg.run.libprep_config is None
        assert cfg.run.workflow_selection is None
        calls.append("work")
        return projects

    monkeypatch.setattr(misc, "enoughFreeSpace", lambda: True)
    monkeypatch.setattr(afterFastq, "analysis_steps", Mock(side_effect=AssertionError("analysis")))
    monkeypatch.setattr(cli, "_run_reporting", check_workflow)
    monkeypatch.setattr(afterFastq, "finalize", check_workflow)
    monkeypatch.setattr(findFlowCells, "markFinished", lambda: projects)
    monkeypatch.setattr(notifications, "send_notification", lambda *_args: None)
    cli._run_state_backed_flowcell(cfg, store, logging.getLogger("test"), prepare=True)
    assert store.read(RUN_ID)["status"] == "completed"
    assert store.read(RUN_ID)["stages"]["analysis"]["metadata"] == {"workflow": "metagenome"}
    assert len(calls) == (2 if stage == "reporting" else 1)
    for name in ("bfq-libprep.config", "bfq-libprep.json"):
        assert (output / name).read_text() == "obsolete: ["


@pytest.mark.parametrize("problem", ["missing", "malformed", "conflicting", "no-workflow"])
def test_legacy_recovery_does_not_guess_from_current_source(workflow_run, problem):
    cfg, source, output, _authoritative = workflow_run
    projects = ["GCF-2026-001", "GCF-2026-002"]
    for project in projects:
        work = cfg.static.paths.analysis_dir / f"{project}_260918"
        work.mkdir(parents=True)
        (work / "config.yaml").write_text("workflow: metagenome\n")
    path = work / "config.yaml"
    if problem == "missing":
        path.unlink()
    elif problem == "malformed":
        path.write_text("workflow: [")
    elif problem == "no-workflow":
        path.write_text("{}")
    else:
        path.write_text("workflow: rnaseq")
    state = new_state(
        RUN_ID, source, output, origin="legacy_rerun", start_stage="reporting", cfg=cfg
    )
    state["projects"] = projects
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(state)
    with (
        store.execution_lease(RUN_ID),
        pytest.raises(LibprepConfigError, match="restart from analysis"),
    ):
        cli._prepare_workflow(cfg, store, state)
    assert store.read(RUN_ID) == state
    assert cfg.run.libprep_config is None


@pytest.mark.parametrize("workflow", [None, "../other", "", 7])
def test_missing_or_invalid_recorded_workflow_fails_closed(workflow_run, workflow):
    cfg, _source, _output, _authoritative = workflow_run
    with pytest.raises(LibprepConfigError, match="flowcell state"):
        workflow_config.prepare_execution(cfg, "finalization", workflow=workflow)
    assert cfg.run.libprep_config is None


@pytest.mark.parametrize("stage", ["demultiplexing", "analysis", "reporting", "finalization"])
def test_restart_invalidates_workflow_only_when_analysis_will_rerun(workflow_run, stage):
    cfg, source, output, _authoritative = workflow_run
    state = completed_state(cfg, source, output)
    state["stages"]["analysis"]["metadata"]["workflow"] = "metagenome"
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.create(state)
    queued = store.queue(RUN_ID, stage)
    expected = {} if stage in ("demultiplexing", "analysis") else {"workflow": "metagenome"}
    assert queued["stages"]["analysis"]["metadata"] == expected
