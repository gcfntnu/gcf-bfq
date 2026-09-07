from pathlib import Path

import pytest

from bcl2fastq_pipeline.config import Paths, PipelineConfig, RunContext


def write_ini(tmp_path: Path, contents: str) -> Path:
    ini_path = tmp_path / "bcl2fastq.ini"
    ini_path.write_text(contents, encoding="utf-8")
    return ini_path


def make_paths(tmp_path: Path) -> Paths:
    return Paths(
        ekista_base_dir=tmp_path / "ekista",
        nova_base_dir=tmp_path / "nova",
        output_dir=tmp_path / "output",
        log_dir=tmp_path / "logs",
        manager_dir=tmp_path / "manager",
    )


def test_loads_static_configuration_and_registers_singleton(tmp_path):
    ini_path = write_ini(
        tmp_path,
        """
[Paths]
ekista_baseDir = /data/ekista
nova_baseDir = /data/nova
outputDir = /data/output
logDir = /data/logs
manager_dir = /data/manager
reportDir = /data/reports
analysisDir = /data/analysis

[System]
sleepTime = 1

[Email]
finishedTo = sequencing@example.org

[Version]
pipeline = 0.3.1

[Commands]
multiqc = /usr/bin/multiqc
""".strip(),
    )

    config = PipelineConfig.load(ini_path)

    assert config is PipelineConfig.get()
    assert config.static.paths.output_dir == Path("/data/output")
    assert config.static.paths.report_dir == Path("/data/reports")
    assert config.static.system == {"sleeptime": "1"}
    assert config.static.email == {"finishedto": "sequencing@example.org"}
    assert config.static.version == {"pipeline": "0.3.1"}
    assert config.static.commands == {"multiqc": "/usr/bin/multiqc"}


def test_load_uses_defaults_for_optional_paths_and_sections(tmp_path):
    ini_path = write_ini(tmp_path, "[Paths]\n")

    config = PipelineConfig.load(ini_path)

    assert config.static.paths.ekista_base_dir == Path("/mnt/seq/ekista")
    assert config.static.paths.nova_base_dir == Path("/mnt/seq/nova")
    assert config.static.paths.output_dir == Path("/mnt/output")
    assert config.static.paths.log_dir == Path("/mnt/logs")
    assert config.static.paths.manager_dir == Path("/mnt/manager")
    assert config.static.paths.report_dir == Path("/mnt/reports")
    assert config.static.paths.analysis_dir == Path("/mnt/analysis")
    assert config.static.system == {}
    assert config.static.email == {}
    assert config.static.version == {}
    assert config.static.commands == {}


def test_get_requires_configuration_to_be_loaded():
    with pytest.raises(RuntimeError, match="PipelineConfig not initialized"):
        PipelineConfig.get()


@pytest.mark.parametrize(
    ("base_name", "expected_source"),
    [("nova", "nova"), ("ekista", "ekista"), ("other", "unknown")],
)
def test_begin_identifies_instrument_source(tmp_path, base_name, expected_source):
    paths = make_paths(tmp_path)
    flowcell_path = tmp_path / base_name / "260901_A01990_0221_AABC123"

    context = RunContext()
    context.begin(flowcell_path, paths)

    assert context.run_id == "260901_A01990_0221_AABC123"
    assert context.flowcell_path == flowcell_path
    assert context.base_dir == flowcell_path.parent
    assert context.instrument_source == expected_source


def test_apply_custom_normalizes_and_promotes_known_options(tmp_path):
    context = RunContext()
    sample_sheet = tmp_path / "SampleSheet.csv"
    submission_form = tmp_path / "SampleSubmissionForm.xlsx"

    context.apply_custom(
        {
            " Libprep ": " RNA-seq ",
            "User": " geir ",
            "Rerun": " YES ",
            "SensitiveData": " 1 ",
            "CustomValue": " retained ",
            "  ": "ignored",
        },
        sample_sheet,
        submission_form,
    )

    assert context.sample_sheet == sample_sheet
    assert context.sample_submission_form == submission_form
    assert context.libprep == "RNA-seq"
    assert context.user == "geir"
    assert context.rerun is True
    assert context.sensitive is True
    assert context.custom == {
        "Libprep": "RNA-seq",
        "User": "geir",
        "Rerun": "YES",
        "SensitiveData": "1",
        "CustomValue": "retained",
    }


def test_reset_clears_run_specific_state(tmp_path):
    context = RunContext(
        run_id="run-id",
        flowcell_path=tmp_path / "run-id",
        base_dir=tmp_path,
        instrument_source="nova",
        sample_sheet=tmp_path / "SampleSheet.csv",
        sample_submission_form=tmp_path / "SampleSubmissionForm.xlsx",
        libprep="RNA-seq",
        pipeline="rnaseq",
        user="geir",
        rerun=True,
        sensitive=True,
        custom={"key": "value"},
    )

    context.reset()

    assert context == RunContext()
