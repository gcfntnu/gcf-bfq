"""Version-specific converter invocation and lossless FASTQ output naming."""

import csv
import io

from types import SimpleNamespace

import pytest

from test_project_execution import PROJECTS, pair
from test_state_integration import configured_bfq

from bcl2fastq_pipeline import makeFastq


@pytest.mark.parametrize("per_sample", [False, True])
def test_legacy_stats_requested_only_for_v1_sheet(tmp_path, monkeypatch, per_sample):
    cfg, source, output = configured_bfq(tmp_path)
    (source / "InterOp").mkdir()
    first = pair(output, PROJECTS[0], settings="AdapterRead1,ACGT\n" if per_sample else "")
    second = pair(output, PROJECTS[1], index="CCCCCCCC")
    cfg.run.input_selection = SimpleNamespace(pairs=(first, second))
    commands = []
    monkeypatch.setattr(makeFastq, "image_version", lambda *_args: "4.2.4")
    monkeypatch.setattr(
        makeFastq.subprocess, "check_call", lambda command, **_kwargs: commands.append(command)
    )
    makeFastq.bcl2fq()
    assert ("--output-legacy-stats" in commands[0]) is not per_sample
    assert cfg.run.demultiplexing_plan["sample_sheet_version"] == (2 if per_sample else 1)
    assert cfg.run.sample_sheet.is_file()


def test_v2_converter_keeps_sample_id_names_and_retains_original_sample_names(
    tmp_path, monkeypatch
):
    cfg, source, output = configured_bfq(tmp_path)
    (source / "InterOp").mkdir()
    first = pair(
        output,
        PROJECTS[0],
        data_settings=",Sample_Name,BarcodeMismatchesIndex1",
        value=",different-display-name,1",
    )
    second = pair(output, PROJECTS[1], index="CCCCCCCC")
    original = first.sample_sheet.read_bytes()
    cfg.run.input_selection = SimpleNamespace(pairs=(first, second))
    commands = []
    monkeypatch.setattr(makeFastq, "image_version", lambda *_args: "4.2.4")
    monkeypatch.setattr(
        makeFastq.subprocess, "check_call", lambda command, **_kwargs: commands.append(command)
    )

    makeFastq.bcl2fq()

    plan = cfg.run.demultiplexing_plan
    rows = list(csv.reader(io.StringIO(plan["sample_sheet_content"])))
    header = rows[rows.index(["[BCLConvert_Data]"]) + 1]
    assert "Sample_Name" not in header
    assert rows[-2][header.index("Sample_ID")] == f"sample-{PROJECTS[0]}"
    assert plan["rows"][0]["effective_parameters"]["Sample_Name"] == "different-display-name"
    assert "--sample-name-column-enabled" not in commands[0]
    assert first.sample_sheet.read_bytes() == original
    retained = (
        output / "provenance" / "demultiplexing" / "sources" / PROJECTS[0] / first.sample_sheet.name
    )
    assert retained.read_bytes() == original


def test_converter_sample_numbers_do_not_erase_embedded_sample_id(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    project = output / PROJECTS[0]
    project.mkdir(parents=True)
    cfg.run.libprep = "Custom"
    (project / "sample_S1_R1_001.fastq.gz").write_bytes(b"first")
    (project / "sample_S1_S2_R1_001.fastq.gz").write_bytes(b"second")
    makeFastq.rename_fastqs()
    assert (project / "sample_R1.fastq.gz").read_bytes() == b"first"
    assert (project / "sample_S1_R1.fastq.gz").read_bytes() == b"second"


def test_rename_conflict_fails_before_moving_any_file(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    project = output / PROJECTS[0]
    project.mkdir(parents=True)
    cfg.run.libprep = "Custom"
    for name in ("alpha_S1_R1_001.fastq.gz", "beta_S2_R1_001.fastq.gz", "beta_R1.fastq.gz"):
        (project / name).write_bytes(name.encode())
    before = {path.name: (path.read_bytes(), path.stat().st_mtime_ns) for path in project.iterdir()}
    with pytest.raises(RuntimeError, match="output-name conflict"):
        makeFastq.rename_fastqs()
    assert {
        path.name: (path.read_bytes(), path.stat().st_mtime_ns) for path in project.iterdir()
    } == before
