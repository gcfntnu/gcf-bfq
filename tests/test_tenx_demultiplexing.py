"""10x conversion contracts; all external execution and mail are doubled."""

import csv
import io
import json
import subprocess

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import flowcell_manager.flowcell_manager as manager
import pytest
import yaml

from bcl2fastq_pipeline.state import FlowcellStateStore
from configmaker.indexes import normalize_indexes
from test_project_execution import PROJECTS, pair
from test_state_integration import configured_bfq, write_inputs

from bcl2fastq_pipeline import (
    demultiplexing,
    index_sequences,
    makeFastq,
    preflight,
    tenx_demultiplexing,
    workflow_config,
)

GEX = "10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3"
MULTIOME = "10X Genomics Chromium Single Cell Multiome Gene Expression"


@pytest.fixture
def profiles(tmp_path, monkeypatch):
    path = tmp_path / "demultiplexing.config"
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "bcl_convert_versions": ["4.2.4"],
                "profiles": {
                    "gex": {
                        "kits": [GEX],
                        "assay": "gex",
                        "minimum_biological_cycles": [28, 25],
                        "index_lengths": [8, 10],
                        "default_backend": "bcl-convert",
                        "legacy_tool": "cellranger-mkfastq",
                    },
                    "multiome-gex": {
                        "kits": [MULTIOME],
                        "assay": "multiome-gex",
                        "minimum_biological_cycles": [28, 50],
                        "index_lengths": [10],
                        "require_dual_index": True,
                        "default_backend": "bcl-convert",
                        "legacy_tool": "cellranger-mkfastq",
                    },
                    "atac": {
                        "kits": ["10X Genomics ATAC"],
                        "assay": "legacy",
                        "default_backend": "mkfastq",
                        "legacy_tool": "cellranger-atac-mkfastq",
                    },
                },
            }
        )
    )
    monkeypatch.setattr(tenx_demultiplexing, "AUTHORITATIVE_PROFILES", path)
    return path


def physical(path, cycles=(28, 10, 10, 90)):
    path.write_text(
        "<RunInfo><Run><Reads>"
        + "".join(
            f'<Read Number="{i + 1}" NumCycles="{count}" IsIndexedRead="{indexed}"/>'
            for i, (count, indexed) in enumerate(
                zip(
                    cycles,
                    ("N", "Y", "Y", "N") if len(cycles) == 4 else ("N", "Y", "N"),
                    strict=True,
                )
            )
        )
        + "</Reads></Run></RunInfo>"
    )
    return path


def tenx_pair(tmp_path, *, index="SI-GA-A1", kit=GEX, custom="", settings="", project=PROJECTS[0]):  # noqa: PLR0913
    result = pair(tmp_path, project, kit=kit, settings=settings, index=index)
    result.sample_sheet.write_text(
        result.sample_sheet.read_text().replace("[Settings]", custom + "[Settings]")
    )
    return result


def test_single_named_index_preserves_all_oligos_and_one_library(profiles, tmp_path):
    source = tenx_pair(tmp_path)
    before = source.sample_sheet.read_bytes()
    result = demultiplexing.plan(
        SimpleNamespace(pairs=(source,)), run_info=physical(tmp_path / "RunInfo.xml")
    )
    assert result["tool"] == "bcl-convert"
    assert len(result["rows"]) == 4
    assert {row["sample_id"] for row in result["rows"]} == {f"sample-{PROJECTS[0]}"}
    assert {row["effective_parameters"]["index"] for row in result["rows"]} == {
        item["index"] for item in normalize_indexes("SI-GA-A1")
    }
    assert {row["effective_parameters"]["OverrideCycles"] for row in result["rows"]} == {
        "Y28;I8N2;N10;Y90"
    }
    assert all(row["index_provenance"]["original_index"] == "SI-GA-A1" for row in result["rows"])
    parsed = list(csv.reader(io.StringIO(result["sample_sheet_content"])))
    header = parsed[parsed.index(["[BCLConvert_Data]"]) + 1]
    assert "TrimUMI" not in result["global_settings"]
    assert "index2" not in header
    assert "BarcodeMismatchesIndex2" not in header
    assert source.sample_sheet.read_bytes() == before
    assert result["sources"][0]["profile"]["sha256"]


@pytest.mark.parametrize(
    "backend,tool", [("bcl-convert", "bcl-convert"), ("mkfastq", "cellranger-mkfastq")]
)
def test_single_sheet_10x_keeps_historical_force_override_precedence(
    profiles, tmp_path, backend, tool
):
    source = tenx_pair(tmp_path, custom=f"Demultiplexer,{backend}\n")
    result = demultiplexing.plan(
        SimpleNamespace(
            pairs=(), sample_sheet=source.sample_sheet, submission_form=source.submission_form
        ),
        run_info=physical(tmp_path / "RunInfo.xml"),
        force_bcl2fastq=True,
    )
    assert result["tool"] == tool


@pytest.mark.parametrize(
    "orientation,expected", [("forward", "AGTGTTACCT"), ("reverse-complement", "AGGTAACACT")]
)
def test_dual_symbolic_sequence_orientation_is_explicit_and_used_once(
    profiles, tmp_path, orientation, expected
):
    source = tenx_pair(
        tmp_path, index="SI-TT-A1", kit=MULTIOME, custom=f"Index2Orientation,{orientation}\n"
    )
    result = demultiplexing.plan(
        SimpleNamespace(pairs=(source,)), run_info=physical(tmp_path / "RunInfo.xml")
    )
    row = result["rows"][0]["effective_parameters"]
    assert row["index2"] == expected
    assert row["OverrideCycles"] == "Y28;I10;I10;Y90"
    assert result["sources"][0]["custom"]["Libprep"] == MULTIOME


@pytest.mark.parametrize(
    "change,match",
    [
        ("unknown", "Unknown/unsupported"),
        ("orientation", "Index2Orientation"),
        ("trim", "trimming"),
        ("trim-umi", "Remove TrimUMI"),
        ("short", "barcode/UMI"),
        ("cycles", "index allocation"),
        ("atac", "retains legacy"),
    ],
)
def test_unsafe_profiles_fail_before_execution(profiles, tmp_path, change, match):
    source = tenx_pair(tmp_path, index="SI-TT-A1", custom="Index2Orientation,forward\n")
    cycles = (28, 10, 10, 90)
    if change == "unknown":
        source.sample_sheet.write_text(
            source.sample_sheet.read_text().replace("SI-TT-A1", "SI-UNKNOWN-A1")
        )
    elif change == "orientation":
        source.sample_sheet.write_text(
            source.sample_sheet.read_text().replace("Index2Orientation,forward\n", "")
        )
    elif change == "trim":
        source.sample_sheet.write_text(
            source.sample_sheet.read_text().replace(
                "[Settings]\n", "[Settings]\nAdapterRead1,ACGT\n"
            )
        )
    elif change == "trim-umi":
        source.sample_sheet.write_text(
            source.sample_sheet.read_text().replace("[Settings]\n", "[Settings]\nTrimUMI,0\n")
        )
    elif change == "short":
        cycles = (26, 10, 10, 90)
    elif change == "cycles":
        source.sample_sheet.write_text(
            source.sample_sheet.read_text().replace(
                "[Settings]\n", "[Settings]\nOverrideCycles,Y28;I8N2;I10;Y90\n"
            )
        )
    else:
        source.sample_sheet.write_text(
            source.sample_sheet.read_text()
            .replace(GEX, "10X Genomics ATAC")
            .replace("[CustomOptions]\n", "[CustomOptions]\nDemultiplexer,bcl-convert\n")
        )
    with pytest.raises(demultiplexing.DemultiplexingPlanError, match=match):
        demultiplexing.plan(
            SimpleNamespace(pairs=(source,)), run_info=physical(tmp_path / "RunInfo.xml", cycles)
        )


def test_mixed_settings_need_explicit_passes_and_keep_independent_profiles(profiles, tmp_path):
    first = tenx_pair(tmp_path)
    second = pair(tmp_path, PROJECTS[1], index="CCCCCCCCCC", settings="MaskShortReads,22\n")
    selection = SimpleNamespace(pairs=(first, second))
    with pytest.raises(demultiplexing.DemultiplexingPlanError, match="different"):
        demultiplexing.plan(selection, run_info=physical(tmp_path / "RunInfo.xml"))
    first.sample_sheet.write_text(
        first.sample_sheet.read_text().replace(
            "[CustomOptions]\n", "[CustomOptions]\nDemultiplexingPass,tenx\n"
        )
    )
    result = demultiplexing.plan(selection, run_info=physical(tmp_path / "RunInfo.xml"))
    assert result["passes"] == 2
    assert {item["global_settings"].get("MaskShortReads") for item in result["pass_plans"]} == {
        "8",
        "22",
    }


def test_deliberate_symbolic_correction_expands_sequences_and_is_reversible(profiles, tmp_path):
    source = tenx_pair(tmp_path)
    original = source.sample_sheet.read_bytes()
    once = index_sequences.transform(original, ("index1",))
    twice = index_sequences.transform(once.content, ("index1",))
    assert b"SI-GA-A1" in once.content  # retained as annotation, never complemented
    assert index_sequences.orientation(once.content, original)["index1"]["status"] == "reversed"
    assert index_sequences.orientation(twice.content, original)["index1"]["status"] == "original"
    assert {item["index"] for item in normalize_indexes("SI-GA-A1")} == {
        line.split(",")[2]
        for line in twice.content.decode().splitlines()
        if line.startswith("sample-")
    }


def test_real_conversion_contract_uses_profile_and_never_retries_after_failure(
    profiles, tmp_path, monkeypatch
):
    cfg, source, output = configured_bfq(tmp_path)
    (source / "InterOp").mkdir()
    physical(source / "RunInfo.xml")
    input_pair = tenx_pair(output)
    cfg.run.input_selection = SimpleNamespace(pairs=(input_pair,))
    monkeypatch.setattr(makeFastq, "image_version", lambda tool: "4.2.4")
    monkeypatch.setattr(makeFastq, "resolve_image", lambda tool: "gcfntnu/bcl-convert:4.2.4")
    call = Mock(side_effect=subprocess.CalledProcessError(2, "bcl-convert"))
    monkeypatch.setattr(makeFastq.subprocess, "check_call", call)
    with pytest.raises(RuntimeError, match="Demultiplexing failed"):
        makeFastq.bcl2fq()
    assert call.call_count == 1
    manifest = json.loads((output / "provenance/demultiplexing/plan.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["tool_image"] == "gcfntnu/bcl-convert:4.2.4"
    assert manifest["sources"][0]["profile"]["name"] == "gex"


def test_separate_pass_publication_and_qc_do_not_double_count(profiles, tmp_path, monkeypatch):
    from bcl2fastq_pipeline import sequencing_qc  # noqa: PLC0415

    cfg, source, output = configured_bfq(tmp_path)
    (source / "InterOp").mkdir()
    physical(source / "RunInfo.xml")
    first = tenx_pair(output, custom="DemultiplexingPass,tenx\n")
    second = pair(output, PROJECTS[1], index="CCCCCCCCCC", settings="MaskShortReads,22\n")
    cfg.run.input_selection = SimpleNamespace(pairs=(first, second))
    monkeypatch.setattr(makeFastq, "image_version", lambda _tool: "4.2.4")
    monkeypatch.setattr(makeFastq, "resolve_image", lambda _tool: "gcfntnu/bcl-convert:4.2.4")
    commands = []

    def convert(command, **_kwargs):
        commands.append(command)
        destination = Path(command[command.index("--output-directory") + 1])
        tenx = destination.name == "tenx"
        project, count = (PROJECTS[0], 100) if tenx else (PROJECTS[1], 80)
        project_dir = destination / project
        project_dir.mkdir()
        (project_dir / f"sample-{project}_S1_R1_001.fastq.gz").write_bytes(
            b"mocked converter output"
        )
        reports = destination / "Reports"
        reports.mkdir()
        (reports / "Demultiplex_Stats.csv").write_text(
            "Lane,SampleID,Sample_Project,# Reads\n"
            f"1,sample-{project},{project},{count}\n1,Undetermined,,{500 - count}\n"
        )

    monkeypatch.setattr(makeFastq.subprocess, "check_call", convert)
    assert makeFastq.bcl2fq()[0] == "multiple conversion passes"
    assert len(commands) == 2
    for project in PROJECTS:
        assert (output / project / f"sample-{project}_S1_R1_001.fastq.gz").is_file()
        manifest = json.loads(
            (output / "provenance/project-demultiplexing" / project / "plan.json").read_text()
        )
        assert manifest["project"] == project
        assert all(row["project"] == project for row in manifest["rows"])
    warnings = []
    module, sources, samples = sequencing_qc._planned_sources(cfg, output, None, warnings)
    lanes, unknown = sequencing_qc._parse_sources(sources, samples)
    assert module == "custom_content" and warnings
    assert lanes["1"]["total_reads"] == 500
    assert lanes["1"]["undetermined_reads"] is None and unknown is None
    assert sorted(sample["assigned_reads"] for sample in samples.values()) == [80, 100]


def test_unverified_converter_version_fails_before_execution(profiles, tmp_path, monkeypatch):
    cfg, source, output = configured_bfq(tmp_path)
    (source / "InterOp").mkdir()
    physical(source / "RunInfo.xml")
    cfg.run.input_selection = SimpleNamespace(pairs=(tenx_pair(output),))
    monkeypatch.setattr(makeFastq, "image_version", lambda _tool: "99.0")
    convert = Mock()
    monkeypatch.setattr(makeFastq.subprocess, "check_call", convert)
    with pytest.raises(
        demultiplexing.DemultiplexingPlanError, match="Profile requires BCL Convert"
    ):
        makeFastq.bcl2fq()
    convert.assert_not_called()


def test_mixed_legacy_pass_cannot_bypass_whole_flowcell_validation(profiles, tmp_path):
    first = tenx_pair(tmp_path, custom="Demultiplexer,mkfastq\nDemultiplexingPass,legacy\n")
    second = pair(tmp_path, PROJECTS[1], index="CCCCCCCCCC")
    with pytest.raises(demultiplexing.DemultiplexingPlanError, match="collision contract"):
        demultiplexing.plan(
            SimpleNamespace(pairs=(first, second)), run_info=physical(tmp_path / "RunInfo.xml")
        )


@pytest.mark.parametrize("multi_pair", [False, True])
def test_analysis_geometry_follows_retained_conversion_after_source_edit(
    profiles, tmp_path, multi_pair
):
    cfg, source, output = configured_bfq(tmp_path)
    input_pair = tenx_pair(output, settings="OverrideCycles,Y28;I8N2;N10;Y80N10\n")
    physical(output / "RunInfo.xml")
    cfg.run.input_selection = SimpleNamespace(
        pairs=(input_pair,) if multi_pair else (),
        sample_sheet=input_pair.sample_sheet,
        submission_form=input_pair.submission_form,
    )
    demultiplexing.prepare(cfg)
    demultiplexing.record(cfg, status="completed")
    demultiplexing.delivery_provenance(cfg, input_pair.project_id, refresh=True)
    # A metadata correction cannot change the read lengths of retained FASTQs.
    input_pair.sample_sheet.write_text(input_pair.sample_sheet.read_text().replace("Y80N10", "Y90"))
    context = workflow_config.prepare_analysis_projects(cfg)[input_pair.project_id]
    assert context.effective_read_geometry == (28, 80)
    assert context.sample_sheet == input_pair.sample_sheet
    retained = output / "provenance/project-demultiplexing" / input_pair.project_id
    (retained / "SampleSheet.project.csv").write_text("corrupted retained sheet")
    with pytest.raises(demultiplexing.DemultiplexingPlanError, match="hash differs"):
        workflow_config.prepare_analysis_projects(cfg)


def test_clean_initialization_and_copied_inputs_retain_physical_geometry(profiles, tmp_path):
    cfg, source, output = configured_bfq(tmp_path)
    write_inputs(source)
    sheet = source / "SampleSheet.csv"
    sheet.write_text(
        sheet.read_text().replace("Illumina DNA Prep", GEX).replace("ACGT", "SI-GA-A1")
    )
    original_info = physical(source / "RunInfo.xml").read_bytes()
    assert not output.exists()
    state = manager.initialize_flowcell(
        flowcell=str(source), from_stage="demultiplexing", force=True
    )
    assert state["status"] == "queued"
    assert (output / "RunInfo.xml").read_bytes() == original_info
    selection = preflight.select_run_inputs(source, output)
    copied = preflight.copy_run_inputs(selection, output / "attempt-inputs")
    assert copied.run_info == output / "attempt-inputs/RunInfo.xml"
    physical(source / "RunInfo.xml", cycles=(30, 10, 10, 90))
    assert copied.run_info.read_bytes() == original_info
    assert preflight.validate_selection(copied).ok
    legacy_selection = preflight.InputSelection(copied.sample_sheet, copied.submission_form)
    retained = preflight.copy_run_inputs(legacy_selection, output / "retained-inputs")
    assert retained.run_info.read_bytes() == original_info
    assert (
        FlowcellStateStore(cfg.static.paths.manager_dir).read(cfg.run.run_id)["status"] == "queued"
    )
