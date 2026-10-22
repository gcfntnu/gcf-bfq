"""10x conversion contracts; all external execution and mail are doubled."""

import json
import subprocess

from types import SimpleNamespace
from unittest.mock import Mock

import pytest
import yaml

from configmaker.indexes import normalize_indexes
from test_project_execution import PROJECTS, pair
from test_state_integration import configured_bfq

from bcl2fastq_pipeline import demultiplexing, index_sequences, makeFastq, tenx_demultiplexing

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
    assert source.sample_sheet.read_bytes() == before
    assert result["sources"][0]["profile"]["sha256"]


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
