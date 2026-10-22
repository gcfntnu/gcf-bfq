"""Conversion plans must preserve validated identities and explicit settings."""

import csv
import io

from types import SimpleNamespace

import pytest

from test_project_preflight import pair

from bcl2fastq_pipeline import demultiplexing, preflight

A, B = "GCF-2099-001", "GCF-2099-002"


def inputs(directory, first="alpha", second="beta"):
    directory.mkdir()
    pair(directory, A, first, "AAAACCCC")
    pair(directory, B, second, "TTTTGGGG")
    return preflight.select_run_inputs(directory, directory)


def add_column(sheet, name, value):
    rows = list(csv.reader(io.StringIO(sheet.read_text())))
    data = rows.index(["[Data]"]) + 1
    rows[data].append(name)
    rows[data + 1].append(value)
    with sheet.open("w", newline="") as handle:
        csv.writer(handle, lineterminator="\n").writerows(rows)


def test_derived_data_preserves_validated_ascii_id_spelling(tmp_path):
    selection = inputs(tmp_path / "source", first="sample ", second="sample")
    result = preflight.validate_selection(selection)
    assert result.ok, result.render_text()
    assert [row["sample_id"] for row in result.assignments] == ["sample ", "sample"]
    generated = demultiplexing.plan(selection)
    assert [row["sample_id"] for row in generated["rows"]] == ["sample ", "sample"]
    data = list(csv.reader(io.StringIO(generated["sample_sheet_content"])))
    header = next(row for row in data if "Sample_ID" in row)
    assert [row[header.index("Sample_ID")] for row in data[-2:]] == ["sample ", "sample"]


@pytest.mark.parametrize(
    "column", ["MaskShortReads", "TrimUMI", "MinimumTrimmedReadLength", "Read1UMIStartFromCycle"]
)
def test_unsupported_populated_conversion_column_blocks_preparation(tmp_path, column):
    source, output = tmp_path / "source", tmp_path / "output"
    selection = inputs(source)
    add_column(selection.pairs[0].sample_sheet, column, "12")
    output.mkdir()
    successful = output / "retained-result.7za"
    successful.write_bytes(b"prior successful delivery")
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in source.iterdir()}
    with pytest.raises(
        preflight.PreflightValidationError, match="unsupported populated demultiplexer data columns"
    ) as error:
        preflight.require_valid_inputs(source, output)
    assert column in str(error.value)
    assert {
        path: (path.read_bytes(), path.stat().st_mtime_ns) for path in source.iterdir()
    } == before
    assert successful.read_bytes() == b"prior successful delivery"


def test_supported_annotations_remain_in_provenance_without_conversion_semantics(tmp_path):
    selection = inputs(tmp_path / "source")
    add_column(selection.pairs[0].sample_sheet, "Description", "  original annotation  ")
    add_column(selection.pairs[0].sample_sheet, "MaskShortReads", "")
    generated = demultiplexing.plan(selection)
    assert generated["rows"][0]["effective_parameters"]["Description"] == "  original annotation  "
    assert "Description" not in generated["sample_sheet_content"]
    assert "MaskShortReads" not in generated["sample_sheet_content"]


def test_forced_bcl2fastq_cannot_validate_a_different_command_line_policy(tmp_path, monkeypatch):
    source = tmp_path / "source"
    selection = inputs(source)
    monkeypatch.setenv("FORCE_BCL2FASTQ", "1")
    result = preflight.validate_selection(selection)
    assert not result.ok
    assert "FORCE_BCL2FASTQ is unsupported for multiple project input pairs" in result.render_text()
    with pytest.raises(preflight.PreflightValidationError, match="Unset FORCE_BCL2FASTQ"):
        preflight.require_valid_inputs(source, tmp_path / "output")
    # The existing legacy backend remains available through its established path.
    legacy = preflight.InputSelection(selection.sample_sheet, selection.submission_form)
    assert preflight.validate_selection(legacy).ok
    assert demultiplexing.plan(SimpleNamespace(pairs=()), force_bcl2fastq=True) is None
