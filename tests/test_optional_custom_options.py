"""Old cleaner options are optional metadata alongside existing scheduling fields."""

from decimal import Decimal

import pytest

from bcl2fastq_pipeline.config import RunContext, parse_custom_options
from test_project_scheduling import PROJECTS, paired_inputs

from bcl2fastq_pipeline import initial_holds, preflight, project_scheduling

RETIRED = {
    "RemoveHumanReads": "true",
    "SingleCell": "False",
    "TrimAdapter": "False",
    "Rerun": "true",
    "ReverseComplementIndexP5": "true",
    "ReverseComplementIndexP7": "False",
}


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("scheduling", [False, True])
def test_reader_and_run_context_accept_old_or_new_producer_options(tmp_path, legacy, scheduling):
    options = {"Libprep": "Synthetic kit", "User": "Synthetic Operator"}
    if legacy:
        options.update(RETIRED)
    if scheduling:
        options.update(Priority="0.5000000000000000000001", InitialHold="true")
    sheet = tmp_path / "SampleSheet.csv"
    original = (
        "[CustomOptions]\n"
        + "".join(f"{key},{value}\n" for key, value in options.items())
        + "[Data]\nSample_ID,Sample_Project\nSAMPLE_001,GCF-2099-001\n"
    ).encode()
    sheet.write_bytes(original)
    parsed, returned = parse_custom_options(sheet)
    context = RunContext()
    context.apply_custom(parsed, returned, tmp_path / "submission.xlsx")
    assert context.custom == options
    assert context.rerun is legacy  # Retain the existing compatibility flag only.
    assert "Index2Orientation" not in context.custom
    priority = project_scheduling.resolve_priority(context.custom)
    assert Decimal(priority["priority"]) == Decimal(options.get("Priority", "1.0"))
    assert initial_holds.resolve(context.custom)["held"] is scheduling
    assert sheet.read_bytes() == original


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("scheduling", [False, True])
def test_complete_preflight_tolerates_missing_retired_options(tmp_path, legacy, scheduling):
    selection = paired_inputs(tmp_path / "inputs", ["0.5" if scheduling else None] * 2)
    options = dict(RETIRED) if legacy else {}
    if scheduling:
        options["InitialHold"] = "true"
    for pair in selection.pairs:
        pair.sample_sheet.write_text(
            pair.sample_sheet.read_text().replace(
                "[Data]", "".join(f"{key},{value}\n" for key, value in options.items()) + "[Data]"
            )
        )
    original = {pair.sample_sheet: pair.sample_sheet.read_bytes() for pair in selection.pairs}
    assert len(original) == 2
    result = preflight.validate_selection(selection)
    assert result.ok, result.render_text()
    schedule = result.summary["project_schedule"]["projects"][PROJECTS[0]]
    assert Decimal(schedule["priority"]) == Decimal("0.5" if scheduling else "1.0")
    # Shared parsing may coerce true to bool; the scheduling consumers support it.
    assert (
        initial_holds.resolve(result.project_results[PROJECTS[0]].custom_options)["held"]
        is scheduling
    )
    assert {sheet: sheet.read_bytes() for sheet in original} == original
