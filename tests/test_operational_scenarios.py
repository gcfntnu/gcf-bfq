"""Small operator journeys; detailed edge cases remain in the existing test suite."""

import gzip
import json
import smtplib
import tarfile

import pytest

from openpyxl import load_workbook
from support.operational import (
    INPUT_NAMES,
    LABEL,
    PROJECT,
    RUN_ID,
    deterministic_inputs,
    make_scenario,
    paired_fastqs,
    smtp_double,
    snapshot,
)

from bcl2fastq_pipeline import afterFastq, preflight


@pytest.fixture
def scenario(tmp_path, monkeypatch):
    return make_scenario(tmp_path / "workspace", monkeypatch)


def initialize(scenario):
    scenario.command("fm", "initialize", RUN_ID, "--from", "analysis", "--force")
    state = scenario.store.read(RUN_ID)
    assert state["status"] == "queued"
    assert state["current_stage"] == "analysis"
    assert state["origin"] == "restored_legacy_fastq"
    assert state["stages"]["demultiplexing"]["status"] == "completed"
    assert state["attempt"] == 0


def protected(scenario):
    paths = [scenario.output / name for name in INPUT_NAMES]
    paths += sorted((scenario.output / PROJECT).glob("*.fastq.gz"))
    paths += [scenario.manifest, scenario.output / "Stats/Stats.json"]
    return {str(path): (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}


def processing(state):
    return {
        key: value
        for key, value in state.items()
        if key not in {"updated_at", "delivery_notifications"}
    }


def test_synthetic_fixture_pairing_metadata_and_installed_validation(scenario):
    reads = []
    for mate in (1, 2):
        path = scenario.output / PROJECT / f"sample_R{mate}.fastq.gz"
        lines = gzip.decompress(path.read_bytes()).decode().splitlines()
        assert len(lines) == 8
        records = [lines[index : index + 4] for index in range(0, len(lines), 4)]
        assert all(record[0].split()[1].startswith(f"{mate}:") for record in records)
        assert all(
            record[2] == "+" and len(record[1]) == len(record[3]) == 12 for record in records
        )
        reads.append([record[0].split()[0] for record in records])
    assert reads[0] == reads[1] == ["@synthetic:0", "@synthetic:1"]
    workbook = load_workbook(scenario.output / INPUT_NAMES[1], read_only=True)
    assert workbook["Sample-Submission-Form"]["A15"].value == "Unique Sample ID"
    assert workbook["Sample-Submission-Form"]["A16"].value == "sample"
    assert workbook["INFO (GCF-lab only)"]["A2"].value == "sample"
    workbook.close()
    result = preflight.validate_selection(
        preflight.select_run_inputs(scenario.source, scenario.output)
    )
    assert result.ok
    assert result.to_dict()["summary"]["planned_sample_count"] == 1
    before = (
        snapshot(scenario.output),
        snapshot(scenario.source),
        snapshot(scenario.cfg.static.paths.manager_dir),
    )
    first = scenario.command("fm", "validate", RUN_ID)
    second = scenario.command("flowcell-manager", "validate", str(scenario.output))
    assert first.stdout == second.stdout
    assert "preflight PASSED" in first.stdout
    assert (
        snapshot(scenario.output),
        snapshot(scenario.source),
        snapshot(scenario.cfg.static.paths.manager_dir),
    ) == before
    assert not scenario.store.exists(RUN_ID)

    duplicate = scenario.root / "fixture-copy"
    deterministic_inputs(duplicate, "-curated")
    paired_fastqs(duplicate / PROJECT)
    for relative in (
        *INPUT_NAMES,
        f"{PROJECT}/sample_R1.fastq.gz",
        f"{PROJECT}/sample_R2.fastq.gz",
    ):
        assert (scenario.output / relative).read_bytes() == (duplicate / relative).read_bytes()


def test_invalid_preflight_blocks_initialization_rerun_and_execution(scenario):
    sheet = scenario.output / "SampleSheet.csv"
    original = sheet.read_text()
    # Valid source inputs must never silently replace malformed curated inputs.
    sheet.write_text(original.replace("sample,GCF-", "missing,GCF-"))
    before = snapshot(scenario.output), snapshot(scenario.source)
    invalid = scenario.command("fm", "validate", RUN_ID, expected=1)
    assert "sample.missing_metadata" in invalid.stderr + invalid.stdout
    scenario.command(
        "flowcell-manager", "initialize", RUN_ID, "--from", "analysis", "--force", expected=1
    )
    assert not scenario.store.exists(RUN_ID)
    assert (snapshot(scenario.output), snapshot(scenario.source)) == before

    sheet.write_text(original)
    initialize(scenario)
    assert scenario.run_attempt()["status"] == "completed"
    sheet.write_text(original.replace("sample,GCF-", "missing,GCF-"))
    before = (
        snapshot(scenario.output),
        snapshot(scenario.source),
        scenario.store.state_path(RUN_ID).read_bytes(),
    )
    scenario.command("fm", "rerun", RUN_ID, "--from", "analysis", "--force", expected=1)
    assert (
        snapshot(scenario.output),
        snapshot(scenario.source),
        scenario.store.state_path(RUN_ID).read_bytes(),
    ) == before

    sheet.write_text(original)
    scenario.command("fm", "rerun", RUN_ID, "--from", "analysis", "--force")
    # Mutation after preparation must be revalidated under the execution lease.
    sheet.write_text(original.replace("sample,GCF-", "missing,GCF-"))
    before_output = snapshot(scenario.output)
    before_doubles = (scenario.root / "logs/external-doubles.jsonl").read_bytes()
    failed = scenario.run_attempt()
    assert failed["status"] == "failed"
    assert "Input preflight failed" in failed["last_error"]["summary"]
    report = failed["stages"]["analysis"]["metadata"]["input_preflight"]
    assert report["validator"]["version"]
    assert report["status"] == "failed"
    assert report == json.loads(
        (scenario.root / "reports" / f"{RUN_ID}.input-preflight.json").read_text()
    )
    assert "sample.missing_metadata" in json.dumps(report)
    assert snapshot(scenario.output) == before_output
    later_doubles = (scenario.root / "logs/external-doubles.jsonl").read_bytes()[
        len(before_doubles) :
    ]
    assert [json.loads(line)["boundary"] for line in later_doubles.splitlines()] == [
        "installed-workflow-provenance"
    ]
    assert "Input preflight failed" in (scenario.root / "reports" / f"{RUN_ID}.error").read_text()


def test_restored_analysis_rerun_preserves_inputs_checksums_and_installed_queries(scenario):
    curated = {name: (scenario.output / name).read_bytes() for name in INPUT_NAMES}
    initialize(scenario)
    completed = scenario.run_attempt()
    assert completed["status"] == "completed"
    assert completed["projects"] == [PROJECT]
    assert completed["stages"]["analysis"]["metadata"]["input_preflight"]["status"] == "passed"
    assert completed["stages"]["analysis"]["metadata"]["workflow"] == "synthetic"
    assert all(entry["status"] == "sent" for entry in completed["delivery_notifications"])
    assert {name: (scenario.output / name).read_bytes() for name in INPUT_NAMES} == curated
    assert "-curated" in (scenario.output / "SampleSheet.csv").read_text()
    for line in scenario.manifest.read_text().splitlines():
        digest, filename = line.split("  ", 1)
        assert digest == afterFastq.file_md5(scenario.output / filename)
    assert len(scenario.manifest.read_text().splitlines()) == 2
    assert LABEL in scenario.archive.read_text()
    with tarfile.open(scenario.retained) as archive:
        assert LABEL.encode() in archive.extractfile("Snakefile").read()
        assert not any(name == "data" or name.startswith("data/") for name in archive.getnames())

    for arguments in (
        ("show", RUN_ID),
        ("status", RUN_ID),
        ("search", PROJECT),
        ("search", "synthetic"),
        ("list", "--status", "completed"),
    ):
        first = scenario.command("fm", *arguments)
        second = scenario.command("flowcell-manager", *arguments)
        assert first.stdout == second.stdout
        assert RUN_ID in first.stdout
    assert json.loads(scenario.command("fm", "show", RUN_ID).stdout)["attempt"] == 1
    protected_before = protected(scenario)
    snapshot_before = scenario.retained.read_bytes()
    before = (
        snapshot(scenario.output),
        snapshot(scenario.root / "scratch"),
        scenario.store.state_path(RUN_ID).read_bytes(),
    )
    preview = scenario.command("fm", "rerun", RUN_ID, "--from", "analysis", "--dry-run")
    assert str(scenario.archive) in preview.stdout
    assert (
        snapshot(scenario.output),
        snapshot(scenario.root / "scratch"),
        scenario.store.state_path(RUN_ID).read_bytes(),
    ) == before
    scenario.command(
        "flowcell-manager", "rerun", str(scenario.output), "--from", "analysis", "--force"
    )
    assert scenario.store.read(RUN_ID)["status"] == "queued"
    assert protected(scenario) == protected_before
    assert not scenario.archive.exists()
    assert not list((scenario.root / "scratch").iterdir())
    assert scenario.retained.read_bytes() == snapshot_before
    assert scenario.run_attempt()["status"] == "completed"
    assert protected(scenario) == protected_before

    # A later restart boundary retains the completed analysis and its snapshot.
    previous = scenario.store.read(RUN_ID)
    previous_work = snapshot(scenario.root / "scratch")
    previous_snapshot = snapshot(scenario.output / "provenance")
    scenario.command("fm", "rerun", RUN_ID, "--from", "finalization", "--force")
    assert scenario.store.read(RUN_ID)["stages"]["analysis"] == previous["stages"]["analysis"]
    assert snapshot(scenario.root / "scratch") == previous_work
    assert scenario.run_attempt()["status"] == "completed"
    assert snapshot(scenario.root / "scratch") == previous_work
    assert snapshot(scenario.output / "provenance") == previous_snapshot
    assert protected(scenario) == protected_before


def test_controlled_processing_failure_snapshot_recovery_and_notification_independence(
    scenario, monkeypatch
):
    initialize(scenario)
    assert scenario.run_attempt()["status"] == "completed"
    first_snapshot = scenario.retained.read_bytes()
    first_identity = scenario.store.read(RUN_ID)["analysis_snapshots"][PROJECT]
    inputs = protected(scenario)
    scenario.command("fm", "rerun", RUN_ID, "--from", "analysis", "--force")
    scenario.fail_workflow = True
    failed = scenario.run_attempt()
    assert failed["status"] == "failed"
    assert failed["current_stage"] == "analysis"
    assert failed["stages"]["analysis"]["status"] == "failed"
    error_report = (scenario.root / "reports" / f"{RUN_ID}.error").read_text()
    assert "controlled synthetic workflow failure" in error_report
    assert "exit status 23" in failed["last_error"]["summary"]
    assert failed["analysis_snapshots"][PROJECT] == first_identity
    assert scenario.retained.read_bytes() == first_snapshot
    assert protected(scenario) == inputs
    # Preserve an explicit failing-attempt diagnostic when recovery updates reports.
    (scenario.root / "logs" / "controlled-failure.error").write_text(error_report)
    scenario.command("flowcell-manager", "status", RUN_ID)
    scenario.command("fm", "rerun", RUN_ID, "--from", "analysis", "--force")
    scenario.fail_workflow = False
    monkeypatch.setattr(smtplib, "SMTP", smtp_double(scenario.root, "fail"))
    completed = scenario.run_attempt()
    assert completed["status"] == "completed"
    assert completed["attempt"] == 3
    assert completed["last_error"] is None
    assert all(stage["status"] == "completed" for stage in completed["stages"].values())
    active_notifications = [
        entry for entry in completed["delivery_notifications"] if entry["status"] != "superseded"
    ]
    assert [entry["kind"] for entry in active_notifications] == ["processed", "finalized"]
    assert all(entry["status"] == "failed" for entry in active_notifications)
    assert all("controlled synthetic SMTP" in entry["last_error"] for entry in active_notifications)
    assert completed["analysis_snapshots"][PROJECT]["analysis_id"] != first_identity["analysis_id"]
    assert protected(scenario) == inputs
    before = snapshot(scenario.output), snapshot(scenario.root / "scratch")
    saved_processing = processing(completed)
    scenario.command("fm", "retry-notifications", RUN_ID, smtp="fail", expected=1)
    assert processing(scenario.store.read(RUN_ID)) == saved_processing
    scenario.command("flowcell-manager", "retry-notifications", RUN_ID)
    retried = scenario.store.read(RUN_ID)
    assert [
        entry["status"]
        for entry in retried["delivery_notifications"]
        if entry["status"] != "superseded"
    ] == ["sent", "sent"]
    assert processing(retried) == saved_processing
    assert (snapshot(scenario.output), snapshot(scenario.root / "scratch")) == before
    before_mail = sorted(path.name for path in (scenario.root / "logs").glob("double-mail-*.eml"))
    scenario.command("fm", "retry-notifications", RUN_ID)
    assert (
        sorted(path.name for path in (scenario.root / "logs").glob("double-mail-*.eml"))
        == before_mail
    )
    with tarfile.open(scenario.retained) as archive:
        assert (
            archive.pax_headers["bfq.analysis_id"]
            == retried["analysis_snapshots"][PROJECT]["analysis_id"]
        )
