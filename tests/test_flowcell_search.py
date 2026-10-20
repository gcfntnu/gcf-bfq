"""Project inventory and the shared installed flowcell-manager/fm interface."""

import copy
import subprocess
import sys

from pathlib import Path

import flowcell_manager.flowcell_manager as manager
import pandas as pd
import pytest

from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.state import FlowcellStateStore, new_state

FAILED = "991003_TESTSEQ_0003_AFAILED"
COMPLETED = "991002_TESTSEQ_0002_ACOMPLETED"
ARCHIVED = "991001_TESTSEQ_0001_AARCHIVED"
LEGACY = "980102_TESTSEQ_0002_ALEGACY"
LEGACY_ARCHIVED = "980101_TESTSEQ_0001_AHISTORY"
PROJECT = "GCF-2099-001"
OTHER_PROJECT = "GCF-2099-002"
LITERAL_PROJECT = "GCF-[literal].+"
DELIVERED_AT = "2099-10-01T12:00:00+00:00"


@pytest.fixture
def search_inventory(tmp_path):
    config = tmp_path / "bfq.ini"
    config.write_text(f"[Paths]\nmanager_dir={tmp_path / 'manager'}\n")
    cfg = PipelineConfig.load(config)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    for run_id, status, projects in (
        (FAILED, "failed", [PROJECT, OTHER_PROJECT, LITERAL_PROJECT, "Other"]),
        (COMPLETED, "completed", [PROJECT, OTHER_PROJECT]),
        (ARCHIVED, "archived", [PROJECT]),
    ):
        state = new_state(
            run_id,
            tmp_path / "source" / run_id,
            tmp_path / "output" / run_id,
            origin="new",
            start_stage="demultiplexing",
        )
        state.update(status=status, current_stage="finalization", projects=projects)
        if status == "failed":
            state["current_stage"] = "analysis"
            state["stages"]["analysis"]["status"] = "failed"
            state["project_progress"] = {
                PROJECT: {
                    "analysis": {
                        "status": "failed",
                        "attempt": 2,
                        "generation": "failed-generation",
                        "scheduling": {"priority": "0.25", "priority_source": "samplesheet"},
                    }
                },
                OTHER_PROJECT: {
                    stage: {
                        "status": "completed",
                        "attempt": 1,
                        "generation": "retained-generation",
                        "completed_at": DELIVERED_AT,
                    }
                    for stage in ("analysis", "finalization")
                },
            }
            state["project_completion"] = {
                OTHER_PROJECT: {
                    "attempt": 1,
                    "generation": "retained-generation",
                    "completed_at": DELIVERED_AT,
                }
            }
        else:
            state["completed_at"] = DELIVERED_AT
            for detail in state["stages"].values():
                detail.update(status="completed", completed_at=DELIVERED_AT)
        if status == "archived":
            state["archive"]["archived_at"] = "2099-10-02T00:00:00+00:00"
        store.create(state)
    rows = [
        # JSON suppresses all historical rows before search or status filtering.
        ("stale-project", str(tmp_path / "output" / FAILED), "2099-10-02"),
        ("stale-project", FAILED, "0"),
        ("stale-project", str(tmp_path / "old-output" / FAILED), "0"),
        (PROJECT, LEGACY, "0"),
        ("Legacy-other", str(tmp_path / "output" / LEGACY), "0"),
        (PROJECT, str(tmp_path / "output" / LEGACY), "0"),
        (PROJECT, LEGACY_ARCHIVED, "2098-01-01"),
        (PROJECT, str(tmp_path / "output" / LEGACY_ARCHIVED), "2098-01-01"),
        # Archiving one legacy project must not archive its neighbour's row.
        (OTHER_PROJECT, LEGACY_ARCHIVED, "0"),
    ]
    pd.DataFrame(
        [(project, path, "2099-09-01", archived) for project, path, archived in rows],
        columns=manager.INVENTORY_COLUMNS,
    ).to_csv(cfg.static.paths.manager_dir / "flowcells.processed", index=False)
    return config, store


@pytest.mark.parametrize(
    ("query", "expected"),
    [
        (
            PROJECT,
            [
                (run_id, PROJECT)
                for run_id in (LEGACY_ARCHIVED, LEGACY, ARCHIVED, COMPLETED, FAILED)
            ],
        ),
        (
            "gCf-2099",
            [
                (LEGACY_ARCHIVED, PROJECT),
                (LEGACY_ARCHIVED, OTHER_PROJECT),
                (LEGACY, PROJECT),
                (ARCHIVED, PROJECT),
                (COMPLETED, PROJECT),
                (COMPLETED, OTHER_PROJECT),
                (FAILED, PROJECT),
                (FAILED, OTHER_PROJECT),
            ],
        ),
        (
            OTHER_PROJECT,
            [(LEGACY_ARCHIVED, OTHER_PROJECT), (COMPLETED, OTHER_PROJECT), (FAILED, OTHER_PROJECT)],
        ),
        (
            "afailed",
            [
                (FAILED, PROJECT),
                (FAILED, OTHER_PROJECT),
                (FAILED, LITERAL_PROJECT),
                (FAILED, "Other"),
            ],
        ),
        ("[literal].+", [(FAILED, LITERAL_PROJECT)]),
        (".*", []),
        ("001,GCF", []),
        ("old-output", []),
        ("stale-project", []),
        ("no-such-run", []),
    ],
)
def test_literal_search_matches_only_project_rows_or_all_rows_of_a_run(
    search_inventory, query, expected
):
    found = manager.combined_list(query=query)
    assert list(found[["run_id", "project"]].itertuples(index=False, name=None)) == expected


def test_one_row_per_project_and_run_has_explicit_parent_context(search_inventory):
    listed = manager.combined_list()
    assert not listed.duplicated(["project", "run_id"]).any()
    assert len(listed) == 11
    assert list(listed.columns) == [
        "project",
        "run_id",
        "status",
        "stage",
        "analysis",
        "delivery",
        "last_delivery",
        "priority",
        "priority_source",
        "parent_status",
        "parent_stage",
        "origin",
        "output_path",
        "archived",
        "analysis_generation",
        "last_delivery_generation",
    ]
    expected = listed.sort_values(["run_id", "project"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(listed.reset_index(drop=True), expected)
    rows = listed.set_index(["run_id", "project"])
    assert rows.loc[(FAILED, PROJECT), "status"] == "failed"
    retained = rows.loc[(FAILED, OTHER_PROJECT)]
    assert retained["status"] == "completed"
    assert retained["analysis"] == retained["delivery"] == "completed"
    assert retained["parent_status"] == "failed" and retained["parent_stage"] == "analysis"
    assert retained["last_delivery"] == DELIVERED_AT
    assert (
        retained["analysis_generation"]
        == retained["last_delivery_generation"]
        == "retained-generation"
    )
    assert rows.loc[(FAILED, PROJECT), "priority"] == "0.25"
    assert rows.loc[(FAILED, PROJECT), "priority_source"] == "samplesheet"
    assert rows.loc[(LEGACY_ARCHIVED, PROJECT), "status"] == "archived"
    assert rows.loc[(LEGACY_ARCHIVED, OTHER_PROJECT), "status"] == "completed"
    assert rows.loc[(LEGACY_ARCHIVED, OTHER_PROJECT), "archived"] == "0"


@pytest.mark.parametrize(
    ("filters", "expected"),
    [
        ({"status": "failed"}, [FAILED]),
        ({"stage": "analysis"}, [FAILED]),
        ({"status": "completed"}, [LEGACY, COMPLETED]),
        ({"status": "archived"}, [LEGACY_ARCHIVED, ARCHIVED]),
        ({"status": "legacy"}, [LEGACY_ARCHIVED, LEGACY]),
        ({"stage": "legacy"}, [LEGACY_ARCHIVED, LEGACY]),
        ({"status": "completed", "stage": "legacy"}, [LEGACY]),
        ({"stage": "finalization"}, [ARCHIVED, COMPLETED]),
        ({"status": "failed", "stage": "finalization"}, []),
        ({"status": "unknown"}, []),
        ({"stage": "unknown"}, []),
    ],
)
def test_project_filters_do_not_resurrect_historical_inventory(search_inventory, filters, expected):
    found = manager.combined_list(query=PROJECT, **filters)
    assert found["run_id"].tolist() == expected
    assert manager.combined_list(query="stale-project", **filters).empty
    # Query filtering composes with the same project-level list filters.
    listed = manager.combined_list(**filters)
    pd.testing.assert_frame_equal(
        found.reset_index(drop=True),
        listed[listed["project"] == PROJECT].reset_index(drop=True),
    )


def test_old_failed_state_preserves_unknown_project_outcomes(search_inventory):
    rows = manager.combined_list(query=FAILED).set_index("project")
    for project in (LITERAL_PROJECT, "Other"):
        assert rows.loc[project, "status"] == "unavailable"
        assert rows.loc[project, "analysis"] == "unavailable"
        assert rows.loc[project, "parent_status"] == "failed"
        assert rows.loc[project, "priority"] == "unavailable"
        assert rows.loc[project, "priority_source"] == "unavailable"
    failed = manager.combined_list(query=FAILED, status="failed")
    assert failed["project"].tolist() == [PROJECT]


def test_search_does_not_change_records_or_require_output_directories(
    search_inventory, monkeypatch
):
    _config, store = search_inventory
    before = {path: path.read_bytes() for path in store.manager_dir.rglob("*") if path.is_file()}
    state_before = copy.deepcopy(store.read(FAILED))

    def forbidden_input_read(*args, **kwargs):
        raise AssertionError("Inventory must use persisted project metadata")

    monkeypatch.setattr(manager.preflight, "validate_selection", forbidden_input_read)
    monkeypatch.setattr(manager.preflight, "select_run_inputs", forbidden_input_read)
    monkeypatch.setattr(FlowcellStateStore, "ensure_writable", forbidden_input_read)
    manager.combined_list(query="GCF-2099")
    after = {path: path.read_bytes() for path in store.manager_dir.rglob("*") if path.is_file()}
    assert after == before
    assert store.read(FAILED) == state_before
    assert not Path(state_before["output_path"]).exists()


def test_unidentified_canonical_run_is_kept_without_inventing_a_project(search_inventory, tmp_path):
    _config, store = search_inventory
    run_id = "991004_TESTSEQ_0004_AUNIDENTIFIED"
    store.create(
        new_state(
            run_id,
            tmp_path / "source" / run_id,
            tmp_path / "output" / run_id,
            origin="new",
            start_stage="demultiplexing",
        )
    )
    found = manager.combined_list(query=run_id)
    assert len(found) == 1 and found.iloc[0]["project"] == "<unidentified>"
    assert found.iloc[0]["status"] == "queued"
    assert found.iloc[0]["stage"] == "demultiplexing"
    assert manager.combined_list(query=PROJECT)["project"].eq(PROJECT).all()


def test_empty_inventory_returns_no_rows_without_creating_manager(tmp_path):
    config = tmp_path / "bfq.ini"
    manager_dir = tmp_path / "manager"
    config.write_text(f"[Paths]\nmanager_dir={manager_dir}\n")
    PipelineConfig.load(config)
    assert manager.combined_list(query="GCF").empty
    assert not manager_dir.exists()


@pytest.mark.parametrize("executable_name", ["flowcell-manager", "fm"])
def test_installed_entry_points_search_help_and_exit_behavior(search_inventory, executable_name):
    config, store = search_inventory
    executable = Path(sys.executable).parent / executable_name
    driver = """
import runpy, smtplib, sys
from bcl2fastq_pipeline.config import PipelineConfig
def forbidden_smtp(*args, **kwargs):
    raise AssertionError("Read-only commands must never send mail")
smtplib.SMTP = forbidden_smtp
smtplib.SMTP_SSL = forbidden_smtp
PipelineConfig.load(sys.argv[1])
sys.argv = sys.argv[2:]
runpy.run_path(sys.argv[0], run_name="__main__")
"""

    def invoke(*args):
        return subprocess.run(
            [sys.executable, "-c", driver, str(config), str(executable), *args],
            check=False,
            cwd="/proc",
            capture_output=True,
            text=True,
        )

    before = {path: path.read_bytes() for path in store.manager_dir.rglob("*") if path.is_file()}
    help_result = invoke("--help")
    assert help_result.returncode == 0, help_result.stderr
    assert "search" in help_result.stdout
    assert "clean-fastqs" in help_result.stdout
    assert "retry-notifications" in help_result.stdout
    search_help = invoke("search", "--help")
    assert search_help.returncode == 0
    assert f"fm search {PROJECT}" in search_help.stdout
    assert "--status" in search_help.stdout and "--stage" in search_help.stdout
    assert "project" in search_help.stdout.casefold()

    result = invoke("search", PROJECT.casefold(), "--status", "failed", "--stage", "analysis")
    assert result.returncode == 0, result.stderr
    assert FAILED in result.stdout and PROJECT in result.stdout
    assert OTHER_PROJECT not in result.stdout and LITERAL_PROJECT not in result.stdout
    listed = invoke("list", "--status", "failed", "--stage", "analysis")
    assert listed.returncode == 0 and listed.stdout == result.stdout
    flowcell = invoke("search", "aFaIlEd", "--status", "failed", "--stage", "analysis")
    assert flowcell.returncode == 0 and flowcell.stdout == result.stdout
    all_projects = invoke("search", FAILED)
    assert all_projects.returncode == 0
    assert OTHER_PROJECT in all_projects.stdout and LITERAL_PROJECT in all_projects.stdout
    assert "parent_status" in all_projects.stdout and "last_delivery" in all_projects.stdout

    no_matches = invoke("search", "not-present")
    assert no_matches.returncode == 0
    assert no_matches.stdout == "No matching flowcells.\n" and not no_matches.stderr
    for query in ("", " \t\n"):
        invalid = invoke("search", query)
        assert invalid.returncode == 2
        assert "query must not be empty or whitespace-only" in invalid.stderr
    missing = invoke("search")
    assert missing.returncode == 2 and "QUERY" in missing.stderr
    state_error = invoke("status", "missing-run")
    assert state_error.returncode == 1 and "No such flowcell" in state_error.stderr
    assert "Traceback" not in state_error.stderr
    after = {path: path.read_bytes() for path in store.manager_dir.rglob("*") if path.is_file()}
    assert after == before
