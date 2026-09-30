"""Structured search and the shared installed flowcell-manager/fm interface."""

import subprocess
import sys

from pathlib import Path

import flowcell_manager.flowcell_manager as manager
import pandas as pd
import pytest

from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.state import FlowcellStateStore, new_state

FAILED = "260925_NB501038_0281_AHL2T7AFXC"
COMPLETED = "260918_MN00686_0026_A000HCMFHF"
ARCHIVED = "260901_NB501038_0279_ARCHIVED"
LEGACY = "250901_NB501038_0200_LEGACY"
LEGACY_ARCHIVED = "240901_NB501038_0100_ARCHIVED"


@pytest.fixture
def search_inventory(tmp_path):
    config = tmp_path / "bfq.ini"
    config.write_text(f"[Paths]\nmanager_dir={tmp_path / 'manager'}\n")
    cfg = PipelineConfig.load(config)
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    for run_id, status, stage, projects in (
        (FAILED, "failed", "analysis", ["GCF-2026-043", "GCF-[literal].+", "Other"]),
        (COMPLETED, "completed", "finalization", ["GCF-2026-043", "GCF-2026-044"]),
        (ARCHIVED, "archived", "finalization", ["GCF-2026-043"]),
    ):
        state = new_state(
            run_id,
            tmp_path / "source" / run_id,
            tmp_path / "output" / run_id,
            origin="new",
            start_stage="demultiplexing",
        )
        state.update(status=status, current_stage=stage, projects=projects)
        if status == "archived":
            state["archive"]["archived_at"] = "2026-09-29T00:00:00+00:00"
        store.create(state)
    rows = [
        # Current JSON takes precedence over every historical path and project.
        ("stale-project", str(tmp_path / "output" / FAILED), "2026-09-29"),
        ("stale-project", FAILED, "0"),
        ("stale-project", str(tmp_path / "old-output" / FAILED), "0"),
        # Merge duplicate projects and path spellings into one legacy run.
        ("GCF-2026-043", LEGACY, "0"),
        ("Legacy-other", str(tmp_path / "output" / LEGACY), "0"),
        ("GCF-2026-043", str(tmp_path / "output" / LEGACY), "0"),
        ("GCF-2026-043", LEGACY_ARCHIVED, "2025-01-01"),
    ]
    pd.DataFrame(
        [(project, path, "2026-09-01", archived) for project, path, archived in rows],
        columns=manager.INVENTORY_COLUMNS,
    ).to_csv(cfg.static.paths.manager_dir / "flowcells.processed", index=False)
    return config, store


@pytest.mark.parametrize(
    ("query", "expected"),
    [
        ("GCF-2026-043", [LEGACY_ARCHIVED, LEGACY, ARCHIVED, COMPLETED, FAILED]),
        ("gCf-2026", [LEGACY_ARCHIVED, LEGACY, ARCHIVED, COMPLETED, FAILED]),
        ("GCF-2026-044", [COMPLETED]),
        (FAILED, [FAILED]),
        ("hl2t7afxc", [FAILED]),
        ("LEGACY", [LEGACY]),
        ("archived", [LEGACY_ARCHIVED, ARCHIVED]),
        ("[literal].+", [FAILED]),
        (".*", []),
        ("043,GCF", []),  # Never search across comma-joined project boundaries.
        ("old-output", []),  # Match IDs and projects, not arbitrary path components.
        ("stale-project", []),
        ("no-such-run", []),
    ],
)
def test_literal_search_matches_structured_identifiers(search_inventory, query, expected):
    assert manager.combined_list(query=query)["run_id"].tolist() == expected


def test_one_row_per_run_preserves_associated_projects_and_list_columns(search_inventory):
    listed = manager.combined_list()
    assert listed["run_id"].tolist() == sorted(
        [FAILED, COMPLETED, ARCHIVED, LEGACY, LEGACY_ARCHIVED]
    )
    found = manager.combined_list(query="GCF-2026-043")
    pd.testing.assert_frame_equal(listed, found)
    assert list(found.columns) == [
        "run_id",
        "status",
        "stage",
        "origin",
        "projects",
        "output_path",
        "archived",
    ]
    rows = found.set_index("run_id")
    assert rows.loc[FAILED, "projects"] == "GCF-2026-043,GCF-[literal].+,Other"
    assert rows.loc[LEGACY, "projects"] == "GCF-2026-043,Legacy-other"
    assert rows.loc[LEGACY_ARCHIVED, "archived"] == "2025-01-01"
    assert rows.loc[ARCHIVED, "status"] == "archived"


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
def test_filters_share_list_semantics_without_resurrecting_inventory(
    search_inventory, filters, expected
):
    for query in (None, "GCF-2026-043"):
        assert manager.combined_list(query=query, **filters)["run_id"].tolist() == expected
    # A run excluded by any filter cannot return via an obsolete project either.
    assert manager.combined_list(query="stale-project", **filters).empty


def test_search_does_not_change_records_or_require_output_directories(search_inventory):
    _config, store = search_inventory
    before = {path: path.read_bytes() for path in store.manager_dir.rglob("*") if path.is_file()}
    manager.combined_list(query="GCF-2026")
    after = {path: path.read_bytes() for path in store.manager_dir.rglob("*") if path.is_file()}
    assert after == before
    assert not Path(store.read(FAILED)["output_path"]).exists()


def test_empty_inventory_returns_no_rows(tmp_path):
    config = tmp_path / "bfq.ini"
    config.write_text(f"[Paths]\nmanager_dir={tmp_path / 'manager'}\n")
    PipelineConfig.load(config)
    assert manager.combined_list(query="GCF").empty


@pytest.mark.parametrize("executable_name", ["flowcell-manager", "fm"])
def test_installed_entry_points_search_help_and_exit_behavior(search_inventory, executable_name):
    config, _store = search_inventory
    executable = Path(sys.executable).parent / executable_name
    # Use the real installed script with isolated configuration and no repository cwd.
    driver = """
import runpy, smtplib, sys
from bcl2fastq_pipeline.config import PipelineConfig
def forbidden_smtp(*args, **kwargs):
    raise AssertionError("Read-only commands must never send mail")
smtplib.SMTP = forbidden_smtp
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

    help_result = invoke("--help")
    assert help_result.returncode == 0, help_result.stderr
    assert "search" in help_result.stdout
    assert "clean-fastqs" in help_result.stdout
    assert "retry-notifications" in help_result.stdout
    search_help = invoke("search", "--help")
    assert search_help.returncode == 0
    assert "fm search GCF-2026-043" in search_help.stdout
    assert "--status" in search_help.stdout and "--stage" in search_help.stdout

    result = invoke("search", "gcf-2026-043", "--status", "failed", "--stage", "analysis")
    assert result.returncode == 0, result.stderr
    assert FAILED in result.stdout and "GCF-[literal].+" in result.stdout
    assert COMPLETED not in result.stdout
    listed = invoke("list", "--status", "failed", "--stage", "analysis")
    assert listed.returncode == 0 and listed.stdout == result.stdout
    flowcell = invoke("search", "hL2t7AfXc")
    assert flowcell.returncode == 0 and flowcell.stdout == result.stdout

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
