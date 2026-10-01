"""Exercise the installed entry point, including script/package shadowing."""

import json
import subprocess
import sys

from pathlib import Path

import pytest

from bcl2fastq_pipeline.state import FlowcellStateStore, new_state


@pytest.mark.parametrize("args", [["--help"], ["retry-notifications", "--help"]])
def test_installed_manager_help_does_not_import_configmaker(args):
    executable = Path(sys.executable).parent / "flowcell-manager"
    result = subprocess.run(
        [str(executable), *args], check=False, cwd="/proc", capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert "retry-notifications" in result.stdout
    assert "Traceback" not in result.stderr


def test_installed_retry_and_status_from_unwritable_directory(tmp_path):
    run_id = "260918_MN00686_0026_A000HCMFHF"
    manager = tmp_path / "manager"
    config = tmp_path / "bfq.ini"
    config.write_text(
        f"[Paths]\nmanager_dir={manager}\n[Email]\nhost=relay.example.test\n"
        "from_address=bfq@example.test\nerror_to=operator@example.test\n"
    )
    store = FlowcellStateStore(manager)
    store.create(
        new_state(
            run_id,
            tmp_path / "absent-source",
            tmp_path / "absent-output",
            origin="restored_legacy_fastq",
            start_stage="finalization",
        )
    )
    store.begin_attempt(run_id)
    store.complete_run(
        run_id,
        ["GCF-2026-001"],
        notification={
            "kind": "finalized",
            "payload": {
                "run_id": run_id,
                "projects": ["GCF-2026-001"],
                "finalize_time": "0:10:00",
                "run_time": "1:00:00",
                "message": "",
            },
        },
    )
    executable = Path(sys.executable).parent / "flowcell-manager"
    # Load an isolated test config, then run the actual installed console script.
    driver = """
import runpy, smtplib, sys
from bcl2fastq_pipeline.config import PipelineConfig
class FakeSMTP:
    def __init__(self, *args, **kwargs): pass
    def send_message(self, message, **kwargs):
        print("TEST_DELIVERY", message["Subject"])
        return {}
    def quit(self): pass
smtplib.SMTP = FakeSMTP
PipelineConfig.load(sys.argv[1])
sys.argv = sys.argv[2:]
runpy.run_path(sys.argv[0], run_name="__main__")
assert "configmaker.configmaker" not in sys.modules
"""

    def command(*args):
        return subprocess.run(
            [sys.executable, "-c", driver, str(config), str(executable), *args],
            check=False,
            cwd="/proc",
            capture_output=True,
            text=True,
        )

    status = command("status", run_id)
    assert status.returncode == 0, status.stderr
    assert "finalized:1: pending" in status.stdout
    sent = command("retry-notifications", run_id, "--kind", "finalized")
    assert sent.returncode == 0, sent.stderr
    assert "TEST_DELIVERY" in sent.stdout
    repeat = command("retry-notifications", run_id)
    assert repeat.returncode == 0, repeat.stderr
    assert "TEST_DELIVERY" not in repeat.stdout
    assert store.read(run_id)["status"] == "completed"
    assert not (tmp_path / "absent-source").exists()
    assert not (tmp_path / "absent-output").exists()

    # A real command failure reports nonzero without a Python traceback.
    before = json.dumps(store.read(run_id), sort_keys=True)
    invalid = command("retry-notifications", run_id, "--kind", "processed")
    assert invalid.returncode == 1
    assert "No retained notifications" in invalid.stderr
    assert "Traceback" not in invalid.stderr
    assert json.dumps(store.read(run_id), sort_keys=True) == before
