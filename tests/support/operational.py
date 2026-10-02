"""Synthetic fixture and declared external doubles for the existing BFQ test suite.

Only tests import this module. No production config switch or daemon simulation
is installed. Real manager commands receive a test singleton at Python startup.
"""

import gzip
import hashlib
import json
import logging
import os
import re
import smtplib
import subprocess
import sys
import zipfile

from configparser import ConfigParser
from dataclasses import dataclass
from pathlib import Path

import yaml

from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.state import FlowcellStateStore

from bcl2fastq_pipeline import afterFastq, cli, state, workflow_config

RUN_ID = "260101_SYNTHETIC_0001_ASYNTHETIC"
PROJECT = "GCF-2026-001"
DATE = RUN_ID.split("_")[0]
INPUT_NAMES = ("SampleSheet.csv", "Sample-Submission-Form.xlsx")
LABEL = "SYNTHETIC EXTERNAL DOUBLE: not a scientific result"


def record(root, boundary, **details):
    with (Path(root) / "logs" / "external-doubles.jsonl").open("a") as handle:
        handle.write(json.dumps({"boundary": boundary, **details}, sort_keys=True) + "\n")


def smtp_double(root, mode="accept"):
    """Keep real message composition/delivery state, replace only the SMTP relay."""
    if mode not in {"accept", "fail"}:
        raise ValueError(f"Unknown scenario SMTP mode: {mode}")

    class SMTPDouble:
        def __init__(self, host, **_kwargs):
            assert host == "smtp.invalid"
            record(root, "SMTP", mode=mode)
            if mode == "fail":
                raise OSError("controlled synthetic SMTP connection failure")

        def send_message(self, message, **kwargs):
            assert all(address.endswith("@example.test") for address in kwargs["to_addrs"])
            name = hashlib.sha256(message["Message-ID"].encode()).hexdigest()[:16]
            (Path(root) / "logs" / f"double-mail-{name}.eml").write_bytes(message.as_bytes())
            return {}

        def quit(self):
            pass

    return SMTPDouble


def scenario_versions(root):
    """Do not inspect a developer server's unrelated /opt workflow checkout."""
    record(root, "installed-workflow-provenance", revision="synthetic-unexecuted")
    return {"bfq": state.package_version(), "gcf_workflows": "synthetic-unexecuted"}


def configure(root):
    """Called by the guarded subprocess startup hook, never by installed BFQ code."""
    root = Path(root).resolve(strict=True)
    assert getattr(sys, "_bfq_external_guard", False), "Missing scenario subprocess guard"
    assert os.environ.get("BFQ_TEST_SCENARIOS") == "1", "Missing strict scenario mode"
    assert (root / "fixture-origin.txt").read_text().startswith("BFQ synthetic fixture v1")
    cfg = PipelineConfig.load(root / "config" / "bcl2fastq.ini")
    for path in vars(cfg.static.paths).values():
        assert path.resolve().is_relative_to(root), f"Scenario path escaped workspace: {path}"
    assert Path(os.environ["TMPDIR"]).resolve() == root / "scratch"
    state.collect_versions = lambda: scenario_versions(root)
    smtplib.SMTP = smtp_double(root, os.environ.get("BFQ_TEST_SMTP_MODE", "accept"))
    smtplib.SMTP_SSL = smtplib.SMTP
    sys._bfq_scenario_configured = str(root)
    return cfg


def deterministic_inputs(directory, suffix):
    """Use the existing metadata fixture; normalize XLSX volatile timestamps."""
    from support.fixtures import write_inputs  # noqa: PLC0415

    write_inputs(directory, suffix)
    path = directory / INPUT_NAMES[1]
    with zipfile.ZipFile(path) as archive:
        members = {name: archive.read(name) for name in archive.namelist()}
    members["docProps/core.xml"] = re.sub(
        rb"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", b"2000-01-01T00:00:00Z", members["docProps/core.xml"]
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in sorted(members.items()):
            item = zipfile.ZipInfo(name, date_time=(2000, 1, 1, 0, 0, 0))
            item.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(item, content)


def paired_fastqs(directory):
    """Two matching read IDs, two 12-base mates; gzip headers contain no paths/time."""
    directory.mkdir(parents=True, exist_ok=True)
    for mate, sequence in ((1, "ACGTTGCAACGT"), (2, "TGCAACGTTGCA")):
        payload = "".join(
            f"@synthetic:{read} {mate}:N:0:ACGT\n{sequence}\n+\n{'I' * len(sequence)}\n"
            for read in range(2)
        )
        with (directory / f"sample_R{mate}.fastq.gz").open("wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as handle:
                handle.write(payload.encode())


def snapshot(directory):
    """Bytes and modification times: detect accidental rewriting as well as removal."""
    return {
        str(path.relative_to(directory)): (path.read_bytes(), path.stat().st_mtime_ns)
        for path in directory.rglob("*")
        if path.is_file()
    }


@dataclass
class Scenario:
    root: Path
    cfg: PipelineConfig
    source: Path
    output: Path
    executables: dict
    fail_workflow: bool = False

    @property
    def store(self):
        return FlowcellStateStore(self.cfg.static.paths.manager_dir)

    @property
    def manifest(self):
        return self.output / f"md5sum_{PROJECT}_fastq.txt"

    @property
    def archive(self):
        return self.output / f"{PROJECT}_{DATE}.7za"

    @property
    def retained(self):
        return self.output / "provenance" / f"{PROJECT}_analysis.tar.gz"

    def command(self, name, *arguments, expected=0, smtp="accept"):
        assert getattr(sys, "_bfq_external_guard", False), "Missing subprocess guard"
        env = os.environ.copy()
        env.update(BFQ_TEST_SCENARIO_ROOT=str(self.root), BFQ_TEST_SMTP_MODE=smtp)
        command = [str(self.executables[name]), *arguments]
        result = subprocess.run(
            command, check=False, cwd=self.root, env=env, capture_output=True, text=True, timeout=30
        )
        with (self.root / "logs" / "commands.jsonl").open("a") as handle:
            handle.write(
                json.dumps(
                    {
                        "argv": command,
                        "returncode": result.returncode,
                        "stdout": result.stdout,
                        "stderr": result.stderr,
                    }
                )
                + "\n"
            )
        assert result.returncode == expected, (
            f"{command} returned {result.returncode}; expected {expected}. "
            f"Inspect {self.root / 'logs/commands.jsonl'}\n{result.stdout}\n{result.stderr}"
        )
        return result

    def run_attempt(self):
        """One real queued attempt, using BFQ's daemon core; never start its scan loop."""
        PipelineConfig._instance = self.cfg
        self.cfg.run.begin(self.source, self.cfg.static.paths)
        cli._run_state_backed_flowcell(
            self.cfg, self.store, logging.getLogger("operational-scenario"), prepare=True
        )
        return self.store.read(RUN_ID)

    def external_command(self, command, cwd=None):
        """Only configmaker's process launch and 7za are replaced here."""
        if command[:2] == ["/opt/conda/bin/python", "/opt/conda/bin/configmaker.py"]:
            record(self.root, "configmaker-process", argv=command)
            work = Path(cwd)
            assert work.is_relative_to(self.root / "scratch")
            (work / "config.yaml").write_text(
                yaml.safe_dump({"workflow": "synthetic", "notice": LABEL})
            )
            (work / "Snakefile").write_text(f"# {LABEL}\n")
            (work / "configmaker.analysis-summary.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "kind": "fastq_discovery",
                        "sample_count": 1,
                        "missing_sample_ids": [],
                        "notice": LABEL,
                    }
                )
            )
        elif command[:2] == ["7za", "a"]:
            target = Path(command[2])
            assert target.parent == self.output
            assert all(Path(value).exists() for value in command[3:])
            record(self.root, "7za", argv=command)
            target.write_text(LABEL + "\n" + json.dumps(command[3:]) + "\n")
        else:
            raise AssertionError(f"Undeclared external command: {command}")
        return 0

    def workflow(self, command, cwd, log_path):
        assert command[0] == "snakemake"
        record(self.root, "snakemake", argv=command, failed=self.fail_workflow)
        Path(log_path).write_text(LABEL + "\n")
        if self.fail_workflow:
            raise subprocess.CalledProcessError(
                23, command, output="controlled synthetic workflow failure; no scientific tools ran"
            )
        work = Path(cwd)
        result = work / "data/tmp/synthetic/bfq"
        result.mkdir(parents=True)
        (result / f"multiqc_{PROJECT}.html").write_text(f"<html><body>{LABEL}</body></html>")
        (result / ".multiqc_config.yaml").write_text("{}\n")
        (work / "data/tmp/sample_info.tsv").write_text("sample\tgroup-curated\n")

    def sequencing_report(self, cfg):
        record(self.root, "sequencing-report-generator")
        (cfg.output_path / f"sequencer_stats_{PROJECT}.html").write_text(
            f"<html><body>{LABEL}</body></html>"
        )
        return [PROJECT]

    def archive_checksums(self, cfg):
        """Replace the external md5sum/worker-pool boundary; hash actual double bytes."""
        record(self.root, "archive-md5sum-process")
        for path in cfg.output_path.glob("*.7za"):
            (path.parent / f"md5sum_{path.stem}_archive.txt").write_text(
                f"{afterFastq.file_md5(path)}  {path.name}\n"
            )


def make_scenario(root, monkeypatch):
    """Build beneath one pytest-owned directory, without acquiring any remote data."""
    from support.fixtures import configured_bfq  # noqa: PLC0415

    executables = {name: Path(sys.executable).parent / name for name in ("fm", "flowcell-manager")}
    assert all(path.is_file() for path in executables.values()), "Run scripts/dev.py setup first"
    monkeypatch.setenv("BFQ_TEST_SCENARIOS", "1")
    monkeypatch.setenv(
        "BFQ_TEST_ALLOWED_EXECUTABLES", json.dumps([str(p) for p in executables.values()])
    )
    assert getattr(sys, "_bfq_external_guard", False), "Missing subprocess guard"
    root = root.resolve()
    cfg, source, output = configured_bfq(root, run_id=RUN_ID)
    cfg.static.system["minspace"] = "0"
    cfg.static.email.update(
        host="smtp.invalid",
        from_address="bfq@example.test",
        finished_to="finished@example.test",
        error_to="operator@example.test",
    )
    (root / "config").mkdir()
    (root / "scratch").mkdir()
    (root / "cache").mkdir()
    (root / "fixture-origin.txt").write_text(
        "BFQ synthetic fixture v1: generated in tests/support/operational.py.\n"
        "Invented identifiers and sequences; no facility data. Two read pairs, one sample.\n"
        + LABEL
        + "\n"
    )
    parser = ConfigParser()
    parser["Paths"] = {
        "ekista_baseDir": str(cfg.static.paths.ekista_base_dir),
        "nova_baseDir": str(cfg.static.paths.nova_base_dir),
        "outputDir": str(cfg.static.paths.output_dir),
        "logDir": str(cfg.static.paths.log_dir),
        "manager_dir": str(cfg.static.paths.manager_dir),
        "reportDir": str(cfg.static.paths.report_dir),
        "analysisDir": str(cfg.static.paths.analysis_dir),
    }
    parser["System"] = cfg.static.system
    parser["Email"] = cfg.static.email
    with (root / "config/bcl2fastq.ini").open("w") as handle:
        parser.write(handle)
    deterministic_inputs(source, "-instrument")
    deterministic_inputs(output, "-curated")
    paired_fastqs(output / PROJECT)
    (output / "Stats").mkdir()
    (output / "Stats/Stats.json").write_text(
        json.dumps(
            {
                "ReadInfosForLanes": [
                    {
                        "ReadInfos": [
                            {"IsIndexedRead": False, "NumCycles": 12, "Number": 1},
                            {"IsIndexedRead": False, "NumCycles": 12, "Number": 2},
                        ]
                    }
                ]
            }
        )
    )
    workflows = root / "config/gcf-workflows"
    workflows.mkdir()
    (workflows / "libprep.config").write_text("Illumina DNA Prep PE:\n  workflow: synthetic\n")
    (workflows / "README.txt").write_text(LABEL + "\n")
    scenario = Scenario(root, cfg, source, output, executables)
    monkeypatch.setenv("BFQ_ENV", "test")
    monkeypatch.setenv("TMPDIR", str(root / "scratch"))
    monkeypatch.setenv("SINGULARITY_CACHEDIR", str(root / "cache"))
    monkeypatch.setenv("APPTAINER_CACHEDIR", str(root / "cache"))
    monkeypatch.setattr(state, "collect_versions", lambda: scenario_versions(root))
    monkeypatch.setattr(workflow_config, "AUTHORITATIVE_CONFIG", workflows / "libprep.config")
    monkeypatch.setattr(afterFastq.subprocess, "check_call", scenario.external_command)
    monkeypatch.setattr(afterFastq, "run_logged_command", scenario.workflow)
    monkeypatch.setattr(afterFastq, "multiqc_stats", scenario.sequencing_report)
    monkeypatch.setattr(afterFastq, "md5sum_archive_worker", scenario.archive_checksums)
    monkeypatch.setattr(smtplib, "SMTP", smtp_double(root))
    monkeypatch.setattr(smtplib, "SMTP_SSL", smtp_double(root))
    return scenario
