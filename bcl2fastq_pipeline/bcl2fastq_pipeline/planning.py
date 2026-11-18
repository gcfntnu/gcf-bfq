"""Versioned, read-only conversion planning for callers outside a run folder."""

import hashlib
import xml.etree.ElementTree as ET

from pathlib import Path
from types import SimpleNamespace

from configmaker.project_pairs import ProjectPair, validate_sheet_set

from bcl2fastq_pipeline import demultiplexing
from bcl2fastq_pipeline.config import parse_custom_options

PLANNING_API_VERSION = 1


def plan_flowcell(project_sheets, *, profiles_path, libprep_config=None, run_info=None):
    """Return a JSON-ready report on explicit project inputs, with no operational state."""
    pairs = [
        ProjectPair(project, str(path), None) for project, path in sorted(project_sheets.items())
    ]
    passes = []
    findings = []
    effective = {}
    legacy = []
    profiles_digest = None
    try:
        if profiles_path is not None:
            profiles_digest = hashlib.sha256(Path(profiles_path).read_bytes()).hexdigest()
        elif any(
            "10X Genomics" in parse_custom_options(Path(p.sample_sheet))[0].get("Libprep", "")
            for p in pairs
        ):
            raise demultiplexing.RunInfoRequired(
                "Supply the workflow-owned demultiplexing.config for 10x planning."
            )
        plan = demultiplexing.plan(
            SimpleNamespace(pairs=pairs), profiles_path=profiles_path, run_info=run_info
        )
        for item in demultiplexing.pass_plans(plan) if plan else []:
            passes.append({key: item[key] for key in ("pass_id", "tool", "projects")})
            for source in item["sources"]:
                effective[source["project"]] = source["settings"]
                if item["tool"] != "bcl-convert":
                    legacy.append(source["project"])
    except demultiplexing.RunInfoRequired as error:
        findings.append(("info", "planning.input_required", str(error)))
    except (ValueError, OSError, ET.ParseError) as error:
        findings.append(("error", "demux.unsupported_plan", str(error)))
    result = validate_sheet_set(
        project_sheets,
        run_info=run_info,
        libprep_config=libprep_config,
        effective_settings=effective,
        allow_separate_passes=True,
        skip_conversion_projects=legacy,
    )
    for severity, code, message in findings:
        result.add(severity, code, message)
    report = result.to_dict()
    report.update(
        kind="flowcell_plan",
        planning_api_version=PLANNING_API_VERSION,
        passes=passes,
        profiles_sha256=profiles_digest,
    )
    if result.ok and any(
        f["code"] in {"planning.geometry_required", "planning.input_required"}
        for f in result.findings
    ):
        report["status"] = "needs_input"
    return report
