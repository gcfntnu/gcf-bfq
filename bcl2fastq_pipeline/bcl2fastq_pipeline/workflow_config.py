"""BFQ source policy and execution-scoped libprep configuration snapshots."""

import hashlib
import json
import logging
import os
import re

from pathlib import Path

import yaml

from configmaker.libprep import LibprepConfig, LibprepConfigError, find_read_geometry

from bcl2fastq_pipeline import analysis_snapshots
from bcl2fastq_pipeline.config import ProjectConfig, parse_custom_options

log = logging.getLogger(__name__)
AUTHORITATIVE_CONFIG = Path("/opt/gcf-workflows/libprep.config")


def physical_run_info(cfg):
    """Prefer the attempt's frozen physical geometry over mutable source files."""
    selection = getattr(cfg.run, "input_selection", None)
    candidates = []
    if selection is not None:
        if getattr(selection, "run_info", None):
            candidates.append(Path(selection.run_info))
        candidates.append(Path(selection.sample_sheet).parent / "RunInfo.xml")
    candidates.append(cfg.output_path / "RunInfo.xml")
    if cfg.run.flowcell_path:
        candidates.append(cfg.run.flowcell_path / "RunInfo.xml")
    return next((path for path in candidates if path.is_file()), None)


def capture_config(cfg):
    """Read the authoritative file once; never consult an environment override."""
    if cfg.run.libprep_config is None:
        if "BFQ_LIBPREP_CONFIG" in os.environ:
            log.warning("BFQ_LIBPREP_CONFIG is retired and ignored; using %s", AUTHORITATIVE_CONFIG)
        cfg.run.libprep_config = LibprepConfig.load(AUTHORITATIVE_CONFIG)
        log.info(
            "Captured libprep configuration source=%s sha256=%s",
            cfg.run.libprep_config.source,
            cfg.run.libprep_config.sha256,
        )
    return cfg.run.libprep_config


def prepare_projects(cfg):
    """Bind each multi-pair project's metadata without changing legacy globals."""
    for name in ("project_configs", "project_workflows", "project_workdirs"):
        if not hasattr(cfg.run, name):
            setattr(cfg.run, name, {})
    pairs = getattr(getattr(cfg.run, "input_selection", None), "pairs", ())
    if not pairs:
        return cfg.run.project_configs
    from configmaker.project_pairs import derive_project_geometry  # noqa: PLC0415

    for pair in pairs:
        project = pair.project_id
        if project in cfg.run.project_configs:
            continue
        custom, _ = parse_custom_options(pair.sample_sheet)
        geometry = derive_project_geometry(pair.sample_sheet, run_info=physical_run_info(cfg))
        cfg.run.project_configs[project] = ProjectConfig(
            project,
            Path(pair.sample_sheet),
            Path(pair.submission_form),
            custom,
            tuple(geometry),
        )
    return cfg.run.project_configs


def project_workflow(cfg, project):
    """Use the completed/selected workflow belonging to this project."""
    return getattr(cfg.run, "project_workflows", {}).get(project) or getattr(
        cfg.run, "pipeline", None
    )


def project_workdir(cfg, project):
    """Corrections never mutate a prior successful analysis working tree."""
    if project in getattr(cfg.run, "project_workdirs", {}):
        return Path(cfg.run.project_workdirs[project])
    name = f"{project}_{cfg.run.run_id.split('_')[0]}"
    if getattr(cfg.run, "output_override", None):
        return Path(cfg.run.output_override).parent / "analysis" / name
    return analysis_snapshots.workdir_path(cfg.run.run_id, project)


def delivery_project_config(cfg, project):
    """Deliver inputs captured by the completed analysis, including later retries."""
    current = getattr(cfg.run, "project_configs", {}).get(project)
    directory = project_workdir(cfg, project) / "inputs"
    path = directory / "project.json"
    if not path.exists():
        return current
    try:
        data = json.loads(path.read_text())
        if data["project"] != project:
            raise ValueError("project identity differs")
        sheet = directory / data["sample_sheet"]
        form = directory / data["submission_form"]
        for source in (sheet, form):
            if (
                source.parent != directory
                or hashlib.sha256(source.read_bytes()).hexdigest() != data["inputs"][source.name]
            ):
                raise ValueError(f"input snapshot does not match its captured hash: {source}")
        return ProjectConfig(project, sheet, form, data["custom"])
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise LibprepConfigError(
            f"Cannot recover completed project inputs from {path}: {error}"
        ) from error


def select_workflow(cfg, project=None):
    """Select once from actual demultiplexer geometry, shared with configmaker."""
    projects = prepare_projects(cfg)
    if project in projects:
        context = projects[project]
        if context.workflow_selection is None:
            context.workflow_selection = capture_config(cfg).select(
                context.libprep, context.effective_read_geometry
            )
        context.pipeline = context.workflow_selection.workflow
        cfg.run.project_workflows[project] = context.pipeline
        cfg.run.pipeline = context.pipeline
        log.info(
            "Project %s libprep selection: %s", project, context.workflow_selection.diagnostics()
        )
        return context.workflow_selection
    if cfg.run.workflow_selection is None:
        snapshot = capture_config(cfg)
        selection = snapshot.select(cfg.run.libprep, find_read_geometry([cfg.output_path]))
        cfg.run.workflow_selection = selection
    cfg.run.pipeline = cfg.run.workflow_selection.workflow
    log.info("Libprep selection: %s", json.dumps(cfg.run.workflow_selection.diagnostics()))
    return cfg.run.workflow_selection


def _checked_workflow(workflow, source):
    if not isinstance(workflow, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", workflow):
        raise LibprepConfigError(
            f"Missing or invalid completed analysis workflow in {source}. "
            "Restore the original analysis configuration or restart from analysis."
        )
    return workflow


def workflow_from_projects(cfg, projects):
    """Recover older state from actual generated configs, never from current /opt."""
    workflows = set()
    run_date = cfg.run.run_id.split("_", 1)[0]
    work_root = Path(os.environ.get("TMPDIR", "/bfq-tmp"))
    for project in sorted(projects):
        path = work_root / f"{project}_{run_date}" / "config.yaml"
        try:
            config = yaml.safe_load(path.read_text(encoding="utf-8"))
            workflow = _checked_workflow(config.get("workflow"), path)
        except (OSError, ValueError, AttributeError, yaml.YAMLError) as error:
            raise LibprepConfigError(
                f"Cannot recover completed analysis workflow from {path}: {error}. "
                "Restore the original project config.yaml or restart from analysis."
            ) from error
        workflows.add(workflow)
    if len(workflows) != 1:
        raise LibprepConfigError(
            f"Cannot recover one completed analysis workflow for {cfg.run.run_id}: "
            f"found {sorted(workflows)}. Restore the original project configurations "
            "or restart from analysis."
        )
    workflow = workflows.pop()
    log.info("Recovered completed analysis workflow=%s from project config.yaml", workflow)
    return workflow


def prepare_execution(cfg, stage, *, workflow=None, project_workflows=None):
    """Capture analysis inputs in memory or restore only the completed workflow name."""
    cfg.run.libprep_config = None
    cfg.run.workflow_selection = None
    cfg.run.pipeline = None
    cfg.run.project_configs.clear()
    cfg.run.project_workflows = dict(project_workflows or {})
    if stage in ("demultiplexing", "analysis"):
        snapshot = capture_config(cfg)
        for context in prepare_projects(cfg).values():
            context.workflow_selection = snapshot.select(
                context.libprep, context.effective_read_geometry
            )
    else:
        if cfg.run.project_workflows:
            for project, value in cfg.run.project_workflows.items():
                _checked_workflow(value, f"flowcell state for {project}")
            unique = set(cfg.run.project_workflows.values())
            cfg.run.pipeline = unique.pop() if len(unique) == 1 else None
        else:
            cfg.run.pipeline = _checked_workflow(workflow, "flowcell state")
        for project in cfg.run.project_workflows:
            context = delivery_project_config(cfg, project)
            if context is not None:
                context.pipeline = cfg.run.project_workflows[project]
                cfg.run.project_configs[project] = context
        log.info("Restored completed analysis workflow=%s from state", cfg.run.pipeline)
