"""BFQ source policy and execution-scoped libprep configuration snapshots."""

import json
import logging
import os
import re

from pathlib import Path

import yaml

from configmaker.libprep import LibprepConfig, LibprepConfigError, find_read_geometry

log = logging.getLogger(__name__)
AUTHORITATIVE_CONFIG = Path("/opt/gcf-workflows/libprep.config")


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


def select_workflow(cfg):
    """Select once from actual demultiplexer geometry, shared with configmaker."""
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


def prepare_execution(cfg, stage, *, workflow=None):
    """Capture analysis inputs in memory or restore only the completed workflow name."""
    cfg.run.libprep_config = None
    cfg.run.workflow_selection = None
    cfg.run.pipeline = None
    if stage in ("demultiplexing", "analysis"):
        capture_config(cfg)
    else:
        cfg.run.pipeline = _checked_workflow(workflow, "flowcell state")
        log.info("Restored completed analysis workflow=%s from state", cfg.run.pipeline)
