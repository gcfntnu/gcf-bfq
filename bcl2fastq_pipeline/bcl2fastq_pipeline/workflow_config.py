"""BFQ source policy and execution-scoped libprep configuration snapshots."""

import json
import logging
import os

from pathlib import Path

from configmaker.libprep import LibprepConfig, LibprepConfigError, find_read_geometry

log = logging.getLogger(__name__)
AUTHORITATIVE_CONFIG = Path("/opt/gcf-workflows/libprep.config")
SNAPSHOT_NAME = "bfq-libprep.config"
SELECTION_NAME = "bfq-libprep.json"


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


def _atomic_write(path, content):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(content)
    temporary.replace(path)


def select_workflow(cfg):
    """Select once from actual demultiplexer geometry, shared with configmaker."""
    if cfg.run.workflow_selection is None:
        snapshot = capture_config(cfg)
        selection = snapshot.select(cfg.run.libprep, find_read_geometry([cfg.output_path]))
        # Retain enough configuration to resume reporting/finalization without
        # reinterpreting completed analysis through a newly edited source file.
        cfg.output_path.mkdir(parents=True, exist_ok=True)
        _atomic_write(cfg.output_path / SNAPSHOT_NAME, snapshot.content)
        _atomic_write(
            cfg.output_path / SELECTION_NAME,
            (json.dumps(selection.diagnostics(), indent=2) + "\n").encode(),
        )
        cfg.run.workflow_selection = selection
    cfg.run.pipeline = cfg.run.workflow_selection.workflow
    log.info("Libprep selection: %s", json.dumps(cfg.run.workflow_selection.diagnostics()))
    return cfg.run.workflow_selection


def prepare_execution(cfg, stage):
    """Start a fresh snapshot under the execution lease, or restore completed analysis."""
    cfg.run.libprep_config = None
    cfg.run.workflow_selection = None
    cfg.run.pipeline = None
    if stage in ("demultiplexing", "analysis"):
        capture_config(cfg)
        return

    snapshot_path = cfg.output_path / SNAPSHOT_NAME
    selection_path = cfg.output_path / SELECTION_NAME
    if not snapshot_path.exists() and not selection_path.exists():
        # Legacy results predate this feature. Make the necessary bootstrap explicit.
        log.warning(
            "No retained libprep selection for %s; resolving legacy results from %s. "
            "Verify the selected workflow matches the existing analysis.",
            cfg.run.run_id,
            AUTHORITATIVE_CONFIG,
        )
        select_workflow(cfg)
        return
    try:
        saved = json.loads(selection_path.read_bytes())
        snapshot = LibprepConfig(saved["source"], snapshot_path.read_bytes())
        selection = snapshot.select(saved["kit"], saved["read_geometry"])
        if selection.diagnostics() != saved:
            raise ValueError("retained configuration and selection/hash disagree")
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise LibprepConfigError(
            f"Cannot restore libprep selection from {selection_path}: {error}. "
            "Restore the matching snapshot files or restart from analysis."
        ) from error
    cfg.run.libprep_config = snapshot
    cfg.run.workflow_selection = selection
    cfg.run.pipeline = selection.workflow
    log.info("Restored libprep selection: %s", json.dumps(selection.diagnostics()))
