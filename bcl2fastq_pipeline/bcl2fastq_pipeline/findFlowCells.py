"""Discovery and completion bookkeeping for BFQ flowcells."""

import csv
import logging
import shutil

from pathlib import Path

import flowcell_manager.flowcell_manager as fm

import bcl2fastq_pipeline.afterFastq as af

from bcl2fastq_pipeline.config import PipelineConfig, parse_custom_options
from bcl2fastq_pipeline.preflight import copy_run_inputs, select_run_inputs
from bcl2fastq_pipeline.state import (
    FlowcellStateStore,
    new_state,
    output_entries,
    resolve_output_path,
    validate_restored_fastqs,
)

log = logging.getLogger(__name__)


def get_sample_sheet(directory):
    """Return the first SampleSheet containing BFQ custom options."""
    for sheet in sorted(Path(directory).glob("SampleSheet*.csv")):
        opts, parsed = parse_custom_options(sheet)
        if opts:
            return opts, parsed
    return None, None


def _submission_form(directory):
    forms = sorted(Path(directory).glob("*Sample-Submission-Form*.xlsx"))
    return forms[0] if forms else None


def _write_discovery_error(cfg, detail):
    report_dir = cfg.static.paths.report_dir
    report_dir.mkdir(parents=True, exist_ok=True)
    path = report_dir / f"{cfg.run.run_id}.error"
    path.write_text(
        f"BFQ refused automatic initialization for {cfg.run.run_id}.\n"
        f"{detail}\n\n"
        "Resolve explicitly with flowcell-manager initialize RUN_ID "
        "--from demultiplexing|analysis.\n",
        encoding="utf-8",
    )
    log.error("%s", path.read_text(encoding="utf-8").strip())
    return path


def flowCellProcessed():
    """Apply #104 discovery precedence and return True when this run must be skipped."""
    cfg = PipelineConfig.get()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    store.ensure_writable()
    run_id = cfg.run.run_id

    if store.exists(run_id):
        state = store.recover_interrupted(run_id)
        return state["status"] != "queued"

    if not fm.list_flowcell_all(str(cfg.output_path)).empty:
        return True

    if cfg.output_path.exists():
        recognized, detail = validate_restored_fastqs(cfg.output_path)
        if recognized:
            store.create(
                new_state(
                    run_id,
                    cfg.run.flowcell_path,
                    cfg.output_path,
                    origin="restored_legacy_fastq",
                    start_stage="analysis",
                    cfg=cfg,
                )
            )
            return False

        # Any pre-existing output without recognized restored FASTQs is deliberately
        # ambiguous. Marker files are ignored, even when they happen to exist.
        _write_discovery_error(cfg, detail)
        return True

    # The canonical state record is BFQ's first durable action for a new run.
    store.create(
        new_state(
            run_id,
            cfg.run.flowcell_path,
            cfg.output_path,
            origin="new",
            start_stage="demultiplexing",
            cfg=cfg,
        )
    )
    return False


def newFlowCell():
    """Load/copy run inputs for a queued state-backed flowcell."""
    cfg = PipelineConfig.get()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    state = store.read(cfg.run.run_id)
    output = resolve_output_path(state["output_path"], cfg.run.run_id, cfg)
    output_entries(output, allow_missing=state["current_stage"] == "demultiplexing")

    selection = select_run_inputs(cfg.run.flowcell_path, output)
    opts = {}
    try:
        opts, _sheet = parse_custom_options(selection.sample_sheet)
    except (OSError, UnicodeError, csv.Error):
        # The shared validator supplies the actionable malformed-input report.
        pass
    curated = selection.sample_sheet.parent == output
    if not opts and not curated and state["origin"] == "new" and not state.get("restart_request"):
        # Automatic instrument discovery remains opt-in via BFQ CustomOptions.
        log.debug("No BFQ [CustomOptions] sample sheet for %s", cfg.run.run_id)
        cfg.run.reset()
        return

    copied = copy_run_inputs(selection, output)
    cfg.run.apply_custom(opts, copied.sample_sheet, copied.submission_form)
    log.info(
        "Prepared %s from state origin=%s stage=%s",
        cfg.run.run_id,
        state["origin"],
        state["current_stage"],
    )


def copy_sample_sub_form(instrument_path, output_path):
    """Compatibility helper retained for callers outside discovery."""
    form = _submission_form(instrument_path)
    if form is None:
        return None
    destination = Path(output_path) / "Sample-Submission-Form.xlsx"
    shutil.copy2(form, destination)
    return destination


def markFinished():
    """Update the compatibility inventory without legacy marker files or duplicates."""
    cfg = PipelineConfig.get()
    project_dirs = af.get_project_dirs(cfg)
    project_names = sorted(af.get_project_names(project_dirs))
    for project in project_names:
        fm.add_flowcell(project=project, path=str(cfg.output_path))
    return project_names
