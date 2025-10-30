"""
This file includes anything involved in finding new flow cells to process.

Note that this includes anything done after a flow cell has been processed,
such as marking it as having been processed and sending emails.
"""

import datetime as dt
import logging
import shutil

from pathlib import Path

import flowcell_manager.flowcell_manager as fm

import bcl2fastq_pipeline.afterFastq as af

from bcl2fastq_pipeline.config import PipelineConfig, parse_custom_options

log = logging.getLogger("bfq")


def modified_time(path: Path):
    return dt.fromtimestamp(path.stat().st_mtime)


# Returns True on processed, False on unprocessed
def flowCellProcessed():
    cfg = PipelineConfig.get()
    flowcells = fm.list_flowcell_all(str(cfg.output_path))
    if not flowcells.empty:
        if rerunFlowcell(cfg):
            return False
        else:
            return True
    return False


# Determine if the flowcell should be rerun
def rerunFlowcell(cfg):
    return False  # deprecated for now
    opts, ss = get_sample_sheet(cfg.run.flowcell_path)
    if not opts:
        return False
    if opts.get("Rerun", False):
        if (cfg.output_path / "SampleSheet.csv").exists():
            prev_start = modified_time(cfg.output_path / "SampleSheet.csv")
        else:
            # The bfq output for the flowcell has been deleted manually. Rerun
            return True
        instr_ss = modified_time(ss)
        if instr_ss > prev_start:
            fm.rerun_flowcell(
                flowcell=cfg.output_path,
                force=True,
            )
            return True
    return False


def get_sample_sheet(d):
    """
    Provide a list of output directories and sample sheets
    """
    ss = list(d.glob("SampleSheet*.csv"))

    if len(ss) == 0:
        return None, None

    for sheet in ss:
        opts, ss_ = parse_custom_options(sheet)
        if ss_ and opts:
            return opts, ss_
    return None, None


"""
Iterate over all folders in config.baseDir from machine SN7001180. For each,
see if it contains an RTAComplete.txt file, as this signifies the completion
of a sequencing run.

If a run has finished, we then check to see if it's already been processed.
Runs are marked as having been processed if they appear in config.finalDir
and have a file called casava.finished or fastq.made. casava.finished is
produced by the old pipeline, while this one creates fastq.made.

This function always returns its configuration. If there's a new flow cell to
process, then the runID is filled in. Otherwise, that's set to None.
"""


def newFlowCell():
    cfg = PipelineConfig.get()

    # EMERGENCY FINNMARK FIX
    if (cfg.output_path).exists():
        opts, ss = get_sample_sheet(cfg.output_path)
        sample_sub_f = cfg.output_path.glob("*Sample-Submission-Form*.xlsx")
        use_bfq_output_ss = True
        log.debug(f"Using samplesheet and submission form from {cfg.output_path}")
    else:
        opts, ss = get_sample_sheet(cfg.run.flowcell_path)
        sample_sub_f = cfg.run.flowcell_path.glob("*Sample-Submission-Form*.xlsx")
        use_bfq_output_ss = False
        log.debug(f"Using samplesheet and submission form from {cfg.run.flowcell_path}")

    if not opts:
        cfg.run.reset()
        log.debug("No custom opts in sample sheet")
        return
    if not list(sample_sub_f):
        cfg.run.reset()
        log.debug("No sample submission form")
        return

    log.info(f"Found a new flow cell: {cfg.run.run_id}")
    if not (cfg.output_path).exists():
        (cfg.output_path).mkdir()

    if not use_bfq_output_ss:
        sample_sub_f = copy_sample_sub_form(cfg.run.flowcell_path, cfg.output_path)
    else:
        sample_sub_f = cfg.output_path / "Sample-Submission-Form.xlsx"

    if ss is not None and use_bfq_output_ss:
        ss = cfg.output_path / "SampleSheet.csv"
        cfg.run.apply_custom(opts, ss, sample_sub_f)
    elif ss is not None and not use_bfq_output_ss:
        shutil.copy2(ss, cfg.output_path / "SampleSheet.csv")
        ss = cfg.output_path / "SampleSheet.csv"
        cfg.run.apply_custom(opts, ss, sample_sub_f)
    else:
        cfg.run.reset()
    return


def copy_sample_sub_form(instrument_path, output_path):
    sample_sub_forms = list(instrument_path.glob("*Sample-Submission-Form*.xlsx"))
    if sample_sub_forms:
        sample_sub_form = sample_sub_forms[0]
        shutil.copy2(sample_sub_form, output_path / "Sample-Submission-Form.xlsx")
        return output_path / "Sample-Submission-Form.xlsx"
    return None


def markFinished():
    cfg = PipelineConfig.get()
    (cfg.output_path / "fastq.made").write_text("")
    project_dirs = af.get_project_dirs(cfg)
    project_names = af.get_project_names(project_dirs)
    now = dt.datetime.now()
    for gcf in project_names:
        fm.add_flowcell(
            project=gcf,
            path=str(cfg.output_path),
            timestamp=now,
        )
