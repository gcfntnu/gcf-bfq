"""
This file includes anything involved in finding new flow cells to process.

Note that this includes anything done after a flow cell has been processed,
such as marking it as having been processed and sending emails.
"""

import datetime as dt
import glob
import os
import syslog
import xml.etree.ElementTree as ET

from shutil import copyfile

import flowcell_manager.flowcell_manager as fm

import bcl2fastq_pipeline.afterFastq as af

CUSTOM_OPTS = [
    "Organism",
    "Libprep",
    "User",
    "Rerun",
    "SingleCell",
    "RemoveHumanReads",
    "SensitiveData",
    "ReverseComplementIndexP5",
    "ReverseComplementIndexP7",
    "TrimAdapter",
    "Adapter",
    "AdapterRead2",
]


# Returns True on processed, False on unprocessed
def flowCellProcessed(config):
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = f"_lanes{lanes}"

    flowcells = fm.list_flowcell_all(
        os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"))
    )
    # path = "%s/%s%s/fastq.made" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes)
    # if os.access(path, os.F_OK) or (not flowcells.empty):
    if not flowcells.empty:
        if rerunFlowcell(config):
            return False
        else:
            return True
    return False


# Determine if the flowcell should be rerun
def rerunFlowcell(config):
    return False #deprecated for now
    # seq_data_path = af.SEQUENCER_OUTPUTFOLDER[config.get("Options","runID").split("_")[-3]]
    ss, opts = getSampleSheets(
        os.path.join(config.get("Paths", "baseDir"), config.get("Options", "runID"))
    )
    if not opts:
        return False
    if opts.get("Rerun", False):
        if os.path.exists(
            os.path.join(
                config.get("Paths", "outputDir"), config.get("Options", "runID"), "SampleSheet.csv"
            )
        ):
            prev_start = os.path.getmtime(
                os.path.join(
                    config.get("Paths", "outputDir"),
                    config.get("Options", "runID"),
                    "SampleSheet.csv",
                )
            )
        else:
            # The bfq output for the flowcell has been deleted manually. Rerun
            return True
        instr_ss = os.path.getmtime(ss)
        if instr_ss > prev_start:
            fm.rerun_flowcell(
                flowcell=os.path.join(
                    config.get("Paths", "outputDir"), config.get("Options", "runID")
                ),
                force=True,
            )
            return True
    return False


# Get the number of lanes in the run. This might not match the number of lanes in the sampleSheet
def getNumLanes(d):
    try:
        tree = ET.parse(f"{d}/RunInfo.xml")
        root = tree.getroot()[0]
        numLanes = root.findall("FlowcellLayout")[0]
        return int(numLanes.get("LaneCount"))
    except:
        return 1


def parseSampleSheet(ss):
    """
    Return a dictionary with keys: (Barcode length 1, Barcode length 2)

    return ss, laneOut, bcLens
    """

    f = open(ss)
    opt_d = None
    opts_data = False
    for line in f:
        if line.startswith("[Data]"):
            return ss, opt_d
        elif line.startswith("[CustomOptions]"):
            opt_d = dict.fromkeys(CUSTOM_OPTS, "")
            opts_data = True
            continue
        elif opts_data:
            key = line.split(",")[0]
            value = line.split(",")[1]
            opt_d[key] = (
                value.rstrip()
                if key in ["Organism", "Libprep", "User", "Adapter", "AdapterRead2"]
                else str2bool(value.rstrip())
            )
    return ss, opt_d


def str2bool(s):
    return s.lower() in ["true", "1"]


def getSampleSheets(d):
    """
    Provide a list of output directories and sample sheets
    """
    ss = glob.glob(f"{d}/SampleSheet*.csv")

    if len(ss) == 0:
        return ([None], None)

    for sheet in ss:
        ss_, opts = parseSampleSheet(sheet)
        if ss_ and opts:
            return ss_, opts
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


def newFlowCell(config):
    # instrument_dir = os.path.dirname(d)
    # instrument_dir = os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),"data",config.get("Options","runID"))
    instrument_dir = os.path.join(config.get("Paths", "baseDir"), config.get("Options", "runID"))
    odir = os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"))

    # EMERGENCY FINNMARK FIX
    if os.path.exists(odir):
        # ss, opts = getSampleSheets(os.path.dirname(d))
        ss, opts = getSampleSheets(odir)
        sample_sub_f = glob.glob(os.path.join(odir, "*Sample-Submission-Form*.xlsx"))
        use_bfq_output_ss = True
    else:
        # ss, opts = getSampleSheets(os.path.dirname(d))
        ss, opts = getSampleSheets(instrument_dir)
        sample_sub_f = glob.glob(os.path.join(instrument_dir, "*Sample-Submission-Form*.xlsx"))
        use_bfq_output_ss = False

    if not opts or not sample_sub_f:
        config.set("Options", "runID", "")
        config.set("Paths", "baseDir", "")
        return config

    syslog.syslog("Found a new flow cell: {}\n".format(config.get("Options", "runID")))
    if not os.path.exists(odir):
        os.makedirs(odir)

    if not use_bfq_output_ss:
        sample_sub_f = copy_sample_sub_form(instrument_dir, odir)
    else:
        sample_sub_f = f"{odir}/Sample-Submission-Form.xlsx"

    if ss is not None and use_bfq_output_ss:
        ss = f"{odir}/SampleSheet.csv"
        config.set("Options", "sampleSheet", ss)
        config.set("Options", "sampleSubForm", sample_sub_f if sample_sub_f else "")
        config = setConfFromOpts(config, opts)
    elif ss is not None and not use_bfq_output_ss:
        copyfile(ss, f"{odir}/SampleSheet.csv")
        ss = f"{odir}/SampleSheet.csv"
        config.set("Options", "sampleSheet", ss)
        config.set("Options", "sampleSubForm", sample_sub_f if sample_sub_f else "")
        config = setConfFromOpts(config, opts)
    else:
        config.set("Options", "runID", "")
        config.set("Paths", "baseDir", "")
        config = setConfFromOpts(config, opts, use_dict_values=False)

    return config


def bool2strint(b):
    return "1" if b else "0"


def copy_sample_sub_form(instrument_path, output_path):
    sample_sub_forms = glob.glob(f"{instrument_path}/*Sample-Submission-Form*.xlsx")
    if bool(sample_sub_forms):
        sample_sub_form = sample_sub_forms[0]
        copyfile(sample_sub_form, f"{output_path}/Sample-Submission-Form.xlsx")
        return f"{output_path}/Sample-Submission-Form.xlsx"
    return None


def setConfFromOpts(config, opts, use_dict_values=True):
    if not opts:
        opts = dict.fromkeys(CUSTOM_OPTS, False)
    for k, v in opts.items():
        if use_dict_values:
            config.set(
                "Options",
                k,
                v
                if k in ["Organism", "Libprep", "User", "Adapter", "AdapterRead2"]
                else bool2strint(v),
            )
        else:
            config.set("Options", k, "")
    return config


def markFinished(config):
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = f"_lanes{lanes}"

    open(
        "{}/{}{}/fastq.made".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes),
        "w",
    ).close()
    project_dirs = af.get_project_dirs(config)
    project_names = af.get_project_names(project_dirs)
    now = dt.datetime.now()
    for gcf in project_names:
        fm.add_flowcell(
            project=gcf,
            path=os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID")),
            timestamp=now,
        )


"""
This function needs to be run after newFlowCell() returns with config.runID
filled in. It creates the output directories.
"""


def MakeTargetDirs(config):
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = f"_lanes{lanes}"

    assert config["Paths"]["runID"] is not None
    os.mkdirs("{}/{}{}".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes))
