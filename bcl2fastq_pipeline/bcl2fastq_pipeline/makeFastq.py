"""
This file contains functions required to actually convert the bcl files to fastq
"""

import glob
import os
import re
import shutil
import subprocess
import syslog

from distutils.dir_util import copy_tree

MKFASTQ_10X = {
    "10X Genomics Visium Spatial Gene Expression Slide & Reagents Kit": "cellranger_spatial_mkfastq",
    "10X Genomics Chromium Next GEM Single Cell ATAC Library & Gel Bead Kit v1.1": "cellranger_atac_mkfastq",
    "10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3": "cellranger_mkfastq",
}


def fixNames(config):
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = f"_lanes{lanes}"

    if "10X Genomics" in config.get("Options", "Libprep"):
        return

    names = glob.glob(
        "{}/{}{}/*/*.fastq.gz".format(
            config.get("Paths", "outputDir"), config.get("Options", "runID"), lanes
        )
    )
    names.extend(
        glob.glob(
            "{}/{}{}/*/*/*.fastq.gz".format(
                config.get("Paths", "outputDir"), config.get("Options", "runID"), lanes
            )
        )
    )
    for fname in names:
        if "_001.fastq.gz" in fname:
            fnew = fname.replace("_001.fastq.gz", ".fastq.gz")
            fnew = re.sub("_S[0-9]+", "", fnew)
            syslog.syslog(f"Moving {fname} to {fnew}\n")
            shutil.move(fname, fnew)


def fix_stats_json(stats_fn):
    stats = list()
    barcodes = False

    with open(stats_fn) as fh:
        for line in fh.readlines():
            if not barcodes:
                stats.append(line)
                if "UnknownBarcodes" in line:
                    barcodes = True
            elif not any(s in line for s in ["{", "}", "Barcodes", "Lane"]):
                stats.append(line.replace("\n", ",\n"))
            elif line.startswith("  }"):
                stats.append(line.replace("\n", ",\n"))
                stats[-2] = stats[-2].replace(",\n", "\n")

    stats[-2] = stats[-2].replace(",\n", "\n")
    stats[-3] = stats[-3].replace(",\n", "\n")
    with open(stats_fn, "w+") as fh:
        fh.writelines(stats)


def bcl2fq(config):
    """
    takes things from /dont_touch_this/solexa_runs/XXX/Data/Intensities/BaseCalls
    and writes most output into config.outputDir/XXX, where XXX is the run ID.
    """
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = f"_lanes{lanes}"

    # Make the output directories
    os.makedirs(
        "{}/{}{}".format(config.get("Paths", "outputDir"), config.get("Options", "runID"), lanes),
        exist_ok=True,
    )
    # Make log directory
    os.makedirs(
        "{}".format(
            os.path.join(
                config.get("Paths", "logDir"), os.path.dirname(config.get("Options", "runID"))
            )
        ),
        exist_ok=True,
    )
    os.makedirs(
        os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"), "InterOp"),
        exist_ok=True,
    )
    copy_tree(
        os.path.join(config.get("Paths", "baseDir"), config.get("Options", "runID"), "InterOp"),
        os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"), "InterOp"),
    )
    old_wd = os.getcwd()
    os.chdir(os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID")))
    force_bcl2fastq = os.environ.get("FORCE_BCL2FASTQ", None)

    if "10X Genomics" in config.get("Options", "Libprep"):
        cmd = "{cellranger_cmd} --output-dir={output_dir} --sample-sheet={sample_sheet} --run={run_dir} {cellranger_options}".format(
            cellranger_cmd=config.get("cellranger", MKFASTQ_10X[config.get("Options", "Libprep")]),
            output_dir="{}/{}".format(
                config.get("Paths", "outputDir"), config.get("Options", "runID")
            ),
            sample_sheet=config.get("Options", "sampleSheet"),
            run_dir="{}/{}/data/{}".format(
                config.get("Paths", "baseDir"),
                config.get("Options", "sequencer"),
                config.get("Options", "runID"),
            ),
            cellranger_options=config.get("cellranger", "cellranger_mkfastq_options"),
        )
        bcl_done = ["cellranger mkfastq", os.environ.get("CR_VERSION")]
    elif force_bcl2fastq:
        cmd = "{} {} --sample-sheet {} -o {}/{}{} -R {}/{} --interop-dir {}/{}/InterOp".format(
            config.get("bcl2fastq", "bcl2fastq"),
            config.get("bcl2fastq", "bcl2fastq_options"),
            config.get("Options", "sampleSheet"),
            config.get("Paths", "outputDir"),
            config.get("Options", "runID"),
            lanes,
            config.get("Paths", "baseDir"),
            config.get("Options", "runID"),
            config.get("Paths", "outputDir"),
            config.get("Options", "runID"),
        )
        bcl_done = ["bcl2fastq", os.environ.get("BCL2FASTQ_VERSION")]
    else:
        cmd = "bcl-convert --force --bcl-input-directory {in_dir} --output-directory {out_dir} --sample-sheet {samplesheet} --bcl-sampleproject-subdirectories true --no-lane-splitting true --output-legacy-stats true".format(
            in_dir=os.path.join(config.get("Paths", "baseDir"), config.get("Options", "runID")),
            out_dir=os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID")),
            samplesheet=os.path.join(
                config.get("Paths", "outputDir"), config.get("Options", "runID"), "SampleSheet.csv"
            ),
        )
        bcl_done = ["bcl-convert", os.environ.get("BCL_CONVERT_VERSION")]
    try:
        syslog.syslog(f"[convert bcl] Running: {cmd}\n")
        with open(
            "{}/{}{}.log".format(
                config.get("Paths", "logDir"), config.get("Options", "runID"), lanes
            ),
            "w",
        ) as logOut:
            subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True)
    except Exception:
        if "10X Genomics" not in config.get("Options", "Libprep") and force_bcl2fastq:
            with open(
                "{}/{}{}.log".format(
                    config.get("Paths", "logDir"), config.get("Options", "runID"), lanes
                )
            ) as logOut:
                log_content = logOut.read()
            if "<bcl2fastq::layout::BarcodeCollisionError>" in log_content:
                cmd += " --barcode-mismatches 0 "
                with open(
                    "{}/{}{}.log".format(
                        config.get("Paths", "logDir"), config.get("Options", "runID"), lanes
                    ),
                    "w",
                ) as logOut:
                    syslog.syslog(f"[bcl2fq] Retrying with --barcode-mismatches 0 : {cmd}\n")
                    subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True)

    if os.path.exists(os.path.join("Reports", "legacy", "Stats")):
        cmd = "ln -sr {} {}".format(
            os.path.join(
                config.get("Paths", "outputDir"),
                config.get("Options", "runID"),
                "Reports",
                "legacy",
                "Stats",
            ),
            os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"), "Stats"),
        )
        subprocess.check_call(cmd, shell=True)
        # fix_stats_json(os.path.join("Stats", "Stats.json"))

    logOut.close()
    os.chdir(old_wd)
    return bcl_done
