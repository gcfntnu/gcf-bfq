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

from bcl2fastq_pipeline.config import PipelineConfig

MKFASTQ_10X = {
    "10X Genomics Visium Spatial Gene Expression Slide & Reagents Kit": "cellranger_spatial_mkfastq",
    "10X Genomics Chromium Next GEM Single Cell ATAC Library & Gel Bead Kit v1.1": "cellranger_atac_mkfastq",
    "10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3": "cellranger_mkfastq",
}

config = None


def fixNames():
    cfg = PipelineConfig.get()

    if "10X Genomics" in cfg.run.libprep:
        return

    names = glob.glob("{cfg.output_path}/*/*.fastq.gz")
    names.extend(glob.glob("{cfg.output_path}/*/*/*.fastq.gz"))
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


def bcl2fq():
    """
    takes things from /dont_touch_this/solexa_runs/XXX/Data/Intensities/BaseCalls
    and writes most output into config.outputDir/XXX, where XXX is the run ID.
    """

    cfg = PipelineConfig.get()
    # Make the output directories
    (cfg.output_path / "InterOp").mkdir(parents=True, exist_ok=True)
    copy_tree(
        cfg.run.flowcell_path / "InterOp",
        cfg.output_path / "InterOp",
    )
    old_wd = os.getcwd()
    os.chdir(cfg.output_path)
    force_bcl2fastq = os.environ.get("FORCE_BCL2FASTQ", None)

    if "10X Genomics" in cfg.run.libprep:
        cellranger_cmd = cfg.static.commands[MKFASTQ_10X[cfg.run.libprep]]
        cellranger_options = cfg.static.commands["cellranger_mkfastq_options"]
        cmd = f"{cellranger_cmd} --output-dir={cfg.output_path} --sample-sheet={cfg.run.sample_sheet} --run={cfg.run.flowcell_path} {cellranger_options}"
        bcl_done = ["cellranger mkfastq", os.environ.get("CR_VERSION")]
    elif force_bcl2fastq:
        bcl2fastq_bin = cfg.static.commands["bcl2fastq"]
        bcl2fastq_opts = cfg.static.commands["bcl2fastq_options"]
        cmd = f"{bcl2fastq_bin} {bcl2fastq_opts} --sample-sheet {cfg.run.sample_sheet} -o {cfg.output_path} -R {cfg.run.flowcell_path} --interop-dir {cfg.output_path}/InterOp"
        bcl_done = ["bcl2fastq", os.environ.get("BCL2FASTQ_VERSION")]
    else:
        cmd = f"bcl-convert --force --bcl-input-directory {cfg.run.flowcell_path} --output-directory {cfg.output_path} --sample-sheet {cfg.run.sample_sheet} --bcl-sampleproject-subdirectories true --no-lane-splitting true --output-legacy-stats true"
        bcl_done = ["bcl-convert", os.environ.get("BCL_CONVERT_VERSION")]
    try:
        syslog.syslog(f"[convert bcl] Running: {cmd}\n")
        with open(cfg.static.paths.log / cfg.run.run_id / ".log", "w") as logOut:
            subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True)
    except Exception:
        if "10X Genomics" not in cfg.run.libprep and force_bcl2fastq:
            with open(cfg.static.paths.log / cfg.run.run_id / ".log") as logOut:
                log_content = logOut.read()
            if "<bcl2fastq::layout::BarcodeCollisionError>" in log_content:
                cmd += " --barcode-mismatches 0 "
                with open(cfg.static.paths.log / cfg.run.run_id / ".log", "w") as logOut:
                    syslog.syslog(f"[bcl2fq] Retrying with --barcode-mismatches 0 : {cmd}\n")
                    subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True)

    if (cfg.output_path / "Reports" / "legacy" / "Stats").exists():
        cmd = f"ln -sr {cfg.output_path}/Reports/legacy/Stats {cfg.output_path}/Stats"
        subprocess.check_call(cmd, shell=True)

    logOut.close()
    os.chdir(old_wd)
    return bcl_done
