"""
This file contains functions required to actually convert the bcl files to fastq
"""

import logging
import os
import re
import shutil
import subprocess

from bcl2fastq_pipeline.config import PipelineConfig

log = logging.getLogger(__name__)

MKFASTQ_10X = {
    "10X Genomics Visium Spatial Gene Expression Slide & Reagents Kit": "cellranger_spatial_mkfastq",
    "10X Genomics Chromium Next GEM Single Cell ATAC Library & Gel Bead Kit v1.1": "cellranger_atac_mkfastq",
    "10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3": "cellranger_mkfastq",
}


def rename_fastqs():
    """
    Find and rename FASTQ files under cfg.output_path:
    - Removes lane suffix `_001`
    - Removes sample numbering `_S<number>`
    """
    cfg = PipelineConfig.get()

    if "10X Genomics" in cfg.run.libprep:
        return

    # Collect FASTQ files 1–2 levels deep
    fastqs = list(cfg.output_path.glob("*/*.fastq.gz"))
    fastqs += list(cfg.output_path.glob("*/*/*.fastq.gz"))

    for fpath in fastqs:
        if fpath.name.endswith("_001.fastq.gz"):
            # Build new filename
            new_name = fpath.name.replace("_001.fastq.gz", ".fastq.gz")
            new_name = re.sub(r"_S[0-9]+", "", new_name)

            fnew = fpath.with_name(new_name)
            log.debug(f"[rename_fastqs] Moving {fpath} → {fnew}")

            # Ensure parent directory exists (should already)
            fnew.parent.mkdir(parents=True, exist_ok=True)

            shutil.move(str(fpath), str(fnew))


def bcl2fq():
    """
    takes things from /dont_touch_this/solexa_runs/XXX/Data/Intensities/BaseCalls
    and writes most output into config.outputDir/XXX, where XXX is the run ID.
    """

    cfg = PipelineConfig.get()
    # Make the output directories
    (cfg.output_path / "InterOp").mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        cfg.run.flowcell_path / "InterOp", cfg.output_path / "InterOp", dirs_exist_ok=True
    )
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

    log_pth = cfg.static.paths.log_dir / "{cfg.run.run_id}.log"
    try:
        log.info(f"[convert bcl] Running: {cmd}\n")
        with log_pth.open("w") as logOut:
            subprocess.check_call(
                cmd, stdout=logOut, stderr=subprocess.STDOUT, shell=True, cwd=cfg.output_path
            )
    except Exception:
        if "10X Genomics" not in cfg.run.libprep and force_bcl2fastq:
            with log_pth.open("w") as logOut:
                log_content = logOut.read()
            if "<bcl2fastq::layout::BarcodeCollisionError>" in log_content:
                cmd += " --barcode-mismatches 0 "
                with log_pth.open("w") as logOut:
                    log.info(f"[bcl2fq] Retrying with --barcode-mismatches 0 : {cmd}\n")
                    subprocess.check_call(
                        cmd,
                        stdout=logOut,
                        stderr=subprocess.STDOUT,
                        shell=True,
                        cwd=cfg.output_path,
                    )

    src = cfg.output_path / "Reports" / "legacy" / "Stats"
    if src.exists():
        dst = cfg.output_path / "Stats"
        if dst.exists() or dst.is_symlink():
            dst.unlink()  # remove old link or directory first
        dst.symlink_to(src, target_is_directory=True)

    return bcl_done
