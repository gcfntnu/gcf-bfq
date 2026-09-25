"""
This file contains functions required to actually convert the bcl files to fastq
"""

import logging
import os
import re
import shlex
import shutil
import subprocess

from bcl2fastq_pipeline.config import PipelineConfig

log = logging.getLogger(__name__)

MKFASTQ_10X = {
    "10X Genomics Visium Spatial Gene Expression Slide & Reagents Kit": "spaceranger",
    "10X Genomics Chromium Next GEM Single Cell ATAC Library & Gel Bead Kit v1.1": "cellranger-atac",
    "10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3": "cellranger",
}

DEMULTIPLEX_LOG_TAIL_LINES = 50


def demultiplexing_error(error, log_path):
    """Build an error that points to the full log and includes useful context."""
    try:
        log_lines = log_path.read_text(errors="replace").splitlines()
    except OSError as log_error:
        log_tail = f"Unable to read demultiplexing log: {log_error}"
    else:
        log_tail = "\n".join(log_lines[-DEMULTIPLEX_LOG_TAIL_LINES:])
        if not log_tail:
            log_tail = "Demultiplexing log is empty."

    return RuntimeError(
        f"Demultiplexing failed with exit code {error.returncode}. "
        f"Full log: {log_path}\n"
        f"Last {DEMULTIPLEX_LOG_TAIL_LINES} log lines:\n{log_tail}"
    )


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
        cellranger_cmd = MKFASTQ_10X[cfg.run.libprep]
        cellranger_options = cfg.static.commands["cellranger_mkfastq_options"]
        cmd = [cellranger_cmd, "mkfastq"]
        cmd.extend(
            [
                f"--output-dir={cfg.output_path}",
                f"--sample-sheet={cfg.run.sample_sheet}",
                f"--run={cfg.run.flowcell_path}",
            ]
        )
        cmd.extend(shlex.split(cellranger_options))
        bcl_done = ["cellranger mkfastq", os.environ.get("CR_VERSION")]
    elif force_bcl2fastq:
        bcl2fastq_opts = cfg.static.commands["bcl2fastq_options"]
        cmd = ["bcl2fastq", *shlex.split(bcl2fastq_opts)]
        cmd.extend(
            [
                "--sample-sheet",
                str(cfg.run.sample_sheet),
                "-o",
                str(cfg.output_path),
                "-R",
                str(cfg.run.flowcell_path),
                "--interop-dir",
                str(cfg.output_path / "InterOp"),
            ]
        )
        bcl_done = ["bcl2fastq", os.environ.get("BCL2FASTQ_VERSION")]
    else:
        cmd = [
            "bcl-convert",
            "--force",
            "--bcl-input-directory",
            str(cfg.run.flowcell_path),
            "--output-directory",
            str(cfg.output_path),
            "--sample-sheet",
            str(cfg.run.sample_sheet),
            "--bcl-sampleproject-subdirectories",
            "true",
            "--no-lane-splitting",
            "true",
            "--output-legacy-stats",
            "true",
        ]
        bcl_done = ["bcl-convert", os.environ.get("BCL_CONVERT_VERSION")]

    log_pth = cfg.static.paths.log_dir / f"{cfg.run.run_id}.log"
    try:
        log.info(f"[convert bcl] Running: {shlex.join(cmd)}\n")
        with log_pth.open("w") as logOut:
            subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, cwd=cfg.output_path)
    except subprocess.CalledProcessError as error:
        if "10X Genomics" not in cfg.run.libprep and force_bcl2fastq:
            log_content = log_pth.read_text(errors="replace")
            if "<bcl2fastq::layout::BarcodeCollisionError>" in log_content:
                cmd.extend(["--barcode-mismatches", "0"])
                try:
                    with log_pth.open("a") as logOut:
                        logOut.write("\nRetrying with --barcode-mismatches 0\n")
                        logOut.flush()
                        log.info(
                            "[bcl2fq] Retrying with --barcode-mismatches 0: %s\n",
                            shlex.join(cmd),
                        )
                        subprocess.check_call(
                            cmd,
                            stdout=logOut,
                            stderr=subprocess.STDOUT,
                            cwd=cfg.output_path,
                        )
                except subprocess.CalledProcessError as retry_error:
                    raise demultiplexing_error(retry_error, log_pth) from retry_error
            else:
                raise demultiplexing_error(error, log_pth) from error
        else:
            raise demultiplexing_error(error, log_pth) from error

    src = cfg.output_path / "Reports" / "legacy" / "Stats"
    if src.exists():
        dst = cfg.output_path / "Stats"
        if dst.exists() or dst.is_symlink():
            dst.unlink()  # remove old link or directory first
        dst.symlink_to(src, target_is_directory=True)

    return bcl_done
