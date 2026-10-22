"""
This file contains functions required to actually convert the bcl files to fastq
"""

import logging
import os
import re
import shlex
import shutil
import subprocess

from bcl2fastq_pipeline import demultiplexing
from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.containers import image_version, resolve_image
from bcl2fastq_pipeline.undetermined_fastqs import is_undetermined

log = logging.getLogger(__name__)

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

    if "10X Genomics" in (cfg.run.libprep or "") and not getattr(cfg.run, "project_configs", {}):
        return

    # Collect FASTQ files 1–2 levels deep
    fastqs = list(cfg.output_path.glob("*/*.fastq.gz"))
    fastqs += list(cfg.output_path.glob("*/*/*.fastq.gz"))

    renames = []
    destinations = set()
    for fpath in fastqs:
        project = next(
            (name for name in fpath.parts if name in getattr(cfg.run, "project_configs", {})), None
        )
        if project and "10X Genomics" in (cfg.run.project_configs[project].libprep or ""):
            continue
        if is_undetermined(fpath):
            # Preserve the reserved S0 identity until successful-run cleanup.
            continue
        if fpath.name.endswith("_001.fastq.gz"):
            # Build new filename
            new_name = fpath.name.replace("_001.fastq.gz", ".fastq.gz")
            new_name = re.sub(r"_S[0-9]+(?=_(?:L[0-9]{3}_)?[RI][0-9]+\.fastq\.gz$)", "", new_name)

            fnew = fpath.with_name(new_name)
            if fnew in destinations or (fnew != fpath and fnew.exists()):
                raise RuntimeError(f"FASTQ output-name conflict: {fpath} would overwrite {fnew}")
            destinations.add(fnew)
            renames.append((fpath, fnew))
    for source, destination in renames:
        log.debug("[rename_fastqs] Moving %s → %s", source, destination)
        shutil.move(str(source), str(destination))


def _pass_command(cfg, item, output):
    tool = item["tool"]
    sheet = item["derived_sample_sheet"]
    if tool.endswith("-mkfastq"):
        return [
            "bfq-container",
            tool,
            "mkfastq",
            f"--output-dir={output}",
            f"--sample-sheet={sheet}",
            f"--run={cfg.run.flowcell_path}",
            *shlex.split(cfg.static.commands["cellranger_mkfastq_options"]),
        ]
    if tool == "bcl2fastq":
        return [
            "bcl2fastq",
            *shlex.split(cfg.static.commands["bcl2fastq_options"]),
            "--sample-sheet",
            sheet,
            "-o",
            str(output),
            "-R",
            str(cfg.run.flowcell_path),
        ]
    command = [
        "bcl-convert",
        "--force",
        "--bcl-input-directory",
        str(cfg.run.flowcell_path),
        "--output-directory",
        str(output),
        "--sample-sheet",
        sheet,
        "--bcl-sampleproject-subdirectories",
        "true",
        "--no-lane-splitting",
        "true",
    ]
    if item["sample_sheet_version"] == 1:
        command.extend(["--output-legacy-stats", "true"])
    return command


def _planned_conversion(cfg):
    """Execute the frozen passes without retry/fallback or shared-setting leakage."""
    result = cfg.run.demultiplexing_plan
    items = demultiplexing.pass_plans(result)
    (cfg.output_path / "InterOp").mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        cfg.run.flowcell_path / "InterOp", cfg.output_path / "InterOp", dirs_exist_ok=True
    )
    identities = []
    for item in items:
        tool = item["tool"]
        version = image_version(tool)
        for source in item["sources"]:
            allowed = source.get("profile", {}).get("bcl_convert_versions")
            if tool == "bcl-convert" and allowed and version not in allowed:
                raise demultiplexing.DemultiplexingPlanError(
                    f"Profile requires BCL Convert {allowed}; selected image reports {version}. Validate a coordinated profile before changing versions."
                )
        output = (
            cfg.output_path
            if len(items) == 1
            else cfg.output_path / "Reports" / "passes" / item["pass_id"]
        )
        output.mkdir(parents=True, exist_ok=True)
        command = _pass_command(cfg, item, output)
        # Ordinary projects keep their established version-only record. The 10x
        # matrix binds both command and full workflow-owned image reference.
        image = (
            resolve_image(tool)
            if any("name" in source.get("profile", {}) for source in item["sources"])
            else None
        )
        log_path = cfg.static.paths.log_dir / (
            f"{cfg.run.run_id}.log"
            if len(items) == 1
            else f"{cfg.run.run_id}.{item['pass_id']}.log"
        )
        item["output_directory"] = str(output)
        demultiplexing.record(
            cfg, item=item, version=version, image=image, command=command, status="running"
        )
        try:
            with log_path.open("w") as handle:
                subprocess.check_call(command, stdout=handle, stderr=subprocess.STDOUT, cwd=output)
        except subprocess.CalledProcessError as error:
            demultiplexing.record(cfg, item=item, status="failed")
            if len(items) > 1:
                demultiplexing.record(cfg, status="failed")
            raise demultiplexing_error(error, log_path) from error
        if len(items) > 1:
            for project in item["projects"]:
                source = output / project
                target = cfg.output_path / project
                if target.exists():
                    raise RuntimeError(f"Separate-pass publication would overwrite {target}")
                if not source.is_dir():
                    raise RuntimeError(
                        f"Conversion did not create expected project directory {source}"
                    )
                shutil.move(str(source), str(target))
        demultiplexing.record(cfg, item=item, status="completed")
        identities.append(f"{tool}:{version}")
    legacy_stats = cfg.output_path / "Reports" / "legacy" / "Stats"
    stats = cfg.output_path / "Stats"
    if legacy_stats.is_dir() and not stats.exists():
        stats.symlink_to(legacy_stats, target_is_directory=True)
    else:
        stats.mkdir(exist_ok=True)
    demultiplexing.record(cfg, status="completed")
    for project in result["projects"]:
        demultiplexing.delivery_provenance(cfg, project, refresh=True)
    if len(items) == 1:
        tool = items[0]["tool"]
        name = tool.removesuffix("-mkfastq") + " mkfastq" if tool.endswith("-mkfastq") else tool
        return [name, items[0]["tool_version"]]
    return ["multiple conversion passes", "; ".join(identities)]


def bcl2fq():
    """
    takes things from /dont_touch_this/solexa_runs/XXX/Data/Intensities/BaseCalls
    and writes most output into config.outputDir/XXX, where XXX is the run ID.
    """

    cfg = PipelineConfig.get()
    multi_pair = bool(getattr(getattr(cfg.run, "input_selection", None), "pairs", ()))
    selection = getattr(cfg.run, "input_selection", None)
    if selection is None and "10X Genomics" in (cfg.run.libprep or ""):
        from types import SimpleNamespace  # noqa: PLC0415

        selection = SimpleNamespace(
            sample_sheet=cfg.run.sample_sheet,
            submission_form=cfg.output_path / "Sample-Submission-Form.xlsx",
            pairs=(),
        )
        cfg.run.input_selection = selection
    if multi_pair or (selection and demultiplexing.needs_plan(selection)):
        demultiplexing.prepare(cfg)
        return _planned_conversion(cfg)
    # Make the output directories
    (cfg.output_path / "InterOp").mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        cfg.run.flowcell_path / "InterOp", cfg.output_path / "InterOp", dirs_exist_ok=True
    )
    force_bcl2fastq = os.environ.get("FORCE_BCL2FASTQ", None)

    if force_bcl2fastq:
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
        bcl_done = ["bcl2fastq", image_version("bcl2fastq")]
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
        ]
        if not multi_pair or cfg.run.demultiplexing_plan["sample_sheet_version"] == 1:
            # 4.2.4 aborts when legacy stats are requested for a v2 sheet.
            cmd.extend(["--output-legacy-stats", "true"])
        bcl_done = ["bcl-convert", image_version("bcl-convert")]

    log_pth = cfg.static.paths.log_dir / f"{cfg.run.run_id}.log"
    if multi_pair:
        demultiplexing.record(cfg, version=bcl_done[1], command=cmd, status="running")
    try:
        log.info(f"[convert bcl] Running: {shlex.join(cmd)}\n")
        with log_pth.open("w") as logOut:
            subprocess.check_call(cmd, stdout=logOut, stderr=subprocess.STDOUT, cwd=cfg.output_path)
    except subprocess.CalledProcessError as error:
        if multi_pair:
            demultiplexing.record(cfg, version=bcl_done[1], command=cmd, status="failed")
        if "10X Genomics" not in (cfg.run.libprep or "") and force_bcl2fastq and not multi_pair:
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
        if dst.is_symlink():
            dst.unlink()
        elif dst.exists():
            shutil.rmtree(dst)
        dst.symlink_to(src, target_is_directory=True)
    elif multi_pair:
        # Native v2 conversion reports stay in Reports. BFQ's additional QC
        # products use Stats without inventing legacy conversion measurements.
        (cfg.output_path / "Stats").mkdir(exist_ok=True)

    if multi_pair:
        demultiplexing.record(cfg, version=bcl_done[1], command=cmd, status="completed")
        for project in cfg.run.demultiplexing_plan["projects"]:
            demultiplexing.delivery_provenance(cfg, project, refresh=True)

    return bcl_done
