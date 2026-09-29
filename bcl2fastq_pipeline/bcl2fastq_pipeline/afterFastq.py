"""
This file includes code that actually runs FastQC and any other tools after the fastq files have actually been made. This uses a pool of workers to process each request.
"""

import codecs
import errno
import hashlib
import json
import logging
import multiprocessing as mp
import os
import pty
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

from collections import deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml

from configmaker.configmaker import SEQUENCERS
from configmaker.validation import VALIDATOR_VERSION

from bcl2fastq_pipeline import workflow_config
from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.interop import prepare_index_metrics, run_interop_csv
from bcl2fastq_pipeline.workflow_config import select_workflow

log = logging.getLogger(__name__)
COMMAND_OUTPUT_TAIL_LINES = 400
ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


def plain_command_output(text):
    """Remove terminal formatting before persisting command output."""
    return ANSI_ESCAPE_RE.sub("", text).replace("\r", "")


def run_logged_command(cmd, cwd, log_path):
    """Run under a PTY, preserving console colour while saving a plain-text log."""
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    output_tail = deque(maxlen=COMMAND_OUTPUT_TAIL_LINES)
    decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
    pending = ""
    master_fd, slave_fd = pty.openpty()

    try:
        with log_path.open("w") as log_fh:
            with subprocess.Popen(
                cmd,
                cwd=cwd,
                stdout=slave_fd,
                stderr=slave_fd,
            ) as process:
                os.close(slave_fd)
                slave_fd = None
                while True:
                    try:
                        data = os.read(master_fd, 4096)
                    except OSError as error:
                        if error.errno == errno.EIO:
                            break
                        raise
                    if not data:
                        break

                    text = decoder.decode(data)
                    sys.stdout.write(text)
                    sys.stdout.flush()
                    pending += text

                    while "\n" in pending:
                        line, pending = pending.split("\n", 1)
                        line = f"{plain_command_output(line)}\n"
                        log_fh.write(line)
                        log_fh.flush()
                        output_tail.append(line)

                pending += decoder.decode(b"", final=True)
                if pending:
                    line = plain_command_output(pending)
                    log_fh.write(line)
                    log_fh.flush()
                    output_tail.append(line)
                returncode = process.wait()
    finally:
        os.close(master_fd)
        if slave_fd is not None:
            os.close(slave_fd)

    if returncode:
        raise subprocess.CalledProcessError(
            returncode,
            cmd,
            output="".join(output_tail),
        )


def command_args(command: str, options: str = "") -> list[str]:
    """Split configured command strings without invoking a shell."""
    return [*shlex.split(command), *shlex.split(options)]


def file_md5(path: Path) -> str:
    """Return the hexadecimal MD5 digest for a file."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as input_fh:
        for chunk in iter(lambda: input_fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def get_read_geometry(run_dir):
    with (run_dir / "Stats" / "Stats.json").open() as stats_file:
        stats_json = json.load(stats_file)
    lane_info = stats_json["ReadInfosForLanes"][0].get("ReadInfos", None)
    if not lane_info:
        return "Read geometry could not be automatically determined."
    R1 = None
    R2 = None
    for read in lane_info:
        if read["IsIndexedRead"]:
            continue
        elif read["Number"] == 1:
            R1 = int(read["NumCycles"])
        elif read["Number"] == 2:
            R2 = int(read["NumCycles"])
    if R1 and R2:
        return f"Paired end - forward read length (R1): {R1}, reverse read length (R2): {R2}"
    elif R1 and not R2:
        return f"Single end - read length (R1): {R1}"
    elif not R1 and not R2:
        return "Read geometry could not be automatically determined."


def to_dirs(files):
    s = set()
    for f in files:
        d = str(f.parent)
        if d.split("/")[-1].startswith("GCF-"):
            s.add(d)
        else:
            s.add(d[: d.rfind("/")])
    return s


def get_sequencer(run_id):
    return SEQUENCERS.get(run_id.split("_")[1], "Sequencer could not be automatically determined.")


def _md5_filename(path):
    """Encode filenames using GNU md5sum's escaped-line convention."""
    name = str(path)
    escaped = any(char in name for char in "\\\n\r")
    name = name.replace("\\", "\\\\").replace("\n", "\\n").replace("\r", "\\r")
    return ("\\" if escaped else "", name)


def _complete_fastq_manifest(manifest, relative_paths):
    """Check syntax and exact file coverage without reading FASTQ contents."""
    expected = {_md5_filename(path) for path in relative_paths}
    seen = set()
    try:
        with manifest.open(encoding="utf-8", newline="") as handle:
            for line in handle:
                match = re.fullmatch(r"(\\?)[0-9a-fA-F]{32} [ *]([^\n]*)\n", line)
                if match is None or match.groups() not in expected or match.groups() in seen:
                    return False
                seen.add(match.groups())
    except (FileNotFoundError, UnicodeError):
        return False
    return seen == expected


def _write_fastq_manifest(manifest, fastqs, output_path):
    """Publish a complete manifest atomically; never expose a partial write."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=manifest.parent,
            prefix=f".{manifest.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            with ThreadPoolExecutor(max_workers=5) as executor:
                checksums = executor.map(file_md5, fastqs)
                for fastq, checksum in zip(fastqs, checksums, strict=True):
                    prefix, name = _md5_filename(fastq.relative_to(output_path))
                    handle.write(f"{prefix}{checksum}  {name}\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, manifest)
        directory_fd = os.open(manifest.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def md5sum_worker(cfg, *, force=False):
    """Generate demultiplexing checksums, or repair legacy manifests once."""
    for project in sorted(get_project_names(get_project_dirs(cfg))):
        manifest = cfg.output_path / f"md5sum_{project}_fastq.txt"
        fastqs = sorted((cfg.output_path / project).rglob("*.fastq.gz"))
        relative_paths = [path.relative_to(cfg.output_path) for path in fastqs]
        if not force and _complete_fastq_manifest(manifest, relative_paths):
            continue
        log.info("[md5sum_worker] Generating FASTQ checksums for %s", project)
        try:
            _write_fastq_manifest(manifest, fastqs, cfg.output_path)
        except Exception as error:
            raise RuntimeError(
                f"FASTQ checksum generation failed for {project}: {error}"
            ) from error


def md5sum_archive(archive_path: Path):
    base = archive_path.with_suffix("")  # remove .7za suffix
    md5_file = archive_path.parent / f"md5sum_{base.name}_archive.txt"

    if not md5_file.exists() or archive_path.stat().st_mtime > md5_file.stat().st_mtime:
        cmd = ["md5sum", archive_path.name]
        log.info(f"[md5sum_worker] Processing {shlex.join(cmd)}")
        with md5_file.open("w") as output_fh:
            subprocess.check_call(cmd, stdout=output_fh, cwd=archive_path.parent)


def md5sum_archive_worker(cfg):
    output_path = Path(cfg.output_path)
    archives = list(output_path.glob("*.7za"))

    with mp.Pool() as pool:
        pool.map(md5sum_archive, archives)


def multiqc_stats(cfg):
    in_confs = sorted(cfg.output_path.glob(".multiqc_config*.yaml"))
    if not in_confs:
        raise RuntimeError(
            "Missing analysis-generated MultiQC configuration in "
            f"{cfg.output_path}; recover with flowcell-manager rerun "
            f"{cfg.output_path.name} --from analysis"
        )

    cwd = cfg.output_path / "Stats"

    shutil.copy2(cfg.run.flowcell_path / "RunInfo.xml", cfg.output_path / "RunInfo.xml")
    # Illumina sequencer update - RunParameters.xml -> runParameters.xml
    shutil.copy2(
        list(cfg.run.flowcell_path.glob("[Rr]unParameters.xml"))[0],
        cfg.output_path / "RunParameters.xml",
    )

    # Illumina interop
    out_f = cfg.output_path / "Stats" / "interop_summary.csv"
    log.info(f"[multiqc_worker] Interop summary on {cfg.output_path}")
    run_interop_csv("interop_summary", cfg.output_path, out_f, cwd)

    prepare_index_metrics(cfg.output_path)
    out_f = cfg.output_path / "Stats" / "interop_index-summary.csv"
    log.info(f"[multiqc_worker] Interop index summary on {cfg.output_path}")
    run_interop_csv("interop_index-summary", cfg.output_path, out_f, cwd)

    samples_custom_data = dict()
    for c in in_confs:
        with c.open() as c_fh:
            mqc_conf = yaml.load(c_fh, Loader=yaml.FullLoader)
        samples_custom_data.update(mqc_conf["custom_data"]["general_statistics"]["data"])

    # use one of the existing multiqc_config.yaml as template
    with in_confs[0].open() as in_conf_fh:
        mqc_conf = yaml.load(in_conf_fh, Loader=yaml.FullLoader)

    projects = sorted(get_project_names(get_project_dirs(cfg)))
    pnames = ", ".join(projects)
    mqc_conf["title"] = pnames
    mqc_conf["intro_text"] = (
        "This report is generated for projects run at Genomics Core Facility, NTNU, Trondheim. The results are reported per sample."
    )
    mqc_conf["custom_data"]["general_statistics"]["data"] = samples_custom_data

    conf_pth = cfg.output_path / "Stats" / ".multiqc_config.yaml"
    with conf_pth.open("w+") as out_conf_fh:
        yaml.dump(mqc_conf, out_conf_fh)

    force_bcl2fastq = os.environ.get("FORCE_BCL2FASTQ", None)
    demultiplexer_module = "bcl2fastq" if force_bcl2fastq else "bclconvert"

    multiqc_opts = cfg.static.commands["multiqc_options"]
    pname = pnames.replace(", ", "_")
    multiqc_out = cfg.output_path / "Stats" / f"sequencer_stats_{pname}.html"

    cmd = command_args("multiqc", multiqc_opts)
    cmd.extend(
        [
            "--config",
            str(conf_pth),
            str(cfg.output_path / "Stats"),
            "--filename",
            str(multiqc_out),
            "-m",
            "interop",
            "-m",
            demultiplexer_module,
        ]
    )
    log.info(f"[multiqc_worker] Processing {cfg.output_path}")

    if os.environ.get("BFQ_TEST", None) and not force_bcl2fastq:
        if not (cfg.output_path / "Stats" / "Demultiplex_Stats.csv").exists():
            log.warning(
                "BFQ-TEST: Testflowcell was generated with bcl2fastq but environment is configured for bcl-convert. Using bcl2fastq paths and mqc modules."
            )
            cmd = [arg.replace("Reports", "Stats") for arg in cmd]
            cmd[cmd.index("bclconvert")] = "bcl2fastq"
            log.info(f"[multiqc_worker] Running: {shlex.join(cmd)}")

    subprocess.check_call(cmd, cwd=cwd)
    return projects


def generate_password(cfg, prefix: str) -> str:
    """
    Generate a one-time archive password and write it to file.

    Parameters
    ----------
    cfg : PipelineConfig
        Global configuration object.
    prefix : str
        Used to name the password file, e.g. 'encryption.<prefix>'.

    Returns
    -------
    str
        The generated password string.
    """
    pw = subprocess.check_output(
        ["xkcdpass", "-n", "5", "-d", "-", "-v", "[a-z]"], text=True
    ).strip("\n")
    pw_file = cfg.output_path / f"encryption.{prefix}"
    pw_file.write_text(f"{pw}\n", encoding="utf-8")
    return pw


def archive_worker(cfg):
    project_dirs = get_project_dirs(cfg)
    pnames = get_project_names(project_dirs)
    run_date = str(cfg.run.run_id).split("_")[0]

    for p in pnames:
        # ------------------------------------------------------------------ #
        # Archive FASTQ
        # ------------------------------------------------------------------ #
        archive_fastq = cfg.output_path / f"{p}_{run_date}.7za"
        if archive_fastq.exists():
            archive_fastq.unlink()

        pw = generate_password(cfg, p) if cfg.run.sensitive else None
        report_dir = cfg.output_path / "Reports"
        archive_inputs = [cfg.output_path / p, cfg.output_path / "Stats"]
        if report_dir.exists():
            archive_inputs.append(report_dir)
        archive_inputs.extend(sorted(cfg.output_path.glob("Undetermined*.fastq.gz")))
        archive_inputs.extend(
            [
                cfg.output_path / f"{p}_samplesheet.tsv",
                cfg.output_path / "SampleSheet.csv",
                cfg.output_path / "Sample-Submission-Form.xlsx",
                cfg.output_path / f"md5sum_{p}_fastq.txt",
            ]
        )

        if cfg.run.libprep and "10X Genomics" in cfg.run.libprep:
            extra = cfg.run.run_id.split("_")[-1][1:]
            archive_inputs.append(cfg.output_path / extra)

        cmd = ["7za", "a"]
        if pw:
            cmd.append(f"-p{pw}")
        cmd.extend([str(archive_fastq), *(str(path) for path in archive_inputs)])

        log.info(f"[archive_worker] Zipping {archive_fastq}")
        subprocess.check_call(cmd)

        # ------------------------------------------------------------------ #
        # Archive pipeline output (QC)
        # ------------------------------------------------------------------ #
        qc_archive = cfg.output_path / f"QC_{p}_{run_date}.7za"
        if qc_archive.exists():
            qc_archive.unlink()

        pw = generate_password(cfg, f"QC_{p}") if cfg.run.sensitive else None
        tmp_dir = Path(os.environ["TMPDIR"])
        qc_dir = tmp_dir / f"{p}_{run_date}" / "data" / "tmp" / cfg.run.pipeline / "bfq"

        # Native 7-Zip follows symlinks by default; -snl would store the links themselves.
        cmd = ["7za", "a"]
        if pw:
            cmd.append(f"-p{pw}")
        cmd.extend([str(qc_archive), str(qc_dir)])

        log.info(f"[archive_worker] Archiving QC output → {qc_archive}\n")
        subprocess.check_call(cmd)


def get_project_names(dirs):
    gcf = set()
    for d in dirs:
        for catalog in d.split("/"):
            if catalog.startswith("GCF-"):
                gcf.add(catalog)
    return gcf


def get_project_dirs(cfg):
    """
    Find project directories under cfg.output_path containing FASTQ files.

    Searches one and two levels deep for *.fastq.gz files,
    then extracts their parent directories via to_dirs().
    """
    fastq_paths = list(cfg.output_path.glob("*/*.fastq.gz"))
    fastq_paths += list(cfg.output_path.glob("*/*/*.fastq.gz"))
    return to_dirs(fastq_paths)


def post_workflow(project_id, base_dir, pipeline):
    run_date = str(base_dir.name).split("_")[0]
    analysis_dir = Path(os.environ["TMPDIR"]) / f"{project_id}_{run_date}"
    bfq_dir = analysis_dir / "data" / "tmp" / pipeline / "bfq"
    odir = base_dir / f"QC_{project_id}" / "bfq"
    odir.mkdir(parents=True, exist_ok=True)

    shutil.copytree(bfq_dir, odir, symlinks=False, dirs_exist_ok=True)

    # Copy sample info
    shutil.copy2(
        analysis_dir / "data" / "tmp" / "sample_info.tsv",
        base_dir / f"{project_id}_samplesheet.tsv",
    )

    return True


def full_align(cfg):
    selection = select_workflow(cfg)
    project_names = get_project_names(get_project_dirs(cfg))
    run_date = str(cfg.output_path.name).split("_")[0]
    for p in sorted(project_names):
        analysis_dir = Path(os.environ["TMPDIR"]) / f"{p}_{run_date}"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        (analysis_dir / "src").mkdir(parents=True, exist_ok=True)
        (analysis_dir / "data").mkdir(parents=True, exist_ok=True)
        log.info(f"Setting up analysis for {analysis_dir}")

        # os.chdir(analysis_dir)

        src = workflow_config.AUTHORITATIVE_CONFIG.parent
        dst = analysis_dir / "src" / "gcf-workflows"

        # copy snakemake pipeline
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
        # copytree sees mutable working-tree files; overwrite the config with the
        # exact bytes captured for this execution, including uncommitted edits.
        selection.config.write(dst / "libprep.config")

        machine = get_sequencer(cfg.run.run_id)
        # create config.yaml
        cmd = [
            "/opt/conda/bin/python",
            "/opt/conda/bin/configmaker.py",
            str(cfg.output_path),
            "-p",
            str(p),
            "--libkit",
            selection.kit,
            "--machine",
            str(machine),
            "--expected-validation-version",
            VALIDATOR_VERSION,
            "--libprep-config",
            str(dst / "libprep.config"),
            "--libprep-sha256",
            selection.config.sha256,
            "--libprep-entry",
            selection.entry,
            "--expected-read-geometry",
            *(str(n) for n in selection.read_geometry),
        ]
        if (analysis_dir / "data/raw/fastq").exists():
            cmd.append("--skip-create-fastq-dir")
        subprocess.check_call(cmd, cwd=analysis_dir)
        discovery_summary = analysis_dir / "configmaker.analysis-summary.json"
        if discovery_summary.is_file():
            shutil.copy2(discovery_summary, cfg.output_path / f"configmaker-analysis-{p}.json")

        # run snakemake pipeline
        cmd = [
            "snakemake",
            "--use-singularity",
            "--singularity-prefix",
            os.environ["SINGULARITY_CACHEDIR"],
            "--cores",
            "32",
            "--scheduler",
            "greedy",
            "-p",
            "multiqc_report",
        ]
        snakemake_log = cfg.static.paths.log_dir / f"{cfg.run.run_id}_{p}_snakemake.log"
        run_logged_command(cmd, cwd=analysis_dir, log_path=snakemake_log)

        # copy report
        shutil.copy2(
            analysis_dir / "data" / "tmp" / cfg.run.pipeline / "bfq" / f"multiqc_{p}.html",
            cfg.output_path / f"multiqc_{p}_{run_date}.html",
        )

        # if additional html reports exists (single cell), copy
        extra_html = (analysis_dir / "data" / "tmp" / cfg.run.pipeline / "bfq" / "summaries").glob(
            "all_samples*.html"
        )
        extra_html = list(extra_html)
        if extra_html:
            shutil.copy2(
                extra_html[0], cfg.output_path / f"all_samples_web_summary_{p}_{run_date}.html"
            )

        # Copy sample info
        shutil.copy2(
            analysis_dir / "data" / "tmp" / "sample_info.tsv",
            cfg.output_path / f"{p}_samplesheet.tsv",
        )

        # copy mqc_config
        shutil.copy2(
            analysis_dir / "data" / "tmp" / cfg.run.pipeline / "bfq" / ".multiqc_config.yaml",
            cfg.output_path / f".multiqc_config_{p}.yaml",
        )

    # os.chdir(old_wd)
    return True


def _disk_usage_message(cfg):
    """Build the operational disk-usage summary used by completion mail."""
    message = ""
    sources = [("output", cfg.output_path.parent)]
    if cfg.run.flowcell_path is not None:
        sources.append(("instruments", cfg.run.flowcell_path.parent))
    for label, path in sources:
        try:
            total, _used, free = shutil.disk_usage(path)
        except OSError:
            message += f"Current free space for {label}: unavailable\n<br>"
            continue
        total /= 1024**3
        free /= 1024**3
        message += (
            f"Current free space for {label}: {free:.0f} of {total:.0f} GiB "
            f"({100 * free / total:5.2f}%)\n<br>"
        )
    return message


def analysis_steps():
    """Run work invalidated by the public analysis restart boundary."""
    cfg = PipelineConfig.get()
    full_align(cfg)


def reporting_steps():
    """Generate reporting products while preserving completed workflow results."""
    cfg = PipelineConfig.get()
    projects = multiqc_stats(cfg)
    cfg.to_file(cfg.output_path / "bcl2fastq.ini")
    return projects


def postMakeSteps():
    """Compatibility wrapper for callers that still expect the combined operation."""
    analysis_steps()
    reporting_steps()
    return _disk_usage_message(PipelineConfig.get())


def finalize():
    cfg = PipelineConfig.get()
    # zip arhive
    archive_worker(cfg)
    # md5sum archive
    md5sum_archive_worker(cfg)
