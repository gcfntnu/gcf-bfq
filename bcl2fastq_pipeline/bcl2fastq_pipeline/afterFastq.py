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
from contextlib import nullcontext
from pathlib import Path

from configmaker.configmaker import SEQUENCERS
from configmaker.validation import VALIDATOR_VERSION

from bcl2fastq_pipeline import (
    analysis_snapshots,
    processing_times,
    project_corrections,
    project_scheduling,
    undetermined_fastqs,
    workflow_config,
)
from bcl2fastq_pipeline.config import PipelineConfig
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
    stats = run_dir / "Stats" / "Stats.json"
    if not stats.exists() and not stats.is_symlink() and (run_dir / "RunInfo.xml").is_file():
        geometry = workflow_config.find_read_geometry([run_dir])
        lane_info = [
            {"Number": number, "NumCycles": cycles, "IsIndexedRead": False}
            for number, cycles in enumerate(geometry, start=1)
        ]
    else:
        with stats.open() as stats_file:
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


def md5sum_worker(cfg, *, force=False, store=None):
    """Generate demultiplexing checksums, or repair legacy manifests once.

    A reused manifest retains its timing; checking coverage is not a new hashing
    execution. When repair is needed, time the enclosing generation once, not
    each project or parallel worker.
    """
    jobs = []
    for project in sorted(get_project_names(get_project_dirs(cfg))):
        manifest = cfg.output_path / f"md5sum_{project}_fastq.txt"
        fastqs = sorted(
            path
            for path in (cfg.output_path / project).rglob("*.fastq.gz")
            if not undetermined_fastqs.is_undetermined(path)
        )
        relative_paths = [path.relative_to(cfg.output_path) for path in fastqs]
        if not force and _complete_fastq_manifest(manifest, relative_paths):
            continue
        jobs.append((project, manifest, fastqs))
    timer = (
        processing_times.measure(store, cfg.run.run_id, "fastq_checksums")
        if store is not None and (jobs or force)
        else nullcontext()
    )
    with timer:
        for project, manifest, fastqs in jobs:
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
    if getattr(getattr(cfg, "run", None), "selected_projects", None) is not None:
        archives = [
            p
            for p in archives
            if any(
                p.name.startswith((f"{project}_", f"QC_{project}_"))
                for project in cfg.run.selected_projects
            )
        ]

    with mp.Pool() as pool:
        pool.map(md5sum_archive, archives)


def multiqc_stats(cfg):
    """Compatibility entry point for metadata-independent sequencing reporting."""
    from bcl2fastq_pipeline import sequencing_qc  # noqa: PLC0415

    return sequencing_qc.generate(cfg)["projects"]


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


def _run_archive_command(cmd, cfg, project, kind):
    log_path = cfg.static.paths.log_dir / f"{cfg.run.run_id}_{project}_archive_{kind}.log"
    try:
        run_logged_command(cmd, cwd=cfg.output_path, log_path=log_path)
    except subprocess.CalledProcessError as error:
        # The password belongs in encryption.<project>, not failure reports/state.
        error.cmd = ["-p<redacted>" if arg.startswith("-p") else arg for arg in cmd]
        raise


def archive_worker(cfg):
    project_dirs = get_project_dirs(cfg)
    pnames = get_project_names(project_dirs)
    run_date = str(cfg.run.run_id).split("_")[0]

    for p in project_scheduling.project_order(cfg, pnames):
        project_scheduling.record_start(cfg, "archiving", p)
        project_config = workflow_config.delivery_project_config(cfg, p)
        sensitive = project_config.sensitive if project_config else cfg.run.sensitive
        # ------------------------------------------------------------------ #
        # Archive FASTQ
        # ------------------------------------------------------------------ #
        archive_fastq = cfg.output_path / f"{p}_{run_date}.7za"
        (cfg.output_path / f"md5sum_{p}_{run_date}_archive.txt").unlink(missing_ok=True)
        if archive_fastq.exists():
            archive_fastq.unlink()

        pw = generate_password(cfg, p) if sensitive else None
        report_dir = cfg.output_path / "Reports"
        archive_inputs = [cfg.output_path / p, cfg.output_path / "Stats"]
        archive_inputs.extend(
            path
            for path in sorted(cfg.output_path.glob("sequencer_stats_*.html"))
            if not project_config or project_corrections.belongs_to_project(path.name, [p])
        )
        if report_dir.exists():
            archive_inputs.append(report_dir)
        archive_inputs.extend(
            path
            for path in sorted(cfg.output_path.glob("Undetermined*.fastq.gz"))
            if not undetermined_fastqs.is_undetermined(path)
        )
        archive_inputs.extend(
            [
                cfg.output_path / f"{p}_samplesheet.tsv",
                project_config.sample_sheet
                if project_config
                else cfg.output_path / "SampleSheet.csv",
                project_config.submission_form
                if project_config
                else cfg.output_path / "Sample-Submission-Form.xlsx",
                cfg.output_path / f"md5sum_{p}_fastq.txt",
            ]
        )
        # A project delivery retains the exact combined conversion inputs and row
        # mappings that contributed its FASTQs, including source hashes.
        if project_config:
            from bcl2fastq_pipeline import demultiplexing  # noqa: PLC0415

            demux_provenance = demultiplexing.delivery_provenance(cfg, p)
            if demux_provenance:
                # 7-Zip flattens absolute input paths to their basename. Keep
                # this relative path so it cannot collide with the FASTQ folder.
                archive_inputs.append(demux_provenance.relative_to(cfg.output_path))

        libprep = project_config.libprep if project_config else cfg.run.libprep
        if libprep and "10X Genomics" in libprep:
            extra = cfg.output_path / cfg.run.run_id.split("_")[-1][1:]
            try:
                extra.lstat()
            except FileNotFoundError:
                # Restored FASTQs need not include the extra mkfastq directory.
                log.info("[archive_worker] Optional 10X archive input absent: %s", extra)
            else:
                # Keep dangling links and unreadable inputs visible to the archiver.
                archive_inputs.append(extra)

        cmd = ["7za", "a"]
        if pw:
            cmd.append(f"-p{pw}")
        cmd.extend(
            [
                str(archive_fastq),
                *undetermined_fastqs.archive_exclusions(),
                *(str(path) for path in archive_inputs),
            ]
        )

        log.info(f"[archive_worker] Zipping {archive_fastq}")
        _run_archive_command(cmd, cfg, p, "fastq")

        # ------------------------------------------------------------------ #
        # Archive pipeline output (QC)
        # ------------------------------------------------------------------ #
        qc_archive = cfg.output_path / f"QC_{p}_{run_date}.7za"
        (cfg.output_path / f"md5sum_QC_{p}_{run_date}_archive.txt").unlink(missing_ok=True)
        if qc_archive.exists():
            qc_archive.unlink()

        pw = generate_password(cfg, f"QC_{p}") if sensitive else None
        qc_dir = (
            workflow_config.project_workdir(cfg, p)
            / "data"
            / "tmp"
            / workflow_config.project_workflow(cfg, p)
            / "bfq"
        )

        # Native 7-Zip follows symlinks by default; -snl would store the links themselves.
        cmd = ["7za", "a"]
        if pw:
            cmd.append(f"-p{pw}")
        cmd.extend([str(qc_archive), *undetermined_fastqs.archive_exclusions(), str(qc_dir)])

        log.info(f"[archive_worker] Archiving QC output → {qc_archive}\n")
        _run_archive_command(cmd, cfg, p, "qc")


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
    dirs = to_dirs(fastq_paths)
    if getattr(getattr(cfg, "run", None), "selected_projects", None) is not None:
        dirs = {d for d in dirs if get_project_names([d]) & cfg.run.selected_projects}
    return dirs


def post_workflow(project_id, base_dir, pipeline, *, analysis_dir=None):
    run_date = str(base_dir.name).split("_")[0]
    analysis_dir = (
        Path(analysis_dir)
        if analysis_dir
        else Path(os.environ["TMPDIR"]) / f"{project_id}_{run_date}"
    )
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


def gpu_enabled():
    """Opt into NVIDIA passthrough for every BFQ Snakemake container job."""
    value = os.environ.get("BFQ_GPU", "0").strip().lower()
    if value in {"1", "true", "yes", "on"}:
        return True
    if value in {"0", "false", "no", "off"}:
        return False
    raise ValueError("BFQ_GPU must be 1 (NVIDIA GPU enabled) or 0 (CPU only)")


def snakemake_command(*, resume=False, target="multiqc_report", cores=32):
    command = [
        "snakemake",
        "--use-singularity",
        "--singularity-prefix",
        os.environ["SINGULARITY_CACHEDIR"],
        "--cores",
        str(cores),
        "--scheduler",
        "greedy",
        "-p",
        target,
    ]
    if gpu_enabled():
        command[1:1] = ["--singularity-args=--nv"]
    if resume:
        # Keep normal failed-job cleanup: Snakemake 9.7.1 can mark a
        # failed job complete when --keep-incomplete retains partial files.
        command[1:1] = ["--rerun-incomplete"]
    return command


def full_align(cfg, *, resume=None):
    workflow_config.prepare_projects(cfg)
    selection = None
    if not resume and not cfg.run.project_configs:
        selection = select_workflow(cfg)
    project_names = get_project_names(get_project_dirs(cfg))
    run_date = str(cfg.run.run_id).split("_")[0]
    workdirs = {}
    if resume and set(project_names) != set(resume["projects"]):
        raise RuntimeError("Resume projects differ from available FASTQs; inspect and requeue")
    for p in project_scheduling.project_order(cfg, project_names):
        project_scheduling.record_start(cfg, "analysis", p)
        if not resume and cfg.run.project_configs:
            selection = select_workflow(cfg, p)
        pipeline = resume["projects"][p]["workflow"] if resume else selection.workflow
        cfg.run.project_workflows[p] = pipeline
        cfg.run.pipeline = pipeline
        analysis_dir = workflow_config.project_workdir(cfg, p)
        if resume:
            entry = resume["projects"][p]
            analysis_dir = Path(entry["identity"]["path"])
            workdirs[p] = entry["identity"]
            log.info("Resuming retained analysis in %s", analysis_dir)
        else:
            analysis_dir.mkdir(parents=True, exist_ok=True)
            (analysis_dir / "src").mkdir(parents=True, exist_ok=True)
            (analysis_dir / "data").mkdir(parents=True, exist_ok=True)
            log.info(f"Setting up analysis for {analysis_dir}")

            workdirs[p] = analysis_snapshots.identify_workdir(analysis_dir, cfg.run.run_id, p)

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
            ]
            project_config = getattr(cfg.run, "project_configs", {}).get(p)
            if project_config:
                retained_inputs = analysis_dir / "inputs"
                retained_inputs.mkdir(exist_ok=True)
                for original in (project_config.sample_sheet, project_config.submission_form):
                    shutil.copy2(original, retained_inputs / original.name)
                (retained_inputs / "project.json").write_text(
                    json.dumps(
                        {
                            "project": p,
                            "sample_sheet": project_config.sample_sheet.name,
                            "submission_form": project_config.submission_form.name,
                            "custom": project_config.custom,
                            "scheduling": (
                                cfg.run.project_schedule["projects"][p]
                                if cfg.run.project_schedule
                                else None
                            ),
                            "libprep_selection": selection.diagnostics(),
                            "inputs": {
                                original.name: hashlib.sha256(original.read_bytes()).hexdigest()
                                for original in (
                                    project_config.sample_sheet,
                                    project_config.submission_form,
                                )
                            },
                        },
                        indent=2,
                    )
                    + "\n",
                    encoding="utf-8",
                )
                cmd.extend(
                    [
                        "--sample-sheet",
                        str(project_config.sample_sheet),
                        "--sample-submission-form",
                        str(project_config.submission_form),
                        "--effective-read-geometry",
                        *(str(n) for n in selection.read_geometry),
                    ]
                )
            cmd.extend(["--expected-read-geometry", *(str(n) for n in selection.read_geometry)])
            if (analysis_dir / "data/raw/fastq").exists():
                cmd.append("--skip-create-fastq-dir")
            subprocess.check_call(cmd, cwd=analysis_dir)
            discovery_summary = analysis_dir / "configmaker.analysis-summary.json"
            if discovery_summary.is_file():
                shutil.copy2(discovery_summary, cfg.output_path / f"configmaker-analysis-{p}.json")

        # run snakemake pipeline
        cmd = snakemake_command(resume=bool(resume))
        snakemake_log = cfg.static.paths.log_dir / f"{cfg.run.run_id}_{p}_snakemake.log"
        run_logged_command(cmd, cwd=analysis_dir, log_path=snakemake_log)

        if resume:
            # Publish current BFQ products only after successful scientific work.
            # Queueing and failed analysis leave prior delivery files untouched.
            qc = cfg.output_path / f"QC_{p}"
            if qc.is_symlink():
                qc.unlink()
            elif qc.exists():
                shutil.rmtree(qc)
            for name in (
                f"all_samples_web_summary_{p}_{run_date}.html",
                f"configmaker-analysis-{p}.json",
            ):
                (cfg.output_path / name).unlink(missing_ok=True)
            discovery = analysis_dir / "configmaker.analysis-summary.json"
            if discovery.is_file():
                shutil.copy2(discovery, cfg.output_path / f"configmaker-analysis-{p}.json")
            post_workflow(p, cfg.output_path, pipeline, analysis_dir=analysis_dir)

        # copy report
        shutil.copy2(
            analysis_dir / "data" / "tmp" / pipeline / "bfq" / f"multiqc_{p}.html",
            cfg.output_path / f"multiqc_{p}_{run_date}.html",
        )

        # if additional html reports exists (single cell), copy
        extra_html = (analysis_dir / "data" / "tmp" / pipeline / "bfq" / "summaries").glob(
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
            analysis_dir / "data" / "tmp" / pipeline / "bfq" / ".multiqc_config.yaml",
            cfg.output_path / f".multiqc_config_{p}.yaml",
        )
        cfg.run.project_workdirs[p] = analysis_dir

    return workdirs


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


def analysis_steps(*, resume=None):
    """Run work invalidated by the public analysis restart boundary."""
    cfg = PipelineConfig.get()
    return full_align(cfg, resume=resume) if resume else full_align(cfg)


def reporting_steps():
    """Generate reporting products while preserving completed workflow results."""
    cfg = PipelineConfig.get()
    # Early reports belong to the conversion execution. Never regenerate them
    # merely because analysis/reporting is rerun. Legacy runs have no early event.
    from bcl2fastq_pipeline.state import FlowcellStateStore  # noqa: PLC0415

    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    state = store.read(cfg.run.run_id) if store.exists(cfg.run.run_id) else {}
    if not state.get("sequencing_qc"):
        multiqc_stats(cfg)
    projects = sorted(get_project_names(get_project_dirs(cfg)))
    cfg.to_file(cfg.output_path / "bcl2fastq.ini")
    return projects


def postMakeSteps():
    """Compatibility wrapper for callers that still expect the combined operation."""
    analysis_steps()
    reporting_steps()
    return _disk_usage_message(PipelineConfig.get())


def finalize(*, store=None):
    cfg = PipelineConfig.get()
    with (
        processing_times.measure(store, cfg.run.run_id, "archiving")
        if store is not None
        else nullcontext()
    ):
        archive_worker(cfg)
    with (
        processing_times.measure(store, cfg.run.run_id, "archive_checksums")
        if store is not None
        else nullcontext()
    ):
        md5sum_archive_worker(cfg)
