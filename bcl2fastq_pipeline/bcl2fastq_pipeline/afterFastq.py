"""
This file includes code that actually runs FastQC and any other tools after the fastq files have actually been made. This uses a pool of workers to process each request.
"""

import json
import logging
import multiprocessing as mp
import os
import shutil
import subprocess

from pathlib import Path

import yaml

from configmaker.configmaker import SEQUENCERS

from bcl2fastq_pipeline.config import PipelineConfig

log = logging.getLogger(__name__)


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


def md5sum_worker(cfg):
    project_dirs = get_project_dirs(cfg)
    pnames = get_project_names(project_dirs)
    for p in pnames:
        md_path = Path(f"md5sum_{p}_fastq.txt")
        if not (cfg.output_path / md_path).exists():
            cmd = (
                f"find {p} -type f -name '*.fastq.gz' | parallel -j 5 md5sum > md5sum_{p}_fastq.txt"
            )
            log.info(f"[md5sum_worker] Processing {cfg.output_path}/{p}")
            subprocess.check_call(cmd, shell=True, cwd=cfg.output_path)


def md5sum_archive(archive_path: Path):
    base = archive_path.with_suffix("")  # remove .7za suffix
    md5_file = archive_path.parent / f"md5sum_{base.name}_archive.txt"

    if not md5_file.exists() or archive_path.stat().st_mtime > md5_file.stat().st_mtime:
        cmd = f"md5sum {archive_path.name} > {md5_file.name}"
        log.info(f"[md5sum_worker] Processing {cmd}")
        subprocess.check_call(cmd, shell=True, cwd=archive_path.parent)


def md5sum_archive_worker(cfg):
    output_path = Path(cfg.output_path)
    archives = list(output_path.glob("*.7za"))

    with mp.Pool() as pool:
        pool.map(md5sum_archive, archives)


def multiqc_stats(cfg):
    cwd = cfg.output_path / "Stats"

    shutil.copy2(cfg.run.flowcell_path / "RunInfo.xml", cfg.output_path / "RunInfo.xml")
    # Illumina sequencer update - RunParameters.xml -> runParameters.xml
    shutil.copy2(
        list(cfg.run.flowcell_path.glob("[Rr]unParameters.xml"))[0],
        cfg.output_path / "RunParameters.xml",
    )

    # Illumina interop
    out_f = cfg.output_path / "Stats" / "interop_summary.csv"
    cmd = f"interop_summary {cfg.output_path} --csv=1 > {out_f}"
    log.info(f"[multiqc_worker] Interop summary on {cfg.output_path}")
    subprocess.check_call(cmd, shell=True, cwd=cwd)

    out_f = cfg.output_path / "Stats" / "interop_index-summary.csv"
    cmd = f"interop_index-summary {cfg.output_path} --csv=1 > {out_f}"
    log.info(f"[multiqc_worker] Interop index summary on {cfg.output_path}")
    subprocess.check_call(cmd, shell=True, cwd=cwd)

    in_confs = list(cfg.output_path.glob(".multiqc_config*.yaml"))
    samples_custom_data = dict()
    for c in in_confs:
        with c.open() as c_fh:
            mqc_conf = yaml.load(c_fh, Loader=yaml.FullLoader)
        samples_custom_data.update(mqc_conf["custom_data"]["general_statistics"]["data"])

    # use one of the existing multiqc_config.yaml as template
    with in_confs[0].open() as in_conf_fh:
        mqc_conf = yaml.load(in_conf_fh, Loader=yaml.FullLoader)

    pnames = get_project_names(get_project_dirs(cfg))
    pnames = ", ".join(pnames)
    mqc_conf["title"] = pnames
    mqc_conf["intro_text"] = (
        "This report is generated for projects run at Genomics Core Facility, NTNU, Trondheim. The results are reported per sample."
    )
    mqc_conf["custom_data"]["general_statistics"]["data"] = samples_custom_data

    conf_pth = cfg.output_path / "Stats" / ".multiqc_config.yaml"
    with conf_pth.open("w+") as out_conf_fh:
        yaml.dump(mqc_conf, out_conf_fh)

    modules = "-m interop "
    FORCE_BCL2FASTQ = os.environ.get("FORCE_BCL2FASTQ", None)
    modules += "-m bclconvert " if not FORCE_BCL2FASTQ else "-m bcl2fastq"

    multiqc_cmd = cfg.static.commands["multiqc_command"]
    multiqc_opts = cfg.static.commands["multiqc_options"]
    pname = pnames.replace(", ", "_")
    multiqc_out = cfg.output_path / "Stats" / f"sequencer_stats_{pname}.html"

    cmd = f"{multiqc_cmd} {multiqc_opts} --config {conf_pth} {cfg.output_path}/Stats --filename {multiqc_out} {modules}"
    log.info(f"[multiqc_worker] Processing {cfg.output_path}")

    if os.environ.get("BFQ_TEST", None) and not FORCE_BCL2FASTQ:
        if not (cfg.output_path / "Stats" / "Demultiplex_Stats.csv").exists():
            log.warning(
                "BFQ-TEST: Testflowcell was generated with bcl2fastq but environment is configured for bcl-convert. Using bcl2fastq paths and mqc modules."
            )
            cmd = cmd.replace("Reports", "Stats")
            cmd = cmd.replace("-m bclconvert", "-m bcl2fastq")
            log.info(f"[multiqc_worker] Running: {cmd}")

    subprocess.check_call(cmd, shell=True, cwd=cwd)


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
    pw = subprocess.check_output("xkcdpass -n 5 -d '-' -v '[a-z]'", shell=True).decode().strip("\n")
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
        opts = f"-p{pw}" if pw else ""

        report_dir = cfg.output_path / "Reports"
        if not report_dir.exists():
            report_dir = ""

        cmd = (
            f"7za a {opts} "
            f"{cfg.output_path}/{p}_{run_date}.7za "
            f"{cfg.output_path}/{p}/ "
            f"{cfg.output_path}/Stats "
            f"{report_dir} "
            f"{cfg.output_path}/Undetermined*.fastq.gz "
            f"{cfg.output_path}/{p}_samplesheet.tsv "
            f"{cfg.output_path}/SampleSheet.csv "
            f"{cfg.output_path}/Sample-Submission-Form.xlsx "
            f"{cfg.output_path}/md5sum_{p}_fastq.txt "
        )

        if cfg.run.libprep and "10X Genomics" in cfg.run.libprep:
            extra = cfg.run.run_id.split("_")[-1][1:]
            cmd += f" {cfg.output_path}/{extra}"

        log.info(f"[archive_worker] Zipping {archive_fastq}")
        subprocess.check_call(cmd, shell=True)

        # ------------------------------------------------------------------ #
        # Archive pipeline output (QC)
        # ------------------------------------------------------------------ #
        qc_archive = cfg.output_path / f"QC_{p}_{run_date}.7za"
        if qc_archive.exists():
            qc_archive.unlink()

        pw = generate_password(cfg, f"QC_{p}") if cfg.run.sensitive else None
        opts = f"-p{pw}" if pw else ""

        tmp_dir = Path(os.environ["TMPDIR"])
        qc_dir = tmp_dir / f"{p}_{run_date}" / "data" / "tmp" / cfg.run.pipeline / "bfq"
        flowdir = cfg.output_path

        cmd = f"7za a -l {opts} {flowdir}/QC_{p}_{run_date}.7za {qc_dir} "

        log.info(f"[archive_worker] Archiving QC output → {qc_archive}\n")
        subprocess.check_call(cmd, shell=True)


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
    # old_wd = Path.cwd()

    # os.chdir(os.environ["TMPDIR"])
    project_names = get_project_names(get_project_dirs(cfg))
    run_date = str(cfg.output_path.name).split("_")[0]
    for p in project_names:
        analysis_dir = Path(os.environ["TMPDIR"]) / f"{p}_{run_date}"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        (analysis_dir / "src").mkdir(parents=True, exist_ok=True)
        (analysis_dir / "data").mkdir(parents=True, exist_ok=True)
        log.info(f"Setting up analysis for {analysis_dir}")

        # os.chdir(analysis_dir)

        src = Path("/opt/gcf-workflows")
        dst = analysis_dir / "src" / "gcf-workflows"

        # copy snakemake pipeline
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)

        create_fastq = " --skip-create-fastq-dir" if Path("data/raw/fastq").exists() else ""
        machine = get_sequencer(cfg.run.run_id)
        # create config.yaml
        cmd = f"/opt/conda/bin/python /opt/conda/bin/configmaker.py {cfg.output_path} -p {p} --libkit '{cfg.run.libprep}' --machine '{machine}' {create_fastq}"
        subprocess.check_call(cmd, shell=True, cwd=analysis_dir)

        # run snakemake pipeline
        cmd = "snakemake --use-singularity --singularity-prefix $SINGULARITY_CACHEDIR --cores 32 --verbose -p multiqc_report"
        subprocess.check_call(cmd, shell=True, cwd=analysis_dir)

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
    (cfg.output_path / "analysis.made").write_text("")
    return True


# All steps that should be run after `make` go here
def postMakeSteps():
    """
    Current steps are:
      1) Run FastQC on each fastq.gz file
      2) Run md5sum on the files in each project directory
    Other steps could easily be added to follow those. Note that this function
    will try to use a pool of threads. The size of the pool is set by config.postMakeThreads
    """
    cfg = PipelineConfig.get()

    cfg.run.set_pipeline_from_yaml(Path("/opt/gcf-workflows/libprep.config"))

    # md5sum fastqs
    md5sum_worker(cfg)

    if not (cfg.output_path / "analysis.made").exists():
        full_align(cfg)
    else:
        log.info("Analysis already made")

    # multiqc_stats
    multiqc_stats(cfg)

    # disk usage
    tot, used, free = shutil.disk_usage(cfg.static.paths.output_dir)
    tot /= 1024**3  # Convert to GiB
    used /= 1024**3
    free /= 1024**3

    message = (
        f"Current free space for output: {free:.0f} of {tot:.0f} GiB "
        f"({100 * free / tot:5.2f}%)\n<br>"
    )

    tot, used, free = shutil.disk_usage(cfg.run.flowcell_path.parent)
    tot /= 1024**3  # Convert to GiB
    used /= 1024**3
    free /= 1024**3

    message += (
        f"Current free space for instruments: {free:.0f} of {tot:.0f} GiB "
        f"({100 * free / tot:5.2f}%)\n<br>\n<br>"
    )
    # save configfile to flowcell
    cfg.to_file(cfg.output_path / "bcl2fastq.ini")

    return message


def finalize():
    cfg = PipelineConfig.get()
    # zip arhive
    archive_worker(cfg)
    # md5sum archive
    md5sum_archive_worker(cfg)
