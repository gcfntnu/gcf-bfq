"""
This file includes code that actually runs FastQC and any other tools after the fastq files have actually been made. This uses a pool of workers to process each request.
"""

import glob
import multiprocessing as mp
import os
import os.path
import shutil
import subprocess
import syslog

import matplotlib

matplotlib.use("Agg")
import json
import re

from pathlib import Path

import yaml

from configmaker.configmaker import SEQUENCERS

from bcl2fastq_pipeline.config import PipelineConfig

config = None
localConfig = None


def set_config(config):
    """Store the configuration globally for this module."""
    global localConfig  # noqa: PLW0603
    localConfig = config


def get_config():
    """Return the current configuration object."""
    return localConfig


def get_gcf_name(fname):
    for name in fname.split("/"):
        if re.match(r"GCF-[0-9]{4}-[0-9]{3,}", name):
            return re.search(r"GCF-[0-9]{4}-[0-9]{3,}", name)[0]
    raise Exception(f"Unable to determine GCF project number for filename {fname}\n")


def to_dirs(files):
    s = set()
    for f in files:
        d = os.path.dirname(f)
        if d.split("/")[-1].startswith("GCF-"):
            s.add(d)
        else:
            s.add(d[: d.rfind("/")])
    return s


def get_sequencer(run_id):
    return SEQUENCERS.get(run_id.split("_")[1], "Sequencer could not be automatically determined.")


def get_read_geometry(run_dir):
    stats_file = open(f"{run_dir}/Stats/Stats.json")
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


def md5sum_worker(config):
    old_wd = os.getcwd()
    os.chdir(os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID")))
    project_dirs = get_project_dirs(config)
    pnames = get_project_names(project_dirs)
    for p in pnames:
        if os.path.exists(f"md5sum_{p}_fastq.txt"):
            continue
        cmd = "find {} -type f -name '*.fastq.gz' | parallel -j 5 md5sum > {}".format(
            p, f"md5sum_{p}_fastq.txt"
        )
        syslog.syslog(
            "[md5sum_worker] Processing {}\n".format(
                os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"), p)
            )
        )
        subprocess.check_call(cmd, shell=True)
    os.chdir(old_wd)


def md5sum_archive(archive_path: Path):
    base = archive_path.with_suffix("")  # remove .7za suffix
    md5_file = archive_path.parent / f"md5sum_{base.name}_archive.txt"

    if not md5_file.exists() or archive_path.stat().st_mtime > md5_file.stat().st_mtime:
        cmd = f"md5sum {archive_path} > {md5_file}"
        subprocess.check_call(cmd, shell=True)


def md5sum_archive_worker(cfg):
    output_path = Path(cfg.output_path)
    archives = list(output_path.glob("*.7za"))

    with mp.Pool() as pool:
        pool.map(md5sum_archive, archives)


def multiqc_stats(cfg):
    oldWd = os.getcwd()
    os.chdir(cfg.ouput_path / "Stats")

    shutil.copy2(cfg.run.flowcell_path / "RunInfo.xml", cfg.output_path / "RunInfo.xml")
    # Illumina sequencer update - RunParameters.xml -> runParameters.xml
    shutil.copy2(
        [f for f in cfg.run.flowcell_path.glob("[Rr]unParameters.xml")][0],
        cfg.output_path / "RunParameters.xml",
    )

    # Illumina interop
    out_f = cfg.output_path / "Stats" / "interop_summary.csv"
    cmd = f"interop_summary {cfg.output_path} --csv=1 > {out_f}"
    syslog.syslog(f"[multiqc_worker] Interop summary on {cfg.output_path}\n")
    subprocess.check_call(cmd, shell=True)

    out_f = cfg.output_path / "Stats" / "interop_index-summary.csv"
    cmd = f"interop_index-summary {cfg.output_path} --csv=1 > {out_f}"
    syslog.syslog(f"[multiqc_worker] Interop index summary on {cfg.output_path}\n")
    subprocess.check_call(cmd, shell=True)

    run_dir = os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"))

    in_confs = glob.glob(os.path.join(run_dir, ".multiqc_config*.yaml"))
    in_confs = [f for f in cfg.output_path.glob(".multiqc_config*.yaml")]
    samples_custom_data = dict()
    for c in in_confs:
        with open(c) as c_fh:
            mqc_conf = yaml.load(c_fh, Loader=yaml.FullLoader)
        samples_custom_data.update(mqc_conf["custom_data"]["general_statistics"]["data"])

    # use one of the existing multiqc_config.yaml as template
    with open(in_confs[0]) as in_conf_fh:
        mqc_conf = yaml.load(in_conf_fh, Loader=yaml.FullLoader)

    pnames = get_project_names(get_project_dirs(cfg))
    pnames = ", ".join(pnames)
    mqc_conf["title"] = pnames
    mqc_conf["intro_text"] = (
        "This report is generated for projects run at Genomics Core Facility, NTNU, Trondheim. The results are reported per sample."
    )
    mqc_conf["custom_data"]["general_statistics"]["data"] = samples_custom_data

    conf_name = cfg.output_path / "Stats" / ".multiqc_config.yaml"
    with open(conf_name, "w+") as out_conf_fh:
        yaml.dump(mqc_conf, out_conf_fh)

    modules = "-m interop "
    FORCE_BCL2FASTQ = os.environ.get("FORCE_BCL2FASTQ", None)
    modules += "-m bclconvert " if not FORCE_BCL2FASTQ else "-m bcl2fastq"
    # stats_dir = "Reports" if not FORCE_BCL2FASTQ else "Stats"

    multiqc_cmd = cfg.static.commands["multiqc_command"]
    multiqc_opts = cfg.static.commands["multiqc_options"]
    pname = pnames.replace(", ", "_")
    multiqc_out = cfg.output_path / "Stats" / f"sequencer_stats_{pname}.html"

    cmd = f"{multiqc_cmd} {multiqc_opts} --config {conf_name} {cfg.output_path}/Stats --filename {multiqc_out} {modules}"
    syslog.syslog(f"[multiqc_worker] Processing {cfg.output_path}\n")

    if os.environ.get("BFQ_TEST", None) and not FORCE_BCL2FASTQ:
        if not (cfg.output_path / "Stats" / "Demultiplex_Stats.csv").exists():
            print(
                "BFQ-TEST: Testflowcell was generated with bcl2fastq but environment is configured for bcl-convert. Using bcl2fastq paths and mqc modules."
            )
            cmd = cmd.replace("Reports", "Stats")
            cmd = cmd.replace("-m bclconvert", "-m bcl2fastq")
            print(cmd)

    subprocess.check_call(cmd, shell=True)

    os.chdir(oldWd)


def archive_worker(cfg):
    project_dirs = get_project_dirs(config)
    pnames = get_project_names(project_dirs)
    run_date = config.get("Options", "runID").split("_")[0]

    for p in pnames:
        """
        Archive fastq
        """
        if os.path.exists(
            os.path.join(
                config.get("Paths", "outputDir"), config.get("Options", "runID"), f"{p}.7za"
            )
        ):
            os.remove(
                os.path.join(
                    config.get("Paths", "outputDir"), config.get("Options", "runID"), f"{p}.7za"
                )
            )
        pw = None
        if config.get("Options", "SensitiveData") == "1":
            pw = (
                subprocess.check_output("xkcdpass -n 5 -d '-' -v '[a-z]'", shell=True)
                .decode()
                .strip("\n")
            )
            with open(
                os.path.join(
                    config.get("Paths", "outputDir"),
                    config.get("Options", "runID"),
                    f"encryption.{p}",
                ),
                "w",
            ) as pwfile:
                pwfile.write(f"{pw}\n")
        opts = f"-p{pw}" if pw else ""
        flowdir = os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"))
        report_dir = (
            os.path.join(flowdir, "Reports")
            if os.path.exists(os.path.join(flowdir, "Reports"))
            else ""
        )

        cmd = "7za a {opts} {flowdir}/{pnr}_{date}.7za {flowdir}/{pnr}/ {flowdir}/Stats {report_dir} {flowdir}/Undetermined*.fastq.gz {flowdir}/{pnr}_samplesheet.tsv {flowdir}/SampleSheet.csv {flowdir}/Sample-Submission-Form.xlsx {flowdir}/md5sum_{pnr}_fastq.txt ".format(
            opts=opts,
            date=run_date,
            flowdir=os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID")),
            report_dir=report_dir,
            pnr=p,
        )
        if "10X Genomics" in config.get("Options", "Libprep"):
            cmd += " {}".format(
                os.path.join(
                    config.get("Paths", "outputDir"),
                    config.get("Options", "runID"),
                    config.get("Options", "runID").split("_")[-1][1:],
                )
            )
        syslog.syslog(
            "[archive_worker] Zipping {}\n".format(
                os.path.join(
                    config.get("Paths", "outputDir"),
                    config.get("Options", "runID"),
                    f"{p}_{run_date}.7za",
                )
            )
        )
        subprocess.check_call(cmd, shell=True)
        """
        Archive pipeline output
        """
        if os.path.exists(
            os.path.join(
                config.get("Paths", "outputDir"), config.get("Options", "runID"), f"QC_{p}.7za"
            )
        ):
            os.remove(
                os.path.join(
                    config.get("Paths", "outputDir"), config.get("Options", "runID"), f"QC_{p}.7za"
                )
            )
        pw = None
        if config.get("Options", "SensitiveData") == "1":
            pw = (
                subprocess.check_output("xkcdpass -n 5 -d '-' -v '[a-z]'", shell=True)
                .decode()
                .strip("\n")
            )
            with open(
                os.path.join(
                    config.get("Paths", "outputDir"),
                    config.get("Options", "runID"),
                    f"encryption.QC_{p}",
                ),
                "w",
            ) as pwfile:
                pwfile.write(f"{pw}\n")

        opts = f"-p{pw}" if pw else ""
        flowdir = os.path.join(config.get("Paths", "outputDir"), config.get("Options", "runID"))
        pipeline = config.get("Options", "pipeline")
        qc_dir = os.path.join(
            os.environ["TMPDIR"], f"{p}_{run_date}", "data", "tmp", pipeline, "bfq"
        )

        cmd = f"7z a -l {opts} {flowdir}/QC_{p}_{run_date}.7za {qc_dir} "

        """
        cmd = "7z a -l {opts} {flowdir}/QC_{pnr}_{date}.7za {qc_dir} ".format(
                opts = opts,
                flowdir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
                pnr = p,
                date = run_date,
                qc_dir = qc_dir
            )
        """
        subprocess.check_call(cmd, shell=True)


def get_project_names(dirs):
    gcf = set()
    for d in dirs:
        for catalog in d.split("/"):
            if catalog.startswith("GCF-"):
                gcf.add(catalog)
    return gcf


def get_project_dirs(cfg):
    project_dirs = glob.glob(f"{cfg.output_path}/*/*.fastq.gz")
    project_dirs.extend(glob.glob(f"{cfg.output_path}/*/*/*.fastq.gz"))
    return to_dirs(project_dirs)


def post_workflow(project_id, base_dir, pipeline):
    analysis_dir = os.path.join(
        os.environ["TMPDIR"], "{}_{}".format(project_id, os.path.basename(base_dir).split("_")[0])
    )
    os.makedirs(os.path.join(base_dir, f"QC_{project_id}", "bfq"), exist_ok=True)
    cmd = "rsync -rvLp {}/ {}".format(
        os.path.join(analysis_dir, "data", "tmp", pipeline, "bfq"),
        os.path.join(base_dir, f"QC_{project_id}", "bfq"),
    )
    subprocess.check_call(cmd, shell=True)

    # Copy sample info
    cmd = "cp {} {}".format(
        os.path.join(analysis_dir, "data", "tmp", "sample_info.tsv"),
        os.path.join(base_dir, f"{project_id}_samplesheet.tsv"),
    )
    subprocess.check_call(cmd, shell=True)

    return True


def full_align(cfg):
    old_wd = os.getcwd()

    os.chdir(os.environ["TMPDIR"])
    project_names = get_project_names(get_project_dirs(cfg))
    run_date = os.path.basename(cfg.output_path).split("_")[0]
    for p in project_names:
        analysis_dir = Path(os.environ["TMPDIR"]) / f"{p}_{run_date}"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        (analysis_dir / "src").mkdir(parents=True, exist_ok=True)
        (analysis_dir / "data").mkdir(parents=True, exist_ok=True)

        os.chdir(analysis_dir)

        # copy snakemake pipeline
        cmd = "rm -rf {dst} && cp -r {src} {dst}".format(
            src="/opt/gcf-workflows", dst=os.path.join(analysis_dir, "src", "gcf-workflows")
        )
        subprocess.check_call(cmd, shell=True)

        # create config.yaml
        cmd = "/opt/conda/bin/python /opt/conda/bin/configmaker.py {runfolder} -p {project} --libkit '{lib}' --machine '{machine}' {create_fastq}".format(
            runfolder=cfg.output_path,
            project=p,
            lib=cfg.run.libprep,
            machine=get_sequencer(cfg.run.run_id),
            create_fastq=" --skip-create-fastq-dir" if Path("data/raw/fastq").exists() else "",
        )
        subprocess.check_call(cmd, shell=True)

        # write template Snakefile for pipeline
        # with open(os.path.join(analysis_dir,"Snakefile"), "w") as sn:
        #    sn.write(SNAKEFILE_TEMPLATE.format(workflow=pipeline))

        # run snakemake pipeline
        cmd = "snakemake --use-singularity --singularity-prefix $SINGULARITY_CACHEDIR --cores 32 --verbose -p multiqc_report"
        subprocess.check_call(cmd, shell=True)

        # copy report
        shutil.copy2(
            Path("data") / "tmp" / cfg.run.pipeline / "bfq" / f"multiqc_{p}.html",
            cfg.output_path / f"multiqc_{p}_{run_date}.html",
        )

        # if additional html reports exists (single cell), copy
        extra_html = glob.glob(
            analysis_dir
            / "data"
            / "tmp"
            / cfg.run.pipeline
            / "bfq"
            / "summaries"
            / "all_samples*.html"
        )
        if extra_html:
            shutil.copy2(
                extra_html[0], cfg.output_path / f"all_samples_web_summary_{p}_{run_date}.html"
            )

        # Copy sample info
        shutil.copy2(
            Path("data") / "tmp" / "sample_info.tsv", cfg.output_path / f"{p}_samplesheet.tsv"
        )

        # copy mqc_config
        shutil.copy2(
            Path("data") / "tmp" / cfg.run.pipeline / "bfq" / ".multiqc_config.yaml",
            cfg.output_path / f".multiqc_config_{p}.yaml",
        )

    os.chdir(old_wd)
    open(cfg.ouput_path / "analysis.made", "w").close()
    return True


# All steps that should be run after `make` go here
def postMakeSteps(config):
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

    tot, used, free = shutil.disk_usage(cfg.static.paths.base_dir)
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
