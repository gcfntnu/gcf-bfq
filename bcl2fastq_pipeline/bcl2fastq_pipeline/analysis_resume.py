"""Validate and bind an explicit reuse of an existing analysis context.

Domain validation stays in gcf-tools. These checks establish BFQ ownership,
input lineage and immutable execution inputs between preview, queue and launch.
"""

import csv
import hashlib
import json
import re

from pathlib import Path

import yaml

from bcl2fastq_pipeline import analysis_snapshots, preflight
from bcl2fastq_pipeline.state import StateConflictError


def _fail(detail):
    raise StateConflictError(
        f"Cannot resume analysis: {detail}. Restore/correct the retained context and curated "
        "inputs, then preview again; use an ordinary analysis rerun to rebuild it."
    )


def _yaml(path):
    try:
        value = yaml.safe_load(path.read_text())
    except (OSError, ValueError, yaml.YAMLError) as error:
        _fail(f"cannot read {path}: {error}")
    if not isinstance(value, dict):
        _fail(f"expected a mapping in {path}")
    return value


def _digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _files(root, paths):
    """Hash effective source bytes, including local edits; do not hash scientific data."""
    result = {}
    for path in paths:
        if path.is_dir():
            children = sorted(path.rglob("*"))
        else:
            children = [path]
        for child in children:
            relative = child.relative_to(root)
            if ".git" in relative.parts or "__pycache__" in relative.parts:
                continue
            if child.is_symlink() and child.is_dir():
                _fail(f"directory symlink in execution inputs: {child}")
            if child.is_file():
                result[str(relative)] = _digest(child)
            elif child.is_symlink():
                _fail(f"broken execution input: {child}")
    return result


def _project(state, project, expected_samples, result):
    output = Path(state["output_path"]).resolve()
    work = analysis_snapshots.workdir_path(state["run_id"], project).absolute()
    if not work.is_dir() or work.is_symlink():
        _fail(f"missing or symlinked workdir {work}")
    marker = work / analysis_snapshots.MARKER
    try:
        identity = json.loads(marker.read_text())
    except (OSError, ValueError) as error:
        _fail(f"missing ownership marker in {work}: {error}")
    if (
        not isinstance(identity, dict)
        or identity.get("run_id") != state["run_id"]
        or identity.get("project") != project
        or not isinstance(identity.get("token"), str)
        or not identity["token"]
    ):
        _fail(f"workdir belongs to a different run/project: {work}")
    recorded = state["stages"]["analysis"]["metadata"].get("workdirs", {}).get(project)
    if recorded and any(recorded.get(key) != identity.get(key) for key in identity):
        _fail(f"workdir ownership changed: {work}")
    config = _yaml(work / "config.yaml")
    if config.get("project_id") not in ([project], project):
        _fail(f"project_id does not match {project} in {work}/config.yaml")
    workflow = config.get("workflow")
    if not isinstance(workflow, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", workflow):
        _fail(f"invalid workflow in {work}/config.yaml")
    samples = config.get("samples")
    if not isinstance(samples, dict) or set(samples) != set(expected_samples):
        _fail(f"sample identities differ from curated inputs for {project}")
    for sample, values in samples.items():
        if not isinstance(values, dict) or values.get("Sample_ID") != sample:
            _fail(f"sample {sample} has inconsistent embedded Sample_ID in {work}/config.yaml")
        # Configmaker repeats flowcell identity once per lane/FASTQ pair.
        for field, expected in (
            ("Flowcell_ID", state["run_id"].split("_")[-1]),
            ("Flowcell_Name", state["run_id"]),
        ):
            value = values.get(field)
            if not isinstance(value, str) or set(value.split(",")) != {expected}:
                _fail(f"sample {sample} has inconsistent {field} in {work}/config.yaml")
    if "wells" in config:
        # Consume the shared parser's effective Cell Multiplexing table, not a
        # second workbook parser or a special exception for invalid identifiers.
        demux = result.forms[0].get("demux", {}).get("data")
        expected = {}
        if demux is not None:
            for row in demux.to_dict(orient="records"):
                expected[str(row["Sample_ID"])] = str(row["Wells"]).replace(" ", "")
        wells = config["wells"]
        if not isinstance(wells, dict) or set(wells) != set(expected):
            _fail(f"wells IDs differ from curated Cell Multiplexing inputs for {project}")
        for key, value in wells.items():
            if (
                not isinstance(value, dict)
                or value.get("Sample_ID") != key
                or value.get("Wells") != expected[key]
            ):
                _fail(f"wells mapping differs from curated inputs for {project}/{key}")
    source = work / "src/gcf-workflows"
    for directory in (work / "src", source, work / "pep"):
        if directory.is_symlink():
            _fail(f"symlinked execution directory {directory}")
    required = [
        work / "Snakefile",
        source / "libprep.config",
        source / workflow / f"{workflow}.smk",
    ]
    pep = work / "pep/pep_config.yaml"
    pep_config = _yaml(pep)
    for key in ("sample_table", "subsample_table"):
        if key not in pep_config and key == "subsample_table":
            continue
        value = pep_config.get(key)
        if not isinstance(value, str):
            _fail(f"missing {key} in {pep}")
        table = (pep.parent / value).resolve()
        if not table.is_relative_to((work / "pep").resolve()):
            _fail(f"PEP table outside retained pep directory: {table}")
        required.append(table)
    for path in required:
        if not path.is_file():
            _fail(f"missing execution input {path}")
    with (pep.parent / pep_config["sample_table"]).open(newline="") as handle:
        names = [row.get("sample_name") for row in csv.DictReader(handle)]
    if len(names) != len(set(names)) or set(names) != set(expected_samples):
        _fail(f"PEP sample identities differ from curated inputs for {project}")
    fastq_dir = config.get("fastq_dir", "data/raw/fastq")
    if not isinstance(fastq_dir, str):
        _fail(f"invalid fastq_dir in {work}/config.yaml")
    raw = work / fastq_dir
    if not raw.resolve().is_relative_to(work.resolve()):
        _fail(f"FASTQ link directory is outside retained workdir: {raw}")
    if not raw.is_dir():
        _fail(f"missing FASTQ directory {raw}")
    links = {}
    for path in sorted(raw.rglob("*.fastq.gz")):
        target = path.resolve()
        if (
            not path.is_symlink()
            or not target.is_file()
            or not target.is_relative_to(output / project)
        ):
            _fail(f"FASTQ must link to this run/project: {path}")
        links[str(path.relative_to(work))] = {
            "target": str(target),
            "size": target.stat().st_size,
            "mtime_ns": target.stat().st_mtime_ns,
        }
    actual = {str(path.resolve()) for path in (output / project).rglob("*.fastq.gz")}
    if not links or {item["target"] for item in links.values()} != actual:
        _fail(f"retained FASTQ links do not cover project {project}")
    return {
        "identity": {**identity, "path": str(work)},
        "workflow": workflow,
        "files": _files(work, [work / "config.yaml", work / "Snakefile", source, work / "pep"]),
        "fastq_links": links,
    }


def inspect(state, *, selection=None, result=None):
    """Read-only validation used before queueing and again under the execution lease."""
    output = Path(state["output_path"])
    if selection is None:
        selection, result = preflight.require_valid_inputs(state["source_path"], output)
    if (
        selection.sample_sheet != output / "SampleSheet.csv"
        or selection.submission_form != output / "Sample-Submission-Form.xlsx"
    ):
        _fail("resume requires both canonical output-side inputs")
    planned = {entry["project_id"]: entry["sample_ids"] for entry in result.summary["projects"]}
    if not planned or (state.get("projects") and set(state["projects"]) != set(planned)):
        _fail("curated project set differs from recorded projects")
    projects = {}
    for project, sample_ids in sorted(planned.items()):
        if Path(project).name != project or not project.startswith("GCF-"):
            _fail(f"invalid project path {project!r}")
        projects[project] = _project(state, project, sample_ids, result)
    workflows = {entry["workflow"] for entry in projects.values()}
    if len(workflows) != 1:
        _fail("projects disagree on retained workflow")
    return {
        "mode": "resume",
        "workflow": workflows.pop(),
        "projects": projects,
        "inputs": {
            str(path): _digest(path) for path in (selection.sample_sheet, selection.submission_form)
        },
    }


def verify(expected, current):
    if expected != current:
        _fail("execution inputs changed since preview/queue; requeue to accept the current context")


def print_plan(context):
    print("Analysis mode: resume (preserve files; Snakemake decides scientific reruns)")
    print(
        "State: analysis queued; reporting/finalization pending; supersede downstream notifications"
    )
    print("Manager preview only; this does not run a Snakemake DAG dry run.")
    for project, entry in context["projects"].items():
        work = entry["identity"]["path"]
        print(f"  {project}: {work}; workflow={entry['workflow']}")
        print(
            f"    Reuse {work}/config.yaml, Snakefile, pep/, src/gcf-workflows/, data/, .snakemake/"
        )
    print("After analysis: refresh BFQ reports, then rebuild delivery archives/checksums.")
