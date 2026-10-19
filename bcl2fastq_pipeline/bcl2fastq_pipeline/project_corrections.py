"""Lease-protected staging and reversible publication of project corrections.

Canonical results stay in place while a correction runs. The publication journal
is state-backed; recovery rolls back an uncommitted installation or finishes an
installation whose processing completion was committed. Prior products are kept
in the correction's history directory, including after successful replacement.
"""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import uuid

from pathlib import Path

from bcl2fastq_pipeline import preflight
from bcl2fastq_pipeline.config import parse_custom_options
from bcl2fastq_pipeline.state import (
    PROVENANCE_DIR,
    StateConflictError,
    apply_cleanup,
    cleanup_plan,
    utcnow,
)

DIRECTORY = "project-corrections"


def project_names(projects):
    """Validate explicit identifiers before they become paths or selectors."""
    result = sorted(set(projects or []))
    if not result:
        raise StateConflictError("Select at least one project with --project")
    for project in result:
        if not isinstance(project, str) or not re.fullmatch(
            r"GCF-[A-Za-z0-9][A-Za-z0-9-]*", project
        ):
            raise StateConflictError(f"Invalid project selector: {project!r}")
    return result


def belongs_to_project(name, projects):
    """Match complete project identifiers in existing delivery product names."""
    return any(
        re.search(rf"(?<![A-Za-z0-9-]){re.escape(project)}(?![A-Za-z0-9-])", name)
        for project in projects
    )


def _pairs(selection):
    pairs = getattr(selection, "pairs", None)
    return pairs if pairs else (selection,)


def validate_scope(state, selection, projects, *, refresh=False):
    projects = project_names(projects)
    pair_projects = {getattr(pair, "project_id", None) for pair in _pairs(selection)} - {None}
    known = pair_projects | set(state.get("projects") or [])
    missing = set(projects) - known
    if missing:
        raise StateConflictError("Unknown selected projects: " + ", ".join(sorted(missing)))
    if refresh and not pair_projects and set(projects) != known:
        raise StateConflictError(
            "Refreshing part of a legacy shared input pair is ambiguous; use project-qualified "
            "input pairs or refresh all projects together"
        )
    return projects


def print_plan(state, selection, projects, from_stage, *, refresh=False):
    """Preview canonical replacements, without creating the transaction directory."""
    projects = validate_scope(state, selection, projects, refresh=refresh)
    print("Selected projects: " + ", ".join(projects))
    print(
        "Retained projects: "
        + (", ".join(sorted(set(state["projects"]) - set(projects))) or "(none)")
    )
    print("Input replacements after successful processing:")
    for pair in _pairs(selection):
        project = getattr(pair, "project_id", None)
        if project is None or project in projects:
            print(
                f"  {project or 'legacy shared pair'}: {pair.sample_sheet}, {pair.submission_form}"
            )
    print("Results are generated in separate staging; publication replaces only selected products.")
    print("Paths replaced or invalidated after successful processing:")
    paths = cleanup_plan(state["output_path"], from_stage, projects=projects)
    for path in paths:
        print(f"  {path}")
    if not paths:
        print("  (none)")
    return paths


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _record(state):
    identifier = (state.get("restart_request") or {}).get("project_correction_id")
    if not identifier:
        return None
    return next(
        (item for item in state.get("project_corrections", []) if item["id"] == identifier), None
    )


def current(state):
    """Return the current correction, if this attempt has an explicit scope."""
    return _record(state)


def prepare(store, run_id, selection, projects, from_stage, *, refresh=False):  # noqa: PLR0913
    """Freeze all effective inputs after validation, without changing canonical files.

    Caller has persisted ``preparing`` and holds the execution lease. Source
    hashes are checked around the copy, so an intervening instrument edit cannot
    silently become this attempt's execution input.
    """
    state = store.read(run_id)
    projects = validate_scope(state, selection, projects, refresh=refresh)
    identifier = uuid.uuid4().hex
    root = Path(state["output_path"]) / PROVENANCE_DIR / DIRECTORY / identifier
    if any(path.is_symlink() for path in (root.parent.parent, root.parent, root)):
        raise StateConflictError("Project correction provenance directories must not be symlinks")
    inputs = root / "inputs"
    inputs.mkdir(parents=True)
    before = [(_hash(pair.sample_sheet), _hash(pair.submission_form)) for pair in _pairs(selection)]
    copied = preflight.copy_run_inputs(selection, inputs)
    copied_pairs = _pairs(copied)
    after = [(_hash(pair.sample_sheet), _hash(pair.submission_form)) for pair in copied_pairs]
    if before != after or before != [
        (_hash(pair.sample_sheet), _hash(pair.submission_form)) for pair in _pairs(selection)
    ]:
        raise StateConflictError(
            "Inputs changed while preparing the correction; validate and retry"
        )
    result = preflight.validate_selection(copied)
    if not result.ok:
        raise preflight.PreflightValidationError(result)
    manifest = [
        {
            "project_id": getattr(pair, "project_id", None),
            "sample_sheet": pair.sample_sheet.name,
            "submission_form": pair.submission_form.name,
            "sample_sheet_sha256": hashes[0],
            "submission_form_sha256": hashes[1],
        }
        for pair, hashes in zip(copied_pairs, after, strict=True)
    ]
    geometry = []
    for name in ("RunInfo.xml", "RunParameters.xml", "runParameters.xml"):
        candidates = (Path(state["output_path"]) / name, Path(state["source_path"]) / name)
        source = next((path for path in candidates if path.is_file()), None)
        if source is not None:
            checksum = _hash(source)
            shutil.copy2(source, inputs / name)
            if _hash(inputs / name) != checksum or _hash(source) != checksum:
                raise StateConflictError(f"Run geometry changed while staging: {source}")
            geometry.append({"name": name, "sha256": checksum})
    record = {
        "id": identifier,
        "projects": projects,
        "from_stage": from_stage,
        "status": "prepared",
        "created_at": utcnow(),
        "root": str(root),
        "inputs": manifest,
        "geometry": geometry,
        "attempt": None,
        "publication": [],
    }

    def update(current_state):
        current_state.setdefault("project_corrections", []).append(record)
        current_state["restart_request"]["projects"] = projects
        current_state["restart_request"]["project_correction_id"] = identifier
        return current_state

    return store.mutate(run_id, update)


def _copy(source, target):
    if source.is_symlink():
        raise StateConflictError(f"Cannot stage a symlink as correction output: {source}")
    if source.is_dir():
        shutil.copytree(source, target, symlinks=True)
    else:
        shutil.copy2(source, target)


def _selection(directory):
    # Resolve only the frozen directory, never fall back to mutable instrument inputs.
    return preflight.select_run_inputs(directory, directory)


def prepare_execution(store, run_id, cfg):
    """Select immutable staged inputs and isolate every mutable correction output."""
    state = store.read(run_id)
    record = _record(state)
    if record is None:
        return None
    if record["status"] in {"publishing", "published"}:
        raise StateConflictError("Recover correction publication before executing another attempt")
    root = Path(record["root"])
    working = root / f"attempt-{state['attempt']}" / state["run_id"]
    if working.exists():
        raise StateConflictError(f"Correction attempt staging already exists: {working}")
    working.mkdir(parents=True)
    inputs = root / "inputs"
    for item in record["inputs"]:
        for role in ("sample_sheet", "submission_form"):
            source = inputs / item[role]
            if _hash(source) != item[f"{role}_sha256"]:
                raise StateConflictError(f"Frozen correction input was changed: {source}")
            if not (working / source.name).exists():
                _copy(source, working / source.name)
    for item in record.get("geometry", []):
        source = inputs / item["name"]
        if _hash(source) != item["sha256"]:
            raise StateConflictError(f"Frozen run geometry was changed: {source}")
        _copy(source, working / source.name)
    canonical = Path(state["output_path"])
    if record["from_stage"] != "demultiplexing":
        # Copying preserves FASTQ bytes/mtime and avoids writable aliases to the
        # successful generation. Hardlinks would allow an external tool to alter it.
        for source in canonical.iterdir():
            if belongs_to_project(source.name, record["projects"]) or source.name in {
                "Stats",
                "Reports",
                "InterOp",
                "RunInfo.xml",
                "RunParameters.xml",
                "runParameters.xml",
            }:
                if not (working / source.name).exists():
                    if source.name == "Stats" and source.is_symlink():
                        resolved = source.resolve(strict=True)
                        if (
                            not resolved.is_relative_to(canonical.resolve())
                            or not resolved.is_dir()
                        ):
                            raise StateConflictError(
                                "Correction Stats must resolve inside the run output"
                            )
                        shutil.copytree(resolved, working / source.name, symlinks=True)
                    else:
                        _copy(source, working / source.name)
        demultiplexing_provenance = canonical / PROVENANCE_DIR / "demultiplexing"
        if demultiplexing_provenance.exists():
            (working / PROVENANCE_DIR).mkdir(exist_ok=True)
            _copy(demultiplexing_provenance, working / PROVENANCE_DIR / "demultiplexing")
        for project in record["projects"]:
            retained = canonical / PROVENANCE_DIR / "project-demultiplexing" / project
            if retained.is_dir():
                destination = working / PROVENANCE_DIR / "project-demultiplexing" / project
                destination.parent.mkdir(parents=True, exist_ok=True)
                _copy(retained, destination)
        apply_cleanup(cleanup_plan(working, record["from_stage"], projects=record["projects"]))
    cfg.run.output_override = working
    cfg.run.selected_projects = set(record["projects"])
    cfg.run.input_selection = _selection(working)
    cfg.run.sample_sheet = cfg.run.input_selection.sample_sheet
    cfg.run.sample_submission_form = cfg.run.input_selection.submission_form
    custom, _section = parse_custom_options(cfg.run.sample_sheet)
    cfg.run.apply_custom(custom, cfg.run.sample_sheet, cfg.run.sample_submission_form)

    def update(current_state):
        current_record = _record(current_state)
        current_record.update(status="running", attempt=state["attempt"], working=str(working))
        return current_state

    store.mutate(run_id, update)
    return working


def _publication_plan(state, record):
    canonical = Path(state["output_path"])
    working = Path(record["working"])
    projects = record["projects"]
    names = {
        path.name for path in canonical.iterdir() if belongs_to_project(path.name, projects)
    } | {path.name for path in working.iterdir() if belongs_to_project(path.name, projects)}
    input_names = set()
    for item in record["inputs"]:
        if item["project_id"] is None or item["project_id"] in projects:
            input_names.update(item[role] for role in ("sample_sheet", "submission_form"))
    names.update(input_names)
    for project in projects:
        relative = Path(PROVENANCE_DIR) / "project-demultiplexing" / project
        if (working / relative).exists() or (canonical / relative).exists():
            names.add(str(relative))
    backup = Path(record["root"]) / "previous"
    operations = []
    for name in sorted(names):
        target = canonical / name
        source = working / name
        # Replacing an unchanged curated input would needlessly change timestamps.
        if (
            name in input_names
            and source.is_file()
            and target.is_file()
            and _hash(source) == _hash(target)
        ):
            continue
        operations.append(
            {
                "target": str(target),
                "source": str(source),
                "backup": str(backup / name),
                "had_original": target.exists() or target.is_symlink(),
                "has_replacement": source.exists(),
                "source_is_directory": source.is_dir() and not source.is_symlink(),
            }
        )
    return operations


def _sync_directory(directory):
    descriptor = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def publish(store, run_id):
    """Journal and install selected products; completion still decides commit."""
    state = store.read(run_id)
    record = _record(state)
    if record is None:
        return
    if record["status"] != "running":
        raise StateConflictError(f"Correction is {record['status']}, not ready for publication")
    operations = _publication_plan(state, record)
    backup = Path(record["root"]) / "previous"
    backup.mkdir()
    for operation in operations:
        for key in ("target", "source", "backup"):
            Path(operation[key]).parent.mkdir(parents=True, exist_ok=True)

    def update(current_state):
        _record(current_state).update(status="publishing", publication=operations)
        return current_state

    store.mutate(run_id, update)
    for operation in operations:
        target, source, previous = (Path(operation[key]) for key in ("target", "source", "backup"))
        if operation["had_original"]:
            os.replace(target, previous)
        if operation["has_replacement"]:
            os.replace(source, target)
            # Retained analysis workdirs and durable notification attachments
            # refer to staged paths. Keep them resolvable after installation.
            source.symlink_to(target, target_is_directory=operation["source_is_directory"])
        for directory in (target.parent, source.parent, previous.parent):
            _sync_directory(directory)
    # Kept separate from 'committed': a crash before complete_run must roll back.
    store.mutate(run_id, lambda current_state: _set_status(current_state, "published"))


def _set_status(state, status):
    _record(state)["status"] = status
    return state


def recover(store, run_id):
    """Repeatable publication recovery. Call under the run execution lease."""
    state = store.read(run_id)
    record = _record(state)
    if record is None or record["status"] not in {"publishing", "published"}:
        return state
    if state["status"] == "completed" and state["attempt"] == record["attempt"]:
        return store.mutate(run_id, lambda current_state: _set_status(current_state, "committed"))
    for operation in reversed(record["publication"]):
        target, source, previous = (Path(operation[key]) for key in ("target", "source", "backup"))
        if previous.exists() or previous.is_symlink():
            if target.exists() or target.is_symlink():
                if source.is_symlink():
                    source.unlink()
                os.replace(target, source)
            os.replace(previous, target)
        elif (
            not operation["had_original"]
            and (not source.exists() or source.is_symlink())
            and target.exists()
        ):
            if source.is_symlink():
                source.unlink()
            os.replace(target, source)
        for directory in (target.parent, source.parent, previous.parent):
            _sync_directory(directory)
    return store.mutate(run_id, lambda current_state: _set_status(current_state, "rolled_back"))


def commit_publication(store, run_id):
    return recover(store, run_id)
