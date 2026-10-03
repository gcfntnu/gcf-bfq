#!/usr/bin/env python

import argparse
import datetime
import json
import os
import socket

from pathlib import Path

import pandas as pd

from bcl2fastq_pipeline.config import PipelineConfig
from bcl2fastq_pipeline.entrypoint import _remove_executable_directory_from_import_path
from bcl2fastq_pipeline.state import (
    FlowcellStateStore,
    StateConflictError,
    StateError,
    apply_cleanup,
    cleanup_plan,
    locate_source_run,
    new_state,
    output_entries,
    resolve_output_path,
    validate_restart_boundary,
    validate_restored_fastqs,
)
from bcl2fastq_pipeline.version import add_version_argument

from bcl2fastq_pipeline import (
    analysis_resume,
    analysis_snapshots,
    fastq_cleanup,
    index_corrections,
    notification_delivery,
    preflight,
    sequencing_delivery,
)

pd.set_option("display.max_rows", 5000)
pd.set_option("display.max_columns", 12)

INVENTORY_COLUMNS = ["project", "flowcell_path", "timestamp", "archived"]


def get_cfg():
    """Return the active PipelineConfig, loading a static one if needed."""
    try:
        return PipelineConfig.get()
    except RuntimeError:
        return PipelineConfig.load("/config/bcl2fastq.ini")


def _inventory_path(cfg):
    return cfg.static.paths.manager_dir / "flowcells.processed"


def _empty_inventory():
    return pd.DataFrame(columns=INVENTORY_COLUMNS)


def _read_inventory(cfg):
    path = _inventory_path(cfg)
    if not path.exists() or path.stat().st_size == 0:
        return _empty_inventory()
    inventory = pd.read_csv(path, dtype=str).fillna("0")
    for column in INVENTORY_COLUMNS:
        if column not in inventory:
            inventory[column] = "0"
    return inventory[INVENTORY_COLUMNS]


def _write_inventory(cfg, inventory):
    path = _inventory_path(cfg)
    path.parent.mkdir(parents=True, exist_ok=True)
    inventory[INVENTORY_COLUMNS].to_csv(path, index=False)


def add_flowcell(**args):
    """Upsert one project/run inventory row without creating duplicates."""
    cfg = get_cfg()
    project = str(args["project"])
    flowcell_path = str(args["path"])
    timestamp = args.get("timestamp") or datetime.datetime.now().isoformat()
    inventory = _read_inventory(cfg)
    match = (inventory["project"] == project) & (inventory["flowcell_path"] == flowcell_path)
    row = {
        "project": project,
        "flowcell_path": flowcell_path,
        "timestamp": str(timestamp),
        "archived": "0",
    }
    if match.any():
        for column, value in row.items():
            inventory.loc[match, column] = value
        # Collapse historical duplicates while preserving the refreshed row.
        inventory = inventory.loc[~match].copy()
        inventory = pd.concat([inventory, pd.DataFrame([row])], ignore_index=True)
    else:
        inventory = pd.concat([inventory, pd.DataFrame([row])], ignore_index=True)
    _write_inventory(cfg, inventory)
    return inventory


def list_processed(**args):
    inventory = _read_inventory(get_cfg())
    return inventory.loc[(inventory["timestamp"] != "0") | (inventory["archived"] != "0")]


def list_all(**args):
    return _read_inventory(get_cfg())


def list_project(project):
    inventory = _read_inventory(get_cfg())
    return inventory.loc[(inventory["project"] == project) & (inventory["timestamp"] != "0")]


def list_flowcell(flowcell):
    inventory = _read_inventory(get_cfg())
    return inventory.loc[
        (inventory["flowcell_path"] == str(flowcell)) & (inventory["timestamp"] != "0")
    ]


def list_flowcell_all(flowcell):
    """Legacy inventory lookup used by BFQ discovery precedence."""
    inventory = _read_inventory(get_cfg())
    return inventory.loc[inventory["flowcell_path"] == str(flowcell)]


def _legacy_row_for_run(cfg, run_id):
    inventory = _read_inventory(cfg)
    by_id = inventory["flowcell_path"].map(lambda value: Path(value).name == run_id)
    return inventory.loc[by_id]


def _run_id(value):
    return Path(str(value)).name


def _resolve_state_or_legacy(value, cfg, store):
    run_id = _run_id(value)
    if store.exists(run_id):
        return run_id, store.read(run_id), _empty_inventory()
    legacy = _legacy_row_for_run(cfg, run_id)
    return run_id, None, legacy


def _legacy_output_path(legacy, run_id, cfg):
    # Validate every historical location, rather than silently choosing the first.
    for path in legacy["flowcell_path"].unique():
        resolve_output_path(path, run_id, cfg)
    return cfg.static.paths.output_dir / run_id


def _confirm(prompt, force):
    if force:
        return True
    return input(f"{prompt} (yes/no): ").strip().lower() == "yes"


def _ensure_not_active(store, run_id, state, force, *, recover=True):
    if state is None:
        return state
    active = store.execution_active(run_id)
    if recover and state["status"] == "running" and not active:
        state = store.recover_interrupted(run_id)
    if active:
        if not force:
            raise StateConflictError(
                f"{run_id} is currently running; repeat with --force only for a deliberate override"
            )

    return state


def _print_plan(  # noqa: PLR0913
    action, run_id, from_stage, paths, refresh_inputs=False, *, output_path=None
):
    print(f"{action}: {run_id}")
    if from_stage:
        print(f"Restart boundary: {from_stage}")
    if output_path is not None:
        print(f"Output directory: {output_path}")
    if refresh_inputs:
        print("Inputs: refresh SampleSheet.csv and Sample-Submission-Form.xlsx from instrument")
    print("Paths to invalidate:")
    if paths:
        for path in paths:
            print(f"  {path}")
    else:
        print("  (none)")


def _cleanup_for_state(state, from_stage):
    output = Path(state["output_path"])
    paths = cleanup_plan(output, from_stage)
    if from_stage == "analysis":
        projects = set(state.get("projects") or [])
        if not projects and output.exists():
            projects = {
                child.name
                for child in output.iterdir()
                if child.is_dir() and child.name.startswith("GCF-")
            }
        work_root = Path(os.environ.get("TMPDIR", "/bfq-tmp"))
        run_date = state["run_id"].split("_", 1)[0]
        for project in projects:
            work = work_root / f"{project}_{run_date}"
            if work.exists():
                paths.append(work)
    return sorted(set(paths), key=str)


def _requested_indexes(args, from_stage):
    indexes = []
    if args.get("reverse_complement_index1", False):
        indexes.append("index1")
    if args.get("reverse_complement_index2", False) or args.get("tom_mode", False):
        indexes.append("index2")
    if indexes and from_stage != "demultiplexing":
        raise StateConflictError(
            "Index reverse-complement options require --from demultiplexing; "
            "no state or inputs were changed"
        )
    return tuple(indexes)


def _require_separate_source(source_path, output_path):
    source = Path(source_path).resolve()
    output = Path(output_path).resolve()
    if source == output or source in output.parents or output in source.parents:
        raise StateConflictError(
            f"Input directory {source_path} overlaps output directory {output_path}; "
            "instrument inputs must remain separate from output"
        )
    if Path(source_path).exists() and Path(output_path).exists():
        if Path(source_path).samefile(output_path):
            raise StateConflictError("Input and output directories refer to the same location")


def _initialize_source(value, cfg):
    raw = str(value)
    explicit = Path(raw).is_absolute() or "/" in raw or raw in {".", ".."}
    if explicit:
        source = Path(os.path.abspath(Path(raw).expanduser()))
        if not source.is_dir():
            raise StateConflictError(f"Input directory is unavailable or not a directory: {source}")
        return source
    matches = []
    for root in (cfg.static.paths.nova_base_dir, cfg.static.paths.ekista_base_dir):
        candidate = Path(root) / raw
        if candidate.is_dir() and not any(candidate.samefile(path) for path in matches):
            matches.append(candidate)
    if not matches:
        raise StateConflictError(
            f"Cannot initialize {raw}: instrument source directory is unavailable"
        )
    if len(matches) > 1:
        raise StateConflictError(
            f"Ambiguous input run {raw}: "
            + ", ".join(map(str, matches))
            + "; supply the exact input directory"
        )
    return matches[0].absolute()


def rerun_flowcell(**args):
    from_stage = args.get("from_stage") or "demultiplexing"
    indexes = _requested_indexes(args, from_stage)
    resume = args.get("resume", False)
    if resume and (from_stage != "analysis" or args.get("refresh_inputs", False)):
        raise StateConflictError("--resume requires --from analysis and forbids --refresh-inputs")
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id, state, legacy = _resolve_state_or_legacy(args["flowcell"], cfg, store)
    force = args.get("force", False)
    dry_run = args.get("dry_run", False)
    refresh_inputs = args.get("refresh_inputs", False)
    reason = args.get("reason")

    if state is None:
        if resume:
            raise StateConflictError("--resume requires existing BFQ state and owned workdirs")
        if legacy.empty:
            raise StateConflictError(f"No state or legacy inventory entry exists for {run_id}")
        output_path = _legacy_output_path(legacy, run_id, cfg)
        source_path = locate_source_run(run_id, cfg)
        if source_path is None:
            raise StateConflictError(
                f"Cannot queue legacy rerun for {run_id}: instrument source directory is unavailable"
            )
        state = new_state(
            run_id,
            source_path,
            output_path,
            origin="legacy_rerun",
            start_stage=from_stage,
            preparing=True,
            cfg=cfg,
        )
        state["restart_request"] = {
            "from": from_stage,
            "reason": reason,
            "hostname": socket.gethostname(),
            "requested_at": datetime.datetime.now(datetime.UTC).isoformat(),
            "refresh_inputs": bool(refresh_inputs),
        }
        creating = True
    else:
        state = _ensure_not_active(store, run_id, state, force, recover=False)
        creating = False

    validate_restart_boundary(state, from_stage)
    output_path = resolve_output_path(state["output_path"], run_id, cfg)
    state = {**state, "output_path": str(output_path)}
    _require_separate_source(state["source_path"], output_path)
    fastq_cleanup.require_restored(state, output_path, from_stage)
    context = analysis_resume.inspect(state) if resume else None
    paths = [] if resume else _cleanup_for_state(state, from_stage)
    if context:
        analysis_resume.print_plan(context)
    _print_plan("Rerun", run_id, from_stage, paths, refresh_inputs, output_path=output_path)
    correction = None
    selected = None
    if not resume and (from_stage in {"demultiplexing", "analysis"} or refresh_inputs):
        selected, _result = preflight.require_valid_inputs(
            state["source_path"], output_path, refresh=refresh_inputs
        )
        if indexes:
            correction = index_corrections.plan(
                store, run_id, selected, state["source_path"], indexes, reason=reason
            )
            index_corrections.print_plan(correction)
    index_corrections.print_orientation_preview(state, selected, correction, refresh=refresh_inputs)
    if dry_run:
        return state
    if not _confirm(f"Queue rerun for {run_id}", force):
        print("Skipped; SampleSheet unchanged.")
        return state

    with store.execution_lease(run_id):
        if not creating and store.read(run_id)["updated_at"] != state["updated_at"]:
            raise StateConflictError("Run changed while preparing the command; inspect and retry")
        if not resume and _cleanup_for_state(state, from_stage) != paths:
            raise StateConflictError("Output paths changed after the preview; inspect and retry")
        _require_separate_source(state["source_path"], output_path)
        fastq_cleanup.require_restored(state, output_path, from_stage)
        selected = None
        if resume:
            analysis_resume.verify(context, analysis_resume.inspect(state))
        if not resume and (from_stage in {"demultiplexing", "analysis"} or refresh_inputs):
            selected, _result = preflight.require_valid_inputs(
                state["source_path"], output_path, refresh=refresh_inputs
            )
        if correction is not None:
            current = index_corrections.plan(
                store, run_id, selected, state["source_path"], indexes, reason=reason
            )
            index_corrections.verify_plan(correction, current)
        if not creating:
            index_corrections.recover(store, run_id)
            analysis_snapshots.recover(store, run_id)
            if state["status"] == "running":
                store.fail_stage(
                    run_id,
                    state["current_stage"],
                    summary="Interrupted before operator rerun; execution lease is inactive",
                    report_path=None,
                    interrupted=True,
                )
        if creating:
            store.create(state)
        else:
            store.set_preparing(
                run_id,
                from_stage,
                reason=reason,
                refresh_inputs=refresh_inputs,
                output_path=output_path,
                analysis_resume=context,
            )

        if selected is not None:
            # Preserve selected legacy filenames before demultiplexing cleanup can
            # remove them, and make copy failures occur before results are removed.
            effective = preflight.copy_run_inputs(selected, output_path)
            result = preflight.validate_selection(effective)
            if not result.ok:
                raise preflight.PreflightValidationError(result)
        if correction is not None:
            index_corrections.apply(store, run_id, correction)
        apply_cleanup(paths)
        return store.queue(run_id, from_stage)


def initialize_flowcell(**args):
    from_stage = args.get("from_stage") or "demultiplexing"
    indexes = _requested_indexes(args, from_stage)
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    source_path = _initialize_source(args["flowcell"], cfg)
    run_id = (
        source_path.resolve().name if str(args["flowcell"]) in {".", ".."} else source_path.name
    )
    force = args.get("force", False)
    dry_run = args.get("dry_run", False)
    refresh_inputs = args.get("refresh_inputs", False)
    reason = args.get("reason")

    if store.exists(run_id):
        raise StateConflictError(f"State already exists for {run_id}; use rerun instead")
    if not _legacy_row_for_run(cfg, run_id).empty:
        raise StateConflictError(
            f"{run_id} is protected by the legacy inventory; use rerun instead"
        )

    output_path = resolve_output_path(run_id, run_id, cfg)
    _require_separate_source(source_path, output_path)
    if from_stage == "analysis":
        recognized, detail = validate_restored_fastqs(output_path)
        if not recognized:
            raise StateConflictError(f"Cannot initialize from analysis: {detail}")

    state = new_state(
        run_id,
        source_path,
        output_path,
        origin="restored_legacy_fastq" if from_stage == "analysis" else "new",
        start_stage=from_stage,
        preparing=True,
        cfg=cfg,
    )
    state["restart_request"] = {
        "from": from_stage,
        "reason": reason,
        "hostname": socket.gethostname(),
        "requested_at": datetime.datetime.now(datetime.UTC).isoformat(),
        "refresh_inputs": bool(refresh_inputs),
    }
    paths = _cleanup_for_state(state, from_stage)
    _print_plan("Initialize", run_id, from_stage, paths, refresh_inputs, output_path=output_path)
    selected, _result = preflight.require_valid_inputs(
        source_path, output_path, refresh=refresh_inputs
    )
    correction = None
    if indexes:
        correction = index_corrections.plan(
            store, run_id, selected, source_path, indexes, reason=reason
        )
        index_corrections.print_plan(correction)
    index_corrections.print_orientation_preview(state, selected, correction, refresh=refresh_inputs)
    if dry_run:
        return state
    if not _confirm(f"Initialize {run_id}", force):
        print("Skipped; SampleSheet unchanged.")
        return state

    with store.execution_lease(run_id):
        if store.exists(run_id):
            raise StateConflictError(f"State was created for {run_id}; inspect and retry")
        _require_separate_source(source_path, output_path)
        if not source_path.is_dir():
            raise StateConflictError(f"Input directory is no longer available: {source_path}")
        if _cleanup_for_state(state, from_stage) != paths:
            raise StateConflictError("Output paths changed after the preview; inspect and retry")
        selected, _result = preflight.require_valid_inputs(
            source_path, output_path, refresh=refresh_inputs
        )
        if correction is not None:
            current = index_corrections.plan(
                store, run_id, selected, source_path, indexes, reason=reason
            )
            index_corrections.verify_plan(correction, current)
        # Persist preparing before copies or destructive work; a partial operation
        # must never become eligible for the daemon to execute automatically.
        store.create(state)
        effective = preflight.copy_run_inputs(selected, output_path)
        result = preflight.validate_selection(effective)
        if not result.ok:
            raise preflight.PreflightValidationError(result)
        if correction is not None:
            index_corrections.apply(store, run_id, correction)
        apply_cleanup(paths)
        return store.queue(run_id, from_stage)


def validate_flowcell(**args):
    """Resolve BFQ's effective input pair without creating state or output paths."""
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id, state, legacy = _resolve_state_or_legacy(args["flowcell"], cfg, store)
    if state is not None:
        output_path = resolve_output_path(state["output_path"], run_id, cfg)
        source_path = Path(state["source_path"])
    else:
        output_path = (
            _legacy_output_path(legacy, run_id, cfg)
            if not legacy.empty
            else resolve_output_path(run_id, run_id, cfg)
        )
        source_path = locate_source_run(run_id, cfg) or cfg.static.paths.nova_base_dir / run_id
    print(f"Input validation: {run_id}")
    _selection, result = preflight.require_valid_inputs(
        source_path, output_path, refresh=args.get("refresh_inputs", False)
    )
    return result


def clean_fastqs(**args):
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id = _run_id(args["flowcell"])
    if not store.exists(run_id):
        raise StateConflictError(
            "clean-fastqs requires state-backed completed finalization; legacy inventory alone "
            "does not establish eligibility"
        )
    state = store.read(run_id)
    output = resolve_output_path(state["output_path"], run_id, cfg)
    resolve_output_path(args["flowcell"], run_id, cfg)
    if store.execution_active(run_id):
        raise StateConflictError(f"{run_id} has an active execution lease")
    files = fastq_cleanup.plan(output)
    retained = fastq_cleanup.required_archives(state, output, files)
    print(f"Clean FASTQs: {run_id}\nOutput directory: {output}")
    print("FASTQs selected for deletion:")
    for name, info in files.items():
        print(f"  {output / name}" + (" (symlink only)" if info["symlink"] else ""))
    if not files:
        print("  (none)")
    print("Required delivery archives/checksums retained:")
    for path in retained:
        print(f"  {path}")
    size = sum(info["size"] for info in files.values() if not info["symlink"])
    print(
        f"Estimated recoverable space: {size:,} bytes ({size / 1024**3:.2f} GiB; excludes symlink targets)"
    )
    if args.get("dry_run", False):
        return state
    if not _confirm(f"Delete the selected FASTQs for {run_id}", args.get("force", False)):
        print("Skipping...")
        return state
    with store.execution_lease(run_id):
        current = store.read(run_id)
        if current["updated_at"] != state["updated_at"]:
            raise StateConflictError("Run changed while preparing the command; inspect and retry")
        resolve_output_path(current["output_path"], run_id, cfg)
        resolve_output_path(args["flowcell"], run_id, cfg)
        if fastq_cleanup.plan(output) != files:
            raise StateConflictError("FASTQs changed after the preview; inspect and retry")
        fastq_cleanup.required_archives(current, output, files)
        result = fastq_cleanup.execute(store, run_id, output, files)
    print(
        f"Removed {len(result['fastq_cleanup']['removed_files'])} FASTQ files/links; processing remains completed."
    )
    return result


def _archive_targets(flowcell, projects):
    targets = [flowcell / project for project in projects]
    targets += list(flowcell.rglob("*.bam"))
    targets += list(flowcell.rglob("*.bam.bai"))
    targets += list(flowcell.rglob("*.fastq.gz"))
    targets += list(flowcell.glob("*.7za"))
    # Avoid deleting descendants twice when a project directory is already included.
    ordered = sorted(set(targets), key=lambda path: (len(path.parts), str(path)))
    result = []
    for path in ordered:
        if flowcell / analysis_snapshots.DIRECTORY in path.parents:
            continue
        if any(parent in path.parents or parent == path for parent in result):
            continue
        if path.exists():
            result.append(path)
    return result


def archive_flowcell(**args):
    force = args.get("force", False)
    dry_run = args.get("dry_run", False)
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id, state, legacy = _resolve_state_or_legacy(args["flowcell"], cfg, store)

    if state is not None:
        state = _ensure_not_active(store, run_id, state, force)
        flowcell = resolve_output_path(state["output_path"], run_id, cfg)
        output_entries(flowcell)
        projects = state.get("projects") or [
            child.name for child in flowcell.glob("GCF-*") if child.is_dir()
        ]
    else:
        if legacy.empty:
            raise StateConflictError(f"No such flowcell: {run_id}")
        flowcell = _legacy_output_path(legacy, run_id, cfg)
        output_entries(flowcell)
        projects = sorted(set(legacy["project"]))

    inventory = _read_inventory(cfg)
    matching_inventory_paths = {
        value
        for value in inventory["flowcell_path"].unique()
        if Path(value).name == run_id and resolve_output_path(value, run_id, cfg) == flowcell
    }
    targets = _archive_targets(flowcell, projects)
    _print_plan("Archive", run_id, None, targets, output_path=flowcell)
    if dry_run:
        return state
    if not _confirm(f"Archive {run_id}", force):
        print("Skipping...")
        return state

    with store.execution_lease(run_id):
        output_entries(flowcell)
        if _archive_targets(flowcell, projects) != targets:
            raise StateConflictError("Output paths changed after the preview; inspect and retry")
        if state is not None:
            if store.read(run_id)["updated_at"] != state["updated_at"]:
                raise StateConflictError(
                    "Run changed while preparing the command; inspect and retry"
                )
            analysis_snapshots.recover(store, run_id)
            store.invalidate_notifications(
                run_id, reason="Delivery outputs archived", output_path=flowcell
            )
        apply_cleanup(targets)
        inventory = _read_inventory(cfg)
        match = inventory["flowcell_path"].isin(matching_inventory_paths)
        if match.any():
            inventory.loc[match, "archived"] = datetime.datetime.now().isoformat()
            _write_inventory(cfg, inventory)
        if state is not None:
            return store.mark_archived(run_id)
        return None


def combined_list(status=None, stage=None, query=None):
    """List or search runs, with canonical state authoritative before filtering."""
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    rows = []
    state_run_ids = set()
    needle = query.casefold() if query is not None else None

    def matches(run_id, projects):
        return needle is None or any(needle in value.casefold() for value in (run_id, *projects))

    for state in store.list_states():
        # Excluded states must not reappear through stale compatibility inventory.
        state_run_ids.add(state["run_id"])
        if status and state["status"] != status:
            continue
        if stage and state["current_stage"] != stage:
            continue
        if not matches(state["run_id"], state["projects"]):
            continue
        output = state["output_path"]
        rows.append(
            {
                "run_id": state["run_id"],
                "status": state["status"],
                "stage": state["current_stage"],
                "origin": state["origin"],
                "projects": ",".join(state["projects"]),
                "output_path": output,
                "archived": state["archive"]["archived_at"] or "0",
            }
        )

    inventory = _read_inventory(cfg)
    inventory = inventory.assign(run_id=inventory["flowcell_path"].map(_run_id))
    for run_id, group in inventory.groupby("run_id", sort=True):
        if run_id in state_run_ids:
            continue
        if stage and stage != "legacy":
            continue
        if status and status not in {"completed", "archived", "legacy"}:
            continue
        archived_values = [value for value in group["archived"] if value != "0"]
        legacy_status = "archived" if archived_values else "completed"
        if status and status not in ("legacy", legacy_status):
            continue
        projects = sorted(set(group["project"]))
        if not matches(run_id, projects):
            continue
        # Historical inventories can contain both a bare ID and an absolute path.
        # Choose a stable representative without requiring archived outputs to exist.
        flowcell_path = sorted(set(group["flowcell_path"]))[0]
        rows.append(
            {
                "run_id": run_id,
                "status": legacy_status,
                "stage": "legacy",
                "origin": "legacy_inventory",
                "projects": ",".join(projects),
                "output_path": flowcell_path,
                "archived": archived_values[-1] if archived_values else "0",
            }
        )

    frame = pd.DataFrame(
        rows,
        columns=["run_id", "status", "stage", "origin", "projects", "output_path", "archived"],
    )
    return frame.sort_values(["run_id", "origin"]) if not frame.empty else frame


def show_flowcell(**args):
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id = _run_id(args["flowcell"])
    if store.exists(run_id):
        state = store.recover_interrupted(run_id)
        try:
            output = resolve_output_path(state["output_path"], run_id, cfg)
            orientation = index_corrections.current_orientation(
                {**state, "output_path": str(output)}
            )
        except StateError as error:
            # Keep show useful for diagnosing a mismatched recorded output path.
            orientation = {"error": str(error)}
        state = {**state, "index_orientation": orientation}
        print(json.dumps(state, indent=2, sort_keys=True))
        return state
    legacy = _legacy_row_for_run(cfg, run_id)
    if legacy.empty:
        raise StateConflictError(f"No such flowcell: {run_id}")
    print(legacy.to_string(index=False))
    return legacy


def status_flowcell(**args):
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id = _run_id(args["flowcell"])
    if store.exists(run_id):
        state = store.recover_interrupted(run_id)
        print(f"{run_id}: {state['status']} ({state['current_stage']})")
        qc = state.get("sequencing_qc")
        if qc:
            print(f"  Sequencing QC (demultiplexing execution {qc['execution']}): {qc['status']}")
            if qc.get("last_error"):
                print(f"    {qc['last_error']}")
            if qc.get("status") == "completed":
                print(f"    {qc['result']['report_path']}")
        cleanup = state.get("fastq_cleanup")
        if cleanup:
            print(f"  FASTQ cleanup: {cleanup['status']} (started {cleanup['started_at']})")
        results = state["stages"]["finalization"]["metadata"].get("analysis_snapshots", {})
        retained = state.get("analysis_snapshots", {})
        for project in sorted(results.keys() | retained.keys()):
            result = results.get(project, {})
            record = retained.get(project, {})
            if result.get("status") == "unavailable":
                print(f"  Analysis snapshot {project}: unavailable ({result['reason']})")
            if record:
                label = "retained"
                if "pending" in record:
                    label = "publication pending; BFQ will retry"
                elif result.get("status") != "available":
                    label = "retained from previous successful finalization"
                print(f"  Analysis snapshot {project}: {label} ({record['archive']})")
        for entry in state.get("delivery_notifications", []):
            if entry["status"] == "superseded":
                continue
            label = entry["status"]
            if label == "sending" and not store.execution_active(run_id):
                label = "uncertain (interrupted delivery)"
            print(f"  Notification {entry['id']}: {label}")
            if entry.get("last_error"):
                print(f"    {entry['last_error']}")
        return state
    legacy = _legacy_row_for_run(cfg, run_id)
    if legacy.empty:
        raise StateConflictError(f"No such flowcell: {run_id}")
    status = "archived" if any(legacy["archived"] != "0") else "completed"
    print(f"{run_id}: {status} (legacy inventory)")
    return legacy


def retry_notifications(**args):
    """Retry saved notification intents without loading inputs or changing processing."""
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id = _run_id(args["flowcell"])
    kind = args.get("kind")
    with store.execution_lease(run_id):
        state = store.read(run_id)
        entries = [
            entry
            for entry in state.get("delivery_notifications", [])
            if (kind is None or entry["kind"] == kind) and entry["status"] != "superseded"
        ]
        if not entries:
            raise StateConflictError(
                "No retained notifications match this request. Legacy notifications are not "
                "reconstructed; superseded outputs cannot be notified."
            )
        if args.get("retry_uncertain"):
            print("Retrying uncertain delivery may duplicate email already accepted by SMTP.")
        successful = notification_delivery.deliver_pending(
            cfg,
            store,
            run_id,
            kind=kind,
            retry=True,
            retry_uncertain=args.get("retry_uncertain", False),
        )
    if not successful:
        raise StateConflictError(
            "Notification delivery remains incomplete. Inspect flowcell-manager status/show; "
            "correct configuration or use --retry-uncertain if duplicate delivery is acceptable."
        )
    print(f"{run_id}: matching notifications sent or already delivered")
    return store.read(run_id)


def retry_sequencing_qc(**args):
    """Recover the report from saved conversion inputs, without BCL or analysis."""
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id = _run_id(args["flowcell"])
    with store.execution_lease(run_id):
        before = store.read(run_id)
        qc = sequencing_delivery.report_and_notify(cfg, store, run_id, retry=True)
    if qc["status"] != "completed":
        raise StateConflictError(f"Sequencing QC remains unavailable: {qc.get('last_error')}")
    print(f"{run_id}: sequencing QC available at {qc['result']['report_path']}")
    print(
        "Existing sent notifications are preserved. Inspect status for delivery failures; "
        "retry with retry-notifications --kind sequencing."
    )
    if before["stages"]["finalization"]["status"] == "completed":
        print(
            "If this recovery changed the report, rerun --from finalization to include it "
            "in delivery archives. Existing archives were not changed."
        )
    return store.read(run_id)


def pretty_print(df):
    if df.empty:
        print("No matching flowcells.")
    else:
        print(df.to_string(index=False))


def _add_index_options(parser):
    parser.add_argument(
        "--reverse-complement-index1",
        action="store_true",
        help="Toggle index1 (index column) before a demultiplexing restart.",
    )
    parser.add_argument(
        "--reverse-complement-index2",
        action="store_true",
        help="Toggle index2 before a demultiplexing restart.",
    )
    parser.add_argument("--tom-mode", action="store_true", help=argparse.SUPPRESS)


def _search_query(value):
    if not value.strip():
        raise argparse.ArgumentTypeError("query must not be empty or whitespace-only")
    return value


def main():
    _remove_executable_directory_from_import_path()
    parser = argparse.ArgumentParser(
        description="Manage flowcells. fm and flowcell-manager provide the same commands.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_version_argument(parser)
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_add = subparsers.add_parser("add", help="Add a project to the compatibility inventory.")
    parser_add.set_defaults(func=add_flowcell)
    parser_add.add_argument("project", type=str)
    parser_add.add_argument("path", type=str)
    parser_add.add_argument("timestamp", nargs="?", type=str)

    parser_archive = subparsers.add_parser("archive", help="Archive flowcell delivery data.")
    parser_archive.set_defaults(func=archive_flowcell)
    parser_archive.add_argument("flowcell", help="Run ID or flowcell path.")
    parser_archive.add_argument("--dry-run", action="store_true")
    parser_archive.add_argument("--force", action="store_true")

    parser_rerun = subparsers.add_parser("rerun", help="Queue a state-backed rerun.")
    parser_rerun.set_defaults(func=rerun_flowcell)
    parser_rerun.add_argument("flowcell", help="Run ID or flowcell path.")
    parser_rerun.add_argument(
        "--from",
        dest="from_stage",
        choices=["demultiplexing", "analysis", "reporting", "finalization"],
        default="demultiplexing",
    )
    parser_rerun.add_argument("--dry-run", action="store_true")
    parser_rerun.add_argument("--force", action="store_true")
    parser_rerun.add_argument("--reason")
    parser_rerun.add_argument("--refresh-inputs", action="store_true")
    parser_rerun.add_argument(
        "--resume",
        action="store_true",
        help="Reuse owned analysis workdirs/configuration; requires --from analysis, forbids --refresh-inputs.",
    )
    _add_index_options(parser_rerun)

    parser_clean = subparsers.add_parser(
        "clean-fastqs", help="Remove delivered FASTQs while retaining archives and other products."
    )
    parser_clean.set_defaults(func=clean_fastqs)
    parser_clean.add_argument("flowcell", help="Run ID or flowcell output path.")
    parser_clean.add_argument("--dry-run", action="store_true")
    parser_clean.add_argument(
        "--force",
        action="store_true",
        help="Confirm deletion non-interactively; retain all eligibility checks.",
    )

    parser_initialize = subparsers.add_parser(
        "initialize", help="Prepare and queue a run from its input directory."
    )
    parser_initialize.set_defaults(func=initialize_flowcell)
    parser_initialize.add_argument(
        "flowcell", help="Input directory or run ID in configured instrument roots."
    )
    parser_initialize.add_argument(
        "--from",
        dest="from_stage",
        default="demultiplexing",
        choices=["demultiplexing", "analysis"],
    )
    parser_initialize.add_argument("--dry-run", action="store_true")
    parser_initialize.add_argument("--force", action="store_true")
    parser_initialize.add_argument("--reason")
    parser_initialize.add_argument("--refresh-inputs", action="store_true")
    _add_index_options(parser_initialize)

    parser_list = subparsers.add_parser("list", help="List state-backed and legacy flowcells.")
    parser_search = subparsers.add_parser(
        "search",
        help="Find flowcells by project name or run/flowcell ID.",
        description=(
            "Search project names and run/flowcell IDs using case-insensitive literal "
            "substrings (not regular expressions). Includes completed and archived runs. "
            "Exit 0 for a successful search, including no matches; 2 for invalid arguments; "
            "1 for state errors."
        ),
        epilog=(
            "Examples: fm search GCF-2026-043; fm search HL2T7AFXC; "
            "fm search GCF-2026 --status failed"
        ),
    )
    parser_search.add_argument("query", metavar="QUERY", type=_search_query)
    for list_parser in (parser_list, parser_search):
        list_parser.set_defaults(
            func=lambda **kwargs: combined_list(
                kwargs.get("status"), kwargs.get("stage"), kwargs.get("query")
            ),
            print_res=True,
        )
        list_parser.add_argument(
            "--status", help="Filter by exact status; legacy selects inventory-only runs."
        )
        list_parser.add_argument(
            "--stage", help="Filter by exact stage; legacy selects inventory-only runs."
        )

    parser_list_processed = subparsers.add_parser(
        "list-processed", help="List compatibility inventory rows."
    )
    parser_list_processed.set_defaults(func=list_processed, print_res=True)

    parser_show = subparsers.add_parser("show", help="Show the full state for one run.")
    parser_show.set_defaults(func=show_flowcell)
    parser_show.add_argument("flowcell")

    parser_status = subparsers.add_parser("status", help="Show concise status for one run.")
    parser_status.set_defaults(func=status_flowcell)
    parser_status.add_argument("flowcell")

    parser_validate = subparsers.add_parser(
        "validate", help="Check the effective input pair read-only, before any processing."
    )
    parser_validate.set_defaults(func=validate_flowcell)
    parser_validate.add_argument("flowcell", help="Run ID or flowcell path.")
    parser_validate.add_argument(
        "--refresh-inputs",
        action="store_true",
        help="Preview validation of both instrument inputs without copying them.",
    )

    parser_retry = subparsers.add_parser(
        "retry-notifications", help="Retry saved completion mail without rerunning processing."
    )
    parser_retry.add_argument("flowcell", help="Run ID or flowcell path.")
    parser_retry.add_argument("--kind", choices=["sequencing", "processed", "finalized"])
    parser_retry.add_argument(
        "--retry-uncertain",
        action="store_true",
        help="Allow resending after uncertain SMTP acceptance; duplicates are possible.",
    )
    parser_retry.set_defaults(func=retry_notifications)

    parser_qc = subparsers.add_parser(
        "retry-sequencing-qc", help="Recover early sequencing QC without rerunning conversion."
    )
    parser_qc.add_argument("flowcell", help="Run ID or flowcell path.")
    parser_qc.set_defaults(func=retry_sequencing_qc)

    args = parser.parse_args()
    values = vars(args)
    try:
        if values.pop("print_res", False):
            pretty_print(values.pop("func")(**values))
        else:
            values.pop("func")(**values)
    except StateError as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
