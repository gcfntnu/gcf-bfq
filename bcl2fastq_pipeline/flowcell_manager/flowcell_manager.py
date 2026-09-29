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
    refresh_run_inputs,
    validate_restart_boundary,
    validate_restored_fastqs,
)

from bcl2fastq_pipeline import notification_delivery

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


def _confirm(prompt, force):
    if force:
        return True
    return input(f"{prompt} (yes/no): ").strip().lower() == "yes"


def _ensure_not_active(store, run_id, state, force):
    if state is None:
        return state
    if state["status"] == "running" and not store.execution_active(run_id):
        state = store.recover_interrupted(run_id)
    if store.execution_active(run_id) or state["status"] == "running":
        if not force:
            raise StateConflictError(
                f"{run_id} is currently running; repeat with --force only for a deliberate override"
            )

    return state


def _print_plan(action, run_id, from_stage, paths, refresh_inputs=False):
    print(f"{action}: {run_id}")
    if from_stage:
        print(f"Restart boundary: {from_stage}")
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


def rerun_flowcell(**args):
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id, state, legacy = _resolve_state_or_legacy(args["flowcell"], cfg, store)
    from_stage = args.get("from_stage") or "demultiplexing"
    force = args.get("force", False)
    dry_run = args.get("dry_run", False)
    refresh_inputs = args.get("refresh_inputs", False)
    reason = args.get("reason")

    if state is None:
        if legacy.empty:
            raise StateConflictError(f"No state or legacy inventory entry exists for {run_id}")
        output_path = Path(legacy.iloc[0]["flowcell_path"])
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
        state = _ensure_not_active(store, run_id, state, force)
        creating = False

    validate_restart_boundary(state, from_stage)
    paths = _cleanup_for_state(state, from_stage)
    _print_plan("Rerun", run_id, from_stage, paths, refresh_inputs)
    if dry_run:
        return state
    if not _confirm(f"Queue rerun for {run_id}", force):
        print("Skipping...")
        return state

    with store.execution_lease(run_id):
        if not creating and store.read(run_id)["updated_at"] != state["updated_at"]:
            raise StateConflictError("Run changed while preparing the command; inspect and retry")
        if creating:
            store.create(state)
        else:
            store.set_preparing(
                run_id,
                from_stage,
                reason=reason,
                refresh_inputs=refresh_inputs,
            )

        apply_cleanup(paths)
        if refresh_inputs:
            refresh_run_inputs(state["source_path"], state["output_path"])
        return store.queue(run_id, from_stage)


def initialize_flowcell(**args):
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    run_id = _run_id(args["flowcell"])
    from_stage = args["from_stage"]
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

    output_path = cfg.static.paths.output_dir / run_id
    source_path = locate_source_run(run_id, cfg)
    if source_path is None:
        raise StateConflictError(
            f"Cannot initialize {run_id}: instrument source directory is unavailable"
        )
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
    _print_plan("Initialize", run_id, from_stage, paths, refresh_inputs)
    if dry_run:
        return state
    if not _confirm(f"Initialize {run_id}", force):
        print("Skipping...")
        return state

    # Persist "preparing" before destructive work. If cleanup fails, the daemon
    # will not run this flowcell until the operator explicitly resolves it.
    store.create(state)
    apply_cleanup(paths)
    if refresh_inputs:
        refresh_run_inputs(source_path, output_path)
    return store.queue(run_id, from_stage)


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
        flowcell = Path(state["output_path"])
        projects = state.get("projects") or [
            child.name for child in flowcell.glob("GCF-*") if child.is_dir()
        ]
    else:
        if legacy.empty:
            raise StateConflictError(f"No such flowcell: {run_id}")
        flowcell = Path(legacy.iloc[0]["flowcell_path"])
        projects = sorted(set(legacy["project"]))

    targets = _archive_targets(flowcell, projects)
    _print_plan("Archive", run_id, None, targets)
    if dry_run:
        return state
    if not _confirm(f"Archive {run_id}", force):
        print("Skipping...")
        return state

    with store.execution_lease(run_id):
        if state is not None:
            if store.read(run_id)["updated_at"] != state["updated_at"]:
                raise StateConflictError(
                    "Run changed while preparing the command; inspect and retry"
                )
            store.invalidate_notifications(run_id, reason="Delivery outputs archived")
        apply_cleanup(targets)
        inventory = _read_inventory(cfg)
        match = inventory["flowcell_path"] == str(flowcell)
        if match.any():
            inventory.loc[match, "archived"] = datetime.datetime.now().isoformat()
            _write_inventory(cfg, inventory)
        if state is not None:
            return store.mark_archived(run_id)
        return None


def combined_list(status=None, stage=None):
    cfg = get_cfg()
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    rows = []
    state_paths = set()

    for state in store.list_states():
        if status and state["status"] != status:
            continue
        if stage and state["current_stage"] != stage:
            continue
        output = state["output_path"]
        state_paths.add(output)
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
    for flowcell_path, group in inventory.groupby("flowcell_path", sort=True):
        if flowcell_path in state_paths:
            continue
        if stage and stage != "legacy":
            continue
        if status and status not in {"completed", "archived", "legacy"}:
            continue
        archived_values = [value for value in group["archived"] if value != "0"]
        legacy_status = "archived" if archived_values else "completed"
        if status and status not in ("legacy", legacy_status):
            continue
        rows.append(
            {
                "run_id": Path(flowcell_path).name,
                "status": legacy_status,
                "stage": "legacy",
                "origin": "legacy_inventory",
                "projects": ",".join(sorted(set(group["project"]))),
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


def pretty_print(df):
    if df.empty:
        print("No matching flowcells.")
    else:
        print(df.to_string(index=False))


def main():
    _remove_executable_directory_from_import_path()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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

    parser_initialize = subparsers.add_parser(
        "initialize", help="Explicitly initialize ambiguous existing output."
    )
    parser_initialize.set_defaults(func=initialize_flowcell)
    parser_initialize.add_argument("flowcell", help="Run ID.")
    parser_initialize.add_argument(
        "--from",
        dest="from_stage",
        required=True,
        choices=["demultiplexing", "analysis"],
    )
    parser_initialize.add_argument("--dry-run", action="store_true")
    parser_initialize.add_argument("--force", action="store_true")
    parser_initialize.add_argument("--reason")
    parser_initialize.add_argument("--refresh-inputs", action="store_true")

    parser_list = subparsers.add_parser("list", help="List state-backed and legacy flowcells.")
    parser_list.set_defaults(
        func=lambda **kwargs: combined_list(kwargs.get("status"), kwargs.get("stage")),
        print_res=True,
    )
    parser_list.add_argument("--status")
    parser_list.add_argument("--stage")

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

    parser_retry = subparsers.add_parser(
        "retry-notifications", help="Retry saved completion mail without rerunning processing."
    )
    parser_retry.add_argument("flowcell", help="Run ID or flowcell path.")
    parser_retry.add_argument("--kind", choices=["processed", "finalized"])
    parser_retry.add_argument(
        "--retry-uncertain",
        action="store_true",
        help="Allow resending after uncertain SMTP acceptance; duplicates are possible.",
    )
    parser_retry.set_defaults(func=retry_notifications)

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
