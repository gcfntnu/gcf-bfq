"""One-time SampleSheet scheduling seeds; persistent manager controls win thereafter."""

import copy
import re

from pathlib import Path

from bcl2fastq_pipeline import project_holds
from bcl2fastq_pipeline.state import StateConflictError, StateValidationError, utcnow


def resolve(custom):
    raw = str(custom.get("InitialHold", "")).strip().lower()
    if raw not in {"", "true", "false"}:
        raise ValueError(f"InitialHold {raw!r} must be true or false (blank defaults to false)")
    return {"held": raw == "true", "source": "samplesheet" if raw else "default"}


def add_preflight_plan(selection, result):
    """Use the shared parser's effective options, including legacy shared sheets."""
    entries = {}
    sources = (
        [
            (pair.project_id, pair.sample_sheet, result.project_results.get(pair.project_id))
            for pair in selection.pairs
        ]
        if selection.multi_pair
        else [(None, selection.sample_sheet, result)]
    )
    for project, sheet, parsed in sources:
        if parsed is None:
            continue
        try:
            value = resolve(parsed.custom_options)
        except ValueError as error:
            result.add(
                "error", "scheduling.initial_hold", str(error), file=str(sheet), project_id=project
            )
            continue
        identity = next(
            (item for item in result.inputs if Path(item["path"]) == Path(sheet).absolute()), {}
        )
        names = (
            [project] if project else [p["project_id"] for p in parsed.summary.get("projects", [])]
        )
        for name in names:
            entries[name] = {
                **value,
                "sample_sheet": Path(sheet).name,
                "sample_sheet_sha256": identity.get("sha256"),
            }
    if result.ok and entries:
        result.summary["project_initial_holds"] = entries
    return result


def _historical(state):
    return bool(
        state.get("started_at")
        or state.get("attempt")
        or state.get("attempts")
        or state.get("completed_at")
        or state.get("project_progress")
        or state.get("project_completion")
        or state["origin"] == "legacy_rerun"
        or state["current_stage"] in {"reporting", "finalization"}
    )


def prepare(state):
    """Classify older work before begin_attempt changes its start/status evidence."""
    if "project_initial_holds" not in state:
        state["project_initial_holds"] = {}
        if _historical(state):
            state["project_initial_holds_legacy"] = True


def _decision(state, project):
    if project in state.get("project_controls", {}):
        return "existing_control"
    if state.get("project_initial_holds_legacy") or (
        "project_initial_holds" not in state and _historical(state)
    ):
        return "historical"
    return "seeded"


def initialize(store, run_id, result):
    """Commit every resolved boolean atomically before batch admission or claims.

    The caller holds the execution lease. The short mutation lock serializes this
    transaction with manager requests; either ordering preserves explicit controls.
    """
    entries = result.summary.get("project_initial_holds", {})
    if not entries:
        return

    def update(state):
        if state["status"] != "running":
            raise StateConflictError("Initial project eligibility requires the running attempt")
        prepare(state)
        state["projects"] = sorted(set(state["projects"]) | set(entries))
        saved = state["project_initial_holds"]
        for project, entry in sorted(entries.items()):
            if project in saved:
                continue
            decision = _decision(state, project)
            record = {
                **copy.deepcopy(entry),
                "decision": decision,
                "at": utcnow(),
                "attempt": state["attempt"],
            }
            if decision == "seeded" and entry["held"]:
                project_holds.request(
                    state,
                    [project],
                    hold=entry["held"],
                    reason="InitialHold from effective SampleSheet",
                    source=entry["source"],
                )
            record["effective_held"] = project_holds.held(state, project)
            saved[project] = record
            state.setdefault("project_control_history", []).append(
                {"project": project, "action": "initial_eligibility", **copy.deepcopy(record)}
            )
        return state

    return store.mutate(run_id, update)


def require_initialized(state, projects):
    if (
        "project_initial_holds" in state
        and not state.get("project_initial_holds_legacy")
        and (
            not state["project_initial_holds"]
            or set(projects) - state["project_initial_holds"].keys()
        )
    ):
        raise StateConflictError(
            "Project eligibility must be initialized before admission or claim"
        )


def validate(state):
    records = state.get("project_initial_holds", {})
    if not isinstance(records, dict):
        raise StateValidationError("Invalid initial project holds")
    if type(state.get("project_initial_holds_legacy", False)) is not bool:
        raise StateValidationError("Invalid historical initial hold policy")
    for project, entry in records.items():
        if (
            project not in state["projects"]
            or not isinstance(entry, dict)
            or type(entry.get("held")) is not bool
            or type(entry.get("effective_held")) is not bool
            or entry.get("source") not in {"samplesheet", "default"}
            or entry.get("decision") not in {"seeded", "existing_control", "historical"}
            or not isinstance(entry.get("sample_sheet"), str)
            or not re.fullmatch(r"[0-9a-f]{64}", str(entry.get("sample_sheet_sha256")))
            or not isinstance(entry.get("at"), str)
            or type(entry.get("attempt")) is not int
        ):
            raise StateValidationError("Invalid initial project hold record")


def print_plan(entries, *, state=None, projects=None):
    if not entries:
        return
    print("Initial project eligibility (applied once; release is manual):")
    for project in sorted(entries):
        if projects is not None and project not in projects:
            continue
        entry = entries[project]
        decision = "seeded"
        effective = entry["held"]
        if state is not None:
            decision = (
                "persisted"
                if project in state.get("project_initial_holds", {})
                else _decision(state, project)
            )
            if decision != "seeded":
                effective = project_holds.held(state, project)
        eligibility = "held" if effective else "released"
        print(
            f"  {project}: InitialHold={str(entry['held']).lower()} ({entry['source']}); "
            f"{eligibility} ({decision})"
        )
