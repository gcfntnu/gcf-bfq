"""Exact decimal project ordering, frozen with each processing attempt."""

import copy
import re

from decimal import Decimal, InvalidOperation
from pathlib import Path

from bcl2fastq_pipeline.state import FlowcellStateStore, StateConflictError, utcnow

_DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")


def resolve_priority(custom):
    """Keep values as decimal strings in JSON; never round through a float."""
    raw = str(custom.get("Priority", "")).strip()
    if not raw:
        return {"priority": "1.0", "priority_source": "default"}
    try:
        value = Decimal(raw)
    except InvalidOperation:
        value = None
    if value is None or not _DECIMAL.fullmatch(raw) or not value.is_finite() or value < 0:
        raise ValueError(f"Priority {raw!r} must be a finite, non-negative decimal value")
    return {"priority": str(value), "priority_source": "samplesheet"}


def ordered(plan, projects=None):
    """Apply the same numeric order independently to each selected stage."""
    entries = plan["projects"]
    selected = set(entries) if projects is None else set(projects)
    unknown = selected - entries.keys()
    if unknown:
        raise StateConflictError(
            f"Projects absent from the frozen scheduling plan: {sorted(unknown)}"
        )
    return sorted(selected, key=lambda project: (Decimal(entries[project]["priority"]), project))


def add_preflight_plan(selection, result):
    """Add BFQ scheduling policy to the shared validator's parsed metadata."""
    entries = {}
    if selection.multi_pair:
        sources = [
            (pair.project_id, pair.sample_sheet, result.project_results.get(pair.project_id))
            for pair in selection.pairs
        ]
    else:
        sources = [(None, selection.sample_sheet, result)]
    for project, sheet, parsed in sources:
        if parsed is None:
            continue
        try:
            priority = resolve_priority(parsed.custom_options)
        except ValueError as error:
            result.add(
                "error", "scheduling.priority", str(error), file=str(sheet), project_id=project
            )
            continue
        identity = next(
            (item for item in result.inputs if Path(item["path"]) == Path(sheet).absolute()), {}
        )
        projects = (
            [project]
            if project
            else [item["project_id"] for item in parsed.summary.get("projects", [])]
        )
        for name in projects:
            entries[name] = {
                **priority,
                "sample_sheet_sha256": identity.get("sha256"),
            }
    if result.ok and entries:
        plan = {"schema_version": 1, "projects": entries}
        plan["order"] = ordered(plan)
        result.summary["project_schedule"] = plan
    return result


def print_plan(plan, *, projects=None):
    if not plan:
        return
    order = ordered(plan, projects)
    label = f"for attempt {plan['attempt']}" if "attempt" in plan else "preview"
    print(f"Project scheduling {label} (eligible analyses before archiving each released batch):")
    for project in order:
        entry = plan["projects"][project]
        print(f"  {project}: priority {entry['priority']} ({entry['priority_source']})")
    print("  Analysis order: " + " -> ".join(order))
    print("  Archiving order: " + " -> ".join(order))


def saved_plan(state):
    """Resume/downstream retries retain the order of the producing analysis."""
    if (state.get("restart_request") or {}).get("analysis_resume"):
        for attempt in reversed(state.get("attempts", [])):
            if attempt.get("project_schedule"):
                return attempt["project_schedule"]
    return state["stages"]["analysis"]["metadata"].get("project_schedule")


def _bind(cfg, store, plan):
    current = store.read(cfg.run.run_id)
    frozen = copy.deepcopy(plan)
    frozen.setdefault("resolved_in_attempt", current["attempt"])
    frozen.update(attempt=current["attempt"], actual_order={"analysis": [], "archiving": []})

    def record(state):
        if state["attempt"] != frozen["attempt"] or state["status"] != "running":
            raise StateConflictError("Attempt changed while freezing project scheduling")
        state["attempts"][-1]["project_schedule"] = copy.deepcopy(frozen)
        state["projects"] = sorted(set(state["projects"]) | set(frozen["projects"]))
        return state

    store.mutate(cfg.run.run_id, record)
    cfg.run.project_schedule = frozen


def restore(cfg, store, state):
    cfg.run.project_schedule = None
    if state["current_stage"] in {"reporting", "finalization"} or (
        state.get("restart_request") or {}
    ).get("analysis_resume"):
        if plan := saved_plan(state):
            _bind(cfg, store, plan)


def freeze(cfg, store, result):
    """Bind once, after effective input copies are verified and before processing."""
    if plan := result.summary.get("project_schedule"):
        if cfg.run.project_schedule is None:
            _bind(cfg, store, plan)
        elif plan["projects"] != cfg.run.project_schedule["projects"]:
            raise StateConflictError(
                "Effective inputs differ from the frozen project scheduling plan"
            )


def project_order(cfg, projects):
    plan = getattr(cfg.run, "project_schedule", None)
    # Historical downstream state has no scheduling plan. Do not infer a new
    # order from today's editable inputs when retrying its completed analysis.
    return ordered(plan, projects) if plan else sorted(projects)


def completed_analysis_plan(cfg, state):
    """Retain unselected projects' producing priorities across scoped corrections."""
    plan = copy.deepcopy(cfg.run.project_schedule)
    selected = cfg.run.selected_projects
    if not plan or selected is None:
        return plan
    retained = {}
    for project in plan["projects"].keys() - selected:
        completion = state.get("project_completion", {}).get(project, {})
        previous = completion.get("stages", {}).get("analysis", {}).get("metadata", {})
        previous = previous.get("project_schedule") or {}
        entry = previous.get("projects", {}).get(project) or state.get("project_progress", {}).get(
            project, {}
        ).get("analysis", {}).get("scheduling")
        plan["projects"][project] = (
            copy.deepcopy(entry) if entry else {**resolve_priority({}), "sample_sheet_sha256": None}
        )
        retained[project] = completion.get("attempt")
    plan["retained_projects"] = retained
    plan["order"] = ordered(plan)
    return plan


def record_start(cfg, stage, project):
    """Persist actual starts before launching a project's external work."""
    plan = getattr(cfg.run, "project_schedule", None)
    if plan is None:
        return
    if getattr(cfg.run, "project_lifecycle", None) is not None:
        # The lifecycle records the start atomically with its eligibility claim.
        return
    store = FlowcellStateStore(cfg.static.paths.manager_dir)
    expected_stage = "analysis" if stage == "analysis" else "finalization"

    def record(state):
        if (
            state["attempt"] != plan["attempt"]
            or state["status"] != "running"
            or state["current_stage"] != expected_stage
        ):
            raise StateConflictError("Project stage no longer belongs to the running attempt")
        current = state["attempts"][-1]["project_schedule"]
        current["actual_order"][stage].append({"project": project, "started_at": utcnow()})
        return state

    saved = store.mutate(cfg.run.run_id, record)
    cfg.run.project_schedule = copy.deepcopy(saved["attempts"][-1]["project_schedule"])
