"""Persistent eligibility and serial project batches under the flowcell lease.

Control requests and stage claims use the state store's short mutation lock.
The execution lease still belongs exclusively to the processing worker.
"""

import copy
import hashlib
import re
import socket

from pathlib import Path

from bcl2fastq_pipeline.state import StateConflictError, StateValidationError, utcnow


def validate(state):
    controls = state.get("project_controls", {})
    if not isinstance(controls, dict):
        raise StateValidationError("Invalid project controls")
    for project, control in controls.items():
        if (
            not isinstance(project, str)
            or not re.fullmatch(r"GCF-[A-Za-z0-9][A-Za-z0-9-]*", project)
            or not isinstance(control, dict)
            or type(control.get("held")) is not bool
            or type(control.get("applied")) is not bool
            or not isinstance(control.get("reason"), str)
            or not isinstance(control.get("requested_at"), str)
        ):
            raise StateValidationError("Invalid project hold")
    if not isinstance(state.get("project_control_history", []), list):
        raise StateValidationError("Invalid project control history")
    if "project_continuation" in state and type(state["project_continuation"]) is not bool:
        raise StateValidationError("Invalid project continuation")
    batch = state.get("project_batch")
    if batch is not None and (
        not isinstance(batch, dict)
        or not isinstance(batch.get("projects"), list)
        or any(not isinstance(p, str) for p in batch["projects"])
        or not set(batch["projects"]).issubset(state.get("projects", []))
        or batch.get("attempt") != state.get("attempt")
    ):
        raise StateValidationError("Invalid project batch")
    active = state.get("active_project_stage")
    if active is not None and (
        not isinstance(active, dict)
        or active.get("project") not in state.get("projects", [])
        or active.get("stage") not in {"analysis", "archiving"}
        or not isinstance(active.get("started_at"), str)
    ):
        raise StateValidationError("Invalid active project stage")


def held(state, project):
    return state.get("project_controls", {}).get(project, {}).get("held", False)


def _event(state, project, action, reason, *, source="manager"):
    state.setdefault("project_control_history", []).append(
        {
            "project": project,
            "action": action,
            "reason": reason,
            "at": utcnow(),
            "attempt": state["attempt"],
            "hostname": socket.gethostname(),
            "source": source,
        }
    )


def apply_requests(state):
    active = state.get("active_project_stage", {}).get("project")
    for project, control in state.get("project_controls", {}).items():
        if project != active and control["held"] != control["applied"]:
            control.update(applied=control["held"], applied_at=utcnow())
            _event(
                state,
                project,
                "held" if control["held"] else "released",
                control["reason"],
                source=control.get("source", "manager"),
            )


def request(state, projects, *, hold, reason, source="manager"):
    """Pure mutation used for both preview and the serialized committed request."""
    projects = set(projects)
    if not projects or projects - set(state["projects"]):
        raise StateConflictError("Select known projects from this flowcell")
    if not isinstance(reason, str) or not reason.strip():
        raise StateConflictError("A nonempty --reason is required")
    for project in sorted(projects):
        controls = state.setdefault("project_controls", {})
        control = controls.get(project)
        if control is not None and control["held"] == hold:
            continue
        controls[project] = {
            "held": hold,
            "applied": control["applied"] if control else False,
            "reason": reason.strip(),
            "requested_at": utcnow(),
            "source": source,
        }
        _event(
            state,
            project,
            "hold_requested" if hold else "release_requested",
            reason.strip(),
            source=source,
        )
    apply_requests(state)
    if state["status"] == "waiting":
        remaining = {
            p
            for p in state["projects"]
            if state.get("project_progress", {}).get(p, {}).get("finalization", {}).get("status")
            != "completed"
        }
        if eligible(state):
            state["status"] = "queued"
            state["stages"][state["current_stage"]]["status"] = "queued"
            state.pop("waiting_reason", None)
        elif any(not held(state, p) for p in remaining):
            continue_pending(state, remaining)
    return state


def unfinished(state):
    progress = state.get("project_progress", {})
    selected = (state.get("restart_request") or {}).get("projects")
    names = set(state["projects"] if selected is None else selected)
    return {
        project
        for project in names
        if progress.get(project, {}).get("finalization", {}).get("status") != "completed"
    }


def retain_pending_contexts(state, projects):
    """A scoped correction must not replace a sibling's suspended input authority."""
    contexts = state.setdefault("project_pending_contexts", {})
    selected = set(state["projects"] if projects is None else projects)
    if projects is not None and state.get("attempts"):
        context = state.get("continuation_context") or state["attempts"][-1]
        for project in set(state["projects"]) - selected:
            progress = state.get("project_progress", {}).get(project, {})
            if progress.get("finalization", {}).get("status") != "completed":
                contexts.setdefault(project, copy.deepcopy(context))
    for project in selected:
        contexts.pop(project, None)
    state.pop("continuation_context", None)


def continue_pending(state, remaining):
    """Finish a scoped attempt, then restore a suspended sibling batch safely."""
    from bcl2fastq_pipeline import project_scheduling  # noqa: PLC0415

    contexts = state.get("project_pending_contexts", {})
    current = state["attempts"][-1].get("project_schedule")
    if current:
        for project in remaining:
            frozen = contexts.get(project, {}).get("project_schedule") or {}
            if entry := frozen.get("projects", {}).get(project):
                current["projects"][project] = copy.deepcopy(entry)
        current["order"] = project_scheduling.ordered(current)
    order = project_scheduling.ordered(current, remaining) if current else sorted(remaining)
    available = [p for p in order if not held(state, p)]
    first = (available or sorted(remaining))[0]
    context = contexts.get(first)
    group = sorted(
        p for p in remaining if contexts.get(p, {}).get("attempt") == (context or {}).get("attempt")
    )
    # Failed work still needs an explicit recovery, even when its siblings deliver.
    progress = state.get("project_progress", {})
    failed = [
        (p, key)
        for p in remaining
        for key in ("analysis", "finalization")
        if progress.get(p, {}).get(key, {}).get("status") == "failed"
    ]
    if failed:
        boundary = "analysis" if any(key == "analysis" for _, key in failed) else "finalization"
        state.update(status="failed", current_stage=boundary, completed_at=None)
        state["stages"][boundary].update(status="failed", completed_at=None)
        state["last_error"] = {
            "summary": "Unfinished projects require explicit recovery: "
            + ", ".join(sorted({p for p, _ in failed})),
            "report_path": None,
        }
        return
    if context:
        state["continuation_context"] = copy.deepcopy(context)
    boundary = "analysis" if any(not analysis_complete(state, p) for p in group) else "finalization"
    state.update(
        status="queued" if available else "waiting",
        current_stage=boundary,
        completed_at=None,
        project_continuation=True,
    )
    state["stages"][boundary].update(status="queued", completed_at=None)
    state["restart_request"].update({"from": boundary, "projects": group})
    if context and context.get("project_correction_id"):
        state["restart_request"]["project_correction_id"] = context["project_correction_id"]
    state["attempts"][-1].update(outcome="running", completed_at=None)
    if not available:
        state["waiting_reason"] = "Unfinished projects are held: " + ", ".join(sorted(remaining))


def eligible(state):
    return {project for project in unfinished(state) if not held(state, project)}


def analysis_complete(state, project):
    progress = state.get("project_progress", {}).get(project, {})
    if progress:
        return progress.get("analysis", {}).get("status") == "completed"
    if state.get("attempts") and state["attempts"][-1].get("project_schedule"):
        request = state.get("restart_request") or {}
        if request.get("from", state["attempts"][-1]["restart_from"]) in {
            "demultiplexing",
            "analysis",
        }:
            return False
    return state["stages"]["analysis"]["status"] == "completed"


def admit(store, run_id):
    """Freeze one batch; a later release is admitted after its delivery phase."""

    def update(state):
        from bcl2fastq_pipeline import initial_holds  # noqa: PLC0415

        initial_holds.require_initialized(state, unfinished(state))
        apply_requests(state)
        if not state.get("project_batch"):
            state["project_batch"] = {
                "projects": sorted(eligible(state)),
                "admitted_at": utcnow(),
                "attempt": state["attempt"],
            }
        return state

    return store.mutate(run_id, update)


def claim(store, run_id, project, stage):
    """Linearize a stage start against controls before any project work begins."""
    claimed = False

    def update(state):
        nonlocal claimed
        from bcl2fastq_pipeline import initial_holds  # noqa: PLC0415

        initial_holds.require_initialized(state, [project])
        apply_requests(state)
        expected = "analysis" if stage == "analysis" else "finalization"
        if state["status"] != "running" or state["current_stage"] != expected:
            raise StateConflictError("Project claim requires its running stage")
        if held(state, project):
            return state
        if state.get("active_project_stage"):
            raise StateConflictError("Another project stage is already claimed")
        state["active_project_stage"] = {"project": project, "stage": stage, "started_at": utcnow()}
        key = "analysis" if stage == "analysis" else "finalization"
        state.setdefault("project_progress", {}).setdefault(project, {})[key] = {
            "status": "running",
            "attempt": state["attempt"],
            "started_at": utcnow(),
        }
        plan = state["attempts"][-1].get("project_schedule")
        if plan:
            plan["actual_order"][stage].append({"project": project, "started_at": utcnow()})
        claimed = True
        return state

    store.mutate(run_id, update)
    return claimed


def finish_claim(state):
    state.pop("active_project_stage", None)
    apply_requests(state)


def advance(store, run_id):
    """Close a batch without manufacturing a completed flowcell or new attempt."""

    def update(state):
        finish_claim(state)
        state.pop("project_batch", None)
        remaining = unfinished(state)
        if not remaining:
            return state
        available = eligible(state)
        boundary = (
            "analysis"
            if any(not analysis_complete(state, p) for p in remaining)
            else "finalization"
        )
        state["current_stage"] = boundary
        state["project_continuation"] = True
        state["status"] = "running" if available else "waiting"
        state["completed_at"] = None
        for name in ("analysis", "reporting", "finalization"):
            detail = state["stages"][name]
            if name == boundary:
                detail.update(status="queued", completed_at=None)
        if not available:
            state["waiting_reason"] = "Unfinished projects are held: " + ", ".join(
                sorted(remaining)
            )
        else:
            state.pop("waiting_reason", None)
        return state

    return store.mutate(run_id, update)


def save_config(cfg, store):
    if cfg.run.libprep_config is None:
        return
    snapshot = cfg.run.libprep_config

    def update(state):
        state["attempts"][-1].setdefault(
            "execution_libprep",
            {
                "source": snapshot.source,
                "content": snapshot.content.decode("utf-8"),
                "sha256": snapshot.sha256,
            },
        )
        return state

    store.mutate(cfg.run.run_id, update)


def restore_execution(cfg, state):
    """Continue the same immutable attempt after a daemon/container restart."""
    from configmaker.libprep import LibprepConfig  # noqa: PLC0415

    from bcl2fastq_pipeline import preflight, workflow_config  # noqa: PLC0415

    attempt = state.get("continuation_context") or state["attempts"][-1]
    record = attempt.get("execution_inputs")
    correction = (state.get("restart_request") or {}).get("project_correction_id")
    resume = (state.get("restart_request") or {}).get("analysis_resume")
    if record:
        for entry in record["inputs"] + record.get("geometry", []):
            if hashlib.sha256(Path(entry["path"]).read_bytes()).hexdigest() != entry["sha256"]:
                raise StateConflictError("Frozen execution inputs changed; use a validated rerun")
        selection = preflight.select_run_inputs(record["directory"], record["directory"])
        cfg.run.input_selection = selection
        custom, _ = workflow_config.parse_custom_options(selection.sample_sheet)
        cfg.run.apply_custom(custom, selection.sample_sheet, selection.submission_form)
    elif resume:
        # Resume retains its curated-path and repaired-workdir authority. The
        # existing resume inspector rechecks every hash under the worker lease.
        selection = preflight.select_run_inputs(state["output_path"], state["output_path"])
        custom, _ = workflow_config.parse_custom_options(selection.sample_sheet)
        cfg.run.apply_custom(custom, selection.sample_sheet, selection.submission_form)
        cfg.run.input_selection = selection
    elif not correction and state["current_stage"] == "analysis":
        raise StateConflictError("Missing frozen execution inputs; use a validated rerun")
    snapshot = attempt.get("execution_libprep")
    if snapshot:
        cfg.run.libprep_config = LibprepConfig(
            snapshot["source"], snapshot["content"].encode("utf-8")
        )
        if cfg.run.libprep_config.sha256 != snapshot["sha256"]:
            raise StateConflictError("Frozen library-prep configuration changed")
    cfg.run.project_schedule = copy.deepcopy(attempt.get("project_schedule"))
    if cfg.run.project_schedule:
        cfg.run.project_schedule["attempt"] = state["attempt"]
    metadata = state["stages"]["analysis"]["metadata"]
    cfg.run.project_workdirs = {
        project: Path(identity["path"])
        for project, identity in metadata.get("workdirs", {}).items()
    }
    cfg.run.project_workflows = dict(metadata.get("project_workflows", {}))
    if state["current_stage"] in {"reporting", "finalization"}:
        workflow_config.prepare_execution(
            cfg,
            state["current_stage"],
            workflow=metadata.get("workflow"),
            project_workflows=metadata.get("project_workflows"),
        )
    else:
        workflow_config.prepare_projects(cfg)


def print_status(state):
    from flowcell_manager.project_listing import project_rows  # noqa: PLC0415

    from bcl2fastq_pipeline import project_scheduling  # noqa: PLC0415

    if state.get("waiting_reason"):
        print("  Waiting: " + state["waiting_reason"])
    if batch := state.get("project_batch"):
        print("  Current batch: " + ", ".join(batch["projects"]))
    attempts = state.get("attempts", [])
    plan = attempts[-1].get("project_schedule") if attempts else None
    pending = eligible(state)
    order = (
        project_scheduling.ordered(plan, pending)
        if plan and pending.issubset(plan["projects"])
        else sorted(pending)
    )
    analyses = [p for p in order if not analysis_complete(state, p)]
    deliveries = [p for p in order if analysis_complete(state, p)]
    print("  Eligible analysis order: " + (" -> ".join(analyses) or "(none)"))
    print("  Eligible delivery order: " + (" -> ".join(deliveries) or "(none)"))
    for row in project_rows(state):
        print(
            f"  {row['project']}: {row['status']}; {row['eligibility']}; next {row['stage']}; "
            f"priority {row['priority']}; {row['waiting_reason']}"
        )
        if seed := state.get("project_initial_holds", {}).get(row["project"]):
            print(
                f"    InitialHold={str(seed['held']).lower()} ({seed['source']}); "
                f"{seed['decision']}, attempt {seed['attempt']}; "
                f"SampleSheet sha256 {seed['sample_sheet_sha256']}"
            )
