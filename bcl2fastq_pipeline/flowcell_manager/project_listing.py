"""Read-only project rows from canonical records, without consulting live inputs.

Project commits take precedence over shared flowcell stages. A project that has
not started waits at its next boundary; only demultiplexing, recorded preflight
failure and reporting are shared work. Older failed/running analysis or archive
stages cannot identify a project's outcome and are labelled unavailable.
"""

UNAVAILABLE = "unavailable"
UNIDENTIFIED = "<unidentified>"
COLUMNS = [
    "project",
    "run_id",
    "status",
    "stage",
    "analysis",
    "delivery",
    "last_delivery",
    "priority",
    "priority_source",
    "parent_status",
    "parent_stage",
    "origin",
    "output_path",
    "archived",
    "analysis_generation",
    "last_delivery_generation",
]
_STAGES = ("demultiplexing", "analysis", "reporting", "finalization")
_WAITING = {"preparing", "queued"}


def _stage(state, name):
    return state.get("stages", {}).get(name, {})


def _selected(state, project):
    projects = (state.get("restart_request") or {}).get("projects")
    return projects is None or project in projects


def _preparing_invalidates(state, project, stage):
    request = state.get("restart_request") or {}
    start = request.get("from")
    return (
        state["status"] == "preparing"
        and _selected(state, project)
        and start in _STAGES
        and _STAGES.index(start) <= _STAGES.index(stage)
    )


def _completed_stage(state, project, name):
    return _stage(state, name).get("status") == "completed" and not _preparing_invalidates(
        state, project, name
    )


def _attempt_plan(state):
    for attempt in reversed(state.get("attempts", [])):
        if attempt.get("attempt") == state.get("attempt"):
            return attempt.get("project_schedule") or {}
    return {}


def _recorded_entry(plan, project):
    return (plan or {}).get("projects", {}).get(project) or {}


def _priority(state, project, progress, completion):
    analysis = progress.get("analysis", {})
    if analysis.get("scheduling"):
        return analysis["scheduling"]
    request = state.get("restart_request") or {}
    retained = not _selected(state, project)
    producing_plan = (
        completion.get("stages", {}).get("analysis", {}).get("metadata", {}).get("project_schedule")
    )
    if retained:
        return _recorded_entry(producing_plan, project)
    # A queued new analysis has no frozen plan yet. Its preceding attempt's
    # priorities are historical, even while preparation still retains metadata.
    reanalysis_waiting = (
        state["status"] in _WAITING
        and request.get("from", state["current_stage"]) in {"demultiplexing", "analysis"}
        and not request.get("analysis_resume")
    )
    if reanalysis_waiting:
        return {}
    saved = _stage(state, "analysis").get("metadata", {}).get("project_schedule")
    # Completed analysis metadata records its producing plan; an active attempt
    # can have a newer plan for other projects or a downstream-only retry.
    if analysis.get("status") == "completed" or _completed_stage(state, project, "analysis"):
        return (
            _recorded_entry(saved, project)
            or _recorded_entry(producing_plan, project)
            or _recorded_entry(_attempt_plan(state), project)
        )
    return _recorded_entry(_attempt_plan(state), project)


def _last_delivery(state, project, progress, completion):
    if completion.get("completed_at"):
        return completion["completed_at"], completion.get("generation") or UNAVAILABLE
    finalization = progress.get("finalization", {})
    if finalization.get("status") == "completed" and finalization.get("completed_at"):
        return finalization["completed_at"], finalization.get("generation") or UNAVAILABLE
    if _stage(state, "finalization").get("status") == "completed":
        # Preparation may retain the previous successful stage temporarily.
        # It remains useful historical evidence, independent of current work.
        completed = state.get("completed_at") or _stage(state, "finalization").get("completed_at")
        return completed or UNAVAILABLE, UNAVAILABLE
    return UNAVAILABLE, UNAVAILABLE


def _outcome(state, detail, stage):
    outcome = detail.get("status", UNAVAILABLE)
    if (
        outcome in {"running", "failed"}
        and state["status"] == "interrupted"
        and state["current_stage"] == stage
        and detail.get("attempt", state.get("attempt")) == state.get("attempt")
    ):
        return "interrupted"
    return outcome


def _completion_is_current(state, project, progress, completion):
    if not completion or progress:
        return False
    if not _selected(state, project):
        return True
    request = state.get("restart_request")
    return (not request and completion.get("attempt") == state.get("attempt")) or _completed_stage(
        state, project, "finalization"
    )


def _pending_analysis(state, project, progress):
    parent = state["status"]
    boundary = state["current_stage"]
    if boundary == "demultiplexing" and _selected(state, project):
        return parent, boundary, "pending", "pending"
    preflight = _stage(state, boundary).get("metadata", {}).get("input_preflight", {})
    if (
        boundary == "analysis"
        and _selected(state, project)
        and preflight.get("status") == "failed"
        and parent in {"failed", "interrupted"}
    ):
        return parent, boundary, "pending", "pending"
    if parent in _WAITING and _selected(state, project):
        return parent, "analysis", "pending", "pending"
    # An explicit reset entry, another project's progress, or a frozen plan
    # establishes the newer project-aware execution model. Absence means this
    # project has not started; it does not inherit its sibling's failure.
    known_pending = (
        project in state.get("project_progress", {})
        or bool(state.get("project_progress"))
        or bool(_recorded_entry(_attempt_plan(state), project))
    )
    if known_pending or progress:
        return "pending", "analysis", "pending", "pending"
    return UNAVAILABLE, "analysis", UNAVAILABLE, "pending"


def _progress_values(state, project, progress, completion):  # noqa: PLR0911
    analysis_detail = progress.get("analysis", {})
    delivery_detail = progress.get("finalization", {})
    analysis = _outcome(state, analysis_detail, "analysis")
    delivery = _outcome(state, delivery_detail, "finalization")
    generation = analysis_detail.get("generation")
    # A delivery of an older analysis must never complete a newer generation.
    if (
        generation
        and delivery_detail.get("generation")
        and delivery_detail["generation"] != generation
    ):
        delivery = UNAVAILABLE
    if delivery == "completed":
        status = "archived" if state.get("archive", {}).get("status") == "archived" else "completed"
        return status, "finalization", "completed", "completed"
    if delivery in {"running", "failed", "interrupted"}:
        return delivery, "finalization", analysis, delivery
    if analysis in {"running", "failed", "interrupted"}:
        return analysis, "analysis", analysis, "pending"
    if _completion_is_current(state, project, progress, completion):
        status = "archived" if state.get("archive", {}).get("status") == "archived" else "completed"
        return status, "finalization", "completed", "completed"
    if not analysis_detail and _completed_stage(state, project, "analysis"):
        analysis = "completed"
    if not progress and _completed_stage(state, project, "finalization"):
        status = "archived" if state.get("archive", {}).get("status") == "archived" else "completed"
        return status, "finalization", analysis, "completed"
    if analysis != "completed":
        return _pending_analysis(state, project, progress)
    # Reporting is shared; finalization has individual commit/start records.
    # An unstarted archive never inherits another project's running/failure.
    boundary = "finalization" if _completed_stage(state, project, "reporting") else "reporting"
    parent = state["status"]
    if parent in _WAITING and _selected(state, project):
        return parent, boundary, "completed", "pending"
    if boundary == "reporting" and state["current_stage"] == "reporting":
        return parent, boundary, "completed", "pending"
    historical_archive = (
        boundary == "finalization"
        and "project_progress" not in state
        and state["current_stage"] == "finalization"
        and parent in {"running", "failed", "interrupted"}
    )
    if historical_archive:
        return UNAVAILABLE, boundary, "completed", UNAVAILABLE
    return "pending", boundary, "completed", "pending"


def project_rows(state):
    """Project one state document without modifying it or inspecting its files.

    Recorded project progress owns current outcomes. Retained completion owns
    last-delivery evidence, and owns current completion only for unselected or
    otherwise still-current work. Completed shared stages support older state;
    missing individual failure/start evidence remains explicitly unavailable.
    """
    progress = state.get("project_progress", {})
    completion = state.get("project_completion", {})
    names = set(state.get("projects", [])) | progress.keys() | completion.keys()
    names.update(_attempt_plan(state).get("projects", {}))
    rows = []
    for project in sorted(names) or [UNIDENTIFIED]:
        current = progress.get(project, {})
        prior = completion.get(project, {})
        status, stage, analysis, delivery = _progress_values(state, project, current, prior)
        priority = _priority(state, project, current, prior)
        last_delivery, last_generation = _last_delivery(state, project, current, prior)
        generation = current.get("analysis", {}).get("generation")
        if not generation and _completion_is_current(state, project, current, prior):
            generation = prior.get("generation")
        rows.append(
            {
                "project": project,
                "run_id": state["run_id"],
                "status": status,
                "stage": stage,
                "analysis": analysis,
                "delivery": delivery,
                "last_delivery": last_delivery,
                "priority": priority.get("priority", UNAVAILABLE),
                "priority_source": priority.get("priority_source", UNAVAILABLE),
                "parent_status": state["status"],
                "parent_stage": state["current_stage"],
                "origin": state["origin"],
                "output_path": state["output_path"],
                "archived": state.get("archive", {}).get("archived_at") or "0",
                "analysis_generation": generation or UNAVAILABLE,
                "last_delivery_generation": last_generation,
            }
        )
    return rows


def legacy_project_row(project, run_id, flowcell_path, archived="0", timestamp="0"):
    """Project historical inventory without inventing analysis or priorities.

    Inventory timestamps record completed processing, so they support historical
    delivery time but contain no generation or separate stage outcomes.
    """
    status = "archived" if archived != "0" else "completed"
    return {
        "project": project or UNIDENTIFIED,
        "run_id": run_id,
        "status": status,
        "stage": "legacy",
        "analysis": UNAVAILABLE,
        "delivery": UNAVAILABLE,
        "last_delivery": timestamp if timestamp != "0" else UNAVAILABLE,
        "priority": UNAVAILABLE,
        "priority_source": UNAVAILABLE,
        "parent_status": status,
        "parent_stage": "legacy",
        "origin": "legacy_inventory",
        "output_path": flowcell_path,
        "archived": archived,
        "analysis_generation": UNAVAILABLE,
        "last_delivery_generation": UNAVAILABLE,
    }
