"""Successful wall-clock processing cost, independent of daemon/SMTP lifetime.

Each step keeps its executions in canonical state. Only its latest, non-invalidated
success is reportable. Reporting also includes the independently retained early QC
execution; the two intervals never overlap.
"""

from __future__ import annotations

import math
import time

from contextlib import contextmanager
from datetime import UTC, datetime, timedelta

STEPS = {
    "demultiplexing": ("Demultiplexing", "demultiplexing"),
    "fastq_checksums": ("FASTQ MD5 checksums", "demultiplexing"),
    "analysis": ("Analysis", "analysis"),
    "reporting": ("Reporting", "reporting"),
    "archiving": ("Archiving", "finalization"),
    "archive_checksums": ("Archive MD5 checksums", "finalization"),
}
STAGES = ("demultiplexing", "analysis", "reporting", "finalization")
OUTCOMES = {"running", "completed", "failed", "interrupted"}


def utcnow():
    return datetime.now(UTC).isoformat()


def _duration(value):
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
        and value >= 0
    )


def validate(timings):
    """Validate optional schema-v1 records without upgrading historical state."""
    if not isinstance(timings, dict) or timings.keys() - STEPS.keys():
        raise ValueError("Invalid processing timing steps")
    for records in timings.values():
        if not isinstance(records, list):
            raise ValueError("Processing timing history must be a list")
        for record in records:
            if (
                not isinstance(record, dict)
                or record.get("outcome") not in OUTCOMES
                or type(record.get("attempt")) is not int
                or record["attempt"] < 1
                or not isinstance(record.get("started_at"), str)
                or "completed_at" not in record
                or "duration_seconds" not in record
            ):
                raise ValueError("Invalid processing timing execution")
            duration = record["duration_seconds"]
            if duration is not None and not _duration(duration):
                raise ValueError("Invalid processing timing duration")
            if record["outcome"] == "completed" and (
                duration is None or not isinstance(record["completed_at"], str)
            ):
                raise ValueError("Successful processing timing requires completion and duration")
            if record["outcome"] == "running" and (
                duration is not None or record["completed_at"] is not None
            ):
                raise ValueError("Running processing timing cannot have a duration")


def invalidate(state, from_stage):
    """Retain history while invalidating precisely the restarted output stages."""
    timings = state.setdefault("processing_timings", {})
    boundary = STAGES.index(from_stage)
    now = utcnow()
    for step, (_label, stage) in STEPS.items():
        if STAGES.index(stage) < boundary:
            continue
        # Empty histories also prevent inferred legacy completion during preparation.
        records = timings.setdefault(step, [])
        if records and not records[-1].get("invalidated_at"):
            records[-1]["invalidated_at"] = now
    return state


def interrupt(state):
    """A hard-crash duration is unknown, never extended until recovery time."""
    for records in state.get("processing_timings", {}).values():
        if records and records[-1]["outcome"] == "running":
            records[-1].update(outcome="interrupted", recovered_at=utcnow())
    return state


def start(state, step):
    if step not in STEPS:
        raise ValueError(f"Unsupported timing step: {step}")
    if state["status"] != "running":
        raise RuntimeError("Processing timing requires an active run attempt")
    records = state.setdefault("processing_timings", {}).setdefault(step, [])
    now = utcnow()
    if records:
        if records[-1]["outcome"] == "running":
            raise RuntimeError(f"Processing timer already running: {step}")
        records[-1].setdefault("invalidated_at", now)
    records.append(
        {
            "attempt": state["attempt"],
            "started_at": now,
            "completed_at": None,
            "duration_seconds": None,
            "outcome": "running",
        }
    )
    return state


def finish(state, step, duration, outcome):
    record = state["processing_timings"][step][-1]
    if record["outcome"] != "running" or record["attempt"] != state["attempt"]:
        raise RuntimeError(f"Processing timer is not running in this attempt: {step}")
    record.update(completed_at=utcnow(), duration_seconds=duration, outcome=outcome)
    return state


@contextmanager
def measure(store, run_id, step):
    """Persist success before the next operation; measure each enclosing interval once."""
    store.mutate(run_id, lambda state: start(state, step))
    started = time.monotonic()
    try:
        yield
    except BaseException as error:
        duration = time.monotonic() - started
        outcome = "failed" if isinstance(error, Exception) else "interrupted"
        store.mutate(run_id, lambda state: finish(state, step, duration, outcome))
        raise
    else:
        duration = time.monotonic() - started
        store.mutate(run_id, lambda state: finish(state, step, duration, "completed"))


def _result(state, step):
    timings = state.get("processing_timings", {})
    records = timings.get(step, [])
    if records:
        record = records[-1]
        if record["outcome"] == "completed" and not record.get("invalidated_at"):
            return {"status": "completed", "duration_seconds": record["duration_seconds"]}
    if step not in timings and state["stages"][STEPS[step][1]]["status"] == "completed":
        return {"status": "unavailable", "duration_seconds": None}
    return {"status": "not_completed", "duration_seconds": None}


def snapshot(state):
    """Capture the six-row breakdown at notification time, including reused successes."""
    steps = {step: _result(state, step) for step in STEPS}
    reporting = steps["reporting"]
    qc = state.get("sequencing_qc")
    if qc:
        if qc["status"] == "completed":
            reporting["duration_seconds"] = (reporting["duration_seconds"] or 0) + qc[
                "duration_seconds"
            ]
        elif reporting["status"] == "completed":
            reporting["status"] = "not_completed"
    return {
        "steps": steps,
        "total_seconds": sum(item["duration_seconds"] or 0 for item in steps.values()),
        "incomplete": any(item["status"] == "not_completed" for item in steps.values()),
        "unavailable": any(item["status"] == "unavailable" for item in steps.values()),
    }


def format_duration(seconds):
    return str(timedelta(seconds=round(seconds)))


def format_summary(timing):
    lines = []
    for step, (label, _stage) in STEPS.items():
        result = timing["steps"][step]
        seconds = result["duration_seconds"]
        if result["status"] == "completed":
            value = format_duration(seconds)
        else:
            value = "Timing unavailable" if result["status"] == "unavailable" else "Not completed"
            if seconds is not None:
                value += f" ({format_duration(seconds)} recorded)"
        lines.append(f"{label}: {value}")
    qualifiers = []
    if timing["incomplete"]:
        qualifiers.append("subtotal of completed steps")
    if timing["unavailable"]:
        qualifiers.append("partial; some timings unavailable")
    qualifier = f" ({'; '.join(qualifiers)})" if qualifiers else ""
    lines.append(f"Total processing time{qualifier}: {format_duration(timing['total_seconds'])}")
    return "\n".join(lines)
