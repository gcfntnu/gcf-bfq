"""Versioned, durable flowcell state for BFQ.

State lives outside flowcell output directories so operator cleanup and archival can
never erase the canonical processing record.  All writes are serialized by a
per-run advisory lock and committed with atomic os.replace().
"""

from __future__ import annotations

import copy
import fcntl
import gzip
import hashlib
import json
import os
import shutil
import socket
import subprocess

from collections.abc import Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path

from bcl2fastq_pipeline.config import parse_custom_options

SCHEMA_VERSION = 1
STAGES = ("demultiplexing", "analysis", "reporting", "finalization")
STAGE_STATUSES = {"pending", "queued", "running", "completed", "failed", "skipped"}
DELIVERY_STATUSES = {"pending", "sending", "sent", "failed", "uncertain", "superseded"}
DELIVERY_ATTEMPT_STATUSES = {"sending", "sent", "failed", "uncertain"}
RUN_STATUSES = {
    "preparing",
    "queued",
    "running",
    "completed",
    "failed",
    "interrupted",
    "archived",
}
ORIGINS = {"new", "legacy_rerun", "restored_legacy_fastq"}
INPUT_FILES = {"SampleSheet.csv", "Sample-Submission-Form.xlsx"}


class StateError(RuntimeError):
    """Base class for flowcell-state failures."""


class StateValidationError(StateError):
    """Raised when a state document does not conform to schema v1."""


class StateConflictError(StateError):
    """Raised for an illegal or conflicting state transition."""


class ExecutionLeaseError(StateConflictError):
    """Raised when another process owns the run execution lease."""


class DeliveryUncertainError(RuntimeError):
    """The relay may have accepted mail, so retrying could deliver a duplicate."""


def utcnow() -> str:
    return datetime.now(UTC).isoformat()


def _safe_run_id(run_id: str) -> str:
    run_id = str(run_id).strip()
    if not run_id or Path(run_id).name != run_id or run_id in {".", ".."}:
        raise StateValidationError(f"Invalid run id: {run_id!r}")
    return run_id


def resolve_output_path(output_path: Path | str, run_id: str, cfg) -> Path:
    """Resolve legacy bare IDs and require agreement with the daemon's output root.

    A bare run ID in historical inventory/state is relative to outputDir, never
    the operator's working directory. Other ambiguous paths require correction.
    Existing aliases (including bind mounts) are accepted only for the same file.
    """
    run_id = _safe_run_id(run_id)
    expected = cfg.static.paths.output_dir / run_id
    if not expected.is_absolute():
        raise StateConflictError("[Paths] outputDir must be absolute; check /config/bcl2fastq.ini")
    recorded = Path(output_path)
    if recorded == Path(run_id):
        return expected
    if not recorded.is_absolute():
        raise StateConflictError(
            f"Ambiguous relative output path {str(recorded)!r} for {run_id}; "
            f"configured output is {expected}. Only a bare run ID can be resolved automatically."
        )
    if recorded == expected:
        return expected
    try:
        same_location = recorded.samefile(expected)
    except OSError:
        same_location = False
    if not same_location:
        raise StateConflictError(
            f"Output path mismatch for {run_id}: recorded {recorded}; configured {expected}. "
            "Check outputDir and mounts before retrying; no output path was changed."
        )
    return expected


def output_entries(output: Path, *, allow_missing: bool = False) -> list[Path]:
    """Inspect a directory without treating a missing/unreadable mount as empty."""
    try:
        return list(output.iterdir())
    except FileNotFoundError as error:
        if allow_missing:
            return []
        raise StateConflictError(
            f"Output directory is unavailable: {output}. Check outputDir and mounts before retrying."
        ) from error
    except OSError as error:
        raise StateConflictError(f"Cannot inspect output directory {output}: {error}") from error


def _stage_map(start_stage: str, completed_before: bool = False) -> dict[str, dict]:
    if start_stage not in STAGES:
        raise StateValidationError(f"Unsupported stage: {start_stage}")
    start_index = STAGES.index(start_stage)
    stages: dict[str, dict] = {}
    for index, stage in enumerate(STAGES):
        if index < start_index:
            status = "completed" if completed_before else "skipped"
        elif index == start_index:
            status = "queued"
        else:
            status = "pending"
        stages[stage] = {
            "status": status,
            "started_at": None,
            "completed_at": utcnow() if status == "completed" else None,
            "failed_at": None,
            "inferred": bool(status == "completed" and completed_before),
            "metadata": {},
        }
    return stages


def collect_versions(cfg=None) -> dict[str, str | None]:
    """Collect version metadata without making it a prerequisite for processing."""
    bfq_version = None
    if cfg is not None:
        bfq_version = cfg.static.version.get("pipeline")
    if not bfq_version:
        try:
            bfq_version = metadata.version("bcl2fastq-pipeline")
        except metadata.PackageNotFoundError:
            bfq_version = None

    workflows_revision = None
    workflows_dir = Path("/opt/gcf-workflows")
    if workflows_dir.exists():
        try:
            workflows_revision = subprocess.check_output(
                ["git", "-C", str(workflows_dir), "rev-parse", "HEAD"],
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        except (OSError, subprocess.CalledProcessError):
            workflows_revision = None

    return {"bfq": bfq_version, "gcf_workflows": workflows_revision}


def new_state(  # noqa: PLR0913
    run_id: str,
    source_path: Path | str,
    output_path: Path | str,
    *,
    origin: str,
    start_stage: str,
    preparing: bool = False,
    cfg=None,
) -> dict:
    run_id = _safe_run_id(run_id)
    if origin not in ORIGINS:
        raise StateValidationError(f"Unsupported origin: {origin}")
    inferred = (
        origin in {"restored_legacy_fastq", "legacy_rerun"} and start_stage != "demultiplexing"
    )
    now = utcnow()
    state = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "source_path": str(Path(source_path)),
        "output_path": str(Path(output_path)),
        "origin": origin,
        "status": "preparing" if preparing else "queued",
        "current_stage": start_stage,
        "attempt": 0,
        "created_at": now,
        "updated_at": now,
        "started_at": None,
        "completed_at": None,
        "failed_at": None,
        "projects": [],
        "stages": _stage_map(start_stage, completed_before=inferred),
        "versions": collect_versions(cfg),
        "demultiplexing": {"tool": None, "version": None},
        "last_error": None,
        "notification": {
            "failure_signature": None,
            "notified": False,
            "notified_at": None,
        },
        "delivery_notifications": [],
        "restart_request": None,
        "archive": {"status": "active", "archived_at": None},
        "attempts": [],
    }
    validate_state(state)
    return state


def validate_state(state: dict) -> None:
    if not isinstance(state, dict):
        raise StateValidationError("State must be a JSON object")
    required = {
        "schema_version",
        "run_id",
        "source_path",
        "output_path",
        "origin",
        "status",
        "current_stage",
        "attempt",
        "created_at",
        "updated_at",
        "started_at",
        "completed_at",
        "failed_at",
        "projects",
        "stages",
        "versions",
        "demultiplexing",
        "last_error",
        "notification",
        "restart_request",
        "archive",
        "attempts",
    }
    missing = sorted(required - state.keys())
    if missing:
        raise StateValidationError(f"Missing state fields: {', '.join(missing)}")
    if state["schema_version"] != SCHEMA_VERSION:
        raise StateValidationError(
            f"Unsupported state schema {state['schema_version']!r}; expected {SCHEMA_VERSION}"
        )
    _safe_run_id(state["run_id"])
    if state["origin"] not in ORIGINS:
        raise StateValidationError(f"Invalid origin: {state['origin']!r}")
    if state["status"] not in RUN_STATUSES:
        raise StateValidationError(f"Invalid run status: {state['status']!r}")
    if state["current_stage"] not in STAGES:
        raise StateValidationError(f"Invalid current stage: {state['current_stage']!r}")
    if not isinstance(state["attempt"], int) or state["attempt"] < 0:
        raise StateValidationError("attempt must be a non-negative integer")
    if not isinstance(state["projects"], list):
        raise StateValidationError("projects must be a list")
    if not isinstance(state["attempts"], list):
        raise StateValidationError("attempts must be a list")
    if set(state["stages"]) != set(STAGES):
        raise StateValidationError("State must contain exactly the supported processing stages")
    for stage, detail in state["stages"].items():
        if not isinstance(detail, dict) or detail.get("status") not in STAGE_STATUSES:
            raise StateValidationError(f"Invalid stage state for {stage}")
    # Optional in v1: loading an older document must not invent mail to send.
    if "delivery_notifications" in state:
        _validate_delivery_notifications(state["delivery_notifications"])


def _validate_delivery_notifications(notifications) -> None:
    if not isinstance(notifications, list):
        raise StateValidationError("delivery_notifications must be a list")
    identifiers = set()
    required = {
        "id",
        "kind",
        "stage",
        "attempt",
        "status",
        "created_at",
        "sent_at",
        "last_error",
        "payload",
        "attempts",
    }
    for entry in notifications:
        if not isinstance(entry, dict) or required - entry.keys():
            raise StateValidationError("Invalid delivery notification fields")
        if (
            not isinstance(entry["kind"], str)
            or not entry["kind"]
            or ":" in entry["kind"]
            or not isinstance(entry["attempt"], int)
            or isinstance(entry["attempt"], bool)
            or entry["attempt"] < 1
            or entry["id"] != f"{entry['kind']}:{entry['attempt']}"
        ):
            raise StateValidationError("Invalid delivery notification identity")
        if entry["id"] in identifiers:
            raise StateValidationError(f"Duplicate delivery notification: {entry['id']}")
        identifiers.add(entry["id"])
        if (
            entry["stage"] not in STAGES
            or not isinstance(entry["status"], str)
            or entry["status"] not in DELIVERY_STATUSES
        ):
            raise StateValidationError("Invalid delivery notification stage or status")
        if not isinstance(entry["payload"], dict) or not isinstance(entry["created_at"], str):
            raise StateValidationError("Invalid delivery notification payload or timestamp")
        for field in ("sent_at", "last_error", "invalidated_at", "invalidation_reason"):
            if entry.get(field) is not None and not isinstance(entry[field], str):
                raise StateValidationError(f"Invalid delivery notification {field}")
        attempts = entry["attempts"]
        if not isinstance(attempts, list):
            raise StateValidationError("Delivery notification attempts must be a list")
        for number, attempt in enumerate(attempts, start=1):
            if (
                not isinstance(attempt, dict)
                or attempt.get("attempt") != number
                or not isinstance(attempt.get("status"), str)
                or attempt["status"] not in DELIVERY_ATTEMPT_STATUSES
                or not isinstance(attempt.get("started_at"), str)
                or "completed_at" not in attempt
                or "error" not in attempt
            ):
                raise StateValidationError("Invalid delivery attempt")
            for field in ("completed_at", "error"):
                if attempt[field] is not None and not isinstance(attempt[field], str):
                    raise StateValidationError(f"Invalid delivery attempt {field}")


def _append_delivery_notification(state: dict, stage: str, notification: dict | None) -> None:
    if notification is None:
        return
    if (
        not isinstance(notification, dict)
        or not isinstance(notification.get("kind"), str)
        or not notification["kind"]
        or not isinstance(notification.get("payload"), dict)
    ):
        raise StateValidationError("Notification requires a kind and payload object")
    state.setdefault("delivery_notifications", []).append(
        {
            "id": f"{notification['kind']}:{state['attempt']}",
            "kind": notification["kind"],
            "stage": stage,
            "attempt": state["attempt"],
            "status": "pending",
            "created_at": utcnow(),
            "sent_at": None,
            "last_error": None,
            "payload": copy.deepcopy(notification["payload"]),
            "attempts": [],
        }
    )


def _invalidate_delivery_notifications(
    state: dict, from_stage: str, reason: str, *, unsent_only: bool = False
) -> None:
    boundary = STAGES.index(from_stage)
    now = utcnow()
    for entry in state.get("delivery_notifications", []):
        if (
            STAGES.index(entry["stage"]) < boundary
            or entry["status"] == "superseded"
            or (unsent_only and entry["status"] == "sent")
        ):
            continue
        entry["status"] = "superseded"
        entry["invalidated_at"] = now
        entry["invalidation_reason"] = reason


def validate_restart_boundary(state: dict, start_stage: str) -> None:
    """Require preserved upstream stages to be trustworthy before cleanup begins."""
    if start_stage not in STAGES:
        raise StateValidationError(f"Unsupported stage: {start_stage}")
    start_index = STAGES.index(start_stage)
    invalid = [
        stage
        for stage in STAGES[:start_index]
        if state["stages"][stage]["status"] not in {"completed", "skipped"}
    ]
    if invalid:
        raise StateConflictError(
            f"Cannot restart from {start_stage}; upstream stages are not complete: "
            + ", ".join(invalid)
        )


class FlowcellStateStore:
    """Read and mutate flowcell state rooted at the configured manager directory."""

    def __init__(self, manager_dir: Path | str):
        self.manager_dir = Path(manager_dir)
        self.states_dir = self.manager_dir / "states"
        self.locks_dir = self.manager_dir / "locks"

    def ensure_writable(self) -> None:
        """Create state infrastructure and prove it is writable before processing."""
        self.states_dir.mkdir(parents=True, exist_ok=True)
        self.locks_dir.mkdir(parents=True, exist_ok=True)
        probe = self.states_dir / f".write-test-{os.getpid()}"
        try:
            with probe.open("w", encoding="utf-8") as handle:
                handle.write("ok")
                handle.flush()
                os.fsync(handle.fileno())
        finally:
            probe.unlink(missing_ok=True)

    def state_path(self, run_id: str) -> Path:
        return self.states_dir / f"{_safe_run_id(run_id)}.json"

    def _lock_path(self, run_id: str) -> Path:
        return self.locks_dir / f"{_safe_run_id(run_id)}.lock"

    def _execution_path(self, run_id: str) -> Path:
        return self.locks_dir / f"{_safe_run_id(run_id)}.execution.lock"

    @contextmanager
    def lock(self, run_id: str) -> Iterator[None]:
        self.ensure_writable()
        with self._lock_path(run_id).open("a+") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    @contextmanager
    def execution_lease(self, run_id: str, *, blocking: bool = False) -> Iterator[None]:
        self.ensure_writable()
        flags = fcntl.LOCK_EX
        if not blocking:
            flags |= fcntl.LOCK_NB
        with self._execution_path(run_id).open("a+") as handle:
            try:
                fcntl.flock(handle.fileno(), flags)
            except BlockingIOError as error:
                raise ExecutionLeaseError(f"Run {run_id} is already active") from error
            try:
                handle.seek(0)
                handle.truncate()
                handle.write(f"pid={os.getpid()} host={socket.gethostname()} acquired={utcnow()}\n")
                handle.flush()
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def execution_active(self, run_id: str) -> bool:
        try:
            with self.execution_lease(run_id):
                return False
        except ExecutionLeaseError:
            return True

    def exists(self, run_id: str) -> bool:
        return self.state_path(run_id).exists()

    def read(self, run_id: str) -> dict:
        path = self.state_path(run_id)
        try:
            state = json.loads(path.read_text(encoding="utf-8"))
        except FileNotFoundError:
            raise StateError(f"No state exists for {run_id}") from None
        except json.JSONDecodeError as error:
            raise StateValidationError(f"Invalid JSON in {path}: {error}") from error
        validate_state(state)
        if state["run_id"] != _safe_run_id(run_id):
            raise StateValidationError(f"State run_id does not match filename: {path}")
        return state

    def _atomic_write_unlocked(self, state: dict) -> None:
        validate_state(state)
        self.ensure_writable()
        path = self.state_path(state["run_id"])
        tmp = self.states_dir / f".{state['run_id']}.{os.getpid()}.tmp"
        payload = json.dumps(state, indent=2, sort_keys=True) + "\n"
        try:
            with tmp.open("w", encoding="utf-8") as handle:
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(tmp, path)
            dir_fd = os.open(self.states_dir, os.O_RDONLY)
            try:
                os.fsync(dir_fd)
            finally:
                os.close(dir_fd)
        finally:
            tmp.unlink(missing_ok=True)

    def create(self, state: dict) -> dict:
        validate_state(state)
        run_id = state["run_id"]
        with self.lock(run_id):
            if self.state_path(run_id).exists():
                raise StateConflictError(f"State already exists for {run_id}")
            self._atomic_write_unlocked(copy.deepcopy(state))
        return state

    def write(self, state: dict) -> dict:
        state = copy.deepcopy(state)
        state["updated_at"] = utcnow()
        with self.lock(state["run_id"]):
            self._atomic_write_unlocked(state)
        return state

    def mutate(self, run_id: str, mutator) -> dict:
        with self.lock(run_id):
            state = self.read(run_id)
            updated = mutator(copy.deepcopy(state))
            updated["updated_at"] = utcnow()
            validate_state(updated)
            self._atomic_write_unlocked(updated)
            return updated

    def list_states(self) -> list[dict]:
        self.ensure_writable()
        states = []
        for path in sorted(self.states_dir.glob("*.json")):
            state = json.loads(path.read_text(encoding="utf-8"))
            validate_state(state)
            states.append(state)
        return states

    def set_preparing(  # noqa: PLR0913
        self,
        run_id: str,
        start_stage: str,
        *,
        reason: str | None,
        refresh_inputs: bool,
        hostname: str | None = None,
        output_path: Path | None = None,
    ) -> dict:
        if start_stage not in STAGES:
            raise StateValidationError(f"Unsupported stage: {start_stage}")

        def update(state: dict) -> dict:
            if output_path is not None:
                state["output_path"] = str(output_path)
            _invalidate_delivery_notifications(state, start_stage, f"Restart from {start_stage}")
            state["status"] = "preparing"
            state["current_stage"] = start_stage
            state["restart_request"] = {
                "from": start_stage,
                "reason": reason,
                "hostname": hostname or socket.gethostname(),
                "requested_at": utcnow(),
                "refresh_inputs": bool(refresh_inputs),
            }
            state["completed_at"] = None
            state["failed_at"] = None
            state["last_error"] = None
            return state

        return self.mutate(run_id, update)

    def queue(self, run_id: str, start_stage: str) -> dict:
        if start_stage not in STAGES:
            raise StateValidationError(f"Unsupported stage: {start_stage}")

        def update(state: dict) -> dict:
            validate_restart_boundary(state, start_stage)
            _invalidate_delivery_notifications(state, start_stage, f"Restart from {start_stage}")
            start_index = STAGES.index(start_stage)
            for index, stage in enumerate(STAGES):
                detail = state["stages"][stage]
                if index < start_index:
                    continue
                detail["status"] = "queued" if index == start_index else "pending"
                detail["started_at"] = None
                detail["completed_at"] = None
                detail["failed_at"] = None
                detail["metadata"] = {}
                detail["inferred"] = False
            state["status"] = "queued"
            state["current_stage"] = start_stage
            state["archive"] = {"status": "active", "archived_at": None}
            return state

        return self.mutate(run_id, update)

    def begin_attempt(self, run_id: str, *, cfg=None) -> dict:
        def update(state: dict) -> dict:
            if state["status"] != "queued":
                raise StateConflictError(
                    f"Run {run_id} is {state['status']!r}, not queued for execution"
                )
            stage = state["current_stage"]
            if state["stages"][stage]["status"] != "queued":
                raise StateConflictError(f"Current stage {stage} is not queued")
            if cfg is not None:
                state["output_path"] = str(resolve_output_path(state["output_path"], run_id, cfg))
            state["attempt"] += 1
            now = utcnow()
            state["status"] = "running"
            state["started_at"] = state["started_at"] or now
            state["failed_at"] = None
            state["versions"].update(collect_versions(cfg))
            state["stages"][stage]["status"] = "running"
            state["stages"][stage]["started_at"] = now
            request = state.get("restart_request") or {}
            state["attempts"].append(
                {
                    "attempt": state["attempt"],
                    "restart_from": stage,
                    "reason": request.get("reason"),
                    "requester_host": request.get("hostname"),
                    "started_at": now,
                    "completed_at": None,
                    "failed_at": None,
                    "outcome": "running",
                    "versions": copy.deepcopy(state["versions"]),
                    "failure": None,
                }
            )
            return state

        return self.mutate(run_id, update)

    def start_stage(self, run_id: str, stage: str) -> dict:
        if stage not in STAGES:
            raise StateValidationError(f"Unsupported stage: {stage}")

        def update(state: dict) -> dict:
            if state["status"] != "running":
                raise StateConflictError(f"Run {run_id} is not running")
            detail = state["stages"][stage]
            if detail["status"] != "queued":
                raise StateConflictError(f"Stage {stage} is not queued")
            state["current_stage"] = stage
            detail["status"] = "running"
            detail["started_at"] = utcnow()
            return state

        return self.mutate(run_id, update)

    def record_input_preflight(self, run_id: str, stage: str, report: dict) -> dict:
        """Attach a complete hash-bound report to its existing restart boundary."""
        if stage not in {"demultiplexing", "analysis"}:
            raise StateValidationError(f"Input preflight is not applicable to {stage}")

        def update(state: dict) -> dict:
            state["stages"][stage]["metadata"]["input_preflight"] = copy.deepcopy(report)
            if state["attempts"] and state["attempts"][-1]["outcome"] == "running":
                state["attempts"][-1].setdefault("input_preflight", {})[stage] = copy.deepcopy(
                    report
                )
            return state

        return self.mutate(run_id, update)

    def complete_stage(
        self,
        run_id: str,
        stage: str,
        metadata_: dict | None = None,
        *,
        notification: dict | None = None,
    ) -> dict:
        """Atomically record completed work and its optional delivery intent."""
        if stage not in STAGES:
            raise StateValidationError(f"Unsupported stage: {stage}")

        def update(state: dict) -> dict:
            detail = state["stages"][stage]
            if state["status"] != "running" or detail["status"] != "running":
                raise StateConflictError(f"Stage {stage} is not running")
            detail["status"] = "completed"
            detail["completed_at"] = utcnow()
            if metadata_:
                detail["metadata"].update(metadata_)
            if stage == "demultiplexing" and metadata_:
                state["demultiplexing"]["tool"] = metadata_.get("tool")
                state["demultiplexing"]["version"] = metadata_.get("version")
            index = STAGES.index(stage)
            if index + 1 < len(STAGES):
                next_stage = STAGES[index + 1]
                state["current_stage"] = next_stage
                state["stages"][next_stage]["status"] = "queued"
            _append_delivery_notification(state, stage, notification)
            return state

        return self.mutate(run_id, update)

    def fail_stage(  # noqa: PLR0913
        self,
        run_id: str,
        stage: str,
        *,
        summary: str,
        report_path: Path | str | None,
        interrupted: bool = False,
        failure_signature: str | None = None,
    ) -> dict:
        if stage not in STAGES:
            raise StateValidationError(f"Unsupported stage: {stage}")

        def update(state: dict) -> dict:
            now = utcnow()
            signature_source = f"{stage}\0{summary}".encode()
            signature = failure_signature or hashlib.sha256(signature_source).hexdigest()
            # Only a just-committed completion may belong to this failure.
            # Queued/preparing runs still refer to the previous successful attempt.
            failed_outcomes = {"running"}
            if state["status"] == "completed":
                failed_outcomes.add("completed")
            # A completion write may have replaced the state file before its
            # durability check failed. Cancel its intent with the failed stage.
            _invalidate_delivery_notifications(state, stage, f"Processing failed at {stage}")
            state["status"] = "interrupted" if interrupted else "failed"
            state["current_stage"] = stage
            state["failed_at"] = now
            state["completed_at"] = None
            detail = state["stages"][stage]
            detail["status"] = "failed"
            detail["failed_at"] = now
            detail["completed_at"] = None
            state["last_error"] = {
                "summary": summary,
                "report_path": str(report_path) if report_path else None,
                "failure_signature": signature,
            }
            if state["notification"].get("failure_signature") != signature:
                state["notification"] = {
                    "failure_signature": signature,
                    "notified": False,
                    "notified_at": None,
                }
            if (
                state["attempts"]
                and state["attempts"][-1].get("attempt") == state["attempt"]
                and state["attempts"][-1].get("outcome") in failed_outcomes
            ):
                attempt = state["attempts"][-1]
                attempt["outcome"] = "interrupted" if interrupted else "failed"
                attempt["failed_at"] = now
                attempt["completed_at"] = None
                attempt["failure"] = copy.deepcopy(state["last_error"])
            return state

        return self.mutate(run_id, update)

    def deliver_failure_notification(self, run_id: str, signature: str, send) -> bool:
        """Make at most one delivery attempt per failure, including across crashes.

        Persist the claim before SMTP and serialize against restart/success updates.
        A crash after claiming may lose mail; automatic retries could duplicate mail
        already accepted by the relay. The on-server report remains authoritative.
        """
        with self.lock(run_id):
            state = self.read(run_id)
            notification = state["notification"]
            if (
                notification.get("failure_signature") != signature
                or notification.get("notified")
                or notification.get("attempted_at")
            ):
                return False
            notification["attempted_at"] = utcnow()
            state["updated_at"] = utcnow()
            self._atomic_write_unlocked(state)
            try:
                send()
            except Exception as error:
                notification["delivery_error"] = f"{type(error).__name__}: {error}"
                state["updated_at"] = utcnow()
                self._atomic_write_unlocked(state)
                raise
            notification["notified"] = True
            notification["notified_at"] = utcnow()
            notification["delivery_error"] = None
            state["updated_at"] = utcnow()
            self._atomic_write_unlocked(state)
            return True

    def deliver_notification(  # noqa: PLR0913
        self,
        run_id: str,
        notification_id: str,
        send,
        *,
        retry: bool = False,
        retry_uncertain: bool = False,
    ) -> bool:
        """Deliver one intent without modifying processing state.

        The caller must hold the run's execution lease. The state lock spans SMTP
        to serialize delivery against cleanup and other retries. ``send(entry)``
        must not call back into state mutations, which would deadlock this lock.

        Pending entries receive one automatic attempt. Failed entries need an
        explicit retry; uncertain outcomes also need duplicate-risk acknowledgement.
        Persisting the claim before SMTP prevents automatic resend after a crash.
        A crash or write failure after SMTP acceptance cannot establish delivery:
        the persisted claim is recovered as uncertain, never assumed sent.

        Return True only after success is persisted. Callback exceptions record a
        failed/uncertain delivery and return False; state I/O failures propagate.
        """
        with self.lock(run_id):
            state = self.read(run_id)
            entry = next(
                (
                    entry
                    for entry in state.get("delivery_notifications", [])
                    if entry["id"] == notification_id
                ),
                None,
            )
            if entry is None:
                raise StateConflictError(f"No notification {notification_id!r} for {run_id}")
            if entry["status"] in {"sent", "superseded"}:
                return False
            if state["stages"][entry["stage"]]["status"] != "completed":
                # Defense for older/inconsistent records: completion mail must
                # never advertise work whose producing stage is not complete.
                _invalidate_delivery_notifications(
                    state, entry["stage"], f"Producing stage {entry['stage']} is not completed"
                )
                state["updated_at"] = utcnow()
                self._atomic_write_unlocked(state)
                return False
            if entry["status"] == "sending":
                message = "Previous delivery was interrupted; SMTP acceptance is unknown"
                self._finish_delivery_unlocked(state, entry, "uncertain", message)
            status = entry["status"]
            if (status == "failed" and not retry) or (
                status == "uncertain" and not (retry and retry_uncertain)
            ):
                return False
            now = utcnow()
            entry["status"] = "sending"
            entry["last_error"] = None
            entry["attempts"].append(
                {
                    "attempt": len(entry["attempts"]) + 1,
                    "status": "sending",
                    "started_at": now,
                    "completed_at": None,
                    "error": None,
                }
            )
            state["updated_at"] = now
            self._atomic_write_unlocked(state)
            try:
                send(copy.deepcopy(entry))
            except Exception as error:
                status = "uncertain" if isinstance(error, DeliveryUncertainError) else "failed"
                self._finish_delivery_unlocked(
                    state, entry, status, f"{type(error).__name__}: {error}"
                )
                return False
            self._finish_delivery_unlocked(state, entry, "sent", None)
            return True

    def _finish_delivery_unlocked(
        self, state: dict, entry: dict, status: str, error: str | None
    ) -> None:
        now = utcnow()
        entry["status"] = status
        entry["last_error"] = error
        if status == "sent":
            entry["sent_at"] = now
        if entry["attempts"]:
            entry["attempts"][-1].update(status=status, completed_at=now, error=error)
        state["updated_at"] = now
        self._atomic_write_unlocked(state)

    def complete_run(
        self, run_id: str, projects: list[str], *, notification: dict | None = None
    ) -> dict:
        """Atomically record finalization and its optional delivery intent."""

        def update(state: dict) -> dict:
            stage = state["stages"]["finalization"]
            if state["status"] != "running" or stage["status"] != "running":
                raise StateConflictError("Finalization must be running before run completion")
            now = utcnow()
            stage["status"] = "completed"
            stage["completed_at"] = now
            state["status"] = "completed"
            state["current_stage"] = "finalization"
            state["completed_at"] = now
            state["failed_at"] = None
            state["last_error"] = None
            state["notification"] = {
                "failure_signature": None,
                "notified": False,
                "notified_at": None,
            }
            state["projects"] = sorted(set(projects))
            if state["attempts"] and state["attempts"][-1].get("outcome") == "running":
                state["attempts"][-1]["outcome"] = "completed"
                state["attempts"][-1]["completed_at"] = now
            _append_delivery_notification(state, "finalization", notification)
            return state

        return self.mutate(run_id, update)

    def mark_archived(self, run_id: str) -> dict:
        def update(state: dict) -> dict:
            now = utcnow()
            _invalidate_delivery_notifications(
                state, "demultiplexing", "Outputs archived", unsent_only=True
            )
            state["archive"] = {"status": "archived", "archived_at": now}
            if state["status"] == "completed":
                state["status"] = "archived"
            return state

        return self.mutate(run_id, update)

    def invalidate_notifications(
        self,
        run_id: str,
        from_stage: str = "demultiplexing",
        *,
        reason: str,
        unsent_only: bool = False,
        output_path: Path | None = None,
    ) -> dict:
        """Persist invalidation before removing outputs, under an execution lease."""
        if from_stage not in STAGES:
            raise StateValidationError(f"Unsupported stage: {from_stage}")

        def update(state: dict) -> dict:
            if output_path is not None:
                state["output_path"] = str(output_path)
            _invalidate_delivery_notifications(state, from_stage, reason, unsent_only=unsent_only)
            return state

        return self.mutate(run_id, update)

    def recover_interrupted(self, run_id: str) -> dict:
        state = self.read(run_id)
        if state["status"] != "running":
            return state
        if self.execution_active(run_id):
            return state
        stage = state["current_stage"]
        return self.fail_stage(
            run_id,
            stage,
            summary=f"Interrupted while {stage} was running; explicit rerun required",
            report_path=None,
            interrupted=True,
        )


def find_fastqs(output_path: Path | str) -> list[Path]:
    output_path = Path(output_path)
    fastqs = list(output_path.glob("GCF-*/*.fastq.gz"))
    fastqs += list(output_path.glob("GCF-*/*/*.fastq.gz"))
    fastqs += list(output_path.glob("Undetermined*.fastq.gz"))
    return sorted(set(fastqs))


def validate_restored_fastqs(output_path: Path | str) -> tuple[bool, str]:  # noqa: PLR0911
    """Recognize an explicitly restored legacy FASTQ tree conservatively."""
    output_path = Path(output_path)
    sample_sheet = output_path / "SampleSheet.csv"
    submission = output_path / "Sample-Submission-Form.xlsx"
    if not sample_sheet.is_file() or not submission.is_file():
        return False, "SampleSheet.csv and Sample-Submission-Form.xlsx are required"

    try:
        options, _ = parse_custom_options(sample_sheet)
    except (OSError, UnicodeError) as error:
        return False, f"Sample sheet is unreadable: {error}"
    if not options or not options.get("Libprep", "").strip():
        return False, "Sample sheet does not contain a usable [CustomOptions] Libprep"

    fastqs = find_fastqs(output_path)
    project_fastqs = [
        path for path in fastqs if any(part.startswith("GCF-") for part in path.parts)
    ]
    if not project_fastqs:
        return False, "No FASTQs were found in recognized GCF project directories"

    for fastq in fastqs:
        try:
            if fastq.stat().st_size == 0:
                return False, f"FASTQ is empty: {fastq}"
            with gzip.open(fastq, "rb") as handle:
                if not handle.read(1):
                    return False, f"FASTQ has no uncompressed content: {fastq}"
        except (OSError, EOFError, gzip.BadGzipFile) as error:
            return False, f"FASTQ is not a readable gzip stream: {fastq}: {error}"
    return True, "recognized restored FASTQs"


def cleanup_plan(output_path: Path | str, from_stage: str) -> list[Path]:
    """Return existing paths invalidated by a restart boundary."""
    output = Path(output_path)
    if from_stage not in STAGES:
        raise StateValidationError(f"Unsupported restart stage: {from_stage}")
    children = output_entries(output, allow_missing=from_stage == "demultiplexing")

    targets: set[Path] = set()
    if from_stage == "demultiplexing":
        targets.update(path for path in children if path.name not in INPUT_FILES)
    else:
        # Finalization products are always downstream of analysis/reporting.
        targets.update(
            path
            for path in children
            if any(
                path.match(pattern) for pattern in ("*.7za", "encryption.*", "md5sum_*_archive.txt")
            )
        )

        if from_stage in {"analysis", "reporting"}:
            targets.update(output.glob("bcl2fastq.ini"))
            stats = output / "Stats"
            if stats.exists():
                targets.update(stats.glob("interop_summary.csv"))
                targets.update(stats.glob("interop_index-summary.csv"))
                targets.update(stats.glob("sequencer_stats_*.html"))
                targets.update(stats.glob(".multiqc_config.yaml"))

        if from_stage == "analysis":
            # full_align creates these inputs/results; reporting consumes them.
            targets.update(output.glob("multiqc_*.html"))
            targets.update(output.glob("all_samples_web_summary_*.html"))
            targets.update(output.glob(".multiqc_config_*.yaml"))
            targets.update(output.glob("QC_*"))
            targets.update(output.glob("GCF-*_samplesheet.tsv"))
            targets.update(output.glob("configmaker-analysis-*.json"))

    # Remove descendants if an ancestor is already scheduled, keeping dry-run output concise.
    ordered = sorted(targets, key=lambda path: (len(path.parts), str(path)))
    collapsed: list[Path] = []
    for path in ordered:
        if any(parent == path or parent in path.parents for parent in collapsed):
            continue
        collapsed.append(path)
    return collapsed


def apply_cleanup(paths: list[Path]) -> None:
    for path in paths:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink(missing_ok=True)


def locate_source_run(run_id: str, cfg) -> Path | None:
    run_id = _safe_run_id(run_id)
    for base in (cfg.static.paths.nova_base_dir, cfg.static.paths.ekista_base_dir):
        candidate = Path(base) / run_id
        if candidate.exists():
            return candidate
    return None


def refresh_run_inputs(source_path: Path | str, output_path: Path | str) -> tuple[Path, Path]:
    """Explicitly replace output-side run inputs from the instrument source."""
    # Local import avoids a module cycle: preflight's operator errors inherit StateError.
    from bcl2fastq_pipeline.preflight import copy_run_inputs, require_valid_inputs  # noqa: PLC0415

    selection, _result = require_valid_inputs(source_path, output_path, refresh=True)
    copied = copy_run_inputs(selection, output_path)
    return copied.sample_sheet, copied.submission_form
