"""Plan, publish and audit explicit SampleSheet index toggles.

Callers hold the flowcell execution lease for apply/recover. A pending history
entry is durable before the atomic sheet replacement. Recovery only reconciles
that entry with the current file; it never applies a transformation itself.
"""

from __future__ import annotations

import copy
import hashlib
import os
import stat
import tempfile
import uuid

from dataclasses import dataclass
from pathlib import Path

from bcl2fastq_pipeline import index_sequences, preflight
from bcl2fastq_pipeline.state import StateConflictError, StateError, utcnow


def _sha256(content):
    return hashlib.sha256(content).hexdigest()


def _unavailable(detail):
    return {index: {"status": "unavailable", "detail": detail} for index in ("index1", "index2")}


def _orientation(content, source_path):
    source_sheet = preflight.select_run_inputs(source_path, source_path, refresh=True).sample_sheet
    report = {"source_sheet": str(source_sheet.absolute()), "source_sha256": None}
    try:
        reference = source_sheet.read_bytes()
        report["source_sha256"] = _sha256(reference)
        report["per_index"] = index_sequences.orientation(content, reference)
    except (OSError, StateError, UnicodeError) as error:
        report["per_index"] = _unavailable(f"Cannot compare against source: {error}")
    return {**report, "effective_sha256": _sha256(content)}


def current_orientation(state):
    """Read current files, not transformation parity or a cached orientation."""
    path = Path(state["output_path"]) / "SampleSheet.csv"
    try:
        content = path.read_bytes()
    except OSError as error:
        return {
            "effective_sheet": str(path),
            "effective_sha256": None,
            "source_directory": state["source_path"],
            "per_index": _unavailable(f"Cannot read effective SampleSheet: {error}"),
        }
    return {"effective_sheet": str(path), **_orientation(content, state["source_path"])}


@dataclass(frozen=True)
class CorrectionPlan:
    operation_id: str
    selected_path: Path
    backup_path: Path
    indexes: tuple[str, ...]
    before: bytes
    after: bytes
    before_sha256: str
    after_sha256: str
    counts: dict
    selected_counts: dict
    orientation: dict
    reason: str | None


def plan(store, run_id, selection, source_path, indexes, *, reason=None):  # noqa: PLR0913
    """Validate a complete toggle without writing state, backups or inputs."""
    # Reuse the store's run-ID validation before constructing a history path.
    store.state_path(run_id)
    selected = Path(selection.sample_sheet).absolute()
    try:
        before = selected.read_bytes()
    except OSError as error:
        raise StateError(
            f"Cannot read SampleSheet for index toggle: {selected}: {error}"
        ) from error
    indexes = tuple(index for index in ("index1", "index2") if index in indexes)
    if not indexes:
        raise StateError("An index toggle must select index1 or index2")
    result = index_sequences.transform(before, indexes)
    operation_id = uuid.uuid4().hex
    backup = (
        store.manager_dir.absolute() / "input-history" / run_id / f"{operation_id}.SampleSheet.csv"
    )
    return CorrectionPlan(
        operation_id=operation_id,
        selected_path=selected,
        backup_path=backup,
        indexes=indexes,
        before=before,
        after=result.content,
        before_sha256=_sha256(before),
        after_sha256=_sha256(result.content),
        counts=result.counts,
        selected_counts=result.selected_counts,
        orientation=_orientation(result.content, source_path),
        reason=reason,
    )


def verify_plan(preview, current):
    """Refuse changed inputs/reference after confirmation; keep the preview ID."""
    for attribute in (
        "selected_path",
        "indexes",
        "before",
        "after",
        "counts",
        "selected_counts",
        "orientation",
        "reason",
    ):
        if getattr(preview, attribute) != getattr(current, attribute):
            raise StateConflictError(
                "Index toggle inputs or source comparison changed after preview; inspect and retry"
            )


def print_plan(correction):
    print(f"Index toggle input: {correction.selected_path}")
    for index in correction.indexes:
        print(
            f"  {index}: reverse-complement {correction.selected_counts[index]} nonempty rows "
            f"({correction.counts[index]} values change)"
        )
        result = correction.orientation["per_index"][index]
        print(f"  Resulting {index} orientation: {result['status']} ({result.get('detail', '')})")
    print(f"Orientation reference: {correction.orientation['source_sheet']}")
    print(f"Input backup: {correction.backup_path}")
    print(f"SampleSheet SHA256: {correction.before_sha256} -> {correction.after_sha256}")


def _sync_directory(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _history_directory(store, run_id, output):
    parent = store.manager_dir.absolute()
    target = parent / "input-history" / run_id
    if target.resolve().is_relative_to(output.resolve()):
        raise StateConflictError("Input history must live outside the flowcell output cleanup tree")
    for directory in (parent / "input-history", target):
        if directory.is_symlink():
            raise StateConflictError(f"Input history directory must not be a symlink: {directory}")
        directory.mkdir(exist_ok=True)
        _sync_directory(directory.parent)
    return target


def _write_backup(path, content):
    with path.open("xb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    _sync_directory(path.parent)


def _replace_sheet(target, content):
    # Atomic replacement also detaches an output symlink/hardlink from its source.
    mode = stat.S_IMODE(target.stat().st_mode)
    fd, name = tempfile.mkstemp(prefix=".SampleSheet.index-", suffix=".tmp", dir=target.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fchmod(handle.fileno(), mode)
            os.fsync(handle.fileno())
        os.replace(temporary, target)
        _sync_directory(target.parent)
    finally:
        temporary.unlink(missing_ok=True)


def _finish(store, run_id, operation_id, status, **details):
    def update(state):
        for entry in state.get("index_corrections", []):
            if entry["id"] == operation_id:
                entry.update(status=status, completed_at=utcnow(), **details)
                if status == "applied":
                    state["restart_request"] = {
                        **(state.get("restart_request") or {}),
                        "index_correction_id": operation_id,
                    }
                return state
        raise StateConflictError(f"Missing index toggle history entry: {operation_id}")

    return store.mutate(run_id, update)


def recover(store, run_id):
    """Reconcile pending records without reapplying or undoing any toggle.

    A new explicit invocation receives a new operation ID and acts on current
    bytes. A plain rerun keeps the current bytes, including a committed toggle.
    """
    state = store.read(run_id)
    for entry in state.get("index_corrections", []):
        if entry["status"] != "pending":
            continue
        target = Path(state["output_path"]) / "SampleSheet.csv"
        try:
            actual = _sha256(target.read_bytes())
        except OSError:
            actual = None
        if actual == entry["after_sha256"]:
            outcome = "applied"
            detail = "Recovered history: effective sheet matches the requested result"
        elif actual == entry["before_sha256"]:
            outcome = "not_applied"
            detail = "Interrupted before sheet replacement; effective sheet is unchanged"
        else:
            outcome = "interrupted_unknown"
            detail = "Effective sheet is unavailable or has changed; prior toggle outcome cannot be established"
        state = _finish(
            store,
            run_id,
            entry["id"],
            outcome,
            recovered_at=utcnow(),
            observed_sha256=actual,
            recovery_detail=detail,
        )
    return state


def apply(store, run_id, correction):
    """Publish one prepared toggle, with its backup and write-ahead audit entry."""
    state = store.read(run_id)
    output = Path(state["output_path"])
    target = output / "SampleSheet.csv"
    for existing in state.get("index_corrections", []):
        if existing["id"] == correction.operation_id:
            entry = existing
            if existing["status"] == "pending":
                state = recover(store, run_id)
                entry = next(
                    item
                    for item in state["index_corrections"]
                    if item["id"] == correction.operation_id
                )
            if (
                entry["status"] == "applied"
                and _sha256(target.read_bytes()) == correction.after_sha256
            ):
                return state
            raise StateConflictError(
                "This toggle operation was interrupted; inspect and issue a new command"
            )
    if state["status"] != "preparing":
        raise StateConflictError("Index toggles require a run in preparing state")
    if any(entry["status"] == "pending" for entry in state.get("index_corrections", [])):
        raise StateConflictError("Recover pending index history before preparing another toggle")
    if target.read_bytes() != correction.before:
        raise StateConflictError(
            "Effective SampleSheet changed after index toggle validation; inspect and retry"
        )
    if output.resolve() == Path(state["source_path"]).resolve():
        raise StateConflictError(
            "Input and output flowcell directories must differ; instrument inputs are read-only"
        )
    history_dir = _history_directory(store, run_id, output)
    if correction.backup_path != history_dir / f"{correction.operation_id}.SampleSheet.csv":
        raise StateConflictError("Index backup destination does not match this run")
    _write_backup(correction.backup_path, correction.before)
    entry = {
        "id": correction.operation_id,
        "operation": "reverse_complement",
        "indexes": list(correction.indexes),
        "status": "pending",
        "requested_at": utcnow(),
        "completed_at": None,
        "reason": correction.reason,
        "attempt_before": state["attempt"],
        "selected_sheet": str(correction.selected_path),
        "effective_sheet": str(target),
        "backup_path": str(correction.backup_path),
        "before_sha256": correction.before_sha256,
        "after_sha256": correction.after_sha256,
        "changed_rows": correction.counts,
        "nonempty_rows": correction.selected_counts,
        "orientation": correction.orientation,
    }

    def record(current):
        current.setdefault("index_corrections", []).append(copy.deepcopy(entry))
        return current

    store.mutate(run_id, record)
    _replace_sheet(target, correction.after)
    return _finish(store, run_id, correction.operation_id, "applied")
