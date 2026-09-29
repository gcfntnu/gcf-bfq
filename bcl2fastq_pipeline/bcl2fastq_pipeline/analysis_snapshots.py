"""Retain the actual analysis working tree after successful delivery finalization.

Archives are staged beside their destination. The durable completion record is
committed before publication, allowing recovery after a crash without publishing
an unsuccessful attempt. Callers hold the flowcell execution lease throughout.
"""

import hashlib
import json
import logging
import os
import tarfile
import tempfile
import uuid

from pathlib import Path

from bcl2fastq_pipeline.state import PROVENANCE_DIR, ExecutionLeaseError

log = logging.getLogger(__name__)
DIRECTORY = PROVENANCE_DIR
MARKER = ".bfq-analysis.json"


def workdir_path(run_id, project):
    return Path(os.environ.get("TMPDIR", "/bfq-tmp")) / f"{project}_{run_id.split('_')[0]}"


def identify_workdir(workdir, run_id, project):
    """Mark ownership before execution; only completed analysis records this token."""
    identity = {"run_id": run_id, "project": project, "token": uuid.uuid4().hex}
    (workdir / MARKER).write_text(json.dumps(identity) + "\n", encoding="utf-8")
    return {**identity, "path": str(workdir.absolute())}


def _analysis_id(state):
    # Completion time survives downstream-only retries and changes on reanalysis.
    completed = state["stages"]["analysis"]["completed_at"]
    if not completed:
        raise RuntimeError("Cannot retain a snapshot without completed analysis")
    return hashlib.sha256(f"{state['run_id']}:{completed}".encode()).hexdigest()


def _sync_directory(directory):
    fd = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _archive_matches(path, analysis_id):
    if not path.is_file() or path.is_symlink():
        return False
    with tarfile.open(path, "r:gz") as archive:
        return archive.pax_headers.get("bfq.analysis_id") == analysis_id


def _source(state, project):
    record = state["stages"]["analysis"]["metadata"].get("workdirs", {}).get(project)
    path = Path(record["path"]) if record else workdir_path(state["run_id"], project)
    if not path.is_dir() or path.is_symlink():
        return None, f"original workdir is unavailable: {path}"
    marker = path / MARKER
    if record:
        if not marker.is_file() or json.loads(marker.read_text()) != {
            key: record[key] for key in ("run_id", "project", "token")
        }:
            return None, f"workdir belongs to a different analysis: {path}"
    elif marker.exists():
        return None, f"legacy analysis cannot be matched to marked workdir: {path}"
    # Older state has no ownership token. Only use its original conventional
    # workdir, never /opt or a reconstructed working tree.
    if not (path / "config.yaml").is_file() or not (path / "Snakefile").is_file():
        return None, f"original config.yaml or Snakefile is unavailable: {path}"
    if not record:
        log.warning("Retaining legacy workdir without an analysis ownership token: %s", path)
    return path, None


def _write_archive(workdir, staging, analysis_id):
    fd, name = tempfile.mkstemp(prefix="analysis-", suffix=".tar.gz", dir=staging)
    path = Path(name)
    try:
        with os.fdopen(fd, "wb") as handle:
            with tarfile.open(
                fileobj=handle,
                mode="w:gz",
                dereference=False,
                format=tarfile.PAX_FORMAT,
                pax_headers={"bfq.analysis_id": analysis_id},
            ) as archive:
                # Enumerating immediate children excludes only top-level data,
                # includes dotfiles, and keeps nested directories named data.
                for child in sorted(workdir.iterdir()):
                    if child.name != "data":
                        archive.add(child, arcname=child.name)
            handle.flush()
            os.fsync(handle.fileno())
        _sync_directory(staging)
        return path
    except BaseException:
        path.unlink(missing_ok=True)
        raise


def prepare(state, projects):
    """Build candidates; existing retained archives are never changed here.

    Return per-project results for the completion record. Missing historical
    sources are explicit unavailable results; write/permission errors fail the
    finalization rather than silently claiming a valid snapshot.
    """
    output = Path(state["output_path"])
    root = output / DIRECTORY
    staging = root / ".staging"
    for directory in (root, staging):
        if directory.is_symlink():
            raise RuntimeError(f"Snapshot directory must not be a symlink: {directory}")
        directory.mkdir(exist_ok=True)
    analysis_id = _analysis_id(state)
    results = {}
    created = []
    try:
        for project in sorted(projects):
            if Path(project).name != project or project in {".", ".."}:
                raise ValueError(f"Invalid snapshot project: {project!r}")
            target = root / f"{project}_analysis.tar.gz"
            retained = state.get("analysis_snapshots", {}).get(project)
            if (
                retained
                and retained["analysis_id"] == analysis_id
                and _archive_matches(target, analysis_id)
            ):
                results[project] = {"status": "available", **retained}
                continue
            source, reason = _source(state, project)
            if source is None:
                log.warning(
                    "Analysis snapshot unavailable for %s/%s: %s", state["run_id"], project, reason
                )
                results[project] = {"status": "unavailable", "reason": reason}
                continue
            candidate = _write_archive(source, staging, analysis_id)
            created.append(candidate)
            results[project] = {
                "status": "available",
                "analysis_id": analysis_id,
                "archive": str(target.relative_to(output)),
                "pending": str(candidate.relative_to(output)),
            }
        return results
    except BaseException:
        for path in created:
            path.unlink(missing_ok=True)
        raise


def recover(store, run_id):
    """Publish only state-committed candidates, then remove abandoned temporary files.

    Safe to repeat after any interruption, including between rename and state
    update. Publication failure leaves the committed candidate for a later retry.
    """
    state = store.read(run_id)
    output = Path(state["output_path"])
    records = state.get("analysis_snapshots", {})
    for project, record in records.items():
        if "pending" not in record:
            continue
        candidate = output / record["pending"]
        target = output / record["archive"]
        if candidate.exists():
            if not _archive_matches(candidate, record["analysis_id"]):
                raise RuntimeError(
                    f"Snapshot candidate does not match completed analysis: {candidate}"
                )
            os.replace(candidate, target)
        elif not _archive_matches(target, record["analysis_id"]):
            raise RuntimeError(f"Committed analysis snapshot unavailable: {target}")
        _sync_directory(target.parent)

        def published(current):
            del current["analysis_snapshots"][project]["pending"]
            return current

        store.mutate(run_id, published)
        log.info("Retained analysis snapshot: %s", target)
    staging = output / DIRECTORY / ".staging"
    if staging.is_dir() and not staging.is_symlink():
        for path in staging.glob("analysis-*.tar.gz"):
            path.unlink()


def recover_pending(store):
    """Recover publication and abandoned writes at startup/each daemon scan."""
    for state in store.list_states():
        output = Path(state["output_path"])
        if (
            not any("pending" in item for item in state.get("analysis_snapshots", {}).values())
            and not (output / DIRECTORY / ".staging").exists()
        ):
            continue
        try:
            with store.execution_lease(state["run_id"]):
                recover(store, state["run_id"])
        except ExecutionLeaseError:
            continue
        except Exception:
            log.exception(
                "Unable to recover analysis snapshots for %s; retained files preserved",
                state["run_id"],
            )
