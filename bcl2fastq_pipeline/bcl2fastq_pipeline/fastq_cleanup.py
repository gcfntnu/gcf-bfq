"""Narrow, restartable removal of delivered FASTQs; never inspect archive contents."""

import os
import stat

from contextlib import contextmanager
from pathlib import Path

from bcl2fastq_pipeline.state import StateConflictError, utcnow

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW


def _raise_walk_error(error):
    raise error


def plan(output):
    """Return relative paths and lstat identity; don't traverse directory symlinks."""
    files = {}
    try:
        if output.is_symlink():
            raise StateConflictError(f"Output directory must not be a symlink: {output}")
        # Opening explicitly also rejects missing/unreadable output roots.
        with os.scandir(output):
            pass
        for root, dirs, names in os.walk(output, followlinks=False, onerror=_raise_walk_error):
            dirs[:] = [name for name in dirs if not (Path(root) / name).is_symlink()]
            for name in sorted(names):
                if not name.endswith(FASTQ_SUFFIXES):
                    continue
                path = Path(root) / name
                info = path.lstat()
                if not (stat.S_ISREG(info.st_mode) or stat.S_ISLNK(info.st_mode)):
                    raise StateConflictError(f"FASTQ is not a regular file or symlink: {path}")
                files[str(path.relative_to(output))] = {
                    "device": info.st_dev,
                    "inode": info.st_ino,
                    "mode": info.st_mode,
                    "size": info.st_size,
                    "mtime_ns": info.st_mtime_ns,
                    "symlink": stat.S_ISLNK(info.st_mode),
                }
    except OSError as error:
        raise StateConflictError(f"Cannot inspect FASTQs in {output}: {error}") from error
    return dict(sorted(files.items()))


def required_archives(state, output, files):
    """Match archive_worker naming, requiring every known project's delivery pair."""
    if state["status"] != "completed" or state["stages"]["finalization"]["status"] != "completed":
        raise StateConflictError("clean-fastqs requires completed finalization")
    projects = set(state.get("projects") or [])
    projects.update(p.name for p in output.iterdir() if p.name.startswith("GCF-") and p.is_dir())
    projects.update(
        part for name in files for part in Path(name).parts[:-1] if part.startswith("GCF-")
    )
    if not projects:
        raise StateConflictError(
            "Cannot establish expected FASTQ delivery archives: no projects recorded"
        )
    retained = []
    run_date = state["run_id"].split("_", 1)[0]
    for project in sorted(projects):
        if Path(project).name != project or project in {".", ".."}:
            raise StateConflictError(f"Invalid project name: {project!r}")
        base = f"{project}_{run_date}"
        retained.extend([output / f"{base}.7za", output / f"md5sum_{base}_archive.txt"])
    missing = [str(path) for path in retained if not path.is_file()]
    if missing:
        raise StateConflictError(
            "Missing required delivery archive/checksum file(s):\n  " + "\n  ".join(missing)
        )
    return retained


@contextmanager
def parent_descriptor(root_fd, relative):
    """Resolve parents without following links, including after preview races."""
    parts = Path(relative).parts
    if not parts or Path(relative).is_absolute() or any(p in {".", ".."} for p in parts):
        raise StateConflictError(f"Invalid cleanup path: {relative!r}")
    descriptor = os.dup(root_fd)
    try:
        for part in parts[:-1]:
            child = os.open(part, DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        yield descriptor, parts[-1]
    finally:
        os.close(descriptor)


def unlink_fastq(root_fd, relative, expected):
    with parent_descriptor(root_fd, relative) as (parent_fd, name):
        info = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        identity = (info.st_dev, info.st_ino, info.st_mode, info.st_size, info.st_mtime_ns)
        if identity != tuple(expected[k] for k in ("device", "inode", "mode", "size", "mtime_ns")):
            raise StateConflictError(f"FASTQ changed after preview: {relative}")
        # unlink never follows the final component, even if replaced by a symlink.
        os.unlink(name, dir_fd=parent_fd)


def execute(store, run_id, output, files):
    """Caller owns execution lease and has revalidated eligibility and preview."""
    previous = store.read(run_id).get("fastq_cleanup", {})
    required = dict(previous.get("required_files", {}))
    required.update(
        {name: None if info["symlink"] else info["size"] for name, info in files.items()}
    )
    record = {
        "status": "in_progress",
        "started_at": utcnow(),
        "finished_at": None,
        "required_files": required,
        "selected_files": list(files),
        "removed_files": [],
        "errors": {},
        "estimated_removed_bytes": 0,
    }

    def save(state):
        state["output_path"] = str(output)
        state["fastq_cleanup"] = record
        return state

    # Persist all potentially removed paths before the first unlink. A hard kill
    # leaves in_progress plus this manifest, so downstream reruns still fail safely.
    store.mutate(run_id, save)
    root_fd = None
    try:
        root_fd = os.open(output, DIRECTORY_FLAGS)
        for relative, info in files.items():
            try:
                unlink_fastq(root_fd, relative, info)
            except (OSError, StateConflictError) as error:
                record["errors"][relative] = str(error)
            else:
                record["removed_files"].append(relative)
                if not info["symlink"]:
                    record["estimated_removed_bytes"] += info["size"]
    except BaseException as error:
        record["status"] = "interrupted"
        record["errors"]["operation"] = str(error) or type(error).__name__
        raise
    else:
        record["status"] = "partial" if record["errors"] else "completed"
    finally:
        if root_fd is not None:
            os.close(root_fd)
        record["finished_at"] = utcnow()
        store.mutate(run_id, save)
    if record["errors"]:
        detail = "\n  ".join(f"{name}: {error}" for name, error in record["errors"].items())
        raise StateConflictError(f"FASTQ cleanup partially failed; safe to retry:\n  {detail}")
    return store.read(run_id)


def require_restored(state, output, from_stage):
    """Check the saved manifest, without hashing or reading FASTQ contents."""
    record = state.get("fastq_cleanup", {})
    if from_stage == "demultiplexing" or record.get("status") == "regenerated":
        return
    required = record.get("required_files", {})
    if not required:
        return
    missing = []
    for relative, size in required.items():
        path = output / relative
        # Treat a dangling restored link as absent. Links to restored data are OK.
        if not path.is_file() or (size is not None and path.stat().st_size != size):
            missing.append(relative)
    if missing:
        raise StateConflictError(
            f"Cannot restart from {from_stage}: FASTQs were cleaned. Restore all cleaned FASTQs "
            "from retained delivery archives to their original output paths before rerunning "
            "(regular files must match their recorded sizes), or rerun from demultiplexing "
            "to regenerate reads. Missing or size-mismatched FASTQs:\n  " + "\n  ".join(missing)
        )
