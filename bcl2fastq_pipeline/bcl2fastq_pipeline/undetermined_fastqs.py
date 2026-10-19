"""Successful-conversion retention policy for the converters' reserved sample zero."""

import logging
import os
import re

from pathlib import Path

log = logging.getLogger(__name__)

# Illumina bcl2fastq/BCL Convert reserve S0 for unassigned reads. mkfastq
# retains this naming in both its published output and Martian staging trees.
# Never match a prefix alone, or the S1+ FASTQs of a sample named Undetermined.
UNDETERMINED_FASTQ_RE = re.compile(
    r"Undetermined_S0(?:_L[0-9]{3})?_(?:R[123]|I[12])_[0-9]{3}\.fastq\.gz"
)


def is_undetermined(path):
    return UNDETERMINED_FASTQ_RE.fullmatch(Path(path).name) is not None


def archive_exclusions():
    """7-Zip recursive basename masks, also covering directory archive inputs.

    7-Zip supports ? (one character), not regex character classes. Keep the
    reserved S0, read identity and fixed-width lane/chunk fields explicit.
    These switches work even when an input is an alias of a mkfastq tree.
    """
    return [
        f"-xr!Undetermined_S0{lane}_{read}_???.fastq.gz"
        for lane in ("", "_L???")
        for read in ("R1", "R2", "R3", "I1", "I2")
    ]


def cleanup(output):
    """Unlink managed files only; caller holds the run's execution lease.

    Demultiplexing must already be committed. A failure belongs to analysis
    preparation, so the normal analysis restart retries without reconversion.
    No retention manifest or second completion marker is needed for these files.
    """
    from bcl2fastq_pipeline import fastq_cleanup  # noqa: PLC0415
    from bcl2fastq_pipeline.state import StateConflictError  # noqa: PLC0415

    files = fastq_cleanup.plan(output, select=is_undetermined)
    removed = 0
    removed_bytes = 0
    errors = {}
    root_fd = os.open(output, fastq_cleanup.DIRECTORY_FLAGS)
    try:
        for relative, info in files.items():
            try:
                fastq_cleanup.unlink_fastq(root_fd, relative, info)
            except FileNotFoundError:
                # A previous/partial cleanup or removal since planning is normal.
                continue
            except (OSError, StateConflictError) as error:
                errors[relative] = str(error)
                log.error("Cannot unlink Undetermined FASTQ %s: %s", relative, error)
            else:
                removed += 1
                size = 0 if info["symlink"] else info["size"]
                removed_bytes += size
                log.info("Removed Undetermined FASTQ %s (%s regular-file bytes)", relative, size)
    finally:
        os.close(root_fd)
        log.info(
            "Undetermined cleanup: removed %s paths, %s regular-file bytes; %s unlink failures",
            removed,
            removed_bytes,
            len(errors),
        )
    if errors:
        detail = "; ".join(f"{name}: {error}" for name, error in errors.items())
        raise StateConflictError(
            "Undetermined cleanup failed; retry from analysis after fixing unlink errors: " + detail
        )
