import sys

from argparse import ArgumentParser
from pathlib import Path

from bcl2fastq_pipeline.version import add_version_argument


def _remove_executable_directory_from_import_path() -> None:
    """Prevent neighboring scripts from shadowing installed Python packages."""
    executable_directory = Path(sys.argv[0]).resolve().parent
    sys.path[:] = [
        entry for entry in sys.path if not entry or Path(entry).resolve() != executable_directory
    ]


def main() -> None:
    _remove_executable_directory_from_import_path()
    parser = ArgumentParser(description="Monitor and process completed Illumina sequencing runs.")
    add_version_argument(parser)
    parser.parse_args()

    # This import must follow sys.path sanitization; configmaker.py otherwise shadows its package.
    from bcl2fastq_pipeline.cli import main as run_pipeline  # noqa: PLC0415

    run_pipeline()
