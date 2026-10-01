"""Read the BFQ generation from installed package metadata."""

from argparse import ArgumentParser
from importlib import metadata


def package_version() -> str | None:
    """Return the installed version, or None for an uninstalled source tree."""
    try:
        return metadata.version("bcl2fastq-pipeline")
    except metadata.PackageNotFoundError:
        return None


def add_version_argument(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--version",
        action="version",
        version=f"BFQ {package_version() or 'unknown (package not installed)'}",
        help="Show the installed BFQ generation and exit.",
    )
