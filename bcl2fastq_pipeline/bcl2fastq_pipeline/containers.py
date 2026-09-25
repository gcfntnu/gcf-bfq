"""Resolve and run external BFQ tools through Apptainer.

Container image references are owned by the checked-out gcf-workflows revision.
Keeping the lookup here means BFQ does not need a second, drifting list of image
versions.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

log = logging.getLogger(__name__)

DEFAULT_DOCKER_CONFIG = Path("/opt/gcf-workflows/docker.config")
DOCKER_CONFIG_ENV = "GCF_WORKFLOWS_DOCKER_CONFIG"
APPTAINER_COMMAND_ENV = "BFQ_APPTAINER_COMMAND"


class ContainerConfigError(RuntimeError):
    """Raised when the workflow container configuration cannot be used."""


@dataclass(frozen=True)
class ToolSpec:
    image_key: str
    executable: str


TOOLS = {
    "bcl-convert": ToolSpec("bcl-convert", "bcl-convert"),
    "bcl2fastq": ToolSpec("bcl2fastq", "bcl2fastq"),
    "cellranger": ToolSpec("cellranger", "cellranger"),
    "cellranger-atac": ToolSpec("cellranger-atac", "cellranger-atac"),
    "multiqc": ToolSpec("multiqc", "multiqc"),
    "spaceranger": ToolSpec("spaceranger", "spaceranger"),
}


def docker_config_path() -> Path:
    """Return the configured gcf-workflows container configuration path."""
    return Path(os.environ.get(DOCKER_CONFIG_ENV, DEFAULT_DOCKER_CONFIG))


def load_images(config_path: Path | None = None) -> dict[str, str]:
    """Load the ``docker`` image mapping from gcf-workflows/docker.config."""
    path = config_path or docker_config_path()
    try:
        with path.open() as config_file:
            config = yaml.safe_load(config_file)
    except OSError as error:
        raise ContainerConfigError(
            f"Cannot read container configuration {path}: {error}"
        ) from error
    except yaml.YAMLError as error:
        raise ContainerConfigError(
            f"Invalid YAML in container configuration {path}: {error}"
        ) from error

    if not isinstance(config, dict) or not isinstance(config.get("docker"), dict):
        raise ContainerConfigError(
            f"Container configuration {path} must contain a 'docker' mapping"
        )

    images = config["docker"]
    invalid = [
        key
        for key, value in images.items()
        if not isinstance(key, str) or not isinstance(value, str)
    ]
    if invalid:
        raise ContainerConfigError(
            f"Container configuration {path} contains non-string image references: {invalid}"
        )
    return images


def resolve_image(tool: str, config_path: Path | None = None) -> str:
    """Resolve a supported BFQ tool to its image reference."""
    try:
        spec = TOOLS[tool]
    except KeyError as error:
        supported = ", ".join(sorted(TOOLS))
        raise ContainerConfigError(
            f"Unknown BFQ container tool {tool!r}; expected: {supported}"
        ) from error

    images = load_images(config_path)
    try:
        image = images[spec.image_key]
    except KeyError as error:
        path = config_path or docker_config_path()
        raise ContainerConfigError(
            f"Container image {spec.image_key!r} required by {tool!r} is missing from {path}"
        ) from error
    return image


def apptainer_image(image: str) -> str:
    """Convert an unqualified Docker image reference for Apptainer."""
    return image if "://" in image else f"docker://{image}"


def build_command(
    tool: str,
    arguments: Sequence[str] = (),
    config_path: Path | None = None,
) -> list[str]:
    """Build an argv-safe Apptainer command for a supported tool."""
    spec = TOOLS.get(tool)
    if spec is None:
        resolve_image(tool, config_path)  # Raise the common, descriptive error.
        raise AssertionError("unreachable")

    image = resolve_image(tool, config_path)
    apptainer = os.environ.get(APPTAINER_COMMAND_ENV, "apptainer")
    return [apptainer, "exec", apptainer_image(image), spec.executable, *arguments]


def run_tool(
    tool: str,
    arguments: Sequence[str] = (),
    config_path: Path | None = None,
) -> int:
    """Run a supported tool while preserving its terminal stdin/stdout/stderr."""
    command = build_command(tool, arguments, config_path)
    log.info("Running %s from %s", tool, command[2])
    print(f"BFQ Apptainer: {tool} -> {command[2]}", file=sys.stderr, flush=True)
    return subprocess.call(command)


def tool_main(tool: str, argv: Sequence[str] | None = None) -> int:
    """Run one configured tool as a transparent console-script wrapper."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    try:
        return run_tool(tool, arguments)
    except ContainerConfigError as error:
        print(f"{tool}: {error}", file=sys.stderr)
        return 2


def bcl_convert_main() -> int:
    return tool_main("bcl-convert")


def bcl2fastq_main() -> int:
    return tool_main("bcl2fastq")


def cellranger_main() -> int:
    return tool_main("cellranger")


def cellranger_atac_main() -> int:
    return tool_main("cellranger-atac")


def multiqc_main() -> int:
    return tool_main("multiqc")


def spaceranger_main() -> int:
    return tool_main("spaceranger")


def main(argv: Sequence[str] | None = None) -> int:
    """Run ``bfq-container TOOL [ARG ...]`` for testing and diagnostics."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments[0] in {"-h", "--help"}:
        tools = ", ".join(sorted(TOOLS))
        print(f"usage: bfq-container TOOL [ARG ...]\n\nsupported tools: {tools}")
        return 0 if arguments else 2

    tool, *tool_arguments = arguments
    return tool_main(tool, tool_arguments)
