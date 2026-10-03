"""Executed by the development runner from outside the checkout, in each venv."""

import importlib
import subprocess
import sys

from importlib import metadata
from pathlib import Path


def main():
    assert getattr(sys, "_bfq_smtp_guard", False), "Test SMTP guard did not load"
    wheel = sys.argv[1] == "wheel"
    prefix = Path(sys.prefix).resolve()
    for module in ("bcl2fastq_pipeline", "flowcell_manager", "configmaker.configmaker"):
        location = Path(importlib.import_module(module).__file__).resolve()
        print(f"{module}: {location}")
        if wheel:
            assert location.is_relative_to(prefix), f"Import escaped wheel environment: {location}"
    entrypoints = {ep.name: ep.value for ep in metadata.entry_points(group="console_scripts")}
    assert entrypoints["fm"] == entrypoints["flowcell-manager"]
    expected_version = f"BFQ {metadata.version('bcl2fastq-pipeline')}"
    for name in ("bfq", "fm", "flowcell-manager"):
        command = prefix / "bin" / name
        for option in ("--version", "--help"):
            result = subprocess.run(
                [str(command), option], check=True, text=True, capture_output=True
            )
            if option == "--version":
                assert result.stdout.strip() == expected_version, result.stdout
            else:
                assert "usage:" in result.stdout, result.stdout
            print(f"{name} {option}: OK")
    print(f"gcf-tools: {metadata.version('gcf-tools')}")


if __name__ == "__main__":
    main()
