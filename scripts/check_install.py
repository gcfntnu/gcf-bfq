"""Executed by the development runner from outside the checkout, in each venv."""

import importlib
import json
import os
import subprocess
import sys
import tempfile

from importlib import metadata
from pathlib import Path


def check_multiqc_wrapper(prefix):
    """Exercise the installed wrapper without a reporting engine or container."""
    candidate = (
        "gcfntnu/gcf-multiqc@sha256:"
        "f3a9873d88241e47f3bd38c0b68ae5b64c746b1d77b63a3182e00bbe3ffcf7b0"
    )
    with tempfile.TemporaryDirectory(prefix="multiqc-wrapper-", dir=Path.cwd()) as directory:
        root = Path(directory)
        config = root / "workflow images.config"
        record = root / "invocation.json"
        fake = root / "fake apptainer"
        fake.write_text(
            f"#!{sys.executable}\n"
            "import json, os, sys\n"
            "from pathlib import Path\n"
            "assert getattr(sys, '_bfq_smtp_guard', False), 'Child SMTP guard missing'\n"
            "Path(os.environ['BFQ_TEST_INVOCATION']).write_text(json.dumps({\n"
            "    'argv': sys.argv[1:], 'cwd': os.getcwd(), 'stdin': sys.stdin.read()}))\n"
            "print('synthetic report stdout')\n"
            "print('synthetic report stderr', file=sys.stderr)\n"
            "sys.exit(int(os.environ['BFQ_TEST_EXIT_CODE']))\n"
        )
        fake.chmod(0o700)
        arguments = [
            "--force",
            "--no-megaqc-upload",
            "--config",
            str(root / "report config.yaml"),
            "--outdir",
            str(root / "report output"),
            "--filename",
            "sequencer stats.html",
            "-m",
            "custom_content",
            "-m",
            "interop",
            "-m",
            "bclconvert",
            str(root / "report inputs"),
        ]
        environment = os.environ | {
            "GCF_WORKFLOWS_DOCKER_CONFIG": str(config),
            "BFQ_APPTAINER_COMMAND": str(fake),
            "BFQ_TEST_INVOCATION": str(record),
        }
        # Replacing the mapping must change selection without changing BFQ.
        for image, exit_code in ((candidate, 0), ("docker://example/multiqc:rollback", 19)):
            config.write_text(f"docker:\n  multiqc: {image}\n")
            environment["BFQ_TEST_EXIT_CODE"] = str(exit_code)
            result = subprocess.run(
                [str(prefix / "bin" / "multiqc"), *arguments],
                env=environment,
                cwd=root,
                text=True,
                input="synthetic stdin\n",
                capture_output=True,
                check=False,
            )
            assert result.returncode == exit_code, result.stderr
            resolved_image = image if "://" in image else f"docker://{image}"
            assert json.loads(record.read_text()) == {
                "argv": ["exec", resolved_image, "multiqc", *arguments],
                "cwd": str(root),
                "stdin": "synthetic stdin\n",
            }
            assert result.stdout == "synthetic report stdout\n", result.stdout
            assert result.stderr.splitlines() == [
                f"BFQ Apptainer: multiqc -> {resolved_image}",
                "synthetic report stderr",
            ], result.stderr

        # A missing mapping fails before external processing is attempted.
        record.unlink()
        config.write_text("docker: {}\n")
        result = subprocess.run(
            [str(prefix / "bin" / "multiqc"), *arguments],
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
        assert result.returncode == 2, result.stderr
        assert "'multiqc' required by 'multiqc' is missing" in result.stderr, result.stderr
        assert str(config) in result.stderr, result.stderr
        assert not record.exists(), "Invalid mapping launched external processing"
        assert not result.stdout, result.stdout
    print(
        "multiqc: installed wrapper, image authority, argv, streams and failures OK (fake Apptainer)"
    )


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
    result = subprocess.run(
        [str(prefix / "bin" / "bfq-gpu-smoke"), "--help"],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "usage:" in result.stdout, result.stdout
    print("bfq-gpu-smoke --help: OK")
    check_multiqc_wrapper(prefix)
    print(f"gcf-tools: {metadata.version('gcf-tools')}")


if __name__ == "__main__":
    main()
