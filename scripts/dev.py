#!/usr/bin/env python3
"""Checkout-owned development setup and offline checks. See docs/development.md."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import tomllib
import venv

from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEV = ROOT / ".dev"
BASELINE = ROOT / "requirements-dev.txt"
GUARD = ROOT / "tests" / "support"


def run(args, *, env, cwd=ROOT):
    print("+ " + " ".join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), cwd=cwd, env=env, check=True)


def identity(path):
    def git(*args):
        result = subprocess.run(
            ["git", "-C", str(path), *args], capture_output=True, text=True, check=False
        )
        return result.stdout.strip() if result.returncode == 0 else "not a Git checkout"

    return {
        "path": str(path),
        "commit": git("rev-parse", "HEAD"),
        "status": git("status", "--short"),
    }


def environment(work, *, guarded=False, python=None):
    env = os.environ.copy()
    for key in ("PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV", "CONDA_PREFIX", "PYTEST_ADDOPTS"):
        env.pop(key, None)
    for key in tuple(env):
        if key.startswith(("PIP_", "UV_")):
            if not guarded and key in {
                "PIP_INDEX_URL",
                "PIP_EXTRA_INDEX_URL",
                "PIP_TRUSTED_HOST",
                "PIP_CERT",
                "PIP_CLIENT_CERT",
                "PIP_TIMEOUT",
                "PIP_RETRIES",
            }:
                continue
            env.pop(key)
    for key, folder in (
        ("TMPDIR", "tmp"),
        ("XDG_CACHE_HOME", "cache"),
        ("PIP_CACHE_DIR", "cache/pip"),
        ("RUFF_CACHE_DIR", "cache/ruff"),
    ):
        directory = work / folder
        directory.mkdir(parents=True, exist_ok=True)
        env[key] = str(directory)
    env.update(
        BFQ_ENV="test",
        PYTHONNOUSERSITE="1",
        PYTHONDONTWRITEBYTECODE="1",
        PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
        PIP_CONFIG_FILE=os.devnull,
        PIP_DISABLE_PIP_VERSION_CHECK="1",
    )
    if python:
        env["PATH"] = str(python.parent) + os.pathsep + os.defpath
    if guarded:
        env["PYTHONPATH"] = str(GUARD)
        env["PIP_NO_INDEX"] = "1"
    return env


def copy_source(source, target):
    shutil.copytree(
        source,
        target,
        ignore=shutil.ignore_patterns(
            ".git",
            ".dev",
            ".venv",
            "venv",
            "env",
            "__pycache__",
            "*.egg-info",
            ".pytest_cache",
            ".ruff_cache",
            "build",
            "dist",
        ),
    )


def digest_tree(path):
    digest = hashlib.sha256()
    for file in sorted(path.rglob("*")):
        if file.is_file():
            digest.update(str(file.relative_to(path)).encode())
            digest.update(file.read_bytes())
    return digest.hexdigest()


@contextmanager
def checkout_lock():
    if DEV.is_symlink():
        raise RuntimeError(f"Refusing symlinked workspace {DEV}")
    DEV.mkdir(exist_ok=True)
    with (DEV / "lock").open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(
                "This checkout already has setup/check running; use another worktree"
            ) from error
        yield


def workspace(kind):
    parent = DEV / kind
    parent.mkdir(exist_ok=True)
    return Path(tempfile.mkdtemp(prefix="run-", dir=parent))


def setup(args):
    companion = Path(args.gcf_tools).expanduser().resolve() if args.gcf_tools else None
    if (
        companion
        and not (companion / "setup.py").is_file()
        and not (companion / "pyproject.toml").is_file()
    ):
        raise RuntimeError(f"No Python project at {companion}; supply a gcf-tools checkout")
    work = workspace("setups")
    print(f"Setup artifacts: {work}", flush=True)
    python = work / "venv" / "bin" / "python"
    venv.EnvBuilder(with_pip=True, symlinks=True).create(python.parent.parent)
    env = environment(work, python=python)
    wheelhouse = work / "wheels"
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    dependencies = [dep for dep in project["dependencies"] if not dep.startswith("gcf-tools")]
    dependencies += project["optional-dependencies"]["dev"]
    dependencies += ["setuptools>=68", "wheel"]
    selection = {"baseline": BASELINE.read_text()}
    if companion:
        snapshot = work / "gcf-tools"
        copy_source(companion, snapshot)
        selection = {"local": identity(companion), "snapshot_sha256": digest_tree(snapshot)}
        dependencies.append(str(snapshot))
    else:
        dependencies += ["-r", str(BASELINE)]
    run([python, "-m", "pip", "wheel", "--wheel-dir", wheelhouse, *dependencies], env=env)
    offline = environment(work, guarded=True, python=python)
    run(
        [
            python,
            "-m",
            "pip",
            "install",
            "--no-index",
            "--find-links",
            wheelhouse,
            "setuptools>=68",
            "wheel",
        ],
        env=offline,
    )
    run(
        [
            python,
            "-m",
            "pip",
            "install",
            "--no-index",
            "--find-links",
            wheelhouse,
            "--no-build-isolation",
            "-e",
            f"{ROOT}[dev]",
        ],
        env=offline,
    )
    run([python, "-m", "pip", "check"], env=offline)
    state = {
        "python": str(python),
        "wheelhouse": str(wheelhouse),
        "companion": selection,
        "python_version": platform.python_version(),
        "source": identity(ROOT),
        "baseline_sha256": hashlib.sha256(BASELINE.read_bytes()).hexdigest(),
        "pyproject_sha256": hashlib.sha256((ROOT / "pyproject.toml").read_bytes()).hexdigest(),
    }
    (work / "setup.json").write_text(json.dumps(state, indent=2) + "\n")
    pending = DEV / "current.json.tmp"
    pending.write_text(json.dumps(state, indent=2) + "\n")
    pending.replace(DEV / "current.json")
    print("Setup ready. Run: python3.11 scripts/dev.py check all")


def smoke(python, mode, env):
    # /tmp is deliberately outside the checkout; delete only our own empty cwd.
    with tempfile.TemporaryDirectory(prefix="bfq-installed-", dir="/tmp") as cwd:
        run([python, "-c", (ROOT / "scripts/check_install.py").read_text(), mode], env=env, cwd=cwd)


def load_setup():
    path = DEV / "current.json"
    if not path.is_file():
        raise RuntimeError("No successful setup. Run: python3.11 scripts/dev.py setup")
    state = json.loads(path.read_text())
    for file, key in ((BASELINE, "baseline_sha256"), (ROOT / "pyproject.toml", "pyproject_sha256")):
        if hashlib.sha256(file.read_bytes()).hexdigest() != state[key]:
            raise RuntimeError(f"{file.name} changed; rerun setup before checks")
    for key in ("python", "wheelhouse"):
        # The executable itself may be a venv symlink to the base interpreter.
        path = Path(state[key])
        if not path.parent.resolve().is_relative_to(DEV.resolve()) or not path.exists():
            raise RuntimeError("Setup moved or is incomplete; rerun setup in this checkout")
    return state


def check(args):
    state = load_setup()
    work = workspace("runs")
    print(f"Check artifacts (retained on success/failure): {work}", flush=True)
    python = Path(state["python"])
    env = environment(work, guarded=True, python=python)
    report = {
        "source": identity(ROOT),
        "setup": state,
        "platform": platform.platform(),
        "profile": args.profile,
    }
    (work / "identity.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    with (work / "packages.txt").open("w") as packages:
        subprocess.run(
            [str(python), "-m", "pip", "freeze", "--all"], env=env, stdout=packages, check=True
        )
    run(
        [
            python,
            "-c",
            "import sys; assert getattr(sys, '_bfq_smtp_guard', False), 'SMTP guard missing'",
        ],
        env=env,
        cwd=work,
    )
    run([python, "-m", "pip", "check"], env=env)
    if args.profile != "wheel":
        run([python, "-m", "ruff", "check", "--no-fix", ROOT], env=env)
        run([python, "-m", "ruff", "format", "--check", ROOT], env=env)
        selection = args.pytest_args or [str(ROOT / "tests")]
        run(
            [
                python,
                "-m",
                "pytest",
                "-o",
                f"cache_dir={work / 'pytest-cache'}",
                "--basetemp",
                work / "pytest",
                *selection,
            ],
            env=env,
        )
        smoke(python, "editable", env)
    if args.profile in ("all", "wheel"):
        source = work / "source"
        source.mkdir()
        for name in ("pyproject.toml", "README.md"):
            shutil.copy2(ROOT / name, source / name)
        copy_source(ROOT / "bcl2fastq_pipeline", source / "bcl2fastq_pipeline")
        run(
            [python, "-m", "build", "--no-isolation", "--outdir", work / "dist", source],
            env=env,
            cwd=work,
        )
        wheel_python = work / "wheel-venv" / "bin" / "python"
        venv.EnvBuilder(with_pip=True, symlinks=True).create(wheel_python.parent.parent)
        wheel_env = environment(work, guarded=True, python=wheel_python)
        wheel = next((work / "dist").glob("*.whl"))
        run(
            [
                wheel_python,
                "-m",
                "pip",
                "install",
                "--no-index",
                "--find-links",
                state["wheelhouse"],
                wheel,
            ],
            env=wheel_env,
            cwd=work,
        )
        run([wheel_python, "-m", "pip", "check"], env=wheel_env, cwd=work)
        smoke(wheel_python, "wheel", wheel_env)
    print(f"Checks passed: {args.profile}; evidence: {work}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    setup_parser = commands.add_parser(
        "setup", help="Acquire dependencies into a new checkout-local venv"
    )
    setup_parser.add_argument(
        "--gcf-tools", metavar="PATH", help="Snapshot an explicit local gcf-tools checkout"
    )
    check_parser = commands.add_parser(
        "check", help="Run offline checks using the last successful setup"
    )
    check_parser.add_argument(
        "profile", nargs="?", default="fast", choices=("fast", "all", "wheel")
    )
    check_parser.add_argument(
        "pytest_args", nargs=argparse.REMAINDER, help="Optional pytest selection after --"
    )
    args = parser.parse_args()
    if sys.platform != "linux" or sys.version_info[:2] != (3, 11):
        parser.error(
            "Supported development baseline: Linux with Python 3.11 and venv/ensurepip. Run with python3.11."
        )
    if shutil.which("git") is None:
        parser.error("Git is required for source identity; install git and retry")
    if args.command == "check" and args.pytest_args[:1] == ["--"]:
        args.pytest_args = args.pytest_args[1:]
    try:
        with checkout_lock():
            setup(args) if args.command == "setup" else check(args)
    except (RuntimeError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(
            1,
            f"Development command failed: {error}\nArtifacts are retained under {DEV}; no other checkout was changed.\n",
        )


if __name__ == "__main__":
    main()
