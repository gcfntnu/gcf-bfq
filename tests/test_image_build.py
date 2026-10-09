"""Exercise image-build command routing with Docker and remote Git doubled."""

import json
import os
import subprocess
import sys

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def build_runner(tmp_path):
    tools = tmp_path / "fake tools"
    tools.mkdir()
    calls = tmp_path / "calls.jsonl"
    executable = (
        f"#!{sys.executable}\n"
        + """
import json
import os
import pathlib
import sys

command = pathlib.Path(sys.argv[0]).name
args = sys.argv[1:]
with open(os.environ["BUILD_CALLS"], "a") as stream:
    stream.write(json.dumps([command, *args]) + "\\n")
if command == "git":
    if args[0] == "check-ref-format":
        sys.exit(1 if " " in args[1] else 0)
    if args[0] == "ls-remote":
        if os.environ.get("FAIL_REMOTE"):
            sys.exit(2)
        sha = "a" * 40 if "gcf-tools" in args[2] else "b" * 40
        print(sha + "\\t" + args[3])
    elif args == ["rev-parse", "HEAD"]:
        print("c" * 40)
elif command == "sudo":
    os.execvp(args[0], args)
elif command == "docker" and args[0] == "build":
    sys.exit(int(os.environ.get("FAIL_BUILD", "0")))
"""
    )
    for name in ("docker", "git", "sudo"):
        path = tools / name
        path.write_text(executable)
        path.chmod(0o755)

    def run(*args, direct=True, **extra_env):
        env = {
            **os.environ,
            "PATH": f"{tools}{os.pathsep}{os.environ['PATH']}",
            "BUILD_CALLS": str(calls),
            **extra_env,
        }
        env.pop("BFQ_DOCKER", None)
        if direct:
            env["BFQ_DOCKER"] = str(tools / "docker")
        result = subprocess.run(
            ["bash", str(ROOT / "build-tag-push.sh"), *args],
            cwd=tmp_path,
            env=env,
            text=True,
            capture_output=True,
            check=False,
        )
        recorded = (
            [json.loads(line) for line in calls.read_text().splitlines()] if calls.exists() else []
        )
        return result, recorded

    return run


def test_test_build_pins_branches_and_uses_selected_base(build_runner):
    # Shell metacharacters in valid refs must remain literal argv values.
    workflow_branch = "issues/148-$(not-a-command)"
    result, calls = build_runner(
        "test", "dev-test", "-b", "gcfntnu/bfq:base-test", "-w", workflow_branch
    )
    assert result.returncode == 0, result.stderr
    docker_calls = [call for call in calls if call[0] == "docker"]
    assert docker_calls == [
        [
            "docker",
            "build",
            "-t",
            "gcfntnu/bfq:dev-test",
            ".",
            "-f",
            "dockerfile-test",
            "--build-arg",
            "BASE_IMAGE=gcfntnu/bfq:base-test",
            "--build-arg",
            f"BFQ_REVISION={'c' * 40}",
            "--build-arg",
            "GCF_TOOLS_BRANCH=bfq-dev",
            "--build-arg",
            f"GCF_TOOLS_REV={'a' * 40}",
            "--build-arg",
            f"GCF_WORKFLOWS_BRANCH={workflow_branch}",
            "--build-arg",
            f"GCF_WORKFLOWS_REV={'b' * 40}",
        ],
        ["docker", "push", "gcfntnu/bfq:dev-test"],
    ]
    assert [
        "git",
        "ls-remote",
        "--exit-code",
        "https://github.com/gcfntnu/gcf-workflows.git",
        f"refs/heads/{workflow_branch}",
    ] in calls


@pytest.mark.parametrize("mode,tag", [("base", "base-test"), ("prod", "prod2-2")])
def test_base_and_production_keep_default_sudo_and_policy(build_runner, mode, tag):
    result, calls = build_runner(mode, tag, direct=False)
    assert result.returncode == 0, result.stderr
    assert calls == [
        ["sudo", "docker", "build", "-t", f"gcfntnu/bfq:{tag}", ".", "-f", f"dockerfile-{mode}"],
        ["docker", "build", "-t", f"gcfntnu/bfq:{tag}", ".", "-f", f"dockerfile-{mode}"],
        ["sudo", "docker", "push", f"gcfntnu/bfq:{tag}"],
        ["docker", "push", f"gcfntnu/bfq:{tag}"],
    ]


def test_failed_build_never_pushes(build_runner):
    result, calls = build_runner("base", "base-test", FAIL_BUILD="7")
    assert result.returncode == 7
    assert len(calls) == 1
    assert calls[0][:2] == ["docker", "build"]


def test_missing_companion_branch_never_builds_or_pushes(build_runner):
    result, calls = build_runner("test", "dev-test", FAIL_REMOTE="1")
    assert result.returncode != 0
    assert not any(call[0] == "docker" for call in calls)


@pytest.mark.parametrize(
    "args",
    [
        (),
        ("other", "tag"),
        ("base", "invalid/tag"),
        ("prod", "prod2-2", "-b", "gcfntnu/bfq:base-test"),
        ("test", "dev-test", "-x"),
        ("test", "dev-test", "-b"),
        ("test", "dev-test", "extra"),
    ],
)
def test_invalid_build_invocation_has_no_side_effects(build_runner, args):
    result, calls = build_runner(*args)
    assert result.returncode != 0
    assert not calls
