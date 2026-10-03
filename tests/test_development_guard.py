"""Prove test-only mail and execution safeguards cross Python subprocess boundaries."""

import json
import os
import subprocess
import sys

import pytest


@pytest.mark.parametrize("smtp_class", ["SMTP", "SMTP_SSL"])
def test_subprocess_smtp_is_blocked_before_socket_on_custom_port(tmp_path, smtp_class):
    driver = f"""
import importlib, smtplib, socket, sys
assert getattr(sys, '_bfq_smtp_guard', False), 'Missing test startup guard'
# Even reloading smtplib must not bypass the audit hook.
importlib.reload(smtplib)
def forbidden_socket(*args, **kwargs):
    raise RuntimeError('SOCKET_REACHED')
socket.create_connection = forbidden_socket
try:
    smtplib.{smtp_class}('relay.example.test', 25252)
except AssertionError as error:
    assert 'BFQ checks must mock SMTP' in str(error)
    print('BLOCKED_BEFORE_SOCKET')
else:
    raise AssertionError('SMTP was not blocked')
"""
    result = subprocess.run(
        [sys.executable, "-c", driver], cwd=tmp_path, capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "BLOCKED_BEFORE_SOCKET"
    assert "SOCKET_REACHED" not in result.stderr


def strict_child(tmp_path, driver, *, allowed=()):
    """Start a guarded Python driver; no scientific executables are launched."""
    env = os.environ.copy()
    env["BFQ_TEST_SCENARIOS"] = "1"
    env["BFQ_TEST_ALLOWED_EXECUTABLES"] = json.dumps(list(allowed))
    return subprocess.run(
        [sys.executable, "-c", driver],
        env=env,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.parametrize(
    "operation",
    [
        "subprocess.run(['snakemake', '--version'], check=True)",
        "subprocess.run([external], check=True)",
        "subprocess.run(external, shell=True, check=True)",
        "subprocess.run([external], executable='/bin/sh', check=True)",
        "os.system(external)",
        "os.execv(external, [external])",
        "os.posix_spawn(external, [external], os.environ)",
    ],
)
def test_scenarios_block_external_execution_before_launch(tmp_path, operation):
    marker = tmp_path / "EXTERNAL_WAS_LAUNCHED"
    external = tmp_path / "external-tool"
    external.write_text(f"#!/bin/sh\ntouch '{marker}'\n")
    external.chmod(0o755)
    driver = f"""
import os, subprocess, sys
assert getattr(sys, '_bfq_external_guard', False), 'Missing external-process guard'
external = {str(external)!r}
try:
    {operation}
except AssertionError as error:
    assert 'BFQ scenarios' in str(error), str(error)
    print('BLOCKED_BEFORE_LAUNCH')
else:
    raise AssertionError('External operation was not blocked')
"""
    result = strict_child(tmp_path, driver)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "BLOCKED_BEFORE_LAUNCH"
    assert not marker.exists()


@pytest.mark.parametrize(
    "child_change",
    [
        "env.pop('BFQ_TEST_SCENARIOS')",
        "env.pop('BFQ_TEST_ALLOWED_EXECUTABLES')",
        "env.pop('PYTHONPATH')",
        "env['PYTHONPATH'] = '/unrelated'",
        "command.insert(1, '-I')",
        "command.insert(1, '-S')",
        "command.insert(1, '-ES')",
    ],
)
def test_scenarios_reject_children_that_drop_startup_guards(tmp_path, child_change):
    driver = f"""
import os, subprocess, sys
command = [sys.executable, '-c', "print('CHILD_LAUNCHED')"]
env = os.environ.copy()
{child_change}
try:
    subprocess.run(command, env=env, check=True)
except AssertionError as error:
    assert 'BFQ scenarios' in str(error), str(error)
    print('GUARD_REMOVAL_BLOCKED')
else:
    raise AssertionError('Child escaped its startup guard')
"""
    result = strict_child(tmp_path, driver, allowed=[sys.executable])
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "GUARD_REMOVAL_BLOCKED"


def test_scenarios_allow_only_declared_child_and_inherit_guard(tmp_path):
    grandchild = """
import subprocess, sys
assert getattr(sys, '_bfq_external_guard', False)
try:
    subprocess.run(['/bin/true'], check=True)
except AssertionError as error:
    assert 'unapproved executable' in str(error)
    print('GRANDCHILD_GUARDED')
else:
    raise AssertionError('Grandchild escaped guard')
"""
    driver = f"""
import subprocess, sys
subprocess.run([sys.executable, '-c', {grandchild!r}], check=True)
"""
    result = strict_child(tmp_path, driver, allowed=[sys.executable])
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "GRANDCHILD_GUARDED"


@pytest.mark.parametrize(
    "operation",
    [
        "socket.getaddrinfo('example.invalid', 80)",
        "socket.socket().connect(('127.0.0.1', 9))",
        "socket.socket(socket.AF_INET, socket.SOCK_DGRAM).sendto(b'blocked', ('127.0.0.1', 9))",
    ],
)
def test_scenarios_block_network_attempts(tmp_path, operation):
    driver = f"""
import socket
try:
    {operation}
except AssertionError as error:
    assert 'BFQ scenarios forbid external operation socket.' in str(error), str(error)
    print('NETWORK_BLOCKED')
else:
    raise AssertionError('Network operation was not blocked')
"""
    result = strict_child(tmp_path, driver)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "NETWORK_BLOCKED"


def test_scenario_config_startup_fails_closed(tmp_path):
    env = os.environ.copy()
    env["BFQ_TEST_SCENARIO_ROOT"] = str(tmp_path)
    # An accidentally missing strict mode must stop before user code starts.
    env.pop("BFQ_TEST_SCENARIOS", None)
    result = subprocess.run(
        [sys.executable, "-c", "print('USER_CODE_REACHED')"],
        env=env,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 86
    assert "Scenario configuration requires BFQ_TEST_SCENARIOS=1" in result.stderr
    assert "USER_CODE_REACHED" not in result.stdout


def test_direct_pytest_scenario_mode_guards_current_process(monkeypatch):
    monkeypatch.setenv("BFQ_TEST_SCENARIOS", "1")
    monkeypatch.setenv("BFQ_TEST_ALLOWED_EXECUTABLES", "[]")
    assert getattr(sys, "_bfq_external_guard", False)
    with pytest.raises(AssertionError, match="blocked unapproved executable"):
        subprocess.run(["/bin/true"], check=True)


def test_missing_scenario_root_fails_before_user_code(tmp_path):
    env = os.environ.copy()
    env.update(
        BFQ_TEST_SCENARIOS="1",
        BFQ_TEST_ALLOWED_EXECUTABLES="[]",
        BFQ_TEST_SCENARIO_ROOT=str(tmp_path / "absent-fixture"),
    )
    result = subprocess.run(
        [sys.executable, "-c", "print('USER_CODE_REACHED')"],
        env=env,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 86
    assert "BFQ scenario startup failed: FileNotFoundError:" in result.stderr
    assert str(tmp_path / "absent-fixture") in result.stderr
    assert "USER_CODE_REACHED" not in result.stdout
