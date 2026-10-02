"""Prove test-only mail protection survives real Python subprocess boundaries."""

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
