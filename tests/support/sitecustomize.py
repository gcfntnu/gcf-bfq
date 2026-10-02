"""Test-only SMTP guard, inherited by Python subprocesses via PYTHONPATH.

Never install this module in a production package. Deliberate SMTP doubles are
allowed; a real smtplib connection is refused before any socket is opened.
"""

import sys


def forbid_smtp(event, _args):
    if event == "smtplib.connect":
        raise AssertionError("BFQ checks must mock SMTP; real connections are forbidden")


if not getattr(sys, "_bfq_smtp_guard", False):
    sys.addaudithook(forbid_smtp)
    sys._bfq_smtp_guard = True
