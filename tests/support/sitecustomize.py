"""Test-only safeguards inherited by Python subprocesses via PYTHONPATH.

SMTP protection always applies. Operational scenarios additionally opt into an
exact executable allowlist and network guard with BFQ_TEST_SCENARIOS=1. These
hooks prevent accidental execution in the normal Python harness, not malicious
code or native libraries: they are not an OS sandbox and are never installed in
production. See docs/development.md for the scenario contract.
"""

import json
import os
import sys

GUARD_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
SCENARIO_ENV = "BFQ_TEST_SCENARIOS"
ALLOWLIST_ENV = "BFQ_TEST_ALLOWED_EXECUTABLES"


def forbid_smtp(event, _args):
    if event == "smtplib.connect":
        raise AssertionError("BFQ checks must mock SMTP; real connections are forbidden")


def scenario_safeguards(event, args):
    # Read the mode dynamically so the direct-pytest scenario fixture can enable
    # the same guard without applying process restrictions to build/pip tooling.
    if os.environ.get(SCENARIO_ENV) != "1":
        return
    if event in {
        "os.system",
        "os.exec",
        "os.posix_spawn",
        "socket.connect",
        "socket.bind",
        "socket.sendto",
        "socket.sendmsg",
        "socket.getaddrinfo",
        "socket.gethostbyname",
        "socket.gethostbyaddr",
    }:
        raise AssertionError(
            f"BFQ scenarios forbid external operation {event}; use a declared double"
        )
    if event != "subprocess.Popen":
        return

    executable, argv, _cwd, environment = args
    executable = os.fsdecode(executable)
    try:
        allowed = json.loads(os.environ.get(ALLOWLIST_ENV, "[]"))
        valid = isinstance(allowed, list) and all(
            isinstance(path, str) and os.path.isabs(path) for path in allowed
        )
    except (ValueError, TypeError):
        valid = False
        allowed = []
    if not valid or not os.path.isabs(executable) or executable not in allowed:
        raise AssertionError(
            f"BFQ scenarios blocked unapproved executable {executable!r}; "
            "only the harness's exact installed command paths are allowed"
        )
    child = os.environ if environment is None else environment
    # A Popen env may use bytes keys/values. Normalize them before comparison.
    child = {os.fsdecode(key): os.fsdecode(value) for key, value in child.items()}
    if (
        child.get(SCENARIO_ENV) != "1"
        or child.get(ALLOWLIST_ENV) != os.environ.get(ALLOWLIST_ENV)
        or child.get("PYTHONPATH", "").split(os.pathsep)[0] != GUARD_DIRECTORY
    ):
        raise AssertionError("BFQ scenarios require inherited subprocess safeguards and PYTHONPATH")
    # Python is not in the normal scenario allowlist. The guard regressions use
    # it deliberately to prove that child startup cannot bypass these hooks.
    if os.path.basename(executable).startswith("python"):
        for raw_argument in argv[1:]:
            argument = os.fsdecode(raw_argument)
            if not argument.startswith("-") or argument in {"-c", "-m", "--"}:
                break
            if not argument.startswith("--") and set(argument[1:]) & set("ISE"):
                raise AssertionError(
                    "BFQ scenarios forbid Python flags that disable startup guards"
                )


if not getattr(sys, "_bfq_smtp_guard", False):
    sys.addaudithook(forbid_smtp)
    sys._bfq_smtp_guard = True
if not getattr(sys, "_bfq_external_guard", False):
    sys.addaudithook(scenario_safeguards)
    sys._bfq_external_guard = True

# Explicit test-only injection for installed CLI subprocesses. Never change BFQ's
# production configuration selection to accommodate the synthetic fixture.
scenario_root = os.environ.get("BFQ_TEST_SCENARIO_ROOT")
if scenario_root:
    try:
        if os.environ.get(SCENARIO_ENV) != "1":
            raise AssertionError("Scenario configuration requires BFQ_TEST_SCENARIOS=1")
        from operational import configure

        configure(scenario_root)
        sys._bfq_scenario_configured = os.path.realpath(scenario_root)
    except BaseException as error:
        # Python normally reports sitecustomize errors and continues executing.
        # A failed config injection must never fall through to /config or SMTP.
        print(f"BFQ scenario startup failed: {type(error).__name__}: {error}", file=sys.stderr)
        sys.stderr.flush()
        os._exit(86)
