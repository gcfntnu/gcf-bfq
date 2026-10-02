# Isolated development and verification

Read [CONTRIBUTING](../CONTRIBUTING.md) for task ownership, compatibility and PR
handoff. These commands are for development; the [README](../README.md) owns
production installation, image builds and operational use.

## Clean checkout

The supported development baseline is **Linux, CPython 3.11, Git and Python's
`venv`/`ensurepip` support**. CI runs that baseline on `ubuntu-latest`; tests use
Linux facilities such as `/proc` and `flock`. This does not narrow the production
package's existing Python `>=3.11` requirement. Other developer platforms/Python
versions are not established by this recipe; use a Linux development environment.

Obtain Python 3.11 with your normal system/environment tooling. On distributions
that package venv separately, install the matching `python3.11-venv` package.
Setup needs ordinary HTTPS access to PyPI and the pinned GitHub archive. It does
not require instrument access, SMTP credentials, containers or a facility server.
The runner reports a missing Python/Git prerequisite; venv or dependency errors
must be resolved in personal setup, not by changing production dependencies.

From the repository root (or use an absolute path to `scripts/dev.py`):

```bash
python3.11 scripts/dev.py setup
python3.11 scripts/dev.py check all
```

No activation is required. Setup creates a fresh venv and dependency wheelhouse
under `.dev/setups/<unique>/`, installs BFQ editable there and runs `pip check`.
Only after success does `.dev/current.json` point to that setup. A failed or
interrupted setup leaves the prior successful environment selected. Repeated
setup creates a new environment; it never upgrades your global Python or another
checkout. It does not edit shell startup files, `/opt`, `/mnt` or site config.

Setup may leave ignored BFQ editable-install metadata in this checkout. Build
checks copy the package into their own run directory, so they do not reuse an old
`build/` or `dist/` tree. Do not move a prepared checkout/venv: Python console
scripts contain absolute paths. Run setup again at the new location.

## Dependency selection

[requirements-dev.txt](../requirements-dev.txt) is the single immutable
`gcf-tools` baseline for this recipe and CI. BFQ's runtime/development requirements
remain in [pyproject.toml](../pyproject.toml). Setup prepares wheels for those
requirements and the build backend. Repeatable checks use only those wheels,
with pip's index disabled and build isolation disabled; they do not download
fixtures or dependencies. This is a fixed companion revision, **not** a lock of
every transitive PyPI version. Each successful setup resolves other requirements
anew; the check report records installed versions for diagnosis.

To intentionally update the companion baseline, change its full commit in
`requirements-dev.txt`, describe the required API/compatibility in the PR, rerun
setup and `check all`, and arrange relevant integration. Do not use a moving branch
in this file. Production Docker builds still deliberately select production
branches; development preparation does not install or execute `gcf-workflows`.

For a coordinated change using an explicit local tools checkout:

```bash
python3.11 scripts/dev.py setup --gcf-tools ../gcf-tools-MY-ISSUE
python3.11 scripts/dev.py check all
```

The runner copies that checkout (including source edits, excluding Git metadata,
common environments/caches and build products), builds its wheel, and records
the source path, Git commit/status and copied-content digest. It does not build
inside or edit the companion checkout. Both editable BFQ checks and clean-wheel
checks use the **same tools snapshot**. Later tools edits require setup again;
they are not picked up silently. Use a separate companion worktree per task and
identify the override in the PR. Omit `--gcf-tools` on the next setup to return to
the declared baseline. Changes to BFQ's `pyproject.toml` or baseline file require
setup again; the runner refuses stale setup metadata.

The runner removes inherited Python import/environment overrides, pip config files
and pytest auto-loaded plugins for predictable checks. Setup honors explicit pip
index, certificate, timeout and retry environment variables; checks disable index
access. Configure private dependency access in your environment, and never paste
authentication into tracked requirements or diagnostics.

## Check profiles

| Command | Evidence | Does not establish |
| --- | --- | --- |
| `python3.11 scripts/dev.py check fast` | `pip check`, non-mutating Ruff lint/format, existing pytest suite (excluding the dedicated scenarios), editable installed-command smoke | Clean wheel packaging or real external processing |
| `python3.11 scripts/dev.py check fast -- tests/test_state.py -q` | Same lint/smoke checks with the explicit pytest selection | Tests outside that selection; report the selection in handoff |
| `python3.11 scripts/dev.py check wheel` | Fresh sdist/wheel build, new venv, offline wheel/dependency install, `pip check`, installed-command smoke | Full pytest suite or scientific workflows |
| `python3.11 scripts/dev.py check all` | Full fast profile, operational scenarios and wheel profile; same command used in CI | Server/container/scientific integration |
| `python3.11 scripts/dev.py check scenarios` | [Synthetic operational scenarios](operational-scenarios.md): real parsing, shared validation, installed manager commands, persistent state, restart/recovery and snapshot retention | Bioinformatics results, real reports, containers or SMTP delivery |
| Relevant [manual integration guides](development-architecture.md#choosing-verification) | Real tools, mounts, reports, workflow outputs, notifications and deployment behavior | Replaced by neither fast nor wheel checks |

Installed-command smoke invokes `bfq`, `fm` and `flowcell-manager` with `--version`
and `--help` using absolute venv executables from a temporary directory outside
the checkout. It checks the shared manager entry point and, for wheel installs,
that BFQ/manager/configmaker imports come from that clean environment. It never
runs bare `bfq` or container wrapper `--help` commands, which can launch work.

Tests named `*_integration.py` in today's suite already test multiple Python
components with temporary data; they remain included in `fast`. The dedicated
`scenarios` profile adds a compact end-to-end operational story;
`all` runs both profiles in separate pytest processes, then checks packaging.
`scenarios` always runs its complete compact set and rejects pytest selections.
A selection supplied to `fast` or `all` changes only the fast test selection.
Neither profile is server integration. See the [fixture and boundary contract](operational-scenarios.md).

## Isolation and mail protection

Each checkout has its own `.dev/`; each check gets a unique `.dev/runs/<unique>/`
for pytest paths, caches, scratch, package builds and the wheel venv. A checkout
lock rejects overlapping setup/check commands **within one checkout**; parallel
independent worktrees are supported. Never copy/symlink `.dev/` between them.
Tests supply temporary config/output/manager paths and mock external commands;
there is no production simulation mode or relaxed production config authority.

The test-only `tests/support/sitecustomize.py` installs a CPython audit hook that
refuses real `smtplib.connect` before sockets open, including SMTP_SSL and custom
ports. The runner loads it through an absolute `PYTHONPATH`, checks its presence,
and propagates it to Python children. `conftest.py` also protects direct pytest
runs and inherited subprocess environments. Tests may replace SMTP with explicit
doubles, as the existing notification tests do. `BFQ_ENV=test` alone is **not** a
mail guard. No guard is shipped in the BFQ wheel or installed in production.

This is an accidental-mail guard for normal Python subprocesses, not an OS
sandbox: `python -I`/`-S`, a replaced environment, non-Python mail tools or custom
socket clients can bypass startup hooks. Do not introduce those execution paths
into local checks. Future subprocess helpers must preserve and verify the guard;
external processing must remain explicitly doubled. The scenario profile adds
a stricter test-only audit guard: only its absolute installed manager executables
may be launched, their child environment must retain the guard, and network
connections and shell/exec/spawn commands are refused. See the scenario guide
for the exact limits. Real notification delivery belongs to separately arranged
server integration.

## Failure evidence and cleanup

Every check prints its owned artifact directory and source/dependency identity.
`identity.json` records current BFQ Git revision/status, setup identity, profile
and platform; `packages.txt` lists installed packages. Each pytest invocation
retains a console transcript, JUnit report and timing/peak-child-memory
summary even on failure. Scenario case directories also retain synthetic inputs,
state, CLI transcripts and external-double outputs. CI uploads identity, test
logs and scenario artifacts on failure for seven days. Include the command, relevant failure and identity in a PR rather than copying all test data.
Local paths and uncommitted status may be meaningful; do not attach sensitive
operational files or credentials.

Successful and failed run/setup artifacts are retained for inspection. They may
be removed once no command is running and their evidence is no longer needed.
Deleting this checkout's `.dev/` discards its disposable environments/caches and
requires setup again; it does not remove source or any other worktree. Do not use
broad cleanup on shared operational paths. No automatic cleanup or retention
policy for production data is introduced here.
