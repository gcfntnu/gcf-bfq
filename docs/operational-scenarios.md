# Small operational scenarios

These scenarios exercise BFQ's operational boundaries using invented data and
explicit external doubles. They help a developer verify a bounded change before
manual server integration. They do not validate biological results or authorize
merge, promotion or deployment.

## Run from a clean checkout

Use Linux and CPython 3.11, as described in the [development guide](development.md):

```bash
python3.11 scripts/dev.py setup
python3.11 scripts/dev.py check scenarios
python3.11 scripts/dev.py check all
```

Setup acquires dependencies, including the immutable `gcf-tools` revision in
`requirements-dev.txt`. A deliberate `setup --gcf-tools PATH` snapshots a local
companion checkout exactly as in the other profiles. Scenario execution itself
needs no downloads, instrument mounts, credentials, bioinformatics applications,
containers or SMTP relay. No companion testdata refurbishment is a prerequisite.

`scenarios` runs the complete compact pytest file, using the installed editable
BFQ and manager commands. It rejects additional pytest arguments to avoid silently
changing the profile. `fast` runs the existing detailed regression suite, lint,
format and installed-command smoke checks; it excludes this new scenario file.
`all` runs fast checks, scenarios in a separate process, and the existing clean
wheel check. CI uses that same `check all` command. The wheel check verifies clean
installation and command entry points; the operational profile uses editable BFQ.

## Runtime and resources

Measured on Linux x86_64, CPython 3.11.16 and the pinned companion on 2026-10-02:
all four scenarios passed in **15.6 seconds**, retaining **1.2 MiB** of case
artifacts. Two concurrent independent runs took **20.6/20.8 seconds** with
about **102 MiB** peak child RSS each; CI scenarios took **21.6 seconds**. The
enclosing full check reported a peak child RSS of **120 MiB**
(the summary reports the maximum across children so far, not total simultaneous
memory). These are observations, not performance assertions.

Budget one CPU, 512 MiB RAM and a few MiB per scenario invocation, plus the
prepared Python environment. This setup used about 536 MiB including acquired
wheels/cache/temp files (252 MiB venv, 48 MiB reusable wheelhouse); allow at least
1 GiB free per checkout for setup and the additional clean-wheel environment.
Retained invocations accumulate. Setup time depends on dependency acquisition;
normal scenarios do not download anything. The shared CI job has a ten-minute
limit for setup, all checks and packaging.

## Fixture contract

`tests/support/operational.py` generates the fixture beneath a unique pytest case
root. It reuses only the small config/input helpers extracted to
`tests/support/fixtures.py`; existing tests keep their helper imports. There is
no second test framework or production simulation mode.

All identifiers and sequences are invented. The synthetic run is
`260101_SYNTHETIC_0001_ASYNTHETIC`, with project `GCF-2026-001` and sample `sample`.
The fixture contains one paired-end sample with two 12-base read pairs, matching
mate IDs and order, valid four-line records, matching sequence/quality lengths,
and compressed R1/R2 files. Gzip headers omit filenames and timestamps; workbook
ZIP entries and core timestamps are normalized. Repeated generation is checked
for byte equality. Generated operational state, logs and reports include normal
runtime timestamps and absolute workspace paths; those are not byte deterministic.

Metadata follows the existing accepted structures: SampleSheet `[CustomOptions]`
and `[Data]` sections with Sample_ID/Sample_Project/index, and an XLSX submission
form with customer headers at row 15 plus the `INFO (GCF-lab only)` sheet. Both
refer to the same sample. Instrument and curated copies deliberately differ in
User/Sample Group so authority and preservation are observable. Tests check mate
structure, shared validation results and expected sample/project identity.

Each case has separate instrument roots (`nova/` and `ekista/`), `output/`,
`manager/`, `scratch/`, `reports/`, `logs/`, `analysis/`, `cache/` and `config/`.
A local INI points only within that case root. The test-only Python startup hook
loads it for the actual installed `fm` and `flowcell-manager` executables, invoked
by absolute path outside the repository working directory. Production continues
to load its normal configuration; no production environment override is added.
The authoritative libprep path is redirected only by pytest monkeypatch in the
harness, to a synthetic local config.

The fixture is intentionally small: it is not an Illumina BCL run or instrument
model, has no meaningful biology, and does not cover every kit, lane layout,
submission column, report version or legacy-state variant. Detailed existing
preflight, restart, checksum, notification and snapshot tests retain that broader
regression coverage. Any future compatible `gcf-tools` fixture can reuse these
conventions without making BFQ depend on another general-purpose testdata system.

## What is real and what is doubled

| Boundary | Executed behavior |
| --- | --- |
| Metadata and preparation | Real BFQ parsing/input selection and shared `gcf-tools` validator; errors retained in preflight JSON and state |
| Installed manager | Real argument parsing, exit codes, initialization, rerun preview/application, show/status/search and notification retry |
| State and recovery | Real file-backed store, attempts, leases, cleanup plans, FASTQ hashes/manifests, completion and delivery intent records |
| Orchestration | One real queued attempt through the daemon's existing core; the daemon scan loop is never started |
| Analysis retention | Real workdir identification and tar snapshot preparation/publication/recovery, containing the synthetic workflow files |
| Installed workflow revision | Declared `synthetic-unexecuted` provenance double avoids probing a facility `/opt/gcf-workflows`; BFQ package version remains real |
| Configmaker process / scientific workflow | Declared doubles write labeled config, Snakefile, tiny report and sample-info output; the scientific programs never run |
| Report generation / delivery archive | Declared report double and 7za double; `.7za` files are labeled text placeholders, not usable 7z archives |
| Archive checksum process | Double hashes the actual placeholder bytes in Python; FASTQ checksum handling remains real BFQ |
| SMTP | Relay double accepts or raises a controlled failure; real message composition and durable delivery/retry logic remain active. Accepted messages are retained as `.eml` files |

Config, workflow, report and archive placeholders carry
`SYNTHETIC EXTERNAL DOUBLE: not a scientific result`. Tiny sample-info and
MultiQC config files keep the compatible structure required by BFQ. `logs/external-doubles.jsonl` identifies invoked boundaries and commands.
Do not use a successful local run as evidence for scientific results, real
MultiQC/report compatibility or email delivery.

## Expected scenario outcomes

The profile passes four tests; failures include the relevant command and retained
case path. The scenario set is deliberately compact:

| Scenario | Expected result |
| --- | --- |
| Fixture and valid preflight | Paired reads and deterministic input bytes verified; shared validator accepts one planned sample; installed manager validation agrees and leaves inputs/state unchanged |
| Invalid preflight | Missing metadata produces `sample.missing_metadata`; initialization and rerun are refused without changing protected files/state. A post-preparation input edit is revalidated, fails the analysis stage, retains diagnostics, and invokes no processing or delivery double (only the synthetic revision lookup) |
| Restored FASTQs and rerun | Initialization queues analysis with demultiplexing complete; one attempt completes. Both installed names agree on show/status/search/list. Preview is read-only; explicit analysis rerun removes downstream products/workdir while preserving curated inputs, FASTQs, manifests, Stats and the successful snapshot. A finalization rerun preserves completed analysis, workdir and snapshot |
| Failure and recovery | A controlled scientific-workflow exit 23 persists failure/error evidence and retains the prior successful snapshot. Explicit rerun succeeds despite SMTP failure; installed notification retry completes delivery without rerunning processing or rewriting products |

Detailed existing tests continue to cover demultiplexing, reporting and
finalization invalidation, explicit input refresh, legacy state, staged snapshot
publication and uncertain notification delivery. See `test_input_preflight.py`,
`test_state_integration.py`, `test_fastq_checksums.py`,
`test_notification_integration.py` and `test_snapshot_integration.py`; `check all`
runs them along with these scenarios instead of duplicating their full matrices.

## Safety and diagnostics

All tests inherit the existing SMTP audit prohibition, including SMTP_SSL,
nonstandard ports and reloaded smtplib. `BFQ_ENV=test` alone is **not** a blocker.
The operational profile additionally uses a strict executable allowlist: only
its absolute installed `fm` and `flowcell-manager` paths are allowed through
subprocess. Shells, external applications, exec/spawn and network connections/DNS
are refused before execution. Children must preserve the guard environment and
startup path. Configuration injection errors exit before the manager can fall
back to production configuration. Negative regression tests verify these failures.

The hooks are test-only Python safeguards, not an OS sandbox against hostile
code, native extensions or a deliberately disabled interpreter environment. Do
not add environment-clearing/native execution paths, allow scientific tools, or
run bare `bfq` to get a test to pass. New doubles must be declared and must not
replace the behavior under test.

The command prints its owned `.dev/runs/run-…/` directory. Both successful and
failed artifacts are retained. Inspect:

- `identity.json` and `packages.txt`: source, setup/companion and package versions.
- `scenarios.log`, `scenarios-junit.xml`, `scenarios-summary.json`: assertions,
  captured failures, elapsed time and peak child RSS observed so far in the check.
- `scenarios/<case>/…`: synthetic inputs, output, persistent state, reports and
  command transcripts (`logs/commands.jsonl`).
- `logs/external-doubles.jsonl`, workflow logs and captured `.eml` messages:
  which doubles ran and what they produced.

An intentional workflow or SMTP failure inside a recovery test is expected;
the profile passes only after asserting the failure and successful recovery.
Unexpected failures return nonzero and point to the retained evidence. CI uploads
check diagnostics and synthetic scenario artifacts on failure for seven days.
No actual facility data or credentials belong in these artifacts.

After the invocation has exited, remove only its printed run directory when no
longer needed. There is no automatic cleanup of other invocation roots. Do not
reuse a pytest `--basetemp` for concurrent tests; the supported runner allocates
it uniquely. Removing a checkout's entire `.dev/` also removes its setup and
requires setup again.

## Concurrent development

The existing checkout lock rejects overlapping setup/check commands in one
checkout. Use two independent worktrees and run setup in each; never share or
symlink `.dev/` or move a prepared venv. From a repository with the changes:

```bash
git worktree add --detach ../bfq-scenario-a HEAD
git worktree add --detach ../bfq-scenario-b HEAD
(cd ../bfq-scenario-a && python3.11 scripts/dev.py setup)
(cd ../bfq-scenario-b && python3.11 scripts/dev.py setup)
(cd ../bfq-scenario-a && python3.11 scripts/dev.py check scenarios > scenario.log 2>&1) &
pid_a=$!
(cd ../bfq-scenario-b && python3.11 scripts/dev.py check scenarios > scenario.log 2>&1) &
pid_b=$!
wait "$pid_a"
wait "$pid_b"
```

Each invocation exercises controlled failure/recovery inside its own case roots.
Inspect each log's printed artifact path and outcome independently. Removal of
one completed invocation's exact run directory must leave the other workspace's
inputs, state, outputs and reports unchanged. This is development isolation, not
permission for concurrent production daemons to share state or scratch.

## Verification record (2026-10-02)

The issue branch started at `7e934f864b64b0af636274aa0d9a52de9fb75543`
(merged #141). Implementation checkpoint
`76ecd27c1587f88a0a5d4400b699d9728c347c32` was verified with the pinned
`gcf-tools` baseline `bdd94a1d944120ac8a096ea8876c0de98dd51ab3`:

- Fresh setup and `check all`: **711 fast/guard tests and four scenarios passed**,
  Ruff lint/format, sdist/wheel build, offline clean-wheel installation and all
  installed command checks passed. [CI run 36975680629](https://github.com/gcfntnu/gcf-bfq/actions/runs/36975680629)
  passed the same complete command from a clean checkout.
- Two detached worktrees at that checkpoint independently ran setup, then
  `check scenarios` concurrently; both passed, using distinct venvs, wheelhouses,
  caches and run roots. Same-checkout overlap was rejected by the existing lock.
- A second concurrent pair appended one deliberately failing test only in
  disposable worktree A: A returned 1 with the injected diagnostic in its log and
  JUnit report; B returned 0. A's failed invocation was archived for inspection
  before its exact printed run directory was removed. All **282 files** in B's
  retained runs kept identical content hashes and modification times; A's earlier
  successful invocation and both setups remained present. The injected test was
  removed, leaving both source worktrees clean.
- Direct pytest with startup `PYTHONPATH` unset and an independent review checked
  harness startup. No production behavior defect was found. The host-dependent
  `/opt/gcf-workflows` revision lookup is explicitly doubled, as described above.

The failure probe tests development diagnostics/isolation; it is not a production
failure or a change to the four committed scenarios. The local failure archive,
transcripts and machine-readable comparison remain under the implementing
checkout's `.dev/isolation-evidence/`; routine scenario artifacts are retained
by the documented runner. PR #142 records subsequent final-head CI evidence.

## Remaining manual integration

Before merging an operational behavior change, select the applicable existing
[server guides](development-architecture.md#choosing-verification) and explicitly
arrange checks with disposable operational data:

- Real bcl-convert/bcl2fastq conversion, instrument metadata, index orientation,
  FASTQ naming/pairing and demultiplexing/checksum restart boundaries.
- Apptainer images, bindings, permissions, scratch/cache space and subprocess
  failures in the deployed environment.
- Actual configmaker and scientific workflows for representative kits/projects,
  effective libprep/workflow snapshot selection, outputs and retained provenance.
- Sequencing and MultiQC report compatibility, attachments and operator-readable
  content; no automatic QC acceptance gate is introduced.
- Deliberately arranged real SMTP recipients/delivery/retry, including outage and
  uncertain-delivery handling. No local command sends that mail.
- Deployment configuration/mount authority, discovery/daemon restart, legacy data,
  resource behavior and production image identity. Promotion and deployment
  remain separate human decisions.
