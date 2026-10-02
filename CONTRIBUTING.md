# Contributing to BFQ

The goal is for a developer and their agent to complete a bounded issue through
local implementation and verification without another developer's chat history.
Operational integration, review and production promotion remain separate steps.

Use these references:

| Need | Reference |
| --- | --- |
| Install an isolated environment and run the checks | [Development guide](docs/development.md) |
| Find code responsibilities and compatibility boundaries | [Architecture map](docs/development-architecture.md) |
| Repository-specific agent instructions | [AGENTS.md](AGENTS.md) |
| Operate BFQ, choose restart boundaries, build or promote production | [README](README.md) |
| Understand the BFQ 2 operational baseline | [BFQ 2](docs/bfq2.md) |

## Start a bounded task

Read the issue and identify its intended behavior, relevant source/tests and
compatibility constraints. Check existing PRs for overlapping files or interfaces.
Record ownership and any dependency on another task in the issue or PR; a branch
name alone is not a coordination record. Several closely related issues may
share a branch when the task explicitly calls for one combined implementation.

Create a feature branch from current `origin/bfq-dev`, with a unique name that
references the issue. Do not implement features directly on either long-lived
branch. For an existing clone, an isolated worktree avoids changing its checkout:

```bash
git fetch origin
git worktree add -b issues/NUMBER-short-description ../bfq-NUMBER origin/bfq-dev
cd ../bfq-NUMBER
python3.11 scripts/dev.py setup
python3.11 scripts/dev.py check fast
```

Replace `NUMBER` and `short-description` with the task's actual values. Separate
clones work equally well. Confirm `git status --short` and `git rev-parse HEAD`
before editing an existing workspace. Do not reset or clean another task's work.

The supported local loop is documented in [development.md](docs/development.md).
It uses a checkout-owned `.dev/` directory. Activation, shell startup changes,
global Python installs, instrument mounts, Docker and `/opt` setup are not part
of that loop. Authentication for cloning/pushing is personal setup: use your
normal Git credentials or the available GitHub connector. Do not commit tokens,
credentials, internal operational configuration or customer data. The workflow
does not require a particular coding agent, hosted session or personal skill.

## Work concurrently

Each independent task gets its own branch, checkout/worktree, `.dev/` and test
paths. Git objects may be shared by worktrees; mutable environments, configuration,
state, outputs and caches may not be shared between tasks. Use `tmp_path` or a
uniquely owned invocation directory for tests. Avoid hard-coded paths even when
they happen to be unused on one developer's machine.

Check especially for overlap in `state.py`, stage orchestration, CLI entry points,
configuration selection, dependency declarations and fixture helpers. If two
tasks need the same interface, agree its contract and implementation order in
the issues/PRs. Rebase or merge the completed prerequisite deliberately, then
repeat the affected checks. Do not silently take over another contributor's file
or rely on an unrecorded chat decision. Setup copies an explicit local
`gcf-tools` override and builds its wheel; checks use that snapshot, not a live
editable companion checkout. Use a dedicated companion checkout for concurrent
coordinated work, rerun setup after its edits and report its revision/dirty state.

A simple rehearsal is two independent worktrees at the same commit, each running
the documented setup and checks. Each should report its own `.dev/` paths and
source/dependency identity, without writing in the other's tree. This verifies
development isolation; it does not authorize two production daemons to share
operational state or scratch space.

## Implement and verify

Preserve the [compatibility contracts](docs/development-architecture.md#compatibility-contracts)
unless the issue explicitly changes them. Prefer small changes with observable
benefits. Reuse BFQ's lifecycle, validation and notification mechanisms instead of
introducing parallel paths. A shared-domain change may need a coordinated
`gcf-tools` issue; a scientific workflow change belongs to `gcf-workflows`.

Run focused existing tests while implementing; add regression coverage where it
proves changed behavior. Use the common check command before handoff:

```bash
python3.11 scripts/dev.py check fast
python3.11 scripts/dev.py check all
```

`fast` covers the existing automated suite, lint/format checks and installed
editable-command smoke checks. `all` additionally runs the compact operational
scenarios and checks the built wheel in a separate environment. Run
`python3.11 scripts/dev.py check scenarios` for the operational profile alone.
The [development guide](docs/development.md) defines the exact profiles,
dependency preparation and diagnostics. Keep lint and formatting
checks non-mutating; make intentional formatting edits as normal source changes.

Tests must not send real email. Use mocked SMTP and the development harness's
subprocess protection; setting `BFQ_ENV=test` alone is insufficient. Local checks
must not run real demultiplexers, Snakemake workflows or containers, or change
production mounts, configuration or state. Do not run bare `bfq` as a smoke test:
it starts the processing service. Installed `bfq --version` and manager help/version
checks are safe entry points exercised by the harness.

Existing tests named `*_integration.py` exercise multiple BFQ components using
temporary files and doubles at external boundaries. They are part of the local
suite, not proof of a server integration pass. The
[operational scenario guide](docs/operational-scenarios.md) describes the
dedicated installed-command/recovery profile, its deterministic synthetic data,
retained artifacts, and external boundaries.

For operational behavior changes, select the applicable manual guide(s) from the
architecture map and specify remaining real-server checks. A maintainer runs
those checks on disposable operational data and inspects actual reports/results
before merge. Real notification delivery belongs to that explicitly arranged
integration test. Documentation/tooling-only work should report the local evidence
and whether any operational behavior changed, rather than invent a scientific
integration result.

## PR and handoff

Commit at logical checkpoints and push the issue branch so work survives a local
session. Open the feature PR against **`bfq-dev`**, explaining the operational
problem and resulting behavior first. A short handoff should include:

- Issue(s), intended outcome, starting base commit and current source commit.
- Compatibility constraints, affected files/interfaces and decisions relevant to
  concurrent work.
- Dependency baseline or explicit override; exact check commands and results.
- Remaining manual integration checks, known limitations and specific unresolved
  questions. Separate observed failures from checks that were not run.

Update documentation in the same PR as behavior it describes. Use `Closes #NUMBER`
for each issue the PR actually completes; use `Refs #NUMBER` for follow-up work.
After merging, verify that the linked issues closed. GitHub closing-keyword
behavior depends on the target/default branch, so an integration/promotion merge
must not be assumed to have handled every issue automatically. If a completed
issue remains open, close it deliberately with the implementing PR and validation
evidence once authorized.

Feature review and relevant manual integration precede merge. Promotion from
`bfq-dev` to **`master`** is a separate deliberate PR, using a regular merge to
preserve the shared history. Follow the existing [production policy](README.md#versions-and-production-updates)
for image tags and builds; a successful local check does not merge, promote,
deploy or send notifications. No new release scheme or per-PR release manifest
is required.
