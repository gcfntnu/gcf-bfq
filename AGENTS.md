# Working on BFQ

Start with [CONTRIBUTING.md](CONTRIBUTING.md), the executable
[development guide](docs/development.md), and the
[architecture and compatibility map](docs/development-architecture.md).
[README.md](README.md) is the operator and production-release reference.

## Task and checkout

- Start bounded issue branches from current `bfq-dev`; feature PRs target
  `bfq-dev`. `master` is production and receives a separate promotion PR.
- Use a separate checkout or worktree for each independent task. Keep its
  environment, cache, generated configuration, state and outputs under that
  checkout's `.dev/` or an invocation-owned temporary directory. Never share a
  mutable working directory between independent agents.
- Read the issue, relevant source/tests and open overlapping work before editing.
  Record scope and affected shared interfaces in the issue or PR. Agree file
  ownership when collaborating; do not overwrite another contributor's work.
- A supplied checkout may already contain work. Inspect its branch, status and
  base before making changes; retain unrelated edits.

## Compatibility and scope

- Preserve current installed commands, input/output formats, configuration
  authority, state/legacy loading and restart behavior unless the issue explicitly
  changes them. `fm` is the same entry point as `flowcell-manager`, not a second CLI.
- Keep curated output-side inputs authoritative. Use the existing state store,
  execution leases and cleanup planning for lifecycle changes; never edit
  operational state JSON or delete outputs to make a test pass.
- Preserve FASTQs and their valid MD5 manifests on downstream reruns, prior
  successful analysis snapshots on failure, and processing completion when
  notification delivery fails.
- BFQ owns orchestration. Shared metadata/libprep rules belong to `gcf-tools`;
  scientific workflows and their tool-image configuration belong to
  `gcf-workflows`. Coordinate interface changes explicitly. Do not extend a BFQ
  issue into companion-repository implementation without task scope to do so.
- Update relevant documentation in the same PR as changed behavior. Keep changes
  bounded; do not bundle unrelated rewrites, formatting sweeps or dependency
  modernization.

## Verification and handoff

- Use `python3.11 scripts/dev.py setup`, then
  `python3.11 scripts/dev.py check fast`; run `check all` for the full local gate.
  Setup acquires dependencies; repeatable checks use the prepared environment.
  See the development guide for a deliberate local `gcf-tools` override.
- Tests must mock SMTP, including subprocesses. `BFQ_ENV=test` does **not**
  prohibit all mail. Never send real mail, start the processing daemon, touch
  production data/configuration or build/push/deploy images as a local check.
- Add meaningful regression coverage for changed behavior and exercise installed
  commands when CLI/packaging changes. Do not mistake mocked external processing
  for validation of scientific outputs or the deployed container.
- State the source/base revisions, changed behavior/interfaces, dependency
  baseline or override, checks and results, limitations and remaining manual
  server checks in the PR. Commit/push at logical checkpoints when authorized.
- Manual operational integration is required before merging affected behavior.
  Merge, production promotion, deployment and issue closure are deliberate steps,
  not effects of the test command. Verify issue closure after an authorized merge.
