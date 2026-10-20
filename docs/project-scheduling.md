# Project priority within a flowcell

Set an optional `Priority` in each effective project SampleSheet's
`[CustomOptions]` section:

```csv
[CustomOptions]
Libprep,Illumina DNA Prep
Priority,0.5
```

Lower values run first. For example, `0`, `0.5`, `1.0`, then `2`. Values must be
finite, non-negative decimals; ordinary decimal and exponent notation are
accepted. Negative values, NaN, infinity and malformed explicit values fail
preflight before cleanup or processing. Missing or blank priority resolves to
`1.0` without rewriting the SampleSheet. Numerically equal values such as `1`,
`1.0` and `1.00` tie, with project ID deciding their order. Values are compared as
exact decimals, independently of floating-point or decimal arithmetic rounding.

Priority applies separately to analysis and archiving. BFQ finishes **all
selected analyses before archiving starts**, then archives in the same priority
order. It continues to run projects sequentially within one flowcell. Priority
does not change workflow selection, scientific parameters or read assignment.
With no priorities, both stages use project-ID order. For a legacy shared input
pair, its one priority applies equally to every project in that pair.

## Preview and frozen order

The existing manager interface shows resolved priorities, whether each is
explicit or defaulted, and the analysis/archive order:

```bash
fm validate 990101_SYNTHETIC_0001_TEST
fm rerun 990101_SYNTHETIC_0001_TEST --from analysis --dry-run
fm status 990101_SYNTHETIC_0001_TEST
```

`fm` and `flowcell-manager` remain the same entry point. Validation previews the
effective inputs. An analysis/demultiplexing rerun previews its proposed new
plan; reporting/finalization retries preview the producing analysis's saved plan.
Status identifies the recorded attempt and actual project starts in each stage.
An attempt's recorded start means processing began, rather than that the project
completed. Existing stage and flowcell outcomes remain authoritative.

Curated output inputs retain their existing precedence. A scoped refresh chooses
instrument inputs only for selected projects, while validating priorities in
**all** effective pairs. Its execution order includes only selected projects.

BFQ freezes the resolved plan with the execution attempt and validated
SampleSheet hashes before launching processing. Editing an original instrument
or output sheet cannot silently reorder that attempt. If execution preflight
finds different effective inputs after a plan was frozen, it fails rather than
changing the plan. Scoped corrections continue using their own frozen inputs.
An explicitly requested new analysis/demultiplexing attempt resolves a new plan.
Analysis resume and reporting/finalization retries retain the producing
analysis's plan. A scoped correction retains the unselected projects' producing
priorities when recording completed analysis, including later whole-run delivery
retries. Historical completed analyses without a recorded plan retain
project-ID ordering on downstream retries; current source edits do not assign
them a new priority retrospectively.

## State, provenance and remaining work

The preflight report adds `summary.project_schedule`, including decimal strings,
explicit/default sources, hashes and planned order. Optional
`attempts[].project_schedule` state records retain the frozen plan, originating
attempt and actual analysis/archive starts. Completed analysis metadata retains
the plan for downstream recovery. The reporting configuration snapshot includes
it; project-paired analysis input provenance also includes each project's
scheduling values. Existing state without these optional records remains usable.
Use the manager commands to inspect state.

This is the priority increment of issue #160. Analysis and finalized notifications
still use the existing flowcell-level timing. Per-project analysis notification,
durable progress and finalization notification are remaining #160 work. Persistent
hold/release controls belong to #161 and depend on those completion boundaries.
Priority alone does not hold a project or allow archiving ahead of other pending
analyses.

## Verification

`tests/test_project_scheduling.py` uses generated metadata, the real validator,
manager/state APIs, analysis/archive loops, checksums and snapshots. Scientific
processing, archive generation and SMTP are explicit doubles. It checks decimal
ties/precision, missing/blank defaults, invalid priority before a selected refresh,
legacy shared pairs, frozen inputs, actual stage order and downstream retry.
The standard `scripts/dev.py check all` gate also runs existing restart,
correction, notification, installed-command and packaging checks.

Before merging operational changes, use disposable integration inputs with
priorities `0`, `0.5`, blank/missing and `2`. Compare manager preview, actual
workflow/archive start order and retained state/provenance. Check that all
analyses finish before the first archive starts. Retry finalization after editing
an original sheet and verify that the completed analysis's order, FASTQs,
manifests and snapshots are retained. Record source/image identities and actual
observations in the implementing PR. Local synthetic checks do not establish
scientific output correctness or actual notification delivery.
