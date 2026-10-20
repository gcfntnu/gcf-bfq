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
selected eligible analyses before archiving starts**, then archives in the same priority
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
`fm list` and `fm search` expose one project per flowcell row with its own
analysis/delivery outcomes and recorded priority/source; see the
[project listing contract](project-listing.md). That projection is the foundation
for [persistent hold/release controls](project-holds.md).
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

## Project progress and notifications

After each project's workflow and required project reports succeed, BFQ commits
its analysis outcome and an Analysis complete notification intent together. It
attempts that notification before starting the next project's workflow. Each
generation retains its own report and effective-input copies under
`provenance/project-analysis/`, so a later failure or replacement cannot change
the metadata or attachments of an earlier email. Analysis complete describes
analysis readiness; delivery follows in the archiving phase.

During archiving BFQ finishes one project's delivery archives, checksums,
successful retained snapshot and publication, commits its completion record,
and attempts its finalized notification before archiving the next project.
Status shows each project's analysis and delivery outcomes while the flowcell
remains incomplete. Mail failures preserve processing outcomes and allow other
projects to continue. Use `fm retry-notifications` for failed delivery;
uncertain acceptance also requires the existing `--retry-uncertain` acknowledgement.
Known delivered intents are not resent automatically.
Flowcell analysis timings exclude these synchronous mail attempts; archiving and
archive-checksum timings accumulate the successful project operations.

If B's analysis fails after A succeeds, a scoped analysis rerun of B preserves
A's workdir, reports, generation and notification. Once B's selected delivery
commits, any other successfully analyzed projects still awaiting delivery are
queued for their remaining finalization. A flowcell with unfinished failed
analysis still requires explicit recovery. [Held projects](project-holds.md) wait
without preventing released projects from reaching delivery.

If B's delivery fails after A commits, a reporting/finalization retry retains A's
archives, checksums, snapshot, completion record and known delivered intent.
For a staged correction with committed analysis, the same command reuses that
frozen analysis/read generation and retries only unfinished delivery products.
The preview identifies staged retries and retained completed projects. An
explicit new analysis/demultiplexing rerun instead deliberately invalidates the
selected projects' downstream outcomes and creates new generations. Omitting
`--project` applies that scientific restart to the whole run. Analysis resume
retains the established repaired-workdir/Snakemake behavior.

## State and provenance

The preflight report adds `summary.project_schedule`, including decimal strings,
explicit/default sources, hashes and planned order. Optional
`attempts[].project_schedule` state records retain the frozen plan, originating
attempt and actual analysis/archive starts. Completed analysis metadata retains
the plan for downstream recovery. The reporting configuration snapshot includes
it; project-paired analysis input provenance also includes each project's
scheduling values. Existing state without these optional records remains usable.
Use the manager commands to inspect state.

Optional `project_progress` records bind each analysis and finalization outcome
to its generation, inputs/workdir and notification identity. New notification
identities include the project (`processed:ATTEMPT:PROJECT` and
`finalized:ATTEMPT:PROJECT`). Existing combined notification records keep their
identities and retry rules; loading historical state creates no new analysis
emails. Legacy execution also captures its effective shared pair without
rewriting authoritative inputs. Priority alone does not hold work or bypass the
analysis phase; persistent [hold/release controls](project-holds.md) set eligibility.

## Verification

`tests/test_project_scheduling.py` uses generated metadata, the real validator,
manager/state APIs, analysis/archive loops, checksums and snapshots. Scientific
processing, archive generation and SMTP are explicit doubles. It checks decimal
ties/precision, missing/blank defaults, invalid priority before a selected refresh,
legacy shared pairs, frozen inputs, actual stage order and downstream retry.
The standard `scripts/dev.py check all` gate also runs existing restart,
correction, notification, installed-command and packaging checks.
`tests/test_project_lifecycle.py` additionally covers actual project boundaries,
later-project failures, worker interruption, before/after-commit publication
faults, partial delivery retries, immutable notification generations and
failed/uncertain SMTP handling. Its external processing and mail use doubles.

Before merging operational changes, use disposable integration inputs with
priorities `0`, `0.5`, blank/missing and `2`. Compare manager preview, actual
workflow/archive start order and retained state/provenance. Check that all
analyses finish before the first archive starts. Retry finalization after editing
an original sheet and verify that the completed analysis's order, FASTQs,
manifests and snapshots are retained. Record source/image identities and actual
observations in the implementing PR. Local synthetic checks do not establish
scientific output correctness or actual notification delivery.
For notification acceptance on an authorized test server, inspect the actual
Analysis complete send before the next workflow starts, and finalized send after
that project's commit before the next archive. Exercise a later-project failure,
staged delivery retry and restart; verify earlier files, completion records and
known delivered intents remain intact. Verify subjects, project scope,
attachments and actual configured test-relay acceptance separately from mailbox
receipt. Never use the local SMTP guard/doubles in an authorized real-mail server
integration or edit state JSON to manufacture success.
