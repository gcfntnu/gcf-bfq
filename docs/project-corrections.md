# Correcting selected projects

`fm` and `flowcell-manager` provide the same interface. A repeated `--project`
selector limits a rerun to those projects. Omission retains the whole-run restart
behavior. The examples below use synthetic project and run identifiers.

For a metadata or workflow correction that retains FASTQs:

```bash
fm rerun 260101_TEST_0001_FLOWCELL --from analysis \
  --project GCF-2026-002 --refresh-inputs --dry-run
fm rerun 260101_TEST_0001_FLOWCELL --from analysis \
  --project GCF-2026-002 --refresh-inputs --reason 'corrected project metadata'
```

Use `--from demultiplexing` when the selected project's FASTQs must be regenerated.
Repeat `--project` to widen the correction scope. The operator chooses this scope;
BFQ does not require proof of historical read-assignment or FASTQ-content
equivalence. Use the dry-run preview to review selected and retained projects,
input replacements, restart stage and affected products before queueing.

Curated output-side inputs remain authoritative. `--refresh-inputs --project`
chooses instrument inputs only for selected projects. All effective pairs,
including retained projects, must pass pairwise and complete-flowcell validation
before queueing and again during processing. Changes to an unselected instrument
pair do not replace its curated counterpart. A legacy shared pair can still be
used without renaming; refreshing only part of that shared pair requires
project-qualified pairs or selection of all its projects.

Scoped corrections require BFQ state. Inventory-only legacy runs must first be
initialized or rerun through the existing whole-run lifecycle. `--resume` and the
index-toggle options cannot be combined with `--project`; edit the selected
input pair and queue a staged correction. `--force` confirms the operation
noninteractively and retains validation and execution-lease checks.

## Generation and publication

Queueing freezes a copy of every effective input pair and available run-geometry
XML under `provenance/project-corrections/<correction-id>/inputs/`. Input hashes,
selected projects and restart boundary are recorded in canonical state. A later
instrument edit cannot change the queued attempt's inputs. Processing uses a
separate directory for each attempt and separate analysis workdirs.

Analysis-only corrections copy the selected project's FASTQs and existing FASTQ
manifest into staging, retaining their bytes and modification times. A combined
demultiplexer may regenerate other projects in staging; publication installs only
selected project products. Unselected canonical FASTQs, manifests, reports,
archives and successful snapshots remain in place throughout.

Publication runs under the existing flowcell execution lease. BFQ journals each
selected replacement, retains the prior products under the correction's
`previous/` directory and then commits processing completion through the state
store. Recovery restores previous products if processing completion was not
committed. After committed completion it retains the replacements. Snapshot
candidates are prepared before this publication and become current only after
the completion commit. Notification delivery remains independent of processing
completion.

Canonical shared conversion statistics remain associated with their original
attempt. Repeated-conversion statistics and generated inputs remain in the
correction attempt directory; BFQ does not add repeated conversions to the
original flowcell totals. State retains per-project successful completion and
workflow metadata, so a retained project's history is distinguishable from the
current correction's stage status. No automatic retention cleanup is introduced.

## Verification before merge

The automated `tests/test_project_corrections.py` scenarios use synthetic files
and mock validation boundaries. They check scoped preview, frozen inputs,
publication and interruption recovery, selected analysis FASTQ/manifest retention,
and preservation of unselected bytes, modification times and file identities.
Existing state, installed-command, notification and snapshot suites cover the
surrounding lifecycle. They do not establish scientific workflow correctness.

On disposable facility data, run two projects with distinct library preparations
and workflows to completion. Record hashes and modification times for the project
that will be retained, including FASTQs, manifests, reports, archives and snapshot.
Correct only the other project's instrument inputs, inspect the scoped dry-run,
and run an analysis correction followed by a demultiplexing correction. Verify the
retained files and analysis workdir remain unchanged, corrected deliveries contain
their own effective inputs/provenance, and shared conversion totals are not
double-counted. Interrupt correction generation and publication separately and
verify prior successful products remain recoverable. Test notification failure
with an arranged mail double or controlled relay and confirm processing remains
completed. Record actual image/source revisions and observed results separately
from the local synthetic checks.
