# Holding and releasing projects

`fm` and `flowcell-manager` expose the same persistent scheduling controls:

```bash
fm hold 990101_SYNTHETIC_0001_TEST --project GCF-2099-002 --reason "Awaiting operator review" --dry-run
fm hold 990101_SYNTHETIC_0001_TEST --project GCF-2099-002 --reason "Awaiting operator review"
fm status 990101_SYNTHETIC_0001_TEST
fm search SYNTHETIC --status held
fm release 990101_SYNTHETIC_0001_TEST --project GCF-2099-002 --reason "Review complete"
```

Repeat `--project` to select several projects. A nonempty reason is required and
retained with the request, application time, attempt and requester host. Repeating
the same requested eligibility is idempotent: it creates no additional transition
or processing attempt. `--dry-run` previews the requested eligibility without
writing files. `fm show RUN` includes the complete control history.

Controls require state-backed, known projects. Before the first daemon scan, use
the existing `fm initialize` command to queue the run, then hold projects before
starting processing. If project identities have not been recorded, the command
validates the complete effective input set read-only to identify them. Legacy
inventory alone cannot safely establish pending work; initialize/rerun it through
the existing interface first.

## Stage boundaries and batches

Shared preflight and demultiplexing still include **every** effective project,
including held projects. Holds never remove reads from the combined conversion
plan or change assignment, priority, scientific configuration or statistics.

After demultiplexing, BFQ admits a batch of released unfinished projects. It
finishes that batch's eligible analyses before reporting and archiving it. The
same frozen decimal priority rules apply independently at each stage: `0`, `0.5`,
`1.0` by default, then `2`, with equal numbers tied by project ID. A held priority
`0` project remains ineligible. With no holds, all selected analyses still precede
the first archive.

A release after batch admission joins the next batch, including releases during
analysis, archiving or finalization. BFQ finishes the admitted batch's delivery
phase safely before admitting more analysis. A new hold is evaluated again at
each project start: a project already analyzed can therefore wait before delivery.

Control updates and project starts serialize under the state store's short lock.
The worker records a stage claim **before** project setup or external work. If
the hold wins that race, the stage cannot start. If the claim wins, status shows
`hold_requested` until that stage finishes; the scientific workflow or complete
archive/checksum/snapshot/publication transaction finishes safely. There is one
worker under the existing flowcell execution lease; manager requests neither
start another worker nor interrupt the active one.

## Waiting, completion and recovery

Released A can finish analysis, its Analysis complete notification, archiving,
snapshot publication and finalized notification while B stays held. A keeps its
own completed project row and compatibility-inventory entry. The parent stays
incomplete with status `waiting`, a next boundary and a waiting reason. All-held
work waits without retrying scientific processing or recording a failure.

`fm list` and `fm search` add `eligibility` and `waiting_reason` to their common
project projection. `--status held` selects held unfinished project rows. Running,
failed, interrupted and completed outcomes remain visible even when their future
eligibility is held. Status includes batch membership, next stage, frozen priority,
requested/applied eligibility and actual starts. A completed project's hold can
remain recorded for a later explicit rerun; releasing it never reopens delivery.

Release queues remaining work for the daemon's next discovery scan. Use the
existing daemon wake-up procedure when an immediate scan is needed. Holds and
release requests survive manager, daemon and container restarts. Continuation
uses the same attempt, verified frozen SampleSheet/workbook/run geometry and
captured library-prep configuration. Editing curated originals cannot refresh or
reprioritize that continuation. Completed sibling files, manifests, generation
records, snapshots and delivered notification intents remain intact.

Analysis and reporting durations accumulate successful released batches. Time
spent held is excluded, and synchronous notification delivery remains excluded
from analysis time.

An explicit whole-run or scoped rerun still performs its normal validated
invalidation and creates a new attempt. **Holds carry forward** until explicitly
released. Scoped correction staging is retained while waiting; continuation does
not recreate it. A correction of another project preserves suspended siblings'
original input context and pending work. Completed deliveries remain protected by
the existing publication journal and snapshot recovery. Analysis resume retains
its existing strict workdir/input checks; a hold does not make a missing or invalid
resume context valid.

A failed stage remains failed until the existing rerun/resume recovery is
requested. Hold/release does not clear errors or fabricate success. Notifications
for committed work remain independent: failed or uncertain delivery uses
`fm retry-notifications` with its existing explicit uncertainty acknowledgement.

## Verification

Run `python3.11 scripts/dev.py check all`. The hold regressions use synthetic
metadata, real state/manager/lifecycle/publication and explicit external-processing
and SMTP doubles, including subprocess mail guards. They exercise waiting/release,
priority exclusion, claim/control serialization, release during delivery, scoped
corrections, immutable continuation and retained completed products.

For authorized server integration, use a disposable flowcell with at least two
projects and the configured real test relay. Initialize it, hold B and process A.
Verify complete shared validation/conversion, A's analysis and finalized sends,
reports, archives, manifests, snapshot and completion record while B stays held
and the parent waits. Restart the daemon, release B, wake discovery and verify
only B's remaining stages and notifications run. Compare A's file hashes,
timestamps/inodes, completion/snapshot records and delivered intents before and
after release. Check manifests from the output directory and archive integrity.

Repeat with all projects held, with a hold during analysis, and with B released
during A's archive transaction. Confirm deferred application, stable batch order
and eventual completion. Exercise a controlled stage/publication failure and
supported recovery, preserving failed-run evidence and earlier successes. Inspect
both installed command names, project filters, priorities and history. Record
actual configured SMTP acceptance separately from mailbox receipt; local doubles
do not establish either. Keep operational identifiers/logs in private evidence.
