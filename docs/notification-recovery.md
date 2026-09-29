# Notification recovery verification

Issue: #122. These checks validate notification recovery without requiring a
failed scientific workflow. No production data changes or real email were used
in automated tests. Use an approved test run/recipient on the server.

## Server checks

1. Build/install the issue branch and verify:

   ```console
   flowcell-manager retry-notifications --help
   ```

2. In the test instance's `/config/bcl2fastq.ini`, keep valid `host`,
   `from_address` and `finished_to`, but rename `error_to` to `errorTo`.
   Start one ordinary test flowcell. Expected: reports, archives and checksums
   are produced; processing completes; finalization notification is `failed`
   with an instruction to use `error_to`. Processing-complete mail can succeed.

   ```console
   flowcell-manager status RUN_ID
   flowcell-manager show RUN_ID
   ```

3. Record the timestamps of representative FASTQs, reports and archives. Correct
   the setting to `error_to`, using an approved test recipient, then run:

   ```console
   flowcell-manager retry-notifications RUN_ID --kind finalized
   ```

   Expected: one finalization message, notification `sent`, unchanged processing
   attempt/stages and output timestamps. No analysis, reporting, archive or
   checksum commands should run. The CLI loads corrected settings immediately;
   an already-running daemon needs its usual restart to reload static config.

4. Repeat the same retry command. Expected: success/already-delivered and no
   additional mail. Restart the daemon and check again: no duplicate completion
   mail or reprocessing.

5. Optional: during a normal finalization-only rerun, use an unavailable SMTP
   hostname in the test configuration. The new finalization notification should
   fail without reverting processing completion. Restore the hostname and retry.
   A finalization-only rerun must retain the existing processed notification.

Uncertain acceptance, interrupted writes and concurrency are covered by automated
fault injection. They do not require deliberately crashing production SMTP.
If a real uncertain outcome occurs, inspect it and use `--retry-uncertain` only
when a possible duplicate is acceptable. Never edit state JSON to fake recovery.

### Legacy output-path recovery

Server verification found an existing legacy state with `output_path` equal to
the bare flowcell ID. BFQ used the configured absolute output root, while cleanup
looked relative to the operator's current directory and incorrectly showed
`(none)`. Bare IDs now resolve consistently under `outputDir`; manager and daemon
reject conflicting or unavailable paths instead of silently selecting another
directory. Preparation and processing share the run execution lease.

After installing the updated branch, preview the affected run:

```console
flowcell-manager rerun 260918_MN00686_0026_A000HCMFHF --from finalization --dry-run
```

Expected: output directory `/mnt/bfq/output/260918_MN00686_0026_A000HCMFHF`,
both project/QC `.7za` archives and their `md5sum_*_archive.txt` files listed.
FASTQ checksum manifests are preserved. Preview and declining confirmation do
not change files or state. If a finalization rerun is wanted, repeat without
`--dry-run` and confirm: the absolute path is recorded atomically with notification
invalidation before cleanup. BFQ then rebuilds only finalization products.
No manual JSON edit is needed. This is a processing rerun, so its new finalization
notification is intentional; `retry-notifications` still performs no cleanup.

The missing-setting, corrected-setting retry and duplicate-suppression scenario
above passed server verification on 2026-09-29. The path-recovery scenario remains
the server check for the follow-up fix. Regression tests include a same-named
directory under the working directory to prove it cannot be selected for cleanup.

## Implementation decisions and difficulties

- Processing completion and notification intent share one atomic JSON mutation.
  There is no crash window between those two records. Composition runs afterward.
- Persist a sending claim before composition/SMTP. A crash can therefore leave an
  uncertain outcome even if no bytes reached the relay. This conservative policy
  avoids automatic duplicate sends; explicit retry remains available.
- State locks serialize delivery records; the execution lease also excludes
  rerun/archive cleanup. Review exposed that the old `--force` live-run override
  could bypass that protection, so it cannot override an active lease anymore.
- Directory fsync may fail after replacing a state file. If recording processing
  completion fails and the stage is failed, its dependent notification is
  superseded. Delivery also verifies that its producing stage is complete.
  Notification-result persistence failures, in contrast, never relabel completed
  processing; they leave an uncertain claim for operator inspection.
- A preparation failure before a new attempt starts leaves prior successful
  attempt history intact. Only a failure after committing the current completion
  can change that current attempt from completed to failed.
- Legacy v1 state has no invented completion intents. Existing error-report
  routing/signature suppression remains separate and unchanged.
- Removed the old catch-all 10x/Parse resend after any email exception: the first
  SMTP attempt might already have been accepted. No implicit second message is
  sent. Partial refusals and interrupted SMTP transactions are uncertain.
- Independent review found the installed CLI importing the neighboring
  `configmaker.py` script instead of the package, plus import-time log creation
  in unwritable directories. Mail dependencies are loaded lazily and the manager
  sanitizes its executable import path. Installed help/status/finalized retry are
  tested from `/proc`; processed composition still uses the existing gcf-tools
  parser, whose logging side effect is addressed by the separate validation work.
- Malformed InterOp CSV with no blank header separator previously looped at EOF.
  It now returns an unavailable-metrics summary instead of hanging notification
  composition. Optional instrument disk statistics also tolerate an absent mount.

## Automated verification scope

Tests cover config errors, malformed recipients, unknown keys, composition and
attachment errors, connection failures, recipient/data rejection, partial
acceptance, send timeout/disconnect, QUIT failure, and stable retry identities.
State tests inject failures before and after atomic replacement, before sending
and after acceptance; exercise interrupted claims, concurrent retry and cleanup,
superseding outputs, legacy loading and retention of genuine processing failures.
Daemon/CLI tests verify saved context, corrected settings, output bytes/mtimes,
completion-before-SMTP, one automatic attempt, restart recovery and meaningful
exit codes. Production integration and actual relay receipt remain server checks.
