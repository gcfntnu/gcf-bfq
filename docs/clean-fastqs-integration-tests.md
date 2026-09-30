# Issue #126: clean-fastqs integration checks

Branch: `126-clean-fastqs`; PR target: `bfq-dev`. Do not merge before operator
integration testing and manual verification.

Use a disposable, completed test flowcell with its normal delivery products, or
a delivered run for which operational copies are confirmed. No delivery tracking,
remote-copy verification or archive integrity scan has been added. The command
requires an existing state record with completed processing/finalization; legacy
inventory alone is rejected. Expected archives use the current `archive_worker`
naming convention, not arbitrary `*.7za` matches or QC archives.

## Preview and prerequisites

1. Run `flowcell-manager clean-fastqs RUN_ID --dry-run` from outside the output
   directory. Repeat with the absolute output path and trailing slash. Confirm
   both previews resolve to the same configured output directory.
2. Confirm the list includes nested project FASTQs and Undetermined reads; verify
   every project delivery archive and archive checksum appears as retained. The
   byte estimate must exclude symlink targets. No checksum/archive command should
   be launched. Check that state JSON and output products remain unchanged.
3. Run without `--dry-run`, answer `no`, and verify the same non-mutating behavior.
4. On the disposable run, temporarily rename one expected delivery archive, then
   its checksum manifest, testing separately. Each must block cleanup even with
   `--force`; restore names afterward. QC archives must not substitute for missing
   FASTQ delivery archives. Test a conflicting absolute path and an unavailable
   output directory as applicable.

## Cleanup and notification independence

1. Include a FASTQ symlink to an external test file and a directory symlink to an
   external tree. Confirm cleanup lists the file link only, skips the directory
   link, and preserves both external targets. Use only disposable test fixtures.
2. Run `flowcell-manager clean-fastqs RUN_ID` and confirm. Verify only selected
   FASTQ files/links disappeared; archives, manifests, reports, input spreadsheets,
   SampleSheet, configurations, passwords, provenance and directories remain.
3. Inspect `flowcell-manager show RUN_ID` and `status RUN_ID`: processing and all
   stages stay completed; `fastq_cleanup` records time/outcome and removed paths.
   Verify this also works with a pending or failed completion notification.
4. If a retained notification needs retrying, use the existing
   `flowcell-manager retry-notifications RUN_ID --kind finalized` procedure in
   the test mail environment. It should use retained inputs without regenerating
   products or requiring FASTQs.
5. Repeat cleanup: zero selected files is harmless, and restoration requirements
   for previously removed files remain recorded.

## Locking, interruption and reruns

1. During an active run/management lease, cleanup (including `--force`) must fail.
   Conversely, management operations and daemon execution must not overlap an
   active cleanup. Automated tests exercise the shared lease and confirmation race.
2. On a disposable copy, simulate one unlink permission failure as the normal
   unprivileged service account. Expect nonzero exit, a `partial` outcome, exact
   per-file errors, and processing still completed. Fix permissions and retry;
   only remaining FASTQs should be removed. Running as root may bypass permissions.
3. If testing process termination, use a disposable run. An uncatchable kill may
   leave `in_progress`; the selected manifest was persisted before any deletion.
   Repeating cleanup reconciles the remaining files. Ordinary interruption records
   `interrupted` and the known deletion outcomes.
4. Before restoring reads, try `rerun RUN_ID --from analysis`, `--from reporting`
   and `--from finalization` (a dry run is sufficient for the refusal). Each must
   explain restoration/regeneration and preserve retained archives/state.
5. Restore every selected FASTQ to its original path using normal archive/password
   handling, including Undetermined reads. A partial or truncated restore should
   still fail. After complete restoration, preview and then test the desired
   downstream rerun. File sizes are checked, not contents/checksums.
6. Separately test `rerun RUN_ID --from demultiplexing` with instrument data present:
   it remains available to regenerate reads. Successful demultiplexing marks the
   old cleanup record `regenerated` and clears its old restoration requirement.

## Implementation decisions and limitations

- Reporting reruns also require restoration because they proceed into finalization
  and rebuild delivery archives.
- Symlink directory traversal is blocked during deletion using directory file
  descriptors and `O_NOFOLLOW`; replacing a parent with an external symlink cannot
  redirect an unlink into that external tree.
- A hard kill cannot provide an exact per-file success count. The durable selected
  manifest and `in_progress` outcome deliberately preserve this uncertainty until
  an operator inspects or repeats cleanup.
- The space estimate is logical regular-file size; filesystem sharing/compression
  can change actual space savings.
- No production/server data or email was used during automated tests.
