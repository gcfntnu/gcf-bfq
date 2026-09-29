# Issue #124: retained analysis snapshot server tests

Use branch `124-retain-analysis-workdir-archive` with the normal BFQ test image
and mounts. No companion repository change or new dependency is required.
The automated suite exercises actual tar creation and state/cleanup transitions;
real Snakemake, 7-Zip, mounted storage permissions and production data volumes
still need this server pass. Keep this PR unmerged until that pass is accepted.

## 1. Run and inspect a recognizable analysis

1. Select a disposable test flowcell with valid inputs and sufficient output
   space. Before a fresh analysis, add a harmless identifying comment to the
   test image's `/opt/gcf-workflows/libprep.config`, and an untracked small file
   such as `/opt/gcf-workflows/snapshot-test-note.txt`. Use an approved parameter
   edit too if you want to verify its effective value in `config.yaml`.
2. Start the flowcell from analysis (or run it fresh from demultiplexing):

   ```bash
   run=260918_MN00686_0026_A000HCMFHF
   project=GCF-2026-044
   flowcell-manager rerun "$run" --from analysis --force
   ```

3. After successful finalization, inspect the retained archive:

   ```bash
   snapshot="/mnt/bfq/output/$run/provenance/${project}_analysis.tar.gz"
   tar -tzf "$snapshot"
   tar -xOzf "$snapshot" src/gcf-workflows/libprep.config
   tar -xOzf "$snapshot" config.yaml
   flowcell-manager status "$run"
   flowcell-manager show "$run"
   ```

   Expect `config.yaml`, `Snakefile`, the copied workflow tree, the untracked
   note, `.bfq-analysis.json` and `.snakemake/log/...`. No member should equal
   `data` or start with `data/`. Confirm the archive contains the actual local
   configuration edit, not only a Git revision. Status should report `retained`;
   top-level `analysis_snapshots` should contain no `pending` field.
4. Extract into a separate directory with `tar -xzf ... -C ...` and inspect
   the logs and configuration. Links should remain symlinks. Do not run the
   extracted workflow until omitted data and external dependencies are restored.

The generic content/symlink cases are covered automatically: external/dangling
links, a top-level `data` symlink, nested `data` directories, hardlinks, dotfiles,
and unrelated newly added files. A local test workdir can be used for additional
inspection without modifying real sequencing data.

## 2. Downstream-only retry after installing different code

1. Save the archive's checksum and modification time (this is a small analysis
   snapshot, not a checksum pass over delivery data):

   ```bash
   sha256sum "$snapshot"
   stat -c '%y' "$snapshot"
   ```

2. Change the test image's installed workflow comment/parameter **after** the
   original analysis has completed. Retain the original workdir/QC data needed
   for normal delivery finalization.
3. Run `flowcell-manager rerun "$run" --from reporting --force`, let BFQ finish,
   and then repeat from `finalization`.
4. The snapshot checksum and timestamp must remain unchanged. The retained
   libprep file and `config.yaml` must still describe the original analysis.
   No workflow should be recopied from `/opt` during these retries.

A matching snapshot is reusable even without its workdir. This is tested
separately in automation because normal QC archive creation still has its own
workdir/data requirements; this feature does not recreate those missing data.

## 3. Failed rerun, followed by successful replacement

1. Record the current snapshot checksum. On the disposable flowcell, introduce
   a reversible analysis failure in the test workflow and queue an analysis
   rerun. Verify the retained archive survives both cleanup and the failure.
2. Restore the workflow, give its local comment/note a new recognizable value,
   and queue a successful analysis rerun. Once finalized, the archive should
   contain the new value and have a new checksum. Exactly one public
   `<project>_analysis.tar.gz` should exist for that project.
3. If practical, repeat with an unavailable test SMTP relay. Processing must
   remain completed and the new snapshot retained; only notification state
   should fail. Do not change production mail routing for this test.

Injected automated failures cover archive-write interruption, disk errors,
pre-commit failure, post-commit state error, failed publication, interruption
between rename and state update, and partial multi-project publication. A hard
kill during snapshot preparation must leave the old public archive intact; a
kill after completion may leave a committed candidate in `.staging/`, which BFQ
publishes on its next scan/startup. Do not manually delete that directory.

## 4. Cleanup retention

1. On a completed disposable run, inspect dry runs at all boundaries:

   ```bash
   flowcell-manager rerun "$run" --from demultiplexing --dry-run
   flowcell-manager rerun "$run" --from analysis --dry-run
   flowcell-manager rerun "$run" --from reporting --dry-run
   flowcell-manager rerun "$run" --from finalization --dry-run
   flowcell-manager archive "$run" --dry-run
   ```

   `provenance/` must never be in the deletion plan. Actual cleanup through all
   four boundaries is also covered automatically.
2. After finishing the retry tests, perform the usual temporary workdir cleanup
   and `flowcell-manager archive "$run" --force` on this disposable/delivered
   data. The snapshot must remain readable and its configuration/logs intact.
   No publication should be pending before removing temporary or staged files.

## Operational decisions and limitations

- Publication is a recoverable two-step operation because the canonical JSON
  state and output archives can be on different filesystems. All candidates are
  complete before finalization is committed; each archive is then atomically
  replaced. Publication errors preserve candidates and appear in logs/status.
  BFQ retries at its next scan/startup; rerun/archive commands recover first and
  stop cleanup if committed publication cannot be completed.
- A state-write error after the atomic completion replacement is treated as a
  committed delivery when the completed record is readable. This avoids marking
  an already published successful snapshot as a failed rerun. The error is logged.
- Workdirs retain the existing `<project>_<date>` naming. An ownership token
  prevents capturing a directory reused by a different run; it does not change
  scheduling or isolate concurrent same-project/same-day analysis directories.
  Avoid such concurrent analyses until that separate layout limitation is
  addressed. The token is checked before and after archive creation.
- Old state has no ownership token. An existing unmarked conventional workdir
  containing `config.yaml` and `Snakefile` can be retained with a warning, but its
  historical ownership cannot be proved. If both original workdir and matching
  retained snapshot are unavailable, finalization records `unavailable` and logs
  the reason rather than substituting installed code. Other finalization
  prerequisites still apply. An older successful archive is kept and labeled
  accordingly by `flowcell-manager status`.
- These archives are internal operational artifacts, with owner-only file
  permissions, and are not added to the FASTQ/QC delivery archives. They are not
  encrypted by `SensitiveData`; retain the same access controls as the workdir
  and output metadata. Size depends on everything outside top-level `data/`,
  including any workflow checkout caches stored there.
