# Issue #115: index orientation integration checks

PR target: `bfq-dev`. Keep this PR unmerged until operator integration testing and
manual verification are satisfactory. No companion repository update is required.

Use disposable fixture flowcells and test output/manager directories. Retain an
untouched copy of each fixture's original SampleSheet for comparison. Select
non-palindromic index values, since an index equal to its own reverse complement
cannot distinguish the two orientations. For example, `AACG` toggles to `CGTT`,
and `GATT` toggles to `AATC`.

Management commands queue work; stop the test daemon while inspecting preparation
and enable it only for the intended conversion pass. This preserves today's
initialize behavior, which can queue a run independently of normal
completion-marker discovery. Development and automated tests must mock SMTP and
must never send real emails. Any server-side notification checks belong to the
operator's separately controlled integration environment.

## 1. Initialize from a selected input directory

1. Place a complete disposable input flowcell under a configured instrument root,
   with valid SampleSheet and submission metadata. Before initialization, run:

   ```bash
   flowcell-manager initialize /instruments/RUN_ID --reverse-complement-index2 --dry-run
   ```

   The preview must identify the exact input directory, configured output
   directory, default `demultiplexing` boundary and proposed index2 correction.
   No output files, backup, queued state or input edit may result. Repeat without
   `--dry-run`, decline confirmation and verify the same lack of mutations.
2. Confirm the operation. Verify the copied output SampleSheet has only index2
   values reverse-complemented, the submission form is copied, and processing is
   queued from demultiplexing. Record the source and output checksums. The source
   file must be byte-for-byte unchanged.
3. Inspect `flowcell-manager show RUN_ID`: source path and output path must refer
   to the selected directories. The `index_corrections` history should record the
   operation and before/after checksums; index2's live comparison should indicate
   reversed relative to source. Inspect the backup under
   `<manager_dir>/input-history/RUN_ID/<operation_id>.SampleSheet.csv`; it should
   contain the effective sheet immediately before the correction.
4. Repeat with another fixture using an exact `/import-kista/RUN_ID` input path.
   With identical run IDs under two configured roots, bare-ID initialization
   must reject the ambiguity and request an explicit path. An explicit path
   must select the intended directory rather than silently use the other root.
5. Verify that a second `initialize` of the same initialized run refuses and
   points the operator towards `rerun`. Separately retain coverage of explicit
   `initialize RUN_ID --from analysis` without index flags, using the existing
   restored-input workflow.

## 2. Toggle either or both indexes

Use the stopped-daemon fixture initialized above. Start from a known copy of its
source sheet so each expected value is clear.

1. Preview and apply index1 alone:

   ```bash
   flowcell-manager rerun RUN_ID --from demultiplexing --reverse-complement-index1 --dry-run
   flowcell-manager rerun RUN_ID --from demultiplexing --reverse-complement-index1
   ```

   Only `index` values change. The preview must still include the normal
   demultiplexing cleanup plan. Check the backup against the pre-operation
   effective sheet and verify a new history record.
2. Repeat index1 with the same explicit flag. Values return to their original
   orientation and a second correction is recorded. There is no repeat-protection
   refusal. Repeat this two-toggle check for index2.
3. Apply both flags together, verify both columns, then apply both again and
   compare the complete output sheet with the starting bytes. Also test index1
   followed by index2 in separate invocations.
4. Perform an ordinary demultiplexing rerun with no correction flag. The effective
   sheet and correction history must remain unchanged. An ordinary downstream
   rerun likewise preserves the sheet. Correction flags with `--from analysis`,
   `reporting` or `finalization` must be rejected before mutation.
5. Test `--refresh-inputs --reverse-complement-index2`: the output should contain
   one index2 reverse complement of the source, regardless of its previous
   orientation. Refresh without a correction should restore source values.
6. Verify the hidden compatibility alias `--tom-mode`: it changes index2 only.
   Combining it with `--reverse-complement-index2` must apply index2 exactly once;
   combining it with `--reverse-complement-index1` must change both columns.
   The alias must be absent from normal command help.

## 3. Format preservation and validation

These cases use small local fixture sheets, not altered instrument data.
Automated tests cover them; repeat representative cases against the installed
CLI to verify packaging and service-account permissions.

1. Include quoted fields containing commas, CRLF line endings, a UTF-8 BOM,
   unrelated settings and project/sample columns. Compare pre/post bytes and
   confirm only the requested index content changes. Repeat for `[Data]` and
   `[BCLConvert_Data]`. Existing `ReverseComplementIndexP5/P7` settings must not
   change.
2. Include mixed-case IUPAC DNA sequences and empty index cells among populated
   rows. Verify the expected case-preserving complements and unchanged empties.
3. Request a missing index column, an entirely empty selected index, and an index
   containing unsupported characters. Each must reject the correction without
   changing the effective sheet, creating correction backups, queueing a new
   restart or invalidating downstream products.
4. Request both columns with a valid index1 and invalid index2. Neither column
   may change. Restore valid data before continuing.
5. Confirm normal rerun cleanup preserves the history and backup files. Inspect
   archive cleanup's plan and test its retention on a disposable finalized run
   if available; backups are outside the output tree and must remain available.

## 4. Source-relative reporting

1. After one index2 toggle, verify the resulting message and `show` identify it
   as reversed relative to source. After the second, verify they report matching
   the source/original orientation.
2. Reorder the rows of a disposable source reference without changing sample/lane
   identity or index values. Comparison should continue to match corresponding
   entries rather than depend on row order.
3. Use separate fixtures with mixed orientations, an unrelated manual index edit
   and ambiguous sample/lane identities. Reporting must describe the uncertainty
   rather than label all rows original or reversed.
4. Temporarily make only the fixture's recorded source SampleSheet unavailable.
   A valid toggle of its existing output sheet must still succeed, with an
   unavailable comparison. Restore the reference afterward.
5. Change the disposable source reference after a correction and inspect `show`.
   Live orientation is evaluated against the current reference, whereas recorded
   history retains the source path/checksum and comparison used at correction
   time. Matching the source never establishes biological/demultiplexing
   correctness. Self-reverse-complementary indexes are also inherently ambiguous.

## 5. Locking and interrupted preparation

1. Hold an active execution/management lease on the fixture. Initialization or
   correction attempts, including `--force`, must refuse to overlap that lease.
   The source, effective sheet, backup history and processing state must remain
   unchanged by the refused operation.
2. Under the normal service account, exercise a backup-write failure on a
   disposable manager directory. A correction must not replace the effective
   sheet unless its backup has been made durable. Restore permissions before
   retrying. Root may bypass permission-based failure fixtures.
3. Use automated fault injection for interruption between backup, replacement and
   history persistence. A completed replacement must remain identifiable from
   its durable operation record, and a temporary partial file must never become
   the effective SampleSheet. Do not approximate this by killing a production
   process.
4. On a disposable fixture, interrupt preparation after the corrected sheet is
   committed but before processing is queued. Inspect the output sheet and
   `show`, then use an ordinary demultiplexing rerun without correction flags.
   It must finish preparation with the existing corrected values. If the flag
   is explicitly supplied again instead, that is a new requested toggle and
   must reverse the current values again.

## 6. Server conversion pass

After the management checks, use a disposable instrument run with a known index
orientation problem, correct it and allow the normal test daemon to process it.
Verify the converter consumes the corrected output-side SampleSheet, expected
reads demultiplex to the intended samples, and normal analysis/reporting and
finalization proceed. Reinspect the instrument SampleSheet to ensure it remains
unchanged. Confirm the correction history and backups remain readable after the
normal completion lifecycle.

Real BCL conversion, network-mounted source availability, service-account storage
permissions and daemon queue pickup need this operator pass. They are not
established by the isolated automated tests.

## Operational decisions

- Explicit initialization retains prepare-and-queue behavior. No new prepared
  state or completion-marker gate has been introduced.
- Repeated explicit flags always toggle current values. Recovery without a new
  toggle uses an ordinary rerun; a new flagged invocation intentionally toggles.
- Orientation is relative to the recorded instrument source, whose content can
  change or become unavailable. It is not a determination of correct index
  orientation for a sequencer or library kit.
- The hidden alias is index2-only; simultaneous aliases do not duplicate the
  operation within one command.

## Implementation findings addressed

- Existing active-run probing rewrote its lease file, and stale-run recovery
  could update state before confirmation. Probing is now read-only; a rerun
  recovers a stale attempt only after confirmation under its execution lease.
- Ordinary rerun preparation replaces the restart request. Attempt provenance
  therefore associates the actual effective SampleSheet checksum with the latest
  applied correction, preserving the link across interrupted preparation while
  avoiding a stale link after input refresh or manual edits.
- Live orientation resolves legacy bare output paths through the configured
  output root. It does not depend on the operator's working directory or rewrite
  the original state merely to display the comparison.
