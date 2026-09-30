# Successful processing times (#118)

Use `118-persist-successful-processing-times` on the integration server. Development
and automated tests mock all email delivery; the server checks below verify the
real tools, persistence across daemon restart, and both completion emails before merge.

## What is measured

| Email row | Work | Invalidated by rerun from |
| --- | --- | --- |
| Demultiplexing | BCL conversion and FASTQ renaming | demultiplexing |
| FASTQ MD5 checksums | FASTQ manifest generation, enclosing all projects/workers | demultiplexing |
| Analysis | Existing analysis/workflow execution | demultiplexing, analysis |
| Reporting | Successful early sequencing QC plus later reporting/QC summary generation | Early part: demultiplexing; later part: demultiplexing, analysis, reporting |
| Archiving | Existing delivery archive generation | any restart boundary |
| Archive MD5 checksums | Existing archive checksum worker, enclosing its parallel workers | any restart boundary |

These are timing categories, not additional restart boundaries. Workflow-generated
reports remain part of Analysis, matching where that work executes today. Input
preflight, daemon setup, manifest coverage checks, inventory/provenance bookkeeping,
notification composition, SMTP, downtime and waiting between attempts are excluded.
Early sequencing QC retains its existing timing and non-blocking failure behavior.

The six rows and **Total processing time** replace the current-attempt runtime in
new processing-complete and finalization emails. No recipients, subjects, attachments
or delivery/retry policies change; the early sequencing email is unchanged.
Processing-complete mail marks archiving/checksums **Not completed** and labels the
total as a subtotal. Unknown legacy durations are **Timing unavailable**, with a
partial total. If only one Reporting component has completed, the row also shows
its recorded duration; that successful component contributes to the subtotal.
Durations are stored in seconds with fractional precision and displayed to whole
seconds. Row/total rounding can differ by a few seconds.

## Inspecting records

```console
flowcell-manager show RUN_ID
```

`processing_timings` contains six step keys, each holding an execution history.
Records include the flowcell attempt, UTC `started_at`/`completed_at`, monotonic
`duration_seconds`, and `outcome` (`running`, `completed`, `failed`, `interrupted`).
Only the latest completed record without `invalidated_at` contributes. Reruns retain
old records but invalidate affected successes before cleanup. The early component
remains in `sequencing_qc` and its `attempts`; its duration is added to Reporting once.

Each completion intent saves a `payload.processing_timing` snapshot after its
relevant timers are committed. Notification retries use that event's saved snapshot,
so a later retry does not change the original email totals. Pre-upgrade intents
retain their original email format. Existing schema-v1 state needs no migration;
missing durations are never reconstructed from timestamps or old runtime strings.

## Server checks

1. **Fresh run:** initialize/process a small flowcell normally. Compare all six
   successful durations with the corresponding tool log intervals. The processing
   email should contain the first four categories and a subtotal; the finalization
   email should contain six durations and their total. Confirm both emails still
   contain the existing QC/operational content and attachments.

2. **Failure, daemon restart and analysis rerun:** induce a controlled analysis
   failure on test data. Inspect the saved successful conversion, checksum and early
   QC timings and the failed analysis duration. Stop/restart BFQ, leave an observable
   idle gap, correct the cause, and run:

   ```console
   flowcell-manager rerun RUN_ID --from analysis --force
   ```

   Before BFQ consumes the request, inspect state: upstream durations remain,
   analysis/reporting/archive durations are invalidated. After success, confirm
   only the new successful analysis/downstream durations appear in the email total.
   The failed interval and idle gap must contribute nothing. FASTQ checksum files
   and their original timing should be unchanged.

3. **Reporting and finalization reruns:** exercise `--from reporting` and then
   `--from finalization`. Each should retain the appropriate upstream timings and
   the original early QC duration; regenerated categories replace, rather than add
   to, their old durations. If archive checksum generation fails, Archiving success
   must already be persisted. A finalization rerun invalidates both categories
   because that boundary regenerates both products.

4. **Demultiplexing rerun:** run `--from demultiplexing`. All six current timings and
   the early QC timing must be invalidated with their outputs. After success, confirm
   the new timings belong to the new attempt and old successes remain in history.

5. **Interrupted execution:** terminate a test worker/daemon while a timed step is
   active. On recovery, the record must be interrupted with no invented duration
   extending to recovery time. An explicit rerun must be required as before. A caught
   interruption may retain its measured interval for diagnostics, never for totals.

6. **Legacy and repair:** rerun an existing/restored flowcell from analysis. Preserved
   results without measured timings must display **Timing unavailable** and a partial
   total. A complete FASTQ manifest must not acquire a tiny replacement duration from
   its coverage check. If a manifest is missing/incomplete, its actual repair is timed;
   a new successful repair replaces that category's previous execution duration.

7. **Early QC/SMTP failure:** let early QC fail while processing continues; Reporting
   must identify incomplete work and exclude the failed interval. Recover it with
   `retry-sequencing-qc` and inspect the successful component. As before, rerun from
   finalization to refresh delivery archives after late report recovery. For an SMTP
   failure, retry the notification and verify identical saved timing values without
   changing processing records or rerunning workers.
