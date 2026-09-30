# Issue #119: server verification before merge

Development tests use mocked SMTP only, with an autouse guard rejecting real
SMTP connections. No email was sent during development. Real relay, container
and instrument verification belongs on the integration servers.

## Expected emails

1. **Demultiplexing complete — sequencing QC**: run/project/user, read geometry,
   sequencing yield and lane quality, assignment/undetermined/zero-read overview,
   disk availability, sequencing MultiQC attachment. Explicitly says analysis is
   separate and there is no automatic QC decision.
2. **Analysis complete — QC summary**: discovered and missing FASTQ samples,
   original submission/sample-group information, optional fastp reads retained and
   after-filter Q30 with report coverage. Project MultiQC/single-cell attachments;
   no duplicate sequencing attachment for runs with an early notification.
3. Existing finalization mail after archive preparation/checksums.

Metrics that do not exist for a workflow are explicitly unavailable. fastp
retention is a read-count ratio; after-filter Q30 is weighted by retained bases.
The input/output totals include both mates, not read pairs. Neither email grades
samples as passing or failing QC.

## Normal run and native input layouts

Build from `119-early-sequencing-qc` and run a representative flowcell normally.
Confirm logs show BCL conversion/rename, early report/email, FASTQ hashing,
analysis, analysis email, archive/checksum finalization, final email in that order.
Verify `flowcell-manager status RUN_ID` and `show` include completed sequencing QC
and one sent `sequencing:<execution>` entry.

Inspect `Stats/sequencing_qc/sequencing_qc.html` and compare values with the source
statistics. Cover bcl-convert `Reports/Demultiplex_Stats.csv` plus `Quality_Metrics.csv`
and `Top_Unknown_Barcodes.csv`, bcl2fastq `Stats/Stats.json`, and supported 10x
mkfastq output. Check actual MultiQC images/wrappers and InterOp binaries. Native
bcl-convert input staging must put RunInfo.xml beside the CSV. A Stats symlink to
Reports/legacy/Stats is supported. The report should not import project analysis
modules or configuration. Check branding/module layout manually in the browser.

Check fastp summaries against `data/tmp/<workflow>/bfq/logs/**/*.fastp.json` in the
original marked workdir. Workflows without those reports must show unavailable,
not zero or an invented universal filtering count. MultiQC remains authoritative
for full workflow-specific metrics.

## Poor assignment and downstream failure

Use a controlled wrong-index run with nearly all reads undetermined, multiple
projects, a known zero-read sample, and preferably no usable project FASTQs.
Confirm planned projects/samples appear from the SampleSheet/statistics and zero
is distinct from unavailable. Check top unknown indexes. Poor values must not
pause or change processing automatically.

Trigger a downstream analysis/configmaker failure after conversion. The early
report/email must already exist and remain readable. Report generation/recovery
must work without analysis YAML/results or a readable submission workbook. Fresh
runs still undergo the existing preflight before conversion; do not disable it
to test this behavior. For the independent report-recovery check, use already
converted output with a failed/missing early report, then make the workbook
unreadable to downstream parsing and run `retry-sequencing-qc`.

## Preservation, retry and new execution

- Restart the daemon: no duplicate sent early mail.
- Run `rerun RUN_ID --from analysis` and complete analysis: early report mtime,
  execution identity and duration stay unchanged; no early resend. Repeat for a
  reporting-only rerun if useful.
- Correct the SampleSheet and rerun from demultiplexing: old early artifacts and
  notification become superseded; new conversion gets a new early report/mail.
- Cause a report-generation failure (e.g. unavailable MultiQC on a controlled
  integration deployment): state records useful diagnostics/log path while
  downstream processing proceeds. Restore the tool, run
  `flowcell-manager retry-sequencing-qc RUN_ID`, and verify no BCL/hash/analysis
  rerun. Repeating the command must not resend sent mail.
- Cause a definite relay rejection: report remains completed and processing
  continues. Correct settings and use
  `flowcell-manager retry-notifications RUN_ID --kind sequencing`.
  Uncertain acceptance requires `--retry-uncertain`, as for other notifications.
- After successful finalization, inspect project `.7za` contents for
  `Stats/sequencing_qc/`. Recovery performed after finalization requires an
  explicit finalization rerun to refresh archives; the command warns about this.

## Timing and compatibility

`sequencing_qc.duration_seconds` records only the latest successful generation for
this conversion. Reporting metadata includes `sequencing_qc_duration_seconds`,
`analysis_reporting_duration_seconds` and their sum in `duration_seconds`.
SMTP time and failed-generation time must not inflate these values. Six-category
email totals remain separate work in #118.

Older records load without the new field. A legacy analysis/reporting rerun does
not invent an early notification; its analysis email includes a sequencing-report
attachment with that explanation. Previously queued v1 combined emails retain
the old template/attachment paths until their own outputs are invalidated.

## Development verification

- Full BFQ pytest suite: **500 passed**, with the SMTP guard enabled.
- Ruff lint, formatting and `git diff --check`: passed.
- Isolated upstream MultiQC **1.18** generated synthetic bcl-convert + InterOp
  and bcl2fastq + InterOp reports. Checked native demultiplexer content and BFQ
  custom sections. SMTP was never used and MegaQC upload was disabled.
- The custom deployed MultiQC image, real InterOp binaries, real instrument
  layouts, archive contents and actual email delivery still need the server
  checks above. Synthetic tests do not replace them.

Compatibility findings resolved during development: RunInfo/CSV colocation for
MultiQC bclconvert; separate custom-content YAML for older MultiQC HTML parsing;
mkfastq's optional flowcell-id directory; duplicate sample IDs disambiguated by
lane/index, with incomplete lane totals marked unavailable; and recovery of a
missing early report without superseding valid later notifications.
