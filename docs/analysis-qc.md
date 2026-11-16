# Analysis filtering QC

BFQ captures optional filtering metrics in the analysis-complete notification.
It reads the recorded workflow's original marked workdir, under
`data/tmp/<workflow>/bfq/logs/`. It does not execute a filtering tool or require
fastp for a cutadapt-only workflow. Full sample-level QC remains in MultiQC.

## Report contract

Each sample delivers one summary named `<sample>.fastp.json` or
`<sample>.cutadapt.json`. The smallRNA cutadapt route in
[gcf-workflows #351](https://github.com/gcfntnu/gcf-workflows/pull/351) merges its
two passes into one native cutadapt report. BFQ reads that delivered summary;
the individual pass reports under `provenance/` are not counted.

| Metric | fastp | cutadapt |
| --- | --- | --- |
| Input reads | `summary.before_filtering.total_reads` | `read_counts.input` |
| Retained reads | `summary.after_filtering.total_reads` | `read_counts.output` |
| Retained bases | `summary.after_filtering.total_bases` | `basepair_counts.output` |
| Q30 bases | `summary.after_filtering.q30_bases`, or `q30_rate * total_bases` | Unavailable |

Cutadapt's required boolean `input.paired` identifies read-pair counts. BFQ
doubles its input/output counts for paired-end reports so totals use individual
reads, including both mates, consistently with fastp. Output bases already include
both mates and are not doubled. See the
[cutadapt 5.2 JSON schema](https://cutadapt.readthedocs.io/en/v5.2/reference.html#json-report).

Retention is total retained reads divided by total input reads. Q30 is weighted
only by retained bases from reports with valid Q30; cutadapt bases never enter
its denominator. A cutadapt-only project shows Q30 as `Unavailable (0 reports)`.
Mixed projects show combined read/base totals and explicit partial Q30 coverage.
Zero counts remain zero; retention for zero input and Q30 for zero covered bases
are unavailable.

The overview displays filtering reports used/found, per-tool coverage, and
separate coverage for retained bases and Q30. Invalid read counts exclude a
report; absent or invalid optional base/quality counts preserve its read metrics.
Warnings disclose partial coverage. Missing reports never imply zero reads.

Delivery filenames determine the tool and sample even when symlinks point to
plain `.json` files. Aliases resolving to the same report are counted once.
If distinct reports name the same sample, BFQ excludes all of those reports and
warns about ambiguous coverage. It does not choose a processing stage or sum the
same sample twice. The workflow must deliver one authoritative filtering summary
per sample.

## Saved notifications

The schema-v1 analysis snapshot gains `projects[].filtering`, with combined
metrics and per-tool summaries in `filtering.tools`. The existing
`projects[].fastp` fields remain fastp-only. Retained base counts and their
coverage are additive fields. These snapshots are notification payload data;
the operational state schema and restart behavior are unchanged.

Notifications replay their saved `summary_text`/`summary_html`; they do not
recollect QC on delivery retries. Older saved fastp-only snapshots remain
deliverable without migration or the original workdir. Newly collected summaries
use tool-neutral labels.

## Server integration before merge

Use the maintained integration harness on disposable test data with this BFQ
branch and the companion workflow route. Local Python checks use synthetic JSON,
mocked SMTP and no scientific workflow/container execution.

1. Complete a smallRNA project with `filter.trim.quantifier: cutadapt`. Confirm
   delivery includes one native cutadapt summary per sample and no fastp report.
   Compare BFQ's read counts, retention and retained bases with the merged reports.
   Confirm Q30 is unavailable and cutadapt coverage is identified in both email
   formats. Inspect the actual project MultiQC attachment separately; MultiQC's
   cutadapt presentation belongs to its companion change.
2. Complete a representative fastp project and compare read/Q30 totals with its
   reports. For mixed-report coverage, use controlled synthetic delivered reports
   in an isolated test workdir; verify Q30 uses only fastp's covered bases and a
   duplicate sample across tools is excluded with a warning.
3. Exercise the existing notification retry path after a controlled delivery
   failure and normal workdir finalization/cleanup. The saved filtering summary
   must replay unchanged without rerunning analysis. Also verify an older pending
   fastp-only notification still delivers. Use the harness's arranged test SMTP.

No changes to scientific parameters, workflow defaults, tool images, companion
repositories or production deployment are included in BFQ issue #185.
