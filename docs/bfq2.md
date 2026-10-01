# BFQ 2: operational modernization

BFQ 2 marks a substantial change in how the facility operates BFQ: runs have
durable state, operators choose explicit recovery boundaries, sequencing QC is
available earlier, and failed notifications can be recovered independently of
completed processing. Inputs, delivery products, and retained analysis records
have clearer ownership and lifetimes.

This document records the meaningful modernization leading to that baseline,
including foundations from late 2024, configuration and logging work in 2025,
and the operational changes in 2026. It is organized around operator and user
effects; it is not a list of every fix or a claim that all changes arrived in one
promotion.

The milestone is recorded by placing the annotated `bfq2-baseline` tag on the
production promotion commit, **after that commit has been merged into `master`**.
The tag is a historical snapshot. A plain clone continues to select `bfq-dev`;
choose `master` explicitly for current stable production code. BFQ remains
version `2` through routine improvements, and production images advance through
`prod2-N`. No recurring release bundles or GitHub Releases are required.

## What changes for operators and recipients

| Area | Operational effect |
| --- | --- |
| Run status and discovery | JSON state records stages, attempts, failures, and recovery requests. `fm status`, `fm show`, and `fm search` replace inference from marker files or manual inventory searches. |
| Input correction | Output-side sample sheets and submission forms are authoritative. Validate before processing; preserve corrections on reruns unless an explicit input refresh is requested. |
| Workflow selection | BFQ and config generation use the same library-preparation configuration. Invalid metadata stops at preflight; unsupported kit/layout selections stop before analysis. |
| Sequencing feedback | A standalone sequencing report and notification are attempted after conversion and FASTQ renaming, before FASTQ hashing and analysis. Operators can inspect yield and index behavior sooner. |
| Recovery | Restart from demultiplexing, analysis, reporting, or finalization with an explicit preview of what will be removed and retained. |
| Email failure | A failed email does not undo successful processing. Retry the notification without rebuilding reports, checksums, or delivery archives. |
| Delivery and storage | Run-dated project archives, portable checksum paths, retained analysis workdirs, and a separate FASTQ cleanup operation make handoff and subsequent recovery more deliberate. |
| Time reporting | Completion messages report persisted timings for the successful processing that produced the current results, including reused work. |

## The data flow and what each notification means

BFQ still watches instrument directories and processes opted-in, completed runs.
Discovery now covers the configured Nova and Ekista roots in the same service;
explicitly queued state is also eligible for processing. The flowcell manager
prepares and queues work; the BFQ daemon executes it.

| Point in the flow | Data and behavior | What an operator or recipient can conclude |
| --- | --- | --- |
| Input selection and preflight | Select the effective `SampleSheet.csv` and `Sample-Submission-Form.xlsx`; validate their relationship and record the result. | Blocking input errors are identified before BCL conversion. A successful check does not establish that index orientation is scientifically correct. |
| Conversion and FASTQ renaming | Produce the sample/project FASTQ layout. Generate early sequencing QC from run metadata, InterOp, demultiplexing statistics, and the sample sheet. Attempt the **sequencing** notification. | Sequencing and demultiplexing metrics are ready for review. FASTQ checksums, analysis, and delivery may still be unfinished. |
| FASTQ checksums | Generate FASTQ MD5 manifests as part of completing demultiplexing. | Failure here fails demultiplexing and blocks analysis. Early QC alone does not imply this step succeeded. |
| Analysis preflight and analysis | Validate again, resolve the library-preparation workflow, and run the project workflows using the selected configuration. | Input corrections are checked at the analysis boundary; analysis results and project QC become available. |
| Reporting | Assemble run/project summaries and reporting products. Persist completion and attempt the **processed** notification. | Processing reports are available. Delivery archives and their checksums may still be unfinished. |
| Finalization | Build delivery archives and archive checksums; retain successful analysis workdir snapshots and record completion. Attempt the **finalized** notification. | BFQ finalization has completed. This does not itself confirm customer delivery, independent backup, or receipt of an email. |

The sequencing and processed notifications go to `finished_to`; the finalized
notification retains the existing `error_to` routing. Early sequencing QC is
informational: it does not automatically pause analysis, apply acceptance
thresholds, or correct indexes. A failure to create that early report is recorded
separately and processing continues; ordinary analysis, reporting, or
archive/checksum failures still fail their respective processing stages.

The early HTML report is stored at the flowcell output root as
`sequencer_stats_<projects>_<run-date>.html`, with supporting artifacts under
`Stats/sequencing_qc`, and is included in project delivery archives. Project
analysis/MultiQC reports remain separate. Processed summaries distinguish
planned samples from samples actually found in FASTQs, making missing or
zero-read samples and coverage findings visible alongside analysis results.

## Inputs are explicit, checked, and preserved

### Correct the effective inputs

Output-side inputs became the normal correction point before the state-based
rewrite. Each existing output-side input is authoritative independently: a
missing partner may be copied from the instrument directory, but a malformed
output file is not silently replaced. Reruns preserve these files by default.
Use `--refresh-inputs` only when intentionally replacing both from the source.

```console
fm validate RUN_ID
fm validate RUN_ID --refresh-inputs
```

Validation is read-only, including the refresh preview. BFQ and the operator
commands use the shared `gcf-tools` validation rules: every sample-sheet sample
must be represented in the effective submission metadata; extra submission
samples are allowed; repeated lane entries must have consistent assignments.
Warnings remain visible without blocking; errors require correction. Validation
does not repair metadata automatically. It runs before destructive restart
preparation and again in the daemon before the relevant processing stage, with
input hashes and findings retained in state.

Changing sample IDs, project assignments, indexes, or lane assignments normally
requires a demultiplexing rerun. Metadata or workflow changes that leave FASTQ
assignments valid can start at analysis. A reporting-only rerun does not apply a
new scientific workflow.

### One library-preparation configuration

`/opt/gcf-workflows/libprep.config` is the shared authority for BFQ and config
generation. The selected single-end or paired-end workflow is checked before
execution, and unknown kits fail explicitly. Generic QC uses the configured
Custom default workflow where appropriate; an unknown kit is not silently
treated as a default.

The configuration used for analysis is captured with its hash and copied into
the analysis working tree. A new analysis attempt selects the current installed
workflow/configuration. Reporting and finalization reuse the recorded analysis
selection and results; they do not silently reinterpret an old analysis using
the current library-preparation file.

### Index orientation has a supported correction path

```console
fm rerun RUN_ID --from demultiplexing --reverse-complement-index2 --dry-run
```

The index1/index2 flags on `initialize` and demultiplexing reruns toggle the
selected index values in the effective output sample sheet. Preview first, then
repeat without `--dry-run` and confirm when the change is wanted. The preview
shows the proposed values; deciding the correct orientation remains an operator
task. The instrument sample sheet is not edited. Original input bytes, hashes,
and change history are retained under the manager directory.

**These flags toggle; they do not mean “ensure reversed.”** Repeating a flag
reverses that index again. With `--refresh-inputs`, refresh happens first and the
toggle applies to the refreshed sheet. If an interrupted request already
committed the edit, finish recovery with an ordinary rerun, without repeating
the toggle. See the [operating instructions](../README.md#initialize-a-run-and-correct-index-orientation)
for preview and recovery details.

## State and recovery replace marker-file operations

Canonical run records live at `<manager_dir>/states/<run-id>.json`, outside the
flowcell output tree. Preserve and persist the manager directory as operational
data. BFQ uses state and execution leases to coordinate processing and cleanup;
the old `bcl.done`, `files.renamed`, `analysis.made`, and `fastq.made` marker files
are no longer authoritative.

Existing installations retain a compatibility path. JSON state takes precedence;
inventory-only legacy completed or archived runs remain protected from automatic
reprocessing. Recognized restored FASTQs with suitable inputs can bootstrap
analysis. Ambiguous existing output is rejected for operator review instead of
being overwritten. The `flowcells.processed` inventory remains readable and is
updated for compatibility; old records are not rewritten into invented attempt
histories.

`fm` is a short alias for `flowcell-manager`, with the same behavior:

```console
fm list
fm search GCF-1234
fm status RUN_ID
fm show RUN_ID
fm rerun RUN_ID --from analysis --dry-run
```

Search accepts project/run text and status/stage filters, including legacy and
archived records whose output directories no longer exist. `status` provides
the operational summary; `show` exposes the saved record and attempt history.

Choose the earliest stage whose products need to change:

| Restart boundary | Work reused | Work invalidated and rebuilt |
| --- | --- | --- |
| `demultiplexing` | Effective run inputs and retained provenance | FASTQs, their checksums, early sequencing QC, analysis, later reports, and delivery products |
| `analysis` | FASTQs, their checksum manifests, effective inputs, and successful early sequencing QC | Analysis workdirs/results, project reports, later reporting, and delivery products |
| `reporting` | FASTQs/checksums, analysis results, project HTML reports/configurations, and successful early sequencing QC | Later reporting products and delivery archives/checksums; sequencing QC is recovered if needed |
| `finalization` | FASTQs/checksums, analysis results, and reports | Delivery archives and their archive checksums |

In particular, reporting does **not** rerun analysis or regenerate the project
MultiQC reports from analysis. Use analysis when those products must change.
Successful early sequencing QC and its notification survive downstream reruns;
returning to demultiplexing invalidates them.

FASTQ checksum manifests now follow the lifetime of their FASTQs, rather than
being recreated by every downstream retry. A missing or incomplete manifest is
repaired before downstream work proceeds. An existing manifest is checked for
syntax and filename coverage, not used to rehash and verify unchanged contents.
Archive MD5s are rebuilt with their archives.

Use `--dry-run` to review the resolved output path and cleanup plan. `--reason`
records operator intent. `--force` skips confirmation only; it cannot override a
live execution lease or failed validation. An interrupted stage does not
silently resume when its lease disappears: inspect and request the appropriate
rerun. Failed cleanup remains `preparing` until explicitly retried, rather than
being accidentally processed as a ready queue entry.

For a new explicit run, use `fm initialize RUN_ID` (or its full input path when
root selection is ambiguous). Initialization intentionally queues work without
requiring the sequencer's completion markers, so it is an operator action to use
only when inputs are ready. Existing state is managed with `rerun`. BFQ does not
automatically extract legacy delivery archives to reconstruct analysis inputs.

## Notifications and useful completion times

Processing completion and notification delivery are recorded independently.
Missing email configuration, attachment/composition errors, or SMTP failures do
not turn completed reporting/finalization into failed processing. Inspect the
notification state and recover the message itself:

```console
fm retry-notifications RUN_ID --kind processed
fm retry-notifications RUN_ID --kind finalized
fm retry-sequencing-qc RUN_ID
```

Pending notifications receive one automatic attempt, including after a daemon
restart. Failed notifications need explicit retry; a notification already marked
sent is not sent again by a retry. When relay acceptance is uncertain after an
interruption or partial delivery, inspect the situation before using
`--retry-uncertain`, which explicitly accepts the risk of a duplicate. This is
not an exactly-once or inbox-delivery guarantee. Historical completed runs do not
receive invented notification intents or a backlog of completion emails.

`retry-sequencing-qc` recovers the report using retained statistics without
conversion, analysis, or FASTQ hashing. It may attempt a sequencing notification
that has never been attempted; failed/uncertain mail still uses the notification
retry commands. If recovering the report after finalization, rerun finalization
before handoff when delivery archives need to contain the recovered report.
Report recovery alone does not rewrite existing archives.

Completion messages show six timing categories: demultiplexing, FASTQ MD5,
analysis, reporting, archiving, and archive MD5. They describe successful work
still contributing to the current products. Downstream retries retain timings
for reused results and replace timings for rebuilt work, rather than summing
all failed or abandoned attempts. These are flowcell wall-clock processing
times, not per-project CPU usage or end-to-end elapsed time; queueing, SMTP, and
setup/preflight time are excluded. Legacy or partial timing information is
identified as such.

Production error notifications provide run/stage context and command/log details,
with suppression of repeated reports for the same failure. `BFQ_ENV=production`
(or `prod`) enables this error-email behavior. The default test environment
suppresses these production **error** messages; it does not disable all ordinary
completion notifications. Test recipients and email handling must still be set
appropriately. Static configuration changes require a daemon restart; `SIGHUP`
only wakes a sleeping scan. An explicit notification retry loads current email
settings.

## Delivery, retained analysis records, and reclaiming space

### Run-scoped deliverables

Earlier modernization brought run dates into archive/report names, copied sample
information into project TSVs, and organized workflow/QC delivery around the
flowcell run folder. Project data and QC archives use
`<project>_<run-date>.7za` and `QC_<project>_<run-date>.7za`, with corresponding
archive checksum files. Relative paths in archive MD5 manifests make them usable
after moving the delivery files. QC archive creation includes symlinked target
data so recipients receive usable content rather than links into the facility
filesystem. Sensitive-data delivery retains its encrypted-archive and password
file convention.

### Keep the actual successful analysis working tree

After successful finalization, each project's retained record is:

```text
<outputDir>/<run-id>/provenance/<project>_analysis.tar.gz
```

This contains the actual analysis workdir: generated configuration and Snakefile,
the copied `gcf-workflows` tree including local edits/untracked files, and working
metadata. Only the top-level `data/` entry is excluded. It provides much stronger
evidence of how an analysis was run than the currently installed workflow tree.

It is one latest successful snapshot per project, not an archive of every
attempt. Failed retries do not replace the previous successful snapshot;
successful reanalysis followed by successful finalization does. Reporting and
finalization reuse a matching retained snapshot. Where an older run has neither
a matching snapshot nor its original workdir, BFQ reports that limitation rather
than substituting today's installed source. These records survive all restart
boundaries and `fm archive`.

If publication of a committed snapshot is interrupted, BFQ retains the staged
candidate and retries publication on a later scan. Leave that pending staging
data intact so recovery can finish.

The snapshot is an internal provenance record, separate from encrypted delivery
archives. It has owner-only file permissions, including for sensitive runs, but
is not itself encrypted. Symlinks are retained as links; referenced datasets,
containers, and the complete external execution environment are not frozen into
it. It is not a self-contained rerun bundle or an automatic restore mechanism.

### FASTQ removal is a separate, deliberate operation

```console
fm clean-fastqs RUN_ID --dry-run
fm clean-fastqs RUN_ID
```

This removes FASTQ files and FASTQ symlinks from a finalized state-backed run,
including nested and Undetermined FASTQs, while retaining delivery archives,
checksum manifests, reports, inputs, passwords, and provenance. It removes links
themselves, not their external targets. Pending email does not block it.

The guard checks completed finalization and the existence of each project's data
archive and archive MD5 file. It does not test archive contents, hashes, delivery,
or backup redundancy. Confirm handoff and retained copies operationally before
cleanup. `fm archive` is the later, broader removal of delivery data from the
output directory; it preserves state and provenance and is distinct from the
pipeline's finalization stage that creates delivery archives.

After FASTQ cleanup, downstream reruns require restoring the recorded FASTQs,
including Undetermined files, to their original paths and sizes. That prerequisite
checks existence/size, not content integrity. Alternatively, rerun demultiplexing
when instrument inputs remain available. **Restore before queueing a downstream
rerun:** queue preparation invalidates delivery archives, so those archives must
remain available until restoration is complete.

## Runtime and maintenance changes supporting these operations

The operational changes build on several less visible modernizations:

- Shared `PipelineConfig` and clearer path handling replaced separate parsing
  approaches in the service and manager. Multiple instrument roots can be used
  without reconfiguring the service for each run. Structured logging and debug
  controls improve diagnosis across processing stages.
- Modern Python packaging declares dependencies and installs `bfq`,
  `flowcell-manager`, and `fm`; the old `.py` executable names are not installed.
  The supported Python baseline is 3.11 or newer. The installed BFQ package is
  the source of `BFQ 2` in CLI output and new attempt metadata; old attempt
  version values remain historical records. State schema versioning is separate.
- The base image and analysis runtime were modernized, including Snakemake and
  Apptainer. BFQ-facing wrappers now run conversion, supported 10x tools, and
  MultiQC through workflow-selected containers. `/opt/gcf-workflows/docker.config`
  supplies those selections; BFQ command configuration supplies options rather
  than a second independent collection of executable/image definitions. Persistent
  input/output paths and required bind mounts remain part of site setup.
- External command construction now handles arguments, working directories, and
  pipeline failures more explicitly. InterOp output parsing was updated for
  current tool behavior. These changes reduce failures caused by quoting, paths,
  or obsolete assumptions and improve the error information operators receive.
- Obsolete QIAseq helper scripts and primer files were removed. This was cleanup
  of unused bundled material, not a new replacement analysis workflow.
- Configuration, state transitions, input validation, notification recovery, and
  cleanup gained automated regression coverage. Ruff checks and CI were added
  and made non-mutating. Real instrument/runtime/SMTP integration remains a
  separate operational verification step.

The scientific analysis version remains the `gcf-workflows` commit recorded as
**Analysis pipeline** in the project report. BFQ generation, per-attempt runtime
metadata, and Docker build number answer different questions. Rebuilding delivery
products with a newer runtime does not relabel the preserved scientific analysis
as though it had been rerun.

## History and further operating detail

The [README](../README.md) is the command/configuration reference. Detailed
[input-preflight guidance](input-preflight.md) and the issue-specific verification
documents in this directory provide further examples. The milestone records
behavior at the tagged commit; consult the selected branch's README for later
changes.

The major lines of work are traceable through these changes:

| Period | Foundations and changes |
| --- | --- |
| Late 2024 | [Output-side sample sheets (#68)](https://github.com/gcfntnu/gcf-bfq/pull/68), expanded instrument discovery including MiniSeq, [run-dated archives and run-folder QC delivery](https://github.com/gcfntnu/gcf-bfq/commit/1f7a27ea4d903ee1b7fa910d96184674cb66a50a), and [project sample metadata](https://github.com/gcfntnu/gcf-bfq/commit/8369a0996bf9ba69953b78ab225bbd24ba1c92e9). |
| 2025 | [Simultaneous instrument-root discovery](https://github.com/gcfntnu/gcf-bfq/commit/1c2754a5122abe4696271c814d60682d73b7a75e), [shared configuration (#79)](https://github.com/gcfntnu/gcf-bfq/pull/79), [path/dependency cleanup (#80)](https://github.com/gcfntnu/gcf-bfq/pull/80), [structured logging (#83)](https://github.com/gcfntnu/gcf-bfq/pull/83), [lint/CI (#77)](https://github.com/gcfntnu/gcf-bfq/pull/77), and [portable archive checksum paths](https://github.com/gcfntnu/gcf-bfq/commit/3015410f25b86a5c4babf7d7dd42154f1e0c5b85). |
| 2026 runtime and configuration | [Configuration tests (#88)](https://github.com/gcfntnu/gcf-bfq/pull/88), [InterOp compatibility (#97)](https://github.com/gcfntnu/gcf-bfq/pull/97), [base image (#98)](https://github.com/gcfntnu/gcf-bfq/pull/98), [packaging (#100)](https://github.com/gcfntnu/gcf-bfq/pull/100), [configuration documentation (#101)](https://github.com/gcfntnu/gcf-bfq/pull/101), [command execution (#105)](https://github.com/gcfntnu/gcf-bfq/pull/105), [obsolete helpers (#107)](https://github.com/gcfntnu/gcf-bfq/pull/107), [non-mutating CI (#108)](https://github.com/gcfntnu/gcf-bfq/pull/108), and [container wrappers (#109)](https://github.com/gcfntnu/gcf-bfq/pull/109). |
| 2026 state and recovery | [State-based management (#111)](https://github.com/gcfntnu/gcf-bfq/pull/111), [production error reporting (#112)](https://github.com/gcfntnu/gcf-bfq/pull/112), [FASTQ checksum lifetime (#120)](https://github.com/gcfntnu/gcf-bfq/pull/120), and [independent notifications (#125)](https://github.com/gcfntnu/gcf-bfq/pull/125). |
| 2026 inputs and visibility | [Library-preparation authority (#127)](https://github.com/gcfntnu/gcf-bfq/pull/127), [input preflight (#128)](https://github.com/gcfntnu/gcf-bfq/pull/128), [early sequencing QC (#131)](https://github.com/gcfntnu/gcf-bfq/pull/131), [index correction (#132)](https://github.com/gcfntnu/gcf-bfq/pull/132), [processing timings (#133)](https://github.com/gcfntnu/gcf-bfq/pull/133), and [manager alias/search (#134)](https://github.com/gcfntnu/gcf-bfq/pull/134). |
| 2026 retention and promotion | [Analysis snapshots (#129)](https://github.com/gcfntnu/gcf-bfq/pull/129), [FASTQ cleanup (#130)](https://github.com/gcfntnu/gcf-bfq/pull/130), [branch/promotion guidance (#113)](https://github.com/gcfntnu/gcf-bfq/pull/113), and [BFQ generation 2 (#135)](https://github.com/gcfntnu/gcf-bfq/pull/135). |

Maintainers: see the separate [one-time tagging procedure](bfq2-tagging.md) and
[annotated-tag message](bfq2-tag-message.txt).
