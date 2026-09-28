# GCF BFQ

BFQ monitors Illumina run directories, demultiplexes completed flowcells, runs
the configured GCF analysis workflow, generates QC reports and archives, and
records completed projects in the flowcell inventory.

The pipeline is developed for the Genomics Core Facility at NTNU. Its paths,
sample-sheet extensions, workflow selection, and delivery conventions are
site-specific. The GCF Docker images are the supported production runtime.

## Runtime overview

BFQ runs as a long-lived process:

1. Load `/config/bcl2fastq.ini`.
2. Search the configured Nova and Ekista roots for completed sequencing runs.
3. Use JSON state to select queued runs; protect inventory-only legacy runs.
4. Require a sample sheet with `[CustomOptions]` and a sample submission form.
5. Demultiplex, run the selected Snakemake workflow, and generate MultiQC output.
6. Create delivery archives and checksums.
7. Record completion in JSON state and update the compatibility inventory.
8. Sleep for the configured interval before scanning again.

Static configuration is loaded once when BFQ starts. Restart the process after
changing `/config/bcl2fastq.ini`. Sending `SIGHUP` only wakes a sleeping process
so that it scans again immediately.

## Installation

BFQ requires Python 3.11 or newer. The complete pipeline additionally depends
on Illumina conversion software, GCF workflows, Apptainer/Singularity, and
several bioinformatics command-line tools supplied by the Docker image.

Install the Python application from the repository root with the in-house
`gcf-tools` dependency:

```console
python -m pip install . \
  "gcf-tools @ https://github.com/gcfntnu/gcf-tools/archive/master.zip"
```

The installation provides two commands:

- `bfq` starts the pipeline service.
- `flowcell-manager` manages the processed-flowcell inventory.

The legacy executable names `bfq.py` and `flowcell_manager.py` are not
installed.

## Development and production releases

`bfq-dev` is the default development branch; `master` is the production branch.

1. Create an issue branch from current `bfq-dev` (including when using an issue's
   **Development** section), and target its feature PR at `bfq-dev`. Link the issue
   with `Closes #NUMBER` in the PR description so merging closes it.
2. Run the relevant automated and integration checks, then merge the feature PR.
3. When the tested development changes are ready for production, open a promotion
   PR from `bfq-dev` to `master` with a short summary of the included changes.
   Review and merge that PR deliberately before building the production image.
4. From the updated `master` checkout, build and push with the existing script,
   manually choosing the next `prod-N` tag. For example, **if the previous release
   was `prod-60`**, the next release is:

   ```bash
   bash build-tag-push.sh prod prod-61
   ```

The script builds and immediately pushes `gcfntnu/bfq:prod-61` in this example;
it does not deploy the image. Release numbers are selected manually. No separate
release manifest or coordinated version increment for supporting tools is required.

Production sources are selected explicitly, independently of repository defaults:

| Repository | Production branch |
| --- | --- |
| `gcf-bfq` | `master` |
| `gcf-tools` | `master` |
| `gcf-workflows` | `main` |

Resolving the current production dependency branches for each new image build is
intentional. The test Dockerfile defaults the tools and workflows to `bfq-dev`
and installs BFQ from the local build context.

### Interpreting analysis versions

The **Analysis pipeline** entry in each project QC report records the
`gcf-workflows` commit that produced that analysis. This is the primary analysis
version for project deliverables. BFQ JSON `versions` describes the current
attempt's execution environment; `attempts[].versions` retains that metadata
for earlier attempts.

`rerun --from analysis` executes the complete analysis workflow using the installed
workflow checkout and regenerates project reports before reporting/finalization.
`rerun --from reporting` regenerates sequencer reporting and downstream delivery
products while preserving project analysis results and reports, including their
original workflow commit. A finalization-only rerun likewise preserves those
reports. Consequently, a later reporting/finalization attempt's runtime revision
can differ from the analysis revision correctly recorded in the preserved report.

## Starting BFQ

Mount the configuration directory at `/config` and the operational storage at
the paths declared below. The container entry point is:

```console
bfq
```

Set `BFQ_DEBUG=1` to enable debug logging. Logs from the demultiplexing command
are written separately to `[Paths] logDir`.

## Static configuration

BFQ reads `/config/bcl2fastq.ini` using Python's `ConfigParser`. Section and
option names are case-insensitive, but underscores are significant. Use the
spellings shown here, especially for email keys. A sanitized configuration is
available at `files/bcl2fastq.ini` and can be copied into the mounted
configuration directory as a starting point.

```ini
[Paths]
ekista_baseDir = /instruments/ekista
nova_baseDir = /instruments/nova
outputDir = /bfq/output
logDir = /bfq/log
manager_dir = /flowcellmanager
reportDir = /bfq/reports
analysisDir = /bfq/analysis

[System]
sleepTime = 0.25
minSpace = 50

[Email]
host = smtp.example.org
from_address = bfq-no-reply@example.org
finished_to = sequencing@example.org
error_to = pipeline-errors@example.org, sequencing-oncall@example.org

[Version]
pipeline = 0.3.1

[Commands]
multiqc_options = -f -q --interactive
bcl2fastq_options = --no-lane-splitting -p 32 -r 12 -w 12 -l WARNING
cellranger_mkfastq_options = --qc --jobmode=local --localcores=32 --localmem=55
```

Do not commit the production configuration if it contains real email addresses,
credentials, internal hostnames, or other sensitive site information.

### `[Paths]`

The `[Paths]` section itself is required. Individual path values have code
defaults, but a production deployment should set them explicitly. Instrument
roots must be readable. The output, log, report, and manager directories must
exist and be writable before BFQ starts. `analysisDir` is currently metadata
only and does not have to exist.

| Option | Default | Purpose |
| --- | --- | --- |
| `ekista_baseDir` | `/mnt/seq/ekista` | Root searched for runs from the Ekista-mounted instrument storage. |
| `nova_baseDir` | `/mnt/seq/nova` | Root searched for runs from Nova-mounted instrument storage. |
| `outputDir` | `/mnt/output` | Parent directory for one BFQ output directory per flowcell. |
| `logDir` | `/mnt/logs` | Destination for demultiplexing command logs. |
| `manager_dir` | `/mnt/manager` | Durable flowcell state root (`states/`, `locks/`) and compatibility `flowcells.processed` inventory. |
| `reportDir` | `/mnt/reports` | Destination for `<run-id>.error` reports. |
| `analysisDir` | `/mnt/analysis` | Reserved analysis path recorded in the run configuration snapshot. Current workflow execution uses `TMPDIR`. |

The flowcell-manager inventory must be a CSV named `flowcells.processed` with
these columns:

```csv
project,flowcell_path,timestamp,archived
```

### `[System]`

These options are required; the current code does not provide runtime defaults
when they are absent.

| Option | Unit | Purpose |
| --- | --- | --- |
| `sleepTime` | hours | Delay between discovery scans. Fractional values are accepted. |
| `minSpace` | GiB | Minimum free space required on `outputDir` before a run starts. |

### `[Email]`

All four options are required when the normal notification paths are used.
Error recipients support comma-separated addresses, including display names; empty
entries and duplicate addresses are ignored. All error recipients receive one message.

| Option | Purpose |
| --- | --- |
| `host` | SMTP relay hostname. |
| `from_address` | Sender used for completion and error messages. |
| `finished_to` | Recipient for the processing-complete email and attached reports. |
| `error_to` | Comma-separated recipients for production error notifications. Also retains the existing final archive-completion notification routing. |

Current SMTP handling does not configure authentication or TLS. Access control
must therefore be provided by the deployment environment or relay.

### Production error notifications

`dockerfile-prod` sets `BFQ_ENV=production`; `dockerfile-test` sets `BFQ_ENV=test`.
Only the exact value `production` enables **error** email. Missing, empty, or
unrecognized values suppress it. This setting does not change the existing
processing-complete and successful-finalization email behavior. In particular,
finalization still uses `error_to`; moving it requires a separate routing decision.
Do not override `BFQ_ENV` to production in ordinary test/development runs.

On a pipeline-stage failure, BFQ first writes `<reportDir>/<run-id>.error`, then
records the failure in canonical JSON state, and only then attempts notification.
`write_error_report()` and `send_error_report()` are separate operations;
`errorEmail()` remains a report-only compatibility alias. The email subject contains
the run ID and stage. Its plain-text body includes the exception, UTC timestamp,
host, and absolute report path, with the saved report attached as `text/plain`.
Captured command diagnostics remain included in the report. SMTP uses the configured
`host` and `from_address`, with a 30-second socket timeout and no added TLS or auth.

Notification protection lives in `<manager_dir>/states/<run-id>.json`, not the
output tree. A signature covers the stage, qualified exception type/message, and
captured command output when present; report timestamps and traceback line numbers
are excluded. Identical failures remain suppressed across daemon restarts and
explicit reruns. A changed failure starts a new notification, and successful run
completion clears the failure/notification record. Attempts retain failure history.

Delivery is **at most once per failure**: BFQ atomically saves `attempted_at` under
the per-run lock before opening SMTP. Success sets `notified` and `notified_at`;
SMTP errors (including partial recipient refusal) are logged and saved as
`delivery_error`, and do not prevent run-context reset or processing other runs.
The original report and pipeline failure remain intact. Relay failure or a crash
after claiming delivery is not automatically retried, since the relay may already
have accepted the mail. `attempted_at` without `notified` indicates an unsuccessful
or interrupted attempt; inspect the report and relay logs. Changed failures or a
successful run followed by another failure enable another notification.

If report writing or durable failure/notification recording fails, no email is
sent. Empty `error_to` suppresses delivery with a log message. Non-production
suppression does not consume a delivery attempt. Discovery-only diagnostic reports
for ambiguous output remain local: they have no initialized pipeline-stage state.

#### Server-side checks for issue #103

Automated tests mock SMTP; they cannot validate the deployed image, relay, or mailbox.
Before rollout:

1. Build the test image from the issue branch and confirm `BFQ_ENV=test` inside it.
   Trigger a controlled stage failure: verify the `.error` report, failed JSON state,
   logged email suppression, and absence of an error-email connection in relay logs.
2. In a controlled production-mode deployment using approved test recipients, verify
   one message reaches every comma-separated `error_to` address, the subject/stage
   and body context are correct, and the attachment matches the saved report.
   `dockerfile-prod` installs `master`; testing unmerged code requires a test build
   containing this branch with an explicit production-mode override.
3. Restart the daemon and explicitly rerun the same failing stage: no duplicate mail.
   Change the failure and verify one new notification. Complete a successful run,
   then reproduce the original failure and verify it is notified again.
4. Simulate an unavailable/refusing relay and verify the report survives, state records
   the original failure plus delivery error, and the daemon can process another run.
   Retrying the same failure must not repeatedly contact the relay.

### `[Version]`

This section is optional. BFQ preserves all entries in the configuration
snapshot written to the run output. It is suitable for site-managed release or
deployment identifiers; the current pipeline does not branch on these values.

### `[Commands]`

| Option | When required | Purpose |
| --- | --- | --- |
| `multiqc_options` | Always | Options passed to the sequencer-level MultiQC invocation. |
| `bcl2fastq_options` | With `FORCE_BCL2FASTQ` | Options passed to legacy bcl2fastq. |
| `cellranger_mkfastq_options` | Any supported 10x run | Shared local execution options for the 10x mkfastq commands. |

BFQ selects the executable internally and runs it through an Apptainer wrapper.
Image references come from the checked-out `gcf-workflows/docker.config`.
Standard runs use `bcl-convert`; `FORCE_BCL2FASTQ` selects `bcl2fastq`. The 10x
wrapper is selected by an exact `Libprep` match in `makeFastq.py`. Legacy INI
executable entries are accepted as extra configuration values but ignored.

## Per-flowcell input

Each candidate flowcell directory must contain:

- The instrument-specific completion marker listed below.
- At least one file matching `SampleSheet*.csv` with a non-empty
  `[CustomOptions]` section.
- At least one file matching `*Sample-Submission-Form*.xlsx`.

For a queued state-backed run, BFQ prefers preserved output-side run inputs. If they are unavailable, it copies the selected files from the instrument run into the output directory as `SampleSheet.csv` and `Sample-Submission-Form.xlsx`. Pre-existing output is first subject to the discovery precedence above; ambiguous output is never guessed into a processing stage.

If several sample sheets exist, BFQ uses the first one returned by the
filesystem that contains non-empty custom options. The ordering is not
guaranteed, so operators should leave only the intended active sheet in a run.

### Completion markers

| Instrument identifier in run name | Required marker |
| --- | --- |
| `SN7001334` | `ImageAnalysis_Netcopy_complete.txt` |
| `NB501038` | `RunCompletionStatus.xml` |
| `M026575`, `M03942`, `M05617`, `M71102` | `ImageAnalysis_Netcopy_complete.txt` |
| `K00251` | `SequencingComplete.txt` |
| `A01990`, `MN00686` | `CopyComplete.txt` |

The run directory name becomes the run ID. Its location under `nova_baseDir` or
`ekista_baseDir` is recorded as the instrument source (`nova` or `ekista`). A
path outside both configured roots is recorded as `unknown`.

## Sample-sheet custom options

BFQ scans the CSV for a case-insensitive `[CustomOptions]` header. It reads the
first two columns until the next section header and retains every non-empty key.
Known keys are promoted into typed per-run state:

```csv
[CustomOptions]
Libprep,Illumina DNA Prep
User,operator-name
Rerun,False
SensitiveData,False
```

| Key | Required | Meaning |
| --- | --- | --- |
| `Libprep` | Yes | Library-preparation name. It selects the workflow through `/opt/gcf-workflows/libprep.config` and selects supported 10x demultiplexing commands. |
| `User` | Operationally expected | Operator or submitter shown in completion reporting. |
| `Rerun` | No | Parsed as a boolean, but automatic rerun detection is currently disabled. Use `flowcell-manager rerun` for an explicit rerun. |
| `SensitiveData` | No | When true, FASTQ and QC 7-Zip archives receive generated passwords written to `encryption.<name>` files in the flowcell output. |

`Rerun` and `SensitiveData` treat `true`, `1`, and `yes` as true, ignoring
letter case. Other custom keys are preserved in the run snapshot for downstream
or future use but are not interpreted by BFQ itself.

Workflow lookup is case-insensitive. BFQ first tries the exact `Libprep` value,
then the same value with ` SE` and ` PE` appended. A missing or unmatched
`libprep.config` results in pipeline value `UNKNOWN` and will normally prevent
successful workflow execution.

## Configuration model

The Python configuration model has three layers:

- `StaticConfig` holds immutable paths plus the `[System]`, `[Email]`,
  `[Version]`, and `[Commands]` mappings loaded at startup.
- `RunContext` holds the current run ID, flowcell path, source, selected input
  files, library preparation, workflow, user, rerun flag, sensitive-data flag,
  and all custom options.
- `PipelineConfig` is the process-wide singleton combining both layers and
  deriving `outputDir/run-id` as the current output path.

At the end of post-processing, BFQ writes a human-readable YAML snapshot to
`<outputDir>/<run-id>/bcl2fastq.ini`. Despite its historical `.ini` suffix,
this generated file is YAML and is not an input configuration file.

## Versioned flowcell state

BFQ stores one canonical schema-v1 JSON record per state-backed run at
`<manager_dir>/states/<run-id>.json`. State is deliberately outside the run
output tree so archival and restart cleanup cannot erase processing history.
Writes are protected by per-run advisory locks and published with atomic file
replacement. BFQ also takes an execution lease while a run is active.

The JSON record is authoritative for state-backed runs and records the current
stage, stage timestamps, attempts, projects, BFQ and workflow revisions when
available, demultiplexer/version, failures and report paths, notification fields,
operator restart requests, and archive state. Previous attempts remain in
history rather than being overwritten.

BFQ no longer reads or writes the legacy `bcl.done`, `files.renamed`,
`analysis.made`, or `fastq.made` files. Existing copies may remain on disk
but have no effect on discovery or restart decisions.

Discovery uses this precedence:

1. If a state JSON exists, it is authoritative even if the compatibility
   inventory also contains the run.
2. Without JSON, an entry in `flowcells.processed` protects a completed or
   archived legacy run from automatic reprocessing.
3. Without either, a complete restored FASTQ tree with BFQ run inputs is
   initialized as `restored_legacy_fastq` and queued from `analysis`.
4. With no state, inventory entry, or output directory, BFQ creates state
   immediately and queues a new run from `demultiplexing`.
5. Other pre-existing output is ambiguous. BFQ writes an actionable
   `<reportDir>/<run-id>.error` report and requires an explicit
   `flowcell-manager initialize` decision.

A stage left `running` without an execution lease is treated as interrupted on
the next inspection. It is not silently resumed; the operator must queue a safe
restart boundary with `flowcell-manager rerun`.

## Waking and restarting

To wake BFQ during its configured sleep interval:

```console
kill -HUP <pid>
```

`SIGHUP` does not reload `/config/bcl2fastq.ini`. Restart the process or
container to apply static configuration changes. Python processing modules are
reloaded between scans, but an installed package update should likewise be
followed by a process restart.

## Flowcell manager

`flowcell-manager` is the operator-facing interface for both state-backed runs
and legacy inventory-only runs. It changes state and performs defined cleanup;
it does not execute pipeline stages. The normal BFQ daemon executes queued work.

Common commands are:

```console
flowcell-manager list
flowcell-manager list --status failed
flowcell-manager list --stage analysis
flowcell-manager show RUN_ID
flowcell-manager status RUN_ID
flowcell-manager rerun RUN_ID --from demultiplexing
flowcell-manager rerun RUN_ID --from analysis --reason "repeat workflow"
flowcell-manager rerun RUN_ID --from reporting
flowcell-manager rerun RUN_ID --from finalization
flowcell-manager initialize RUN_ID --from demultiplexing
flowcell-manager initialize RUN_ID --from analysis
flowcell-manager archive RUN_ID
flowcell-manager list-processed
```

Destructive operations show their cleanup plan and prompt by default. Use
`--dry-run` to preview without changing files or state, and `--force` only
for deliberate non-interactive operation. `--reason` is retained in state.
`--refresh-inputs` explicitly recopies the sample sheet and submission form
from the instrument source; without it, output-side run inputs are preserved.

The supported restart boundaries invalidate these products:

| Restart boundary | Preserved | Invalidated |
| --- | --- | --- |
| `demultiplexing` | `SampleSheet.csv`, `Sample-Submission-Form.xlsx` | FASTQs and all downstream products |
| `analysis` | FASTQs, FASTQ checksums and run inputs | workflow/QC output, reports, archives, archive checksums, matching workflow work directories |
| `reporting` | FASTQs, FASTQ checksums, workflow results, project HTML reports and project MultiQC configurations | sequencer reports/metrics, aggregate MultiQC configuration, archives, completion products |
| `finalization` | FASTQs, FASTQ checksums, workflow results, reports | delivery archives and archive checksums |

FASTQ manifests (`md5sum_<project>_fastq.txt`) belong to demultiplexing.
BFQ generates them after FASTQ renaming, before marking demultiplexing complete.
A checksum failure therefore fails demultiplexing and prevents analysis from
starting. Each manifest is written to a temporary file and atomically replaced
only after every checksum has been written and flushed to disk.

Analysis, reporting and finalization reruns preserve complete FASTQ manifests
without re-reading FASTQ contents. A demultiplexing rerun invalidates them and
regenerates checksums for the newly produced FASTQs. Archive checksums remain
part of finalization and are invalidated whenever their archives are rebuilt.

For older runs, restored FASTQs or an interrupted legacy checksum write, queue
the desired downstream boundary as usual, for example:

```bash
flowcell-manager rerun RUN_ID --from analysis
```

Before the requested stage, BFQ checks each project's manifest syntax and exact
coverage of the current FASTQ filenames. Missing, malformed, duplicate or
incomplete entries cause that project's manifest to be rebuilt once; complete
manifests are reused unchanged. This compatibility repair does not rerun BCL
conversion and runs under the flowcell execution lease. If repair fails, BFQ
stops at the requested restart boundary and reports a FASTQ checksum error;
retry that boundary after correcting the underlying problem. Temporary files
left by an abrupt interruption are never treated as completed manifests.
This check validates manifest completeness, not file integrity: it does not
detect manual edits to FASTQ contents with unchanged filenames.

A state record is put into `preparing` before restart cleanup begins and is
queued only after cleanup succeeds. This prevents partial destructive work from
being mistaken for a runnable state. Active runs are refused unless `--force`
is explicitly supplied.

`archive` remains separate from pipeline finalization. It removes delivery data
from the output tree while preserving canonical JSON state. The compatibility
`flowcells.processed` inventory is still maintained for legacy protection,
project-to-flowcell search, and external consumers; successful state-backed
completion updates it without duplicate project/run rows.

BFQ never extracts legacy `.7za` archives automatically. To resume a legacy
run from restored FASTQs, extract them explicitly and retain
`SampleSheet.csv` plus `Sample-Submission-Form.xlsx`; readable FASTQs must
exist in recognized `GCF-*` project directories.

## Environment controls

| Variable | Effect |
| --- | --- |
| `BFQ_DEBUG` | Enables debug logging when set. |
| `BFQ_TEST` | Enables compatibility handling for test flowcells generated with bcl2fastq while the image defaults to bcl-convert. |
| `FORCE_BCL2FASTQ` | Uses legacy bcl2fastq instead of bcl-convert for non-10x runs. |
| `GCF_WORKFLOWS_DOCKER_CONFIG` | Overrides the default `/opt/gcf-workflows/docker.config` image mapping. |
| `BFQ_APPTAINER_COMMAND` | Overrides the default `apptainer` executable, for example with `singularity`. |
| `APPTAINER_WRITABLE_TMPFS` | Makes container filesystems temporarily writable. Set to `true` by the base image. |
| `SINGULARITY_WRITABLE_TMPFS` | Backwards-compatible equivalent of `APPTAINER_WRITABLE_TMPFS`. |
| `TMPDIR` | Root for per-project workflow work directories and QC archive sources. Set by the base image. |

Apptainer/Singularity cache, temporary-directory, and bind-path variables are
also provided by the base image for the downstream Snakemake workflows. The
writable tmpfs overlay is discarded after each container command; output that
must persist still needs to be written to a bind-mounted path.

The demultiplexer name and image tag (or digest) are recorded in the flowcell state JSON. They are resolved from the active `gcf-workflows/docker.config`; BFQ no longer maintains separate version environment variables for these tools.

## Output overview

Each run is written below `<outputDir>/<run-id>`. Important products include:

- Demultiplexed project FASTQs and `Undetermined` FASTQs.
- `Stats`, `InterOp`, `RunInfo.xml`, and `RunParameters.xml`.
- Per-project sample information, MultiQC reports, archives, and archive MD5s.
- `QC_<project>` workflow outputs and `QC_<project>_<date>.7za` archives.
- A static/run configuration snapshot. Durable processing state is stored under `manager_dir`, not in this output tree.
- `encryption.*` password files when `SensitiveData` is true.

The exact project analysis content is defined by the selected workflow in
`gcf-workflows` rather than by BFQ itself.

## Development and testing

Create an editable environment from the repository root:

```console
python -m pip install -e ".[dev]" \
  "gcf-tools @ https://github.com/gcfntnu/gcf-tools/archive/master.zip"
```

Run the local checks:

```console
python -m pip check
python -m pytest
ruff check .
ruff format --check .
python -m build
```

Unit tests use temporary files and do not require sequencing data, mounted
instrument storage, external services, or bioinformatics applications.

### Server-side integration testing for state management

Issue #104 changes filesystem authority and destructive restart behavior, so the
automated suite is necessary but not sufficient for production rollout. Before
deployment, exercise the branch on the BFQ test server with real mounts,
Apptainer images, workflow work directories, SMTP/report paths, and the
production-style `manager_dir`.

The integration pass should verify at least:

- a completely new standard Illumina run creates JSON state before output work,
  completes every stage, and never creates legacy restart markers;
- a state-backed run already present in `flowcells.processed` follows JSON
  authority when an explicit rerun is queued;
- inventory-only historical runs remain skipped and searchable;
- explicitly restored legacy FASTQs bootstrap at `analysis`, while incomplete
  or ambiguous restored output is refused;
- `rerun --dry-run` and each of the four restart boundaries invalidate exactly
  the displayed products, preserve curated inputs/FASTQs as documented, and
  `--refresh-inputs` replaces inputs only when requested;
- killing BFQ during each stage leaves a recoverable interrupted state rather
  than a false completion, and two daemon processes cannot execute the same run;
- an unavailable or unwritable `manager_dir` prevents any untracked
  processing;
- completion updates `flowcells.processed` without duplicate project/run rows,
  and archive/rerun state survives output cleanup;
- standard bcl-convert, forced bcl2fastq, supported 10x demultiplexers, and every
  configured downstream workflow complete successfully from a clean start.

Run the full supported-workflow matrix as an overnight integration test before
production deployment. These checks intentionally are not simulated by the unit
test suite because they depend on server mounts, real external tools, and
deployment configuration.

For a version supporting Illumina bcl2fastq v1, see the historical
[`bcl2fastqV1`](https://github.com/maxplanck-ie/bcl2fastq_pipeline/tree/bcl2fastqV1)
branch of the upstream project.
