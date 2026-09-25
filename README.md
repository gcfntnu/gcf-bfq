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
3. Ignore flowcells already present in the flowcell-manager inventory.
4. Require a sample sheet with `[CustomOptions]` and a sample submission form.
5. Demultiplex, run the selected Snakemake workflow, and generate MultiQC output.
6. Create delivery archives and checksums.
7. Write completion markers and add each project to the inventory.
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
error_to = pipeline-errors@example.org

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
| `manager_dir` | `/mnt/manager` | Directory containing `flowcells.processed`. |
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
Recipient values are passed to the local SMTP server as configured strings.

| Option | Purpose |
| --- | --- |
| `host` | SMTP relay hostname. |
| `from_address` | Sender used for completion and error messages. |
| `finished_to` | Recipient for the processing-complete email and attached reports. |
| `error_to` | Recipient for the final archive-completion message. Runtime error details are currently written below `reportDir`. |

Current SMTP handling does not configure authentication or TLS. Access control
must therefore be provided by the deployment environment or relay.

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

If an output directory already exists, BFQ first looks there for the sample
sheet and submission form. Otherwise it copies the selected files from the
instrument run into the new output directory as `SampleSheet.csv` and
`Sample-Submission-Form.xlsx`.

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

## Processing and restart markers

BFQ uses empty marker files to skip completed stages when an interrupted run is
encountered again:

| Marker | Meaning |
| --- | --- |
| `bcl.done` | Demultiplexing completed. The file records the converter and version. |
| `files.renamed` | FASTQ renaming completed. |
| `analysis.made` | Project-level Snakemake analysis completed. |
| `fastq.made` | The full BFQ process completed. |

The authoritative discovery-time processed check is currently the
flowcell-manager inventory, not `fastq.made`. On successful completion BFQ
writes `fastq.made` and adds one inventory row per detected GCF project.

For an interrupted run that has not yet entered the inventory, deleting an
individual marker causes its corresponding stage to run again. Once the run is
present in the inventory, discovery skips it before inspecting these markers.
For a full, explicit rerun, use the flowcell manager as described below.

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

Common inventory operations are:

```console
flowcell-manager list
flowcell-manager list-processed
flowcell-manager add GCF-2026-001 /bfq/output/RUN_ID 2026-09-23T12:00:00
flowcell-manager archive /bfq/output/RUN_ID
flowcell-manager rerun /bfq/output/RUN_ID
```

`archive` removes project directories, BAM files, FASTQ files, and 7-Zip
archives while retaining the flowcell directory and marking the inventory row
as archived. `rerun` removes the entire flowcell output directory and its
inventory rows. Both commands prompt unless `--force` is supplied. Review the
displayed target carefully before confirming either destructive operation.

## Environment controls

| Variable | Effect |
| --- | --- |
| `BFQ_DEBUG` | Enables debug logging when set. |
| `BFQ_TEST` | Enables compatibility handling for test flowcells generated with bcl2fastq while the image defaults to bcl-convert. |
| `FORCE_BCL2FASTQ` | Uses legacy bcl2fastq instead of bcl-convert for non-10x runs. |
| `GCF_WORKFLOWS_DOCKER_CONFIG` | Overrides the default `/opt/gcf-workflows/docker.config` image mapping. |
| `BFQ_APPTAINER_COMMAND` | Overrides the default `apptainer` executable, for example with `singularity`. |
| `TMPDIR` | Root for per-project workflow work directories and QC archive sources. Set by the base image. |
| `BCL_CONVERT_VERSION`, `BCL2FASTQ_VERSION`, `CR_VERSION` | Version strings recorded in `bcl.done`. Set by the image. |

Apptainer/Singularity cache, temporary-directory, and bind-path variables are
also provided by the base image for the downstream Snakemake workflows.

## Output overview

Each run is written below `<outputDir>/<run-id>`. Important products include:

- Demultiplexed project FASTQs and `Undetermined` FASTQs.
- `Stats`, `InterOp`, `RunInfo.xml`, and `RunParameters.xml`.
- Per-project sample information, MultiQC reports, archives, and archive MD5s.
- `QC_<project>` workflow outputs and `QC_<project>_<date>.7za` archives.
- A static/run configuration snapshot and the restart markers listed above.
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

For a version supporting Illumina bcl2fastq v1, see the historical
[`bcl2fastqV1`](https://github.com/maxplanck-ie/bcl2fastq_pipeline/tree/bcl2fastqV1)
branch of the upstream project.
