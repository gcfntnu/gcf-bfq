# BFQ architecture and compatibility map

This is a source-navigation guide for bounded changes. The [README](../README.md)
remains the operator reference; [CONTRIBUTING](../CONTRIBUTING.md) describes the
development and review agreement. Paths below are relative to
`bcl2fastq_pipeline/bcl2fastq_pipeline/` unless stated otherwise.

## Responsibilities and data flow

The installed `bfq` command enters `entrypoint.py`, then `cli.py`. Configuration
loads once from `/config/bcl2fastq.ini`. Discovery selects eligible sequencing
runs; canonical state controls which attempts are queued. Under a per-run
execution lease, orchestration progresses through demultiplexing, analysis,
reporting and finalization. Completion creates durable notification intents;
delivery failures remain separate from processing failures.

| Area | Source | Responsibility / first tests to read |
| --- | --- | --- |
| Configuration | `config.py` | `StaticConfig` plus mutable `RunContext`, held by a process-wide `PipelineConfig` singleton. `tests/test_config.py`. |
| Discovery and inventory | `findFlowCells.py`, `cli.py: candidate_flowcells` | Run selection, curated input preparation, explicit restored-FASTQ/legacy handling, compatibility inventory updates. `tests/test_state_integration.py`. |
| Persistent state and locking | `state.py` | Schema validation, atomic state updates, per-run execution lease, attempts, restart state and notification records. `tests/test_state.py`, `tests/test_notification_state.py`. |
| Input validation | `preflight.py` | Resolve effective inputs and call shared `gcf-tools` validation before demultiplexing/analysis. `tests/test_input_preflight.py`. |
| Stage orchestration | `cli.py: _run_state_backed_flowcell` | Run the four restart boundaries; commit processing outcomes separately from notification delivery. `tests/test_notification_integration.py`, `tests/test_early_qc_integration.py`. |
| Demultiplexing and corrections | `makeFastq.py`, `index_corrections.py`, `index_sequences.py` | Conversion selection, FASTQ naming and explicit index-orientation corrections. `tests/test_index_cli.py`, `tests/test_index_corrections.py`. |
| Analysis and finalization | `afterFastq.py`, `workflow_config.py` | Capture library-prep configuration, invoke configmaker/Snakemake, copy reports, create archives/checksums. `tests/test_after_fastq.py`, `tests/test_workflow_config.py`, `tests/test_fastq_checksums.py`. |
| Operator CLI and cleanup | `../flowcell_manager/flowcell_manager.py`, `state.py: cleanup_plan`, `fastq_cleanup.py` | Shared `fm`/`flowcell-manager`, validation and rerun previews, explicit initialization, search/status, cleanup and archive operations. `tests/test_state_integration.py`, `tests/test_flowcell_search.py`, `tests/test_fastq_cleanup.py`. |
| Analysis retention | `analysis_snapshots.py` | Retain the actual successful analysis workdir, stage/commit/recover publication, preserve the previous successful archive. `tests/test_analysis_snapshots.py`, `tests/test_snapshot_integration.py`. |
| QC and notifications | `sequencing_qc.py`, `sequencing_delivery.py`, `analysis_qc.py`, `notifications.py`, `notification_delivery.py`, `misc.py` | Early sequencing report, analysis summary, durable delivery/retry, production error reporting. `tests/test_sequencing_qc.py`, `tests/test_analysis_qc.py`, `tests/test_notification_mail.py`, `tests/test_error_reporting.py`. |
| Timing | `processing_times.py` | Successful processing measurements and accumulated reporting. `tests/test_processing_times.py`, `tests/test_processing_time_integration.py`. |
| External command wrappers | `containers.py` | Resolve images from workflow-owned `docker.config` and execute via Apptainer/Singularity. `tests/test_containers.py`. |
| Packaging and entry points | Root `pyproject.toml`, `entrypoint.py`, `version.py` | Distribution metadata and installed commands; avoid neighboring-script import shadowing. `tests/test_packaging.py`, `tests/test_versions.py`. |

FASTQ checksums belong to demultiplexing. Early sequencing QC/email happens after
conversion and rename, before FASTQ checksums and analysis. Analysis invokes the
scientific workflow and produces project reports; the subsequent reporting stage
does not rerun scientific analysis. Finalization creates delivery archives and
their checksums, records completion, and publishes retained analysis snapshots.

## Companion repositories and runtime paths

| Owner | Contract |
| --- | --- |
| `gcf-bfq` | Scheduling, lifecycle/state, input selection, operational CLI, delivery and recovery. |
| `gcf-tools` | Shared SampleSheet/workbook validation and domain parsing, configmaker, shared libprep selection API. BFQ must not independently reimplement those rules. |
| `gcf-workflows` | Scientific workflows, authoritative `libprep.config`, tool-image `docker.config`; copied into the actual analysis workdir before execution. |
| Production Docker environment | Illumina/bioinformatics applications, `/opt/conda` configmaker, workflow checkout, Apptainer/Singularity and operational mounts. Local Python checks do not reproduce this environment. |

Production intentionally selects BFQ `master`, tools `master` and workflows
`main` when building a new image. Development checks instead use an immutable
companion baseline, with an explicit local override for coordinated changes;
see [development.md](development.md). Updating that baseline does not alter the
production branch policy. The installed package currently requires
`gcf-tools>=0.3.0`; older feature documents describe their original companion PRs
and must not be used as current dependency pins.

These are operational paths, **not local development targets**:

| Location | Meaning |
| --- | --- |
| `/config/bcl2fastq.ini` | Runtime site configuration; static values load once at process start. |
| Configured instrument roots | Sequencer inputs; ordinary discovery requires completed copying and BFQ custom options. |
| `<outputDir>/<run-id>/` | Effective curated inputs, project FASTQs, Stats, reports, delivery archives and checksums. |
| `<manager_dir>/states/`, `<manager_dir>/locks/` | Canonical JSON state and lifecycle locks; `flowcells.processed` remains the compatibility CSV inventory. |
| `<TMPDIR>/<project>_<run-date>/` | Actual analysis working tree, including generated config, copied workflows and `data/`. Same-project/date workdirs can collide; developer isolation must not reuse a server's TMPDIR. |
| `<outputDir>/<run-id>/provenance/` | Retained successful analysis snapshots and staged publication recovery. |
| Configured `logDir`, `reportDir` | External-command logs and failure/preflight reports. |

Tests construct temporary paths and inject configuration. The CLI deliberately
retains production configuration authority; no production simulation mode is
required for local testing. Container-config overrides and test mocks are not
permission to weaken the authoritative libprep path.

## Compatibility contracts

Treat changes to these behaviors as explicit interface changes needing suitable
regression coverage and review, even when a refactor appears internal:

- **Commands:** `bfq`, `flowcell-manager` and `fm` remain installed. `fm` is the
  exact same function as `flowcell-manager`. Existing subcommands, flags, exit
  behavior and supported run-ID/path forms must remain compatible. The old
  `bfq.py` and `flowcell_manager.py` names are not installed aliases.
- **Inputs and outputs:** retain accepted SampleSheet/workbook structures, IDs,
  FASTQ/project layouts, report/archive naming, inventory columns and configmaker
  interfaces. Do not silently repair metadata or change a scientific parameter
  as part of tooling/maintenance work.
- **Input ownership:** each canonical output-side input wins independently,
  including a malformed one. Only explicit `--refresh-inputs` chooses both
  instrument copies. Validation happens before destructive rerun preparation and
  again at execution; `--force` bypasses confirmation, not validation or leases.
- **State and legacy runs:** state is authoritative; inventory-only legacy runs
  remain protected. Legacy marker files do not drive new execution. Preserve
  loading of existing schema-v1 state and missing optional fields; extend the
  existing state mechanisms rather than introducing competing completion markers.
- **Restart boundaries:** retain `demultiplexing`, `analysis`, `reporting` and
  `finalization`. An analysis rerun preserves FASTQs and their valid MD5
  manifests, then rebuilds downstream analysis/results. Reporting/finalization
  preserve completed project analysis and its recorded workflow identity. Use
  the actual cleanup plan and CLI preview rather than reproducing a glob list in
  new code. `clean-fastqs` and `archive` are distinct explicit operator actions.
- **Configuration:** `/opt/gcf-workflows/libprep.config` is authoritative; retired
  `BFQ_LIBPREP_CONFIG` does not override it. New demultiplexing/analysis attempts
  capture bytes once and pass the same selection to each project. Downstream
  retries use the completed analysis workflow rather than current installed
  libprep settings. A wake-up signal does not reload static site configuration.
- **Snapshots:** retain the actual successful workdir with only top-level
  `data/` excluded. Failed attempts and cleanup preserve the previous successful
  snapshot. Respect staged publication and recovery; never reconstruct an older
  analysis from today's installed workflow checkout.
- **Notifications:** SMTP failure does not undo successful processing. Retry
  existing durable intents through the established state mechanisms; uncertain
  delivery needs explicit handling. Local tests mock mail. `BFQ_ENV=test` is not
  a general mail blocker; error mail specifically uses `BFQ_ENV=production`.
- **Operator decisions:** early QC informs manual decisions; it is not an
  automatic acceptance threshold or pause. FASTQ cleanup relies on the documented
  delivered-data/archive-presence policy, not a new mandatory checksum reread.

## Choosing verification

Start with the tests listed above and the common
[local checks](development.md). Fixtures in `tests/test_state_integration.py`
already supply temporary configuration, metadata and tiny FASTQs; several tests
reuse them. Keep any extraction of those helpers bounded. Dedicated reusable
operational scenarios are follow-up [#140](https://github.com/gcfntnu/gcf-bfq/issues/140),
not a claim made by the current local suite.

Use the relevant existing guide to specify the remaining manual integration:

| Changed boundary | Operational reference / manual checks |
| --- | --- |
| Input selection and validation | [Input preflight](input-preflight.md), [coordinated preflight checks](preflight-integration-tests.md) |
| Libprep selection and configmaker propagation | [Library-prep authority](libprep-configuration.md); use current dependency/build policy, not its historical companion branch |
| Restart/state and notification recovery | [README](../README.md), [notification recovery](notification-recovery.md) |
| Index corrections and conversion | [Index orientation checks](index-orientation-integration-tests.md) |
| Early QC and notifications | [Early sequencing QC checks](early-sequencing-qc-integration-tests.md) |
| Snapshot lifecycle | [Analysis snapshot checks](analysis-snapshot-integration-tests.md) |
| FASTQ removal/archive prerequisites | [Cleanup checks](clean-fastqs-integration-tests.md) |
| Search and CLI inventory | [Search checks](flowcell-search-integration-tests.md) |
| Processing-time reporting | [Timing checks](processing-time-integration-tests.md) |

Real BCL conversion, container bindings, scientific workflow outputs, actual
MultiQC compatibility, SMTP delivery and deployed mounts need the relevant server
checks. Tests that replace these operations prove BFQ's surrounding logic only.
